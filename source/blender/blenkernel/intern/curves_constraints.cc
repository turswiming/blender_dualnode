/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "BKE_curves_constraints.hh"

#include "BLI_index_mask_ops.hh"

#include "BKE_bvhutils.h"
#include "BKE_mesh_runtime.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "PIL_time.h"

namespace blender::bke::curves {

using blender::bke::CurvesGeometry;
using threading::EnumerableThreadSpecific;

static const int curves_grain_size = 64;

const ConstraintSolver::Params &ConstraintSolver::params() const
{
  return params_;
}

Span<float> ConstraintSolver::segment_lengths() const
{
  return segment_lengths_cu_;
}

const ConstraintSolver::Result &ConstraintSolver::result() const
{
  return result_;
}

void ConstraintSolver::clear_result()
{
  result_ = Result{};
}

void ConstraintSolver::initialize(const Params &params,
                                  const CurvesGeometry &curves,
                                  IndexMask curve_selection)
{
  params_ = params;

  if (params_.use_length_constraints) {
    segment_lengths_cu_.reinitialize(curves.points_num());

    const Span<float3> positions_cu = curves.positions();
    threading::parallel_for(curve_selection.index_range(), 256, [&](const IndexRange range) {
      for (const int curve_i : curve_selection.slice(range)) {
        const IndexRange points = curves.points_for_curve(curve_i).drop_back(1);

        float length_cu = 0.0f, prev_length_cu;
        for (const int point_i : points) {
          const float3 &p1_cu = positions_cu[point_i];
          const float3 &p2_cu = positions_cu[point_i + 1];
          prev_length_cu = length_cu;
          length_cu = math::distance(p1_cu, p2_cu);
          segment_lengths_cu_[point_i] = length_cu;
        }
      }
    });
  }
  else {
    segment_lengths_cu_.reinitialize(0);
  }

  if (params_.use_collision_constraints) {
    contacts_num_.reinitialize(curves.points_num());
    contacts_.reinitialize(params_.max_contacts_per_point * curves.points_num());
  }
  else {
    contacts_num_.reinitialize(0);
    contacts_.reinitialize(0);
  }
}

void ConstraintSolver::step_curves(CurvesGeometry &curves,
                                   const Mesh *surface,
                                   const CurvesSurfaceTransforms &transforms,
                                   const Span<float3> start_positions,
                                   VArray<int> changed_curves)
{
  double step_start = PIL_check_seconds_timer();

  MutableSpan<float3> positions = curves.positions_for_write();

  const float max_substep_travel_distance = params_.max_travel_distance / params_.substep_count;
  const float max_collision_distance = 2.0f * max_substep_travel_distance;
  const float max_travel_distance_sq = params_.max_travel_distance * params_.max_travel_distance;

  const bool clamp_travel = true;

  /* Compute position delta per substep ("velocity") */
  Array<float3> delta_substep(curves.points_num());
  threading::parallel_for(
      changed_curves.index_range(), curves_grain_size, [&](const IndexRange range) {
        for (const int i : range) {
          const int curve_i = changed_curves[i];
          const IndexRange points = curves.points_for_curve(curve_i);
          for (const int point_i : points) {
            const float3 delta_step = positions[point_i] - start_positions[point_i];
            positions[point_i] = start_positions[point_i];

            const float travel_sq = len_squared_v3(delta_step);
            if (travel_sq <= max_travel_distance_sq) {
              delta_substep[point_i] = delta_step / params_.substep_count;
              continue;
            }

            result_.max_travel_exceeded = true;
            if (clamp_travel) {
              delta_substep[point_i] = math::normalize(delta_step) * max_substep_travel_distance;
            }
            else {
              delta_substep[point_i] = delta_step / params_.substep_count;
            }
          }
        }
      });

  for (const int substep : IndexRange(params_.substep_count)) {
    /* Set unconstrained position: x <- x + v*dt */
    threading::parallel_for(
        changed_curves.index_range(), curves_grain_size, [&](const IndexRange range) {
          for (const int i : range) {
            const int curve_i = changed_curves[i];
            const IndexRange points = curves.points_for_curve(curve_i);
            for (const int point_i : points) {
              positions[point_i] += delta_substep[point_i];
            }
          }
        });

    if (params_.use_collision_constraints) {
      find_contact_points(curves, surface, transforms, max_collision_distance, changed_curves);
    }
    solve_constraints(curves, changed_curves);
  }

  result_.timing.step_total += PIL_check_seconds_timer() - step_start;
}

void ConstraintSolver::find_contact_points(const CurvesGeometry &curves,
                                           const Mesh *surface,
                                           const CurvesSurfaceTransforms &transforms,
                                           const float max_dist,
                                           VArray<int> changed_curves)
{
  /* Should be set when initializing constraints */
  BLI_assert(contacts_num_.size() == curves.points_num());
  BLI_assert(contacts_.size() == curves.points_num() * params_.max_contacts_per_point);

  contacts_num_.fill(0);

  if (surface == nullptr) {
    return;
  }

  const float curves_to_surface_scale = mat4_to_scale(transforms.curves_to_surface.ptr());
  const float surface_to_curves_scale = mat4_to_scale(transforms.curves_to_surface.ptr());
  const float max_dist_su = curves_to_surface_scale * max_dist;
  const float max_dist_sq_su = max_dist_su * max_dist_su;

  Span<MLoopTri> surface_looptris = {BKE_mesh_runtime_looptri_ensure(surface),
                                     BKE_mesh_runtime_looptri_len(surface)};

  /** BVH tree of the surface mesh for finding collisions. */
  double build_bvh_start = PIL_check_seconds_timer();
  BVHTreeFromMesh surface_bvh;
  BKE_bvhtree_from_mesh_get(&surface_bvh, surface, BVHTREE_FROM_LOOPTRI, 2);
  BLI_SCOPED_DEFER([&]() { free_bvhtree_from_mesh(&surface_bvh); });
  if (!surface_bvh.cached) {
    result_.timing.build_bvh += PIL_check_seconds_timer() - build_bvh_start;
  }

  double find_contacts_start = PIL_check_seconds_timer();
  threading::parallel_for(
      changed_curves.index_range(), curves_grain_size, [&](const IndexRange range) {
        for (const int i : range) {
          const int curve_i = changed_curves[i];
          /* First point is anchored to the surface, ignore collisions. */
          const IndexRange points = params_.use_root_constraints ?
                                        curves.points_for_curve(curve_i).drop_front(1) :
                                        curves.points_for_curve(curve_i);
          for (const int point_i : points) {
            const float3 &new_p = curves.positions()[point_i];

            MutableSpan<Contact> contacts(contacts_.begin() +
                                              params_.max_contacts_per_point * point_i,
                                          params_.max_contacts_per_point);

            float3 p_su = transforms.curves_to_surface * new_p;
            BLI_bvhtree_range_query_cpp(
                *surface_bvh.tree,
                p_su,
                max_dist_su,
                [&](const int triangle_i, const float3 &co_su, float dist_sq_su) {
                  const MLoopTri &looptri = surface_looptris[triangle_i];
                  const float3 v0_su = surface->mvert[surface->mloop[looptri.tri[0]].v].co;
                  const float3 v1_su = surface->mvert[surface->mloop[looptri.tri[1]].v].co;
                  const float3 v2_su = surface->mvert[surface->mloop[looptri.tri[2]].v].co;
                  float3 closest_su;
                  closest_on_tri_to_point_v3(closest_su, co_su, v0_su, v1_su, v2_su);
                  dist_sq_su = len_squared_v3v3(co_su, closest_su);
                  if (dist_sq_su <= max_dist_sq_su) {
                    const int contacts_num = contacts_num_[point_i];
                    int insert_i;
                    if (contacts_num < params_.max_contacts_per_point) {
                      insert_i = contacts_num;
                      ++contacts_num_[point_i];
                    }
                    else {
                      /* Replace the contact with the largest distance. */
                      insert_i = -1;
                      float max_dist_sq_su = dist_sq_su;
                      for (int contact_i : IndexRange(4)) {
                        if (contacts[contact_i].dist_sq_ > max_dist_sq_su) {
                          insert_i = contact_i;
                          max_dist_sq_su = contacts[contact_i].dist_sq_;
                        }
                      }
                    }
                    if (insert_i >= 0) {
                      float3 normal_su;
                      normal_tri_v3(normal_su, v0_su, v1_su, v2_su);
                      contacts[insert_i] = Contact{
                          max_dist_sq_su,
                          transforms.surface_to_curves_normal * float3{normal_su},
                          transforms.surface_to_curves * float3{closest_su}};
                    }
                  }
                });
          }
        }
      });
  result_.timing.find_contacts += PIL_check_seconds_timer() - find_contacts_start;
}

void ConstraintSolver::solve_constraints(CurvesGeometry &curves, VArray<int> changed_curves) const
{
  /* Gauss-Seidel method for solving length and contact constraints.
   * See for example "Position-Based Simulation Methods in Computer Graphics"
   * by Mueller et. al. for an in-depth description.
   */

  double solve_constraints_start = PIL_check_seconds_timer();
  MutableSpan<float3> positions_cu = curves.positions_for_write();
  VArray<float> radius = curves.attributes().lookup_or_default<float>(
      "radius", ATTR_DOMAIN_POINT, 0.0f);

  /* Compliance (inverse stiffness)
   * Alpha is used in physical simulation to control the softness of a constraint:
   * For alpha == 0 the constraint is stiff and the maximum correction factor is applied.
   * For values > 0 the constraint becomes squishy, and some violation is
   * permitted, and the constraint gets corrected over multiple time steps.
   */
  const float alpha = 0.0f;

  threading::parallel_for(
      changed_curves.index_range(), curves_grain_size, [&](const IndexRange range) {
        for (const int idx_curve : range) {
          const int curve_i = changed_curves[idx_curve];
          const IndexRange points = curves.points_for_curve(curve_i);

          /* Solve constraints */
          for (const int solver_i : IndexRange(params_.max_solver_iterations)) {
            /* Distance constraints */
            if (params_.use_length_constraints) {
              for (const int segment_i : IndexRange(points.size() - 1)) {
                const int point_a = points[segment_i];
                const int point_b = points[segment_i + 1];
                float3 &pa = positions_cu[point_a];
                float3 &pb = positions_cu[point_b];

                const float w0 = (params_.use_root_constraints && segment_i == 0 ? 0.0f : 1.0f);
                const float w1 = 1.0f;
                const float w = w0 + w1;

                float length_cu;
                const float3 direction = math::normalize_and_get_length(pb - pa, length_cu);
                const float expected_length_cu = segment_lengths_cu_[point_a];
                /* Constraint function */
                const float C = length_cu - expected_length_cu;
                const float3 gradient = direction;
                /* Lagrange multiplier */
                const float lambda = -C / (w + alpha);
                pa -= gradient * lambda * w0;
                pb += gradient * lambda * w1;
              }
            }

            /* Contact constraints */
            if (params_.use_collision_constraints) {
              const IndexRange points = params_.use_root_constraints ?
                                            curves.points_for_curve(curve_i).drop_front(1) :
                                            curves.points_for_curve(curve_i);
              for (const int point_i : points) {
                float3 &p = positions_cu[point_i];
                const float radius_p = radius[point_i];
                const int contacts_num = contacts_num_[point_i];
                Span<Contact> contacts(
                    contacts_.begin() + params_.max_contacts_per_point * point_i,
                                       contacts_num);
                for (const Contact &c : contacts) {
                  /* Constraint function */
                  const float C = math::dot(p - c.point_, c.normal_) - radius_p;
                  const float3 gradient = c.normal_;
                  /* Lagrange multiplier */
                  const float lambda = -C / (1.0f + alpha);
                  if (lambda > 0.0f) {
                    p += gradient * lambda;
                  }
                }
              }
            }
          }

          /* Compute residuals */

          /* Distance constraints */
          if (params_.use_length_constraints) {
            result_.constraint_count += points.size() - 1;
            for (const int segment_i : IndexRange(points.size() - 1)) {
              const int point_a = points[segment_i];
              const int point_b = points[segment_i + 1];
              float3 &pa = positions_cu[point_a];
              float3 &pb = positions_cu[point_b];

              const float w0 = (params_.use_root_constraints && segment_i == 0 ? 0.0f : 1.0f);
              const float w1 = 1.0f;
              const float w = w0 + w1;

              float length_cu = math::length(pb - pa);
              const float expected_length_cu = segment_lengths_cu_[point_a];
              /* Constraint function */
              const float C = length_cu - expected_length_cu;
              /* Lagrange multiplier */
              const float lambda = -C / (w + alpha);
              const double error_sq = lambda * lambda;
              result_.error_squared_sum += error_sq;
              result_.max_error_squared = std::max(result_.max_error_squared, error_sq);
            }
          }

          /* Contact constraints */
          if (params_.use_collision_constraints) {
            const IndexRange points = params_.use_root_constraints ?
                                          curves.points_for_curve(curve_i).drop_front(1) :
                                          curves.points_for_curve(curve_i);
            for (const int point_i : points) {
              float3 &p = positions_cu[point_i];
              const float radius_p = radius[point_i];
              const int contacts_num = contacts_num_[point_i];
              result_.constraint_count += contacts_num;
              Span<Contact> contacts(contacts_.begin() + params_.max_contacts_per_point * point_i,
                                     contacts_num);
              for (const Contact &c : contacts) {
                /* Constraint function */
                const float C = math::dot(p - c.point_, c.normal_) - radius_p;
                /* Lagrange multiplier */
                const float lambda = -C / (1.0f + alpha);
                if (lambda > 0.0f) {
                  const double error_sq = lambda * lambda;
                  result_.error_squared_sum += error_sq;
                  result_.max_error_squared = std::max(result_.max_error_squared, error_sq);
                }
              }
            }
          }
        }
      });

  if (result_.constraint_count > 0) {
    result_.rms_residual = sqrt(result_.error_squared_sum / result_.constraint_count);
  }
  else {
    result_.rms_residual = 0.0;
  }
  if (result_.max_error_squared > params_.error_threshold * params_.error_threshold) {
    result_.error = ErrorType::NotConverged;
  }

  result_.timing.solve_constraints += PIL_check_seconds_timer() - solve_constraints_start;
}

}  // namespace blender::bke::curves
