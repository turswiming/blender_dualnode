/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "GEO_constraint_solver.hh"

#include "BLI_index_mask_ops.hh"

#include "BKE_bvhutils.h"
#include "BKE_mesh_runtime.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_types.h"

#include "PIL_time.h"

namespace blender::geometry {

using bke::CurvesGeometry;
using bke::CurvesSurfaceTransforms;
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
                                   VArray<int> changed_curves,
                                   bool update_error)
{
  const double step_start = PIL_check_seconds_timer();

  clear_result();

  const MutableSpan<float3> positions = curves.positions_for_write();

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

  if (update_error) {
    compute_error(curves, changed_curves);
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

  const Span<MLoopTri> surface_looptris = {BKE_mesh_runtime_looptri_ensure(surface),
                                           BKE_mesh_runtime_looptri_len(surface)};

  const double build_bvh_start = PIL_check_seconds_timer();

  /** BVH tree of the surface mesh for finding collisions. */
  BVHTreeFromMesh surface_bvh;
  BKE_bvhtree_from_mesh_get(
      &surface_bvh, surface, BVHTREE_FROM_LOOPTRI, params_.bvh_branching_factor);
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

            const MutableSpan<Contact> contacts = contacts_.as_mutable_span().slice(
                params_.max_contacts_per_point * point_i, params_.max_contacts_per_point);

            const float3 p_su = transforms.curves_to_surface * new_p;
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

void ConstraintSolver::apply_distance_constraint(float3 &point_a,
                                                 float3 &point_b,
                                                 const float segment_length,
                                                 const float weight_a,
                                                 const float weight_b) const
{
  float length_cu;
  const float3 direction = math::normalize_and_get_length(point_b - point_a, length_cu);
  /* Constraint function */
  const float C = length_cu - segment_length;
  const float3 gradient = direction;

  /* Lagrange multiplier */
  const float lambda = -C / (weight_a + weight_b + params_.alpha);
  point_a -= gradient * lambda * weight_a;
  point_b += gradient * lambda * weight_b;
}

float ConstraintSolver::get_distance_constraint_error(const float3 &point_a,
                                                      const float3 &point_b,
                                                      const float segment_length) const
{
  float length_cu = math::length(point_b - point_a);
  /* Constraint function */
  const float C = length_cu - segment_length;
  return C;
}

void ConstraintSolver::apply_contact_constraint(float3 &point,
                                                const float radius,
                                                const ConstraintSolver::Contact &contact) const
{
  /* Constraint function */
  const float C = math::dot(point - contact.point_, contact.normal_) - radius;
  const float3 gradient = contact.normal_;

  /* Lagrange multiplier */
  const float lambda = -C / (1.0f + params_.alpha);
  if (lambda > 0.0f) {
    point += gradient * lambda;
  }
}

float ConstraintSolver::get_contact_constraint_error(
    const float3 &point, const float radius, const ConstraintSolver::Contact &contact) const
{
  /* Constraint function */
  const float C = math::dot(point - contact.point_, contact.normal_) - radius;
  return math::min(C, 0.0f);
}

void ConstraintSolver::solve_constraints(CurvesGeometry &curves, VArray<int> changed_curves) const
{
  const int solver_max_iterations = [&]() {
    switch (params_.solver_type) {
      case SolverType::Sequential:
        return 1;
      case SolverType::PositionBasedDynamics:
        return params_.max_solver_iterations;
    }
    return 0;
  }();

  const double solve_constraints_start = PIL_check_seconds_timer();

  VArray<float> radius = curves.attributes().lookup_or_default<float>(
      "radius", ATTR_DOMAIN_POINT, 0.0f);

  threading::parallel_for(
      changed_curves.index_range(), curves_grain_size, [&](const IndexRange range) {
        for (const int idx_curve : range) {
          const int curve_i = changed_curves[idx_curve];
          const IndexRange points = curves.points_for_curve(curve_i);

          /* Solve constraints */
          for (const int solver_i : IndexRange(solver_max_iterations)) {
            solve_curve_constraints(curves, radius, points);
          }
        }
      });

  result_.timing.solve_constraints += PIL_check_seconds_timer() - solve_constraints_start;
}

void ConstraintSolver::solve_curve_constraints(CurvesGeometry &curves,
                                               const VArray<float> radius,
                                               const IndexRange points) const
{
  const int skipped_root_points = [&]() {
    switch (params_.solver_type) {
      /* The sequential solver only moves the 2nd point of each segment. */
      case SolverType::Sequential:
        return (int)points.size();
      /* The PBD solver moves both points, except for the pinned root. */
      case SolverType::PositionBasedDynamics:
        return params_.use_root_constraints ? 1 : 0;
    }
    return 0;
  }();

  MutableSpan<float3> positions_cu = curves.positions_for_write();

  result_.constraint_count = 0;

  /* Distance constraints */
  if (params_.use_length_constraints) {
    result_.constraint_count += points.size() - 1;
    for (const int segment_i : IndexRange(points.size() - 1)) {
      const int point_a = points[segment_i];
      const int point_b = points[segment_i + 1];
      const float segment_length = segment_lengths_cu_[point_a];
      float3 &pa = positions_cu[point_a];
      float3 &pb = positions_cu[point_b];

      const float w0 = (segment_i < skipped_root_points ? 0.0f : 1.0f);
      const float w1 = 1.0f;
      apply_distance_constraint(pa, pb, segment_length, w0, w1);
    }
  }

  /* Contact constraints */
  if (params_.use_collision_constraints) {
    for (const int point_i : points.drop_front(skipped_root_points)) {
      float3 &p = positions_cu[point_i];
      const float radius_p = radius[point_i];
      const int contacts_num = contacts_num_[point_i];
      result_.constraint_count += contacts_num;

      const Span<Contact> contacts = contacts_.as_span().slice(
          params_.max_contacts_per_point * point_i, contacts_num);
      for (const Contact &c : contacts) {
        apply_contact_constraint(p, radius_p, c);
      }
    }
  }
}

void ConstraintSolver::compute_error(const CurvesGeometry &curves, VArray<int> changed_curves) const
{
  VArray<float> radius = curves.attributes().lookup_or_default<float>(
      "radius", ATTR_DOMAIN_POINT, 0.0f);

  /* Accumulate error (no threading) */
  for (int i : changed_curves.index_range()) {
    const int curve_i = changed_curves[i];
    const IndexRange points = curves.points_for_curve(curve_i);
    compute_curve_error(curves, radius, points);
  }

  if (result_.constraint_count > 0) {
    result_.residual.rms_error = sqrt(result_.residual.error_squared_sum /
                                      result_.constraint_count);
  }
  else {
    result_.residual.rms_error = 0.0;
  }

  if (result_.residual.max_error_squared > params_.error_threshold * params_.error_threshold) {
    result_.status = Result::Status::ErrorNoConvergence;
  }
}

void ConstraintSolver::compute_curve_error(const CurvesGeometry &curves,
                                           const VArray<float> radius,
                                           const IndexRange points) const
{
  Span<float3> positions_cu = curves.positions();

  /* Distance constraints */
  if (params_.use_length_constraints) {
    for (const int segment_i : IndexRange(points.size() - 1)) {
      const int point_a = points[segment_i];
      const int point_b = points[segment_i + 1];
      const float segment_length = segment_lengths_cu_[point_a];
      const float3 &pa = positions_cu[point_a];
      const float3 &pb = positions_cu[point_b];

      const float error = get_distance_constraint_error(pa, pb, segment_length);
      const double error_sq = error * error;
      result_.residual.error_squared_sum += error_sq;
      result_.residual.max_error_squared = std::max(result_.residual.max_error_squared, error_sq);
    }
  }

  /* Contact constraints */
  if (params_.use_collision_constraints) {
    for (const int point_i : points) {
      const float3 &p = positions_cu[point_i];
      const float radius_p = radius[point_i];
      const int contacts_num = contacts_num_[point_i];
      const Span<Contact> contacts = contacts_.as_span().slice(
          params_.max_contacts_per_point * point_i, contacts_num);
      for (const Contact &c : contacts) {
        const float error = get_contact_constraint_error(p, radius_p, c);
        const double error_sq = error * error;
        result_.residual.error_squared_sum += error_sq;
        result_.residual.max_error_squared = std::max(result_.residual.max_error_squared,
                                                      error_sq);
      }
    }
  }
}

}  // namespace blender::geometry
