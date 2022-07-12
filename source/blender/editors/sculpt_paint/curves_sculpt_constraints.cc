/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "curves_sculpt_intern.hh"

#include "BLI_index_mask_ops.hh"

#include "BKE_bvhutils.h"
#include "BKE_mesh_runtime.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

namespace blender::ed::sculpt_paint {

using blender::bke::CurvesGeometry;
using threading::EnumerableThreadSpecific;

static const int MAX_CONTACTS = 4;

void CurvesConstraintSolver::initialize(const CurvesGeometry *curves)
{
  Span<float3> positions_cu = curves->positions();

  segment_lengths_cu_.reinitialize(curves->points_num());
  threading::parallel_for(curves->curves_range(), 128, [&](const IndexRange range) {
    for (const int curve_i : range) {
      const IndexRange points = curves->points_for_curve(curve_i);

      float length_cu = 0.0f, prev_length_cu;
      for (const int point_i : points.drop_back(1)) {
        const float3 &p1_cu = positions_cu[point_i];
        const float3 &p2_cu = positions_cu[point_i + 1];
        prev_length_cu = length_cu;
        length_cu = math::distance(p1_cu, p2_cu);
        segment_lengths_cu_[point_i] = length_cu;
      }
    }
  });

  contacts_num_.reinitialize(curves->points_num());
  contacts_.reinitialize(MAX_CONTACTS * curves->points_num());
}

void CurvesConstraintSolver::find_contact_points(
    const Depsgraph *depsgraph,
    Object *object,
    const CurvesGeometry *curves,
    const Object *surface_ob,
    Span<float3> orig_positions,
    threading::EnumerableThreadSpecific<Vector<int>> &changed_curves)
{
  if (surface_ob == nullptr || surface_ob->type != OB_MESH) {
    contacts_.reinitialize(0);
    return;
  }

  const Mesh *surface = static_cast<Mesh *>(surface_ob->data);
  Span<MLoopTri> surface_looptris = {BKE_mesh_runtime_looptri_ensure(surface),
                                     BKE_mesh_runtime_looptri_len(surface)};

  BKE_bvhtree_from_mesh_get(&surface_bvh_, surface, BVHTREE_FROM_LOOPTRI, 2);
  BLI_SCOPED_DEFER([&]() { free_bvhtree_from_mesh(&surface_bvh_); });

  contacts_num_.fill(0);

  //std::cout << "FIND contacts:" << std::endl;
  threading::parallel_for_each(changed_curves, [&](const Vector<int> &changed_curves) {
    threading::parallel_for(changed_curves.index_range(), 256, [&](const IndexRange range) {
      for (const int curve_i : changed_curves.as_span().slice(range)) {
        const IndexRange points = curves->points_for_curve(curve_i);
        for (const int point_i : points) {
          const float3 &old_p = orig_positions[point_i];
          const float3 &new_p = curves->positions()[point_i];
          const float margin = 0.01f; /* XXX TODO get curve radius */
          const float margin_sq = margin * margin;

          MutableSpan<Contact> contacts(contacts_.begin() + MAX_CONTACTS * point_i, MAX_CONTACTS);

          float3 dir = new_p - old_p;
          float hit_dist = normalize_v3(dir) + margin;
          //std::cout << "  p " << point_i << " move distance " << hit_dist << std::endl;
          BLI_bvhtree_ray_cast_all_cpp(
              *surface_bvh_.tree,
              old_p,
              dir,
              margin,
              hit_dist,
              [&](const int triangle_i, const BVHTreeRay &ray, BVHTreeRayHit &hit) {
                surface_bvh_.raycast_callback(&surface_bvh_, triangle_i, &ray, &hit);

                if (hit.index >= 0) {
                  const float dist = hit.dist;

                  const int contacts_num = contacts_num_[point_i];
                  int insert_i;
                  if (contacts_num < MAX_CONTACTS) {
                    insert_i = contacts_num;
                    ++contacts_num_[point_i];
                  }
                  else {
                    /* Replace the contact with the largest distance. */
                    /* XXX this is ugly, can be optimized a good deal (simd?) */
                    float max_dist = dist;
                    insert_i = -1;

                    if (contacts[0].dist_ > max_dist) {
                      max_dist = contacts[0].dist_;
                      insert_i = 0;

                      if (contacts[1].dist_ > max_dist) {
                        max_dist = contacts[1].dist_;
                        insert_i = 1;

                        if (contacts[2].dist_ > max_dist) {
                          max_dist = contacts[2].dist_;
                          insert_i = 2;

                          if (contacts[3].dist_ > max_dist) {
                            max_dist = contacts[3].dist_;
                            insert_i = 3;
                          }
                        }
                      }
                    }
                  }

                  if (insert_i >= 0) {
                    contacts[insert_i] = Contact{dist, hit.no, hit.co};
                    // std::cout << "  point " << point_i << " normal " << normal_su << " distance
                    // "
                    // << sqrtf(dist_sq) << std::endl;
                  }
                }
              });
        }
      }
    });
  });
}

void CurvesConstraintSolver::solve_constraints(
    CurvesGeometry *curves, EnumerableThreadSpecific<Vector<int>> &changed_curves) const
{
  /* Gauss-Seidel method for solving length and contact constraints.
   * See for example "Position-Based Simulation Methods in Computer Graphics"
   * by Mueller et. al. for an in-depth description.
   */
  const int solver_iterations = 5;

  const Span<float> expected_lengths_cu = segment_lengths_cu_;
  MutableSpan<float3> positions_cu = curves->positions_for_write();

  threading::parallel_for_each(changed_curves, [&](const Vector<int> &changed_curves) {
    threading::parallel_for(changed_curves.index_range(), 256, [&](const IndexRange range) {
      for (const int curve_i : changed_curves.as_span().slice(range)) {
        const IndexRange points = curves->points_for_curve(curve_i);

        for (int solver_i : IndexRange(solver_iterations)) {
          /* Solve contact constraints */
          //if (solver_i == 0) {
          //  std::cout << "SOLVE:" << std::endl;
          //}
          for (const int point_i : points) {
            const float3 &p = positions_cu[point_i];
            const float margin = 0.01f; /* XXX TODO get curve radius */
            const float margin_sq = margin * margin;

            const int contacts_num = contacts_num_[point_i];
            //if (solver_i == 0 && point_i == 7 && contacts_num > 0) {
            //  std::cout << "  point " << point_i << " contacts " << contacts_num << std::endl;
            //}
            Span<Contact> contacts(contacts_.begin() + MAX_CONTACTS * point_i, contacts_num);
            int contact_i = 0;
            for (const Contact &c : contacts) {
              //float dist_sq = len_squared_v3v3(p, c.point_);
              //if (dist_sq > margin_sq) {
              //  continue;
              //}

              /* Lagrange multiplier for solving a single contact constraint. */
              float lambda = dot_v3v3(p - c.point_, c.normal_) - margin;
              //if (solver_i == 0 && point_i == 7) {
              //  std::cout << "    contact " << contact_i << " lambda=" << lambda << std::endl;
              //}
              if (lambda < 0.0f) {
                positions_cu[point_i] = p - lambda * c.normal_;
              }
              //++contact_i;
            }
          }

          /* Solve distance constraints */
          for (const int segment_i : points.drop_back(1)) {
            const float3 &p1_cu = positions_cu[segment_i];
            float3 &p2_cu = positions_cu[segment_i + 1];
            const float3 direction = math::normalize(p2_cu - p1_cu);
            const float expected_length_cu = expected_lengths_cu[segment_i];
            p2_cu = p1_cu + direction * expected_length_cu;
          }
        }
      }
    });
  });
}

}  // namespace blender::ed::sculpt_paint
