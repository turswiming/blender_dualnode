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

void CurvesConstraintSolver::find_contact_points(const Depsgraph *depsgraph,
                                                 Object *object,
                                                 const CurvesGeometry *curves,
                                                 const Object *surface_ob,
                                                 const CurvesSurfaceTransforms &transforms,
                                                 Span<float3> orig_positions,
                                                 Span<int> changed_curves)
{
  contacts_num_.fill(0);

  if (surface_ob == nullptr || surface_ob->type != OB_MESH) {
    contacts_.reinitialize(0);
    return;
  }

  const float curves_to_surface_scale = mat4_to_scale(transforms.curves_to_surface.ptr());
  const float surface_to_curves_scale = mat4_to_scale(transforms.curves_to_surface.ptr());

  const Mesh *surface = static_cast<Mesh *>(surface_ob->data);
  Span<MLoopTri> surface_looptris = {BKE_mesh_runtime_looptri_ensure(surface),
                                     BKE_mesh_runtime_looptri_len(surface)};

  /** BVH tree of the surface mesh for finding collisions. */
  BVHTreeFromMesh surface_bvh;
  BKE_bvhtree_from_mesh_get(&surface_bvh, surface, BVHTREE_FROM_LOOPTRI, 2);
  BLI_SCOPED_DEFER([&]() { free_bvhtree_from_mesh(&surface_bvh); });

  VArray<float> radius = curves->attributes().lookup_or_default<float>(
      "radius", ATTR_DOMAIN_POINT, 0.0f);

  threading::parallel_for(changed_curves.index_range(), 256, [&](const IndexRange range) {
    for (const int curve_i : changed_curves.slice(range)) {
      const IndexRange points = curves->points_for_curve(curve_i);
      /* First point is anchored to the surface, ignore collisions. */
      for (const int point_i : points.drop_front(1)) {
        const float3 &old_p = orig_positions[point_i];
        const float3 &new_p = curves->positions()[point_i];
        const float margin = radius[point_i];

        MutableSpan<Contact> contacts(contacts_.begin() + MAX_CONTACTS * point_i, MAX_CONTACTS);

        float3 start_su = transforms.curves_to_surface * old_p;
        float3 dir_su = transforms.curves_to_surface.ref_3x3() * (new_p - old_p);
        float margin_su = curves_to_surface_scale * margin;
        float hit_dist_su = normalize_v3(dir_su) + margin_su;
        BLI_bvhtree_ray_cast_all_cpp(
            *surface_bvh.tree,
            start_su,
            dir_su,
            margin_su,
            hit_dist_su,
            [&](const int triangle_i, const BVHTreeRay &ray, BVHTreeRayHit &hit) {
              surface_bvh.raycast_callback(&surface_bvh, triangle_i, &ray, &hit);

              if (hit.index >= 0) {
                const float dist_cu = surface_to_curves_scale * hit.dist;

                const int contacts_num = contacts_num_[point_i];
                int insert_i;
                if (contacts_num < MAX_CONTACTS) {
                  insert_i = contacts_num;
                  ++contacts_num_[point_i];
                }
                else {
                  /* Replace the contact with the largest distance. */
                  insert_i = -1;
                  float max_dist_cu = dist_cu;
                  for (int contact_i : IndexRange(4)) {
                    if (contacts[contact_i].dist_ > max_dist_cu) {
                      insert_i = contact_i;
                      max_dist_cu = contacts[contact_i].dist_;
                    }
                  }
                }
                if (insert_i >= 0) {
                  contacts[insert_i] = Contact{dist_cu,
                                                transforms.surface_to_curves_normal *
                                                    float3{hit.no},
                                                transforms.surface_to_curves * float3{hit.co}};
                }
              }
            });
      }
    }
  });
}

void CurvesConstraintSolver::solve_constraints(
    CurvesGeometry *curves, Span<int> changed_curves) const
{
  /* Gauss-Seidel method for solving length and contact constraints.
   * See for example "Position-Based Simulation Methods in Computer Graphics"
   * by Mueller et. al. for an in-depth description.
   */
  const int solver_iterations = 5;

  const Span<float> expected_lengths_cu = segment_lengths_cu_;
  MutableSpan<float3> positions_cu = curves->positions_for_write();
  VArray<float> radius = curves->attributes().lookup_or_default<float>(
      "radius", ATTR_DOMAIN_POINT, 0.0f);

  threading::parallel_for(changed_curves.index_range(), 256, [&](const IndexRange range) {
    for (const int curve_i : changed_curves.slice(range)) {
      const IndexRange points = curves->points_for_curve(curve_i);

      /* First point is anchored to surface, contact and length constraints to no apply. */
      for (const int point_i : points.drop_front(1)) {
        float3 &p = positions_cu[point_i];
        const int contacts_num = contacts_num_[point_i];
        for (int solver_i : IndexRange(solver_iterations)) {
          /* Solve contact constraints */
          Span<Contact> contacts(contacts_.begin() + MAX_CONTACTS * point_i, contacts_num);
          for (const Contact &c : contacts) {
            /* Lagrange multiplier for solving a single contact constraint.
              * Note: The contact point is already offset from the surface by the radius due to the raycast callback,
              *       no need to subtract the radius from lambda. */
            const float lambda = dot_v3v3(p - c.point_, c.normal_);
            if (lambda < 0.0f) {
              p -= lambda * c.normal_;
            }
          }

          /* Solve distance constraint */
          {
            const float3 &p_prev = positions_cu[point_i - 1];
            const float3 direction = math::normalize(p - p_prev);
            const float expected_length_cu = expected_lengths_cu[point_i - 1];
            p = p_prev + direction * expected_length_cu;
          }
        }
      }
    }
  });
}

}  // namespace blender::ed::sculpt_paint
