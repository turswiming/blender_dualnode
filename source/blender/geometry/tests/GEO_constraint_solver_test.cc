/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "testing/testing.h"

#include "BKE_curves.hh"

#include "GEO_constraint_solver.hh"

namespace blender::geometry::tests {

using bke::CurvesGeometry;
using bke::CurvesSurfaceTransforms;

inline CurvesSurfaceTransforms create_curves_surface_transforms()
{
  CurvesSurfaceTransforms transforms;
  transforms.curves_to_world = float4x4::identity();
  transforms.curves_to_surface = float4x4::identity();
  transforms.world_to_curves = float4x4::identity();
  transforms.world_to_surface = float4x4::identity();
  transforms.surface_to_world = float4x4::identity();
  transforms.surface_to_curves = float4x4::identity();
  transforms.surface_to_curves_normal = float4x4::identity();
  return transforms;
}

TEST(curves_constraints, LengthAndRootConstraint)
{
  CurvesGeometry curves(9, 3);
  curves.fill_curve_types(CURVE_TYPE_POLY);
  const MutableSpan<int> point_offsets = curves.offsets_for_write();
  point_offsets[0] = 0;
  point_offsets[1] = 2;
  point_offsets[2] = 6;
  point_offsets[3] = 9;
  const MutableSpan<float3> positions = curves.positions_for_write();
  positions[0] = float3(.0f, 0, 0);
  positions[1] = float3(.1f, 0, 0);

  positions[2] = float3(.0f, 0, 1);
  positions[3] = float3(.1f, 0, 1);
  positions[4] = float3(.15f, 0, 1);

  positions[5] = float3(.0f, 0, 2);
  positions[6] = float3(.1f, 0, 2);
  positions[7] = float3(.15f, 0, 2);
  positions[8] = float3(.25f, 0, 2);

  const CurvesSurfaceTransforms transforms = create_curves_surface_transforms();

  ConstraintSolver solver;
  ConstraintSolver::Params params;
  params.use_collision_constraints = false;

  solver.initialize(params, curves, curves.curves_range());

  const VArray<int> changed_curves = VArray<int>::ForFunc(curves.curves_num(),
                                                    [](int64_t i) { return (int)i; });
  const Array<float3> orig_positions = curves.positions();

  positions[1].y += 0.1f;
  positions[3].y += 0.1f;
  positions[6].y += 0.1f;

  solver.step_curves(curves, nullptr, transforms, orig_positions, changed_curves);

  EXPECT_EQ(solver.result().status, ConstraintSolver::Result::Status::Ok);
  EXPECT_LE(solver.result().residual.max_error_squared, params.error_threshold * params.error_threshold);
  EXPECT_LE(solver.result().residual.rms_error, params.error_threshold);

  const float eps = 0.005f;
  EXPECT_V3_NEAR(positions[0], float3(0.0f, 0.0f, 0.0f), eps);
  EXPECT_V3_NEAR(positions[1], float3(0.065f, 0.075f, 0.0f), eps);
  EXPECT_V3_NEAR(positions[2], float3(0.0f, 0.0f, 1.0f), eps);
  EXPECT_V3_NEAR(positions[3], float3(0.0824425220, 0.057f, 1.0f), eps);
  EXPECT_V3_NEAR(positions[4], float3(0.118486673, 0.022f, 1.0f), eps);
  EXPECT_V3_NEAR(positions[5], float3(0.0f, 0.0f, 2.0f), eps);
  EXPECT_V3_NEAR(positions[6], float3(0.1f, 0.1f, 2.0f), eps);
  EXPECT_V3_NEAR(positions[7], float3(0.125f, 0.057f, 2.0f), eps);
  EXPECT_V3_NEAR(positions[8], float3(0.213f, 0.009f, 2.0f), eps);
}

}  // namespace blender::geometry::tests
