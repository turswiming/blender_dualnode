/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "testing/testing.h"

#include "CLG_log.h"

#include "BLI_noise.hh"
#include "BLI_rand.hh"

#include "BKE_curves.hh"
#include "BKE_curves_constraints.hh"
#include "BKE_idtype.h"
#include "BKE_lib_id.h"
#include "BKE_mesh.h"

#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

namespace blender::bke::curves::tests {

using curves::ConstraintSolver;

class CurveConstraintSolverPerfTestSuite : public testing::Test {
 public:
  /* Sets up Blender just enough to not crash on creating a mesh. */
  static void SetUpTestSuite()
  {
    testing::Test::SetUpTestSuite();

    CLG_init();
    BLI_threadapi_init();

    //DNA_sdna_current_init();
    //BKE_blender_globals_init();

    BKE_idtype_init();
  }

  static void TearDownTestSuite()
  {
    BLI_threadapi_exit();

    //BKE_blender_atexit();

    //BKE_tempdir_session_purge();
    //BKE_appdir_exit();
    CLG_exit();

    testing::Test::TearDownTestSuite();
  }
};

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

  EXPECT_EQ(solver.result().error, ConstraintSolver::ErrorType::Ok);
  EXPECT_LE(solver.result().max_error_squared, params.error_threshold * params.error_threshold);
  EXPECT_LE(solver.result().rms_residual, params.error_threshold);

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

#define DO_PERF_TESTS 1

#if DO_PERF_TESTS

const uint32_t randomized_curves_seed = 12345;
const uint32_t randomized_offset_seed = 23451;

static float3 random_point_inside_unit_sphere(RandomNumberGenerator &rng)
{
  const float theta = (float)(M_PI * 2.0) * rng.get_float();
  const float phi = acosf(2.0f * rng.get_float() - 1.0f);
  const float r = sqrt3f(rng.get_float());
  return r * float3{sinf(phi) * cosf(theta), sinf(phi) * sinf(theta), cosf(phi)};
}

static Mesh *create_noise_grid(const float noise = 0.1f, const int resolution = 100, const float size = 1.0f)
{
  const int tot_verts = resolution * resolution;
  const int tot_polys = (resolution - 1) * (resolution - 1);
  const int tot_loops = tot_polys * 4;
  Mesh *mesh = BKE_mesh_new_nomain(tot_verts, 0, 0, tot_loops, tot_polys);
  const MutableSpan<MVert> mverts(mesh->mvert, mesh->totvert);
  const MutableSpan<MPoly> mpolys(mesh->mpoly, mesh->totpoly);
  const MutableSpan<MLoop> mloops(mesh->mloop, mesh->totloop);

  for (const int i : IndexRange(resolution)) {
    for (const int j : IndexRange(resolution)) {
      float3 co = float3((float)i - 0.5f, (float)j - 0.5f, 0) * size;
      const float3 offset = noise::perlin_float3_fractal_distorted(co, 2.0f, 0.5f, 0.0f);
      co += noise * (offset - float3(0.5f));

      const int vert_i = i * resolution + j;
      copy_v3_v3(mverts[vert_i].co, co);
    }
  }

  for (const int i : IndexRange(resolution - 1)) {
    for (const int j : IndexRange(resolution - 1)) {
      const int vert_i = i * resolution + j;
      const int poly_i = i * (resolution - 1) + j;
      const int loop_i = poly_i * 4;

      mpolys[poly_i].loopstart = loop_i;
      mpolys[poly_i].totloop = 4;

      mloops[loop_i + 0].v = vert_i;
      mloops[loop_i + 1].v = vert_i + 1;
      mloops[loop_i + 2].v = vert_i + resolution + 1;
      mloops[loop_i + 3].v = vert_i + resolution;
    }
  }

  BKE_mesh_calc_edges(mesh, false, false);

  return mesh;
}

static CurvesGeometry create_randomized_curves(const int curve_num,
                                               const int point_num_min,
                                               const int point_num_max,
                                               const float dist_min,
                                               const float dist_max)
{
  RandomNumberGenerator rng(randomized_curves_seed);

  Array<int> point_offsets(curve_num + 1);
  const int point_num_range = std::max(point_num_max - point_num_min, 0);
  int point_num = 0;
  for (int curve_i : IndexRange(curve_num)) {
    point_offsets[curve_i] = point_num;
    point_num += point_num_min + (rng.get_int32() % point_num_range);
  }
  point_offsets.last() = point_num;

  CurvesGeometry curves(point_num, curve_num);
  curves.fill_curve_types(CURVE_TYPE_POLY);
  curves.offsets_for_write().copy_from(point_offsets);

  const MutableSpan<float3> positions = curves.positions_for_write();
  for (int curve_i : curves.curves_range()) {
    float3 pos = random_point_inside_unit_sphere(rng);
    for (int point_i : curves.points_for_curve(curve_i)) {
      positions[point_i] = pos;
      pos += rng.get_unit_float3() * interpf(dist_max, dist_min, rng.get_float());
    }
  }

  return curves;
}

static void randomized_point_offset(CurvesGeometry &curves,
                                    VArray<int> changed_curves,
                                    float dist_min,
                                    float dist_max,
                                    float3 direction = float3(1, 0, 0),
                                    float cone_angle = M_PI)
{
  RandomNumberGenerator rng(randomized_offset_seed);

  math::normalize(direction);
  float3 n1, n2;
  ortho_basis_v3v3_v3(n1, n2, direction);

  const MutableSpan<float3> positions = curves.positions_for_write();
  for (int i : changed_curves.index_range()) {
    const int curve_i = changed_curves[i];
    for (int point_i : curves.points_for_curve(curve_i)) {
      const float theta = rng.get_float() * cone_angle;
      const float phi = rng.get_float() * (float)(2.0 * M_PI);
      const float3 dir_p = direction * cosf(theta) + n1 * sinf(theta) * cosf(phi) +
                      n2 * sinf(theta) * sinf(phi);
      const float dist_p = interpf(dist_max, dist_min, rng.get_float());
      positions[point_i] += dir_p * dist_p;
    }
  }
}

static void print_test_stats(const CurvesGeometry &curves, const ConstraintSolver::Result &result)
{
  const testing::TestInfo *const test_info = testing::UnitTest::GetInstance()->current_test_info();

  std::cout << "Curves: " << curves.curves_num() << " Points: " << curves.points_num() << std::endl;
  std::cout << "RMS=" << result.rms_residual << " max error=" << sqrt(result.max_error_squared)
            << std::endl;
  std::cout << "total: " << result.timing.step_total * 1000.0f
            << " ms, build bvh: " << result.timing.build_bvh * 1000.0f
            << " ms, find contacts: " << result.timing.find_contacts * 1000.0f
            << " ms, solve: " << result.timing.solve_constraints * 1000.0f << " ms" << std::endl;
}

template <typename ParamType>
class SolverPerfTestSuite : public CurveConstraintSolverPerfTestSuite,
                            public testing::WithParamInterface<ParamType> {
 public:
  struct TestResult {
    ParamType param;
    ConstraintSolver::Result solver_result;
  };

  static Vector<TestResult> &results()
  {
    static Vector<TestResult> results_;
    return results_;
  }

  static void SetUpTestSuite()
  {
    CurveConstraintSolverPerfTestSuite::SetUpTestSuite();

    results().clear();
  }

  static void TearDownTestSuite()
  {
    /* CSV printout for simple data import */
    std::cout << "Parameter,RMS Error,Max Error,Collision Detection (ms),Solve Time (ms)" << std::endl;
    for (const TestResult &result : results()) {
      std::cout << result.param << "," << result.solver_result.rms_residual << ","
                << sqrtf(result.solver_result.max_error_squared) << ","
                << result.solver_result.timing.find_contacts << ","
                << result.solver_result.timing.solve_constraints << std::endl;
    }

    CurveConstraintSolverPerfTestSuite::TearDownTestSuite();
  }

  void randomized_test(const ParamType &test_param, const ConstraintSolver::Params &solver_params)
  {
    CurvesGeometry curves = create_randomized_curves(10000, 4, 50, 0.1f, 0.2f);
    const CurvesSurfaceTransforms transforms = create_curves_surface_transforms();

    Mesh *surface = create_noise_grid();

    ConstraintSolver solver;
    solver.initialize(solver_params, curves, curves.curves_range());

    const VArray<int> changed_curves = VArray<int>::ForFunc(curves.curves_num(),
                                                            [](int64_t i) { return (int)i; });
    const Array<float3> orig_positions = curves.positions();
    randomized_point_offset(curves, changed_curves, 0.0f, 1.0f);
    solver.step_curves(curves, surface, transforms, orig_positions, changed_curves, true);

    BKE_id_free(nullptr, surface);

    results().append(TestResult{test_param, solver.result()});
  }
};

using SolverIterationsTestSuite = SolverPerfTestSuite<int>;

TEST_P(SolverIterationsTestSuite, RandomizedTest)
{
  ConstraintSolver::Params params;
  params.max_solver_iterations = GetParam();

  randomized_test(GetParam(), params);
}

INSTANTIATE_TEST_SUITE_P(SolverIterationsTests,
                         SolverIterationsTestSuite,
                         testing::Values(1, 5, 10, 20, 50, 100),
                         [](const testing::TestParamInfo<int> info) {
                           std::stringstream ss;
                           ss << "SolverIterationsTest_" << info.param;
                           return ss.str();
                         });

using BVHBranchingFactorTestSuite = SolverPerfTestSuite<int>;

TEST_P(BVHBranchingFactorTestSuite, RandomizedTest)
{
  ConstraintSolver::Params params;
  params.bvh_branching_factor = GetParam();

  randomized_test(GetParam(), params);
}

INSTANTIATE_TEST_SUITE_P(BVHBranchingFactorTests,
                         BVHBranchingFactorTestSuite,
                         testing::Values(2, 4, 6, 8, 16, 32),
                         [](const testing::TestParamInfo<int> info) {
                           std::stringstream ss;
                           ss << "BVHBranchingFactorTest_" << info.param;
                           return ss.str();
                         });

#endif  // DO_PERF_TESTS

}  // namespace blender::bke::curves::tests
