/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "BLI_math_matrix.h"
#include "BLI_math_matrix.hh"

TEST(math_matrix, interp_m4_m4m4_regular)
{
  /* Test 4x4 matrix interpolation without singularity, i.e. without axis flip. */

  /* Transposed matrix, so that the code here is written in the same way as print_m4() outputs. */
  /* This matrix represents T=(0.1, 0.2, 0.3), R=(40, 50, 60) degrees, S=(0.7, 0.8, 0.9) */
  float matrix_a[4][4] = {
      {0.224976f, -0.333770f, 0.765074f, 0.100000f},
      {0.389669f, 0.647565f, 0.168130f, 0.200000f},
      {-0.536231f, 0.330541f, 0.443163f, 0.300000f},
      {0.000000f, 0.000000f, 0.000000f, 1.000000f},
  };
  transpose_m4(matrix_a);

  float matrix_i[4][4];
  unit_m4(matrix_i);

  float result[4][4];
  const float epsilon = 1e-6;
  interp_m4_m4m4(result, matrix_i, matrix_a, 0.0f);
  EXPECT_M4_NEAR(result, matrix_i, epsilon);

  interp_m4_m4m4(result, matrix_i, matrix_a, 1.0f);
  EXPECT_M4_NEAR(result, matrix_a, epsilon);

  /* This matrix is based on the current implementation of the code, and isn't guaranteed to be
   * correct. It's just consistent with the current implementation. */
  float matrix_halfway[4][4] = {
      {0.690643f, -0.253244f, 0.484996f, 0.050000f},
      {0.271924f, 0.852623f, 0.012348f, 0.100000f},
      {-0.414209f, 0.137484f, 0.816778f, 0.150000f},
      {0.000000f, 0.000000f, 0.000000f, 1.000000f},
  };

  transpose_m4(matrix_halfway);
  interp_m4_m4m4(result, matrix_i, matrix_a, 0.5f);
  EXPECT_M4_NEAR(result, matrix_halfway, epsilon);
}

TEST(math_matrix, interp_m3_m3m3_singularity)
{
  /* A singularity means that there is an axis mirror in the rotation component of the matrix.
   * This is reflected in its negative determinant.
   *
   * The interpolation of 4x4 matrices performs linear interpolation on the translation component,
   * and then uses the 3x3 interpolation function to handle rotation and scale. As a result, this
   * test for a singularity in the rotation matrix only needs to test the 3x3 case. */

  /* Transposed matrix, so that the code here is written in the same way as print_m4() outputs. */
  /* This matrix represents R=(4, 5, 6) degrees, S=(-1, 1, 1) */
  float matrix_a[3][3] = {
      {-0.990737f, -0.098227f, 0.093759f},
      {-0.104131f, 0.992735f, -0.060286f},
      {0.087156f, 0.069491f, 0.993768f},
  };
  transpose_m3(matrix_a);
  EXPECT_NEAR(-1.0f, determinant_m3_array(matrix_a), 1e-6);

  /* This matrix represents R=(0, 0, 0), S=(-1, 0, 0) */
  float matrix_b[3][3] = {
      {-1.0f, 0.0f, 0.0f},
      {0.0f, 1.0f, 0.0f},
      {0.0f, 0.0f, 1.0f},
  };
  transpose_m3(matrix_b);

  float result[3][3];
  interp_m3_m3m3(result, matrix_a, matrix_b, 0.0f);
  EXPECT_M3_NEAR(result, matrix_a, 1e-5);

  interp_m3_m3m3(result, matrix_a, matrix_b, 1.0f);
  EXPECT_M3_NEAR(result, matrix_b, 1e-5);

  interp_m3_m3m3(result, matrix_a, matrix_b, 0.5f);
  float expect[3][3] = {
      {-0.997681f, -0.049995f, 0.046186f},
      {-0.051473f, 0.998181f, -0.031385f},
      {0.044533f, 0.033689f, 0.998440f},
  };
  transpose_m3(expect);
  EXPECT_M3_NEAR(result, expect, 1e-5);

  /* Interpolating between a matrix with and without axis flip can cause it to go through a zero
   * point. The determinant det(A) of a matrix represents the change in volume; interpolating
   * between matrices with det(A)=-1 and det(B)=1 will have to go through a point where
   * det(result)=0, so where the volume becomes zero. */
  float matrix_i[3][3];
  unit_m3(matrix_i);
  zero_m3(expect);
  interp_m3_m3m3(result, matrix_a, matrix_i, 0.5f);
  EXPECT_NEAR(0.0f, determinant_m3_array(result), 1e-5);
  EXPECT_M3_NEAR(result, expect, 1e-5);
}

namespace blender::tests {

using namespace blender::math;

TEST(math_matrix, MatrixInverse)
{
  float3x3 mat(2);
  float3x3 inv = inverse(mat);
  EXPECT_NEAR(inv[0][0], 0.5f, 1e-8f);
  EXPECT_NEAR(inv[0][1], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[0][2], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[1][0], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[1][1], 0.5f, 1e-8f);
  EXPECT_NEAR(inv[1][2], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[2][0], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[2][1], 0.0f, 1e-8f);
  EXPECT_NEAR(inv[2][2], 0.5f, 1e-8f);
}

TEST(math_matrix, MatrixDeterminant)
{
  float2x2 m2({1, 2}, {3, 4});
  float3x3 m3({1, 2, 3}, {-3, 4, -5}, {5, -6, 7});
  float4x4 m4({1, 2, -3, 3}, {3, 4, -5, 3}, {5, 6, 7, -3}, {5, 6, 7, 1});
  EXPECT_NEAR(determinant(m2), -2.0f, 1e-8f);
  EXPECT_NEAR(determinant(m3), -16.0f, 1e-8f);
  EXPECT_NEAR(determinant(m4), -112.0f, 1e-8f);
  EXPECT_NEAR(determinant(double2x2(m2)), -2.0f, 1e-8f);
  EXPECT_NEAR(determinant(double3x3(m3)), -16.0f, 1e-8f);
  EXPECT_NEAR(determinant(double4x4(m4)), -112.0f, 1e-8f);
}

TEST(math_matrix, MatrixAccess)
{
  float4x4 m({1, 2, 3, 4}, {5, 6, 7, 8}, {9, 1, 2, 3}, {4, 5, 6, 7});
  /** Access helpers. */
  EXPECT_EQ(m.forward(), float3(1, 2, 3));
  EXPECT_EQ(m.right(), float3(5, 6, 7));
  EXPECT_EQ(m.up(), float3(9, 1, 2));
  EXPECT_EQ(m.location(), float3(4, 5, 6));
}

TEST(math_matrix, MatrixInit)
{
  float4x4 expect;

  float4x4 m = mat4x4::from_location<float>({1, 2, 3});
  expect = float4x4({1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {1, 2, 3, 1});
  EXPECT_TRUE(compare(m, expect, 0.00001f));

  m = mat4x4::from_rotation<float>({1, 2, 3});
  expect = float4x4({0.411982, -0.0587266, -0.909297, 0},
                    {-0.833738, -0.426918, -0.350175, 0},
                    {-0.36763, 0.902382, -0.224845, 0},
                    {0, 0, 0, 1});
  EXPECT_TRUE(compare(m, expect, 0.00001f));

  m = mat4x4::from_scale<float>({1, 2, 3});
  expect = float4x4({1, 0, 0, 0}, {0, 2, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 1});
  EXPECT_TRUE(compare(m, expect, 0.00001f));

  m = mat4x4::from_loc_rot<float>({1, 2, 3}, {1, 2, 3});
  expect = float4x4({0.411982, -0.0587266, -0.909297, 0},
                    {-0.833738, -0.426918, -0.350175, 0},
                    {-0.36763, 0.902382, -0.224845, 0},
                    {1, 2, 3, 1});
  EXPECT_TRUE(compare(m, expect, 0.00001f));

  m = mat4x4::from_loc_rot_scale<float>({1, 2, 3}, {1, 2, 3}, {1, 2, 3});
  expect = float4x4({0.411982, -0.0587266, -0.909297, 0},
                    {-1.66748, -0.853835, -0.700351, 0},
                    {-1.10289, 2.70714, -0.674535, 0},
                    {1, 2, 3, 1});
  EXPECT_TRUE(compare(m, expect, 0.00001f));
}

TEST(math_matrix, MatrixMethods)
{
  float4x4 m = float4x4({0, 3, 0, 0}, {2, 0, 0, 0}, {0, 0, 2, 0}, {0, 0, 0, 1});
  auto expect_eul = rotation::EulerXYZ<float>(0, 0, M_PI_2);
  EXPECT_V3_NEAR(to_euler(m), expect_eul, 0.0002f);
  auto expect_qt = rotation::Quaternion<float>(0, -M_SQRT1_2, M_SQRT1_2, 0);
  EXPECT_V4_NEAR(to_quaternion(m), expect_qt, 0.0002f);
  EXPECT_EQ(to_scale(m), float3(3, 2, 2));
  EXPECT_TRUE(is_negative(m));
  EXPECT_FALSE(is_unit_scale(m));
  m = normalize(m);
  EXPECT_TRUE(is_unit_scale(m));
}

TEST(math_matrix, MatrixTranspose)
{
  float4x4 m({1, 2, 3, 4}, {5, 6, 7, 8}, {9, 1, 2, 3}, {2, 5, 6, 7});
  float4x4 expect({1, 5, 9, 2}, {2, 6, 1, 5}, {3, 7, 2, 6}, {4, 8, 3, 7});
  EXPECT_EQ(transpose(m), expect);
}

TEST(math_matrix, MatrixInterpolationRegular)
{
  /* Test 4x4 matrix interpolation without singularity, i.e. without axis flip. */

  /* Transposed matrix, so that the code here is written in the same way as print_m4() outputs. */
  /* This matrix represents T=(0.1, 0.2, 0.3), R=(40, 50, 60) degrees, S=(0.7, 0.8, 0.9) */
  float4x4 m2 = transpose(float4x4({0.224976f, -0.333770f, 0.765074f, 0.100000f},
                                   {0.389669f, 0.647565f, 0.168130f, 0.200000f},
                                   {-0.536231f, 0.330541f, 0.443163f, 0.300000f},
                                   {0.000000f, 0.000000f, 0.000000f, 1.000000f}));
  float4x4 m1 = float4x4::identity();
  float4x4 result;
  const float epsilon = 1e-6;
  result = interpolate(m1, m2, 0.0f);
  EXPECT_M4_NEAR(result, m1, epsilon);
  result = interpolate(m1, m2, 1.0f);
  EXPECT_M4_NEAR(result, m2, epsilon);

  /* This matrix is based on the current implementation of the code, and isn't guaranteed to be
   * correct. It's just consistent with the current implementation. */
  float4x4 expect = transpose(float4x4({0.690643f, -0.253244f, 0.484996f, 0.050000f},
                                       {0.271924f, 0.852623f, 0.012348f, 0.100000f},
                                       {-0.414209f, 0.137484f, 0.816778f, 0.150000f},
                                       {0.000000f, 0.000000f, 0.000000f, 1.000000f}));
  result = interpolate(m1, m2, 0.5f);
  EXPECT_M4_NEAR(result, expect, epsilon);
}

TEST(math_matrix, MatrixInterpolationSingularity)
{
  /* A singularity means that there is an axis mirror in the rotation component of the matrix.
   * This is reflected in its negative determinant.
   *
   * The interpolation of 4x4 matrices performs linear interpolation on the translation component,
   * and then uses the 3x3 interpolation function to handle rotation and scale. As a result, this
   * test for a singularity in the rotation matrix only needs to test the 3x3 case. */

  /* Transposed matrix, so that the code here is written in the same way as print_m4() outputs. */
  /* This matrix represents R=(4, 5, 6) degrees, S=(-1, 1, 1) */
  float3x3 matrix_a = transpose(float3x3({-0.990737f, -0.098227f, 0.093759f},
                                         {-0.104131f, 0.992735f, -0.060286f},
                                         {0.087156f, 0.069491f, 0.993768f}));
  EXPECT_NEAR(-1.0f, determinant(matrix_a), 1e-6);

  /* This matrix represents R=(0, 0, 0), S=(-1, 0, 0) */
  float3x3 matrix_b = transpose(
      float3x3({-1.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f}, {0.0f, 0.0f, 1.0f}));

  float3x3 result = interpolate(matrix_a, matrix_b, 0.0f);
  EXPECT_M3_NEAR(result, matrix_a, 1e-5);

  result = interpolate(matrix_a, matrix_b, 1.0f);
  EXPECT_M3_NEAR(result, matrix_b, 1e-5);

  result = interpolate(matrix_a, matrix_b, 0.5f);

  float3x3 expect = transpose(float3x3({-0.997681f, -0.049995f, 0.046186f},
                                       {-0.051473f, 0.998181f, -0.031385f},
                                       {0.044533f, 0.033689f, 0.998440f}));
  EXPECT_M3_NEAR(result, expect, 1e-5);

  /* Interpolating between a matrix with and without axis flip can cause it to go through a zero
   * point. The determinant det(A) of a matrix represents the change in volume; interpolating
   * between matrices with det(A)=-1 and det(B)=1 will have to go through a point where
   * det(result)=0, so where the volume becomes zero. */
  float3x3 matrix_i = float3x3::identity();
  expect = float3x3::zero();
  result = interpolate(matrix_a, matrix_i, 0.5f);
  EXPECT_NEAR(0.0f, determinant(result), 1e-5);
  EXPECT_M3_NEAR(result, expect, 1e-5);
}

}  // namespace blender::tests
