/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "BLI_math_matrix.hh"
#include "BLI_math_matrix_types.hh"

namespace blender::tests {

using float2x2 = mat_base<float, 2, 2>;
using float3x2 = mat_base<float, 3, 2>;
using double3x2 = mat_base<double, 3, 2>;

using namespace blender::math;

TEST(math_matrix_types, ScalarConstructor)
{
  float2x2 m(5.0f);
  EXPECT_EQ(m[0][0], 5.0f);
  EXPECT_EQ(m[1][1], 5.0f);
  EXPECT_EQ(m[0][1], 0.0f);
  EXPECT_EQ(m[1][0], 0.0f);
}

TEST(math_matrix_types, VectorConstructor)
{
  float3x2 m({1.0f, 2.0f}, {3.0f, 4.0f}, {5.0f, 6.0f});
  EXPECT_EQ(m[0][0], 1.0f);
  EXPECT_EQ(m[0][1], 2.0f);
  EXPECT_EQ(m[1][0], 3.0f);
  EXPECT_EQ(m[1][1], 4.0f);
  EXPECT_EQ(m[2][0], 5.0f);
  EXPECT_EQ(m[2][1], 6.0f);
}

TEST(math_matrix_types, SmallerMatrixConstructor)
{
  float2x2 m2({1.0f, 2.0f}, {3.0f, 4.0f});
  float3x3 m3(m2);
  EXPECT_EQ(m3[0][0], 1.0f);
  EXPECT_EQ(m3[0][1], 2.0f);
  EXPECT_EQ(m3[0][2], 0.0f);
  EXPECT_EQ(m3[1][0], 3.0f);
  EXPECT_EQ(m3[1][1], 4.0f);
  EXPECT_EQ(m3[1][2], 0.0f);
  EXPECT_EQ(m3[2][0], 0.0f);
  EXPECT_EQ(m3[2][1], 0.0f);
  EXPECT_EQ(m3[2][2], 1.0f);
}

TEST(math_matrix_types, ComponentMasking)
{
  float3x3 m3({1.1f, 1.2f, 1.3f}, {2.1f, 2.2f, 2.3f}, {3.1f, 3.2f, 3.3f});
  float2x2 m2(m3);
  EXPECT_EQ(m2[0][0], 1.1f);
  EXPECT_EQ(m2[0][1], 1.2f);
  EXPECT_EQ(m2[1][0], 2.1f);
  EXPECT_EQ(m2[1][1], 2.2f);
}

TEST(math_matrix_types, PointerConversion)
{
  float array[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  float2x2 m2(array);
  EXPECT_EQ(m2[0][0], 1.0f);
  EXPECT_EQ(m2[0][1], 2.0f);
  EXPECT_EQ(m2[1][0], 3.0f);
  EXPECT_EQ(m2[1][1], 4.0f);
}

TEST(math_matrix_types, TypeConversion)
{
  float3x2 m(double3x2({1.0f, 2.0f}, {3.0f, 4.0f}, {5.0f, 6.0f}));
  EXPECT_EQ(m[0][0], 1.0f);
  EXPECT_EQ(m[0][1], 2.0f);
  EXPECT_EQ(m[1][0], 3.0f);
  EXPECT_EQ(m[1][1], 4.0f);
  EXPECT_EQ(m[2][0], 5.0f);
  EXPECT_EQ(m[2][1], 6.0f);

  double3x2 d(m);
  EXPECT_EQ(d[0][0], 1.0f);
  EXPECT_EQ(d[0][1], 2.0f);
  EXPECT_EQ(d[1][0], 3.0f);
  EXPECT_EQ(d[1][1], 4.0f);
  EXPECT_EQ(d[2][0], 5.0f);
  EXPECT_EQ(d[2][1], 6.0f);
}

TEST(math_matrix_types, PointerArrayConversion)
{
  float array[2][2] = {{1.0f, 2.0f}, {3.0f, 4.0f}};
  float(*ptr)[2] = array;
  float2x2 m2(ptr);
  EXPECT_EQ(m2[0][0], 1.0f);
  EXPECT_EQ(m2[0][1], 2.0f);
  EXPECT_EQ(m2[1][0], 3.0f);
  EXPECT_EQ(m2[1][1], 4.0f);
}

TEST(math_matrix_types, ComponentAccess)
{
  float3x3 m3({1.1f, 1.2f, 1.3f}, {2.1f, 2.2f, 2.3f}, {3.1f, 3.2f, 3.3f});
  EXPECT_EQ(m3.x.x, 1.1f);
  EXPECT_EQ(m3.x.y, 1.2f);
  EXPECT_EQ(m3.y.x, 2.1f);
  EXPECT_EQ(m3.y.y, 2.2f);
}

TEST(math_matrix_types, AddOperator)
{
  float3x3 m3({1.1f, 1.2f, 1.3f}, {2.1f, 2.2f, 2.3f}, {3.1f, 3.2f, 3.3f});

  m3 = m3 + float3x3(2);
  EXPECT_EQ(m3[0][0], 3.1f);
  EXPECT_EQ(m3[0][2], 1.3f);
  EXPECT_EQ(m3[2][0], 3.1f);
  EXPECT_EQ(m3[2][2], 5.3f);

  m3 += float3x3(-1.0f);
  EXPECT_EQ(m3[0][0], 2.1f);
  EXPECT_EQ(m3[0][2], 1.3f);
  EXPECT_EQ(m3[2][0], 3.1f);
  EXPECT_EQ(m3[2][2], 4.3f);

  m3 += 1.0f;
  EXPECT_EQ(m3[0][0], 3.1f);
  EXPECT_EQ(m3[0][2], 2.3f);
  EXPECT_EQ(m3[2][0], 4.1f);
  EXPECT_EQ(m3[2][2], 5.3f);

  m3 = m3 + 1.0f;
  EXPECT_EQ(m3[0][0], 4.1f);
  EXPECT_EQ(m3[0][2], 3.3f);
  EXPECT_EQ(m3[2][0], 5.1f);
  EXPECT_EQ(m3[2][2], 6.3f);

  m3 = 1.0f + m3;
  EXPECT_EQ(m3[0][0], 5.1f);
  EXPECT_EQ(m3[0][2], 4.3f);
  EXPECT_EQ(m3[2][0], 6.1f);
  EXPECT_EQ(m3[2][2], 7.3f);
}

TEST(math_matrix_types, SubtractOperator)
{
  float3x3 m3({10.0f, 10.2f, 10.3f}, {20.1f, 20.2f, 20.3f}, {30.1f, 30.2f, 30.3f});

  m3 = m3 - float3x3(2);
  EXPECT_EQ(m3[0][0], 8.0f);
  EXPECT_EQ(m3[0][2], 10.3f);
  EXPECT_EQ(m3[2][0], 30.1f);
  EXPECT_EQ(m3[2][2], 28.3f);

  m3 -= float3x3(-1.0f);
  EXPECT_EQ(m3[0][0], 9.0f);
  EXPECT_EQ(m3[0][2], 10.3f);
  EXPECT_EQ(m3[2][0], 30.1f);
  EXPECT_EQ(m3[2][2], 29.3f);

  m3 -= 1.0f;
  EXPECT_EQ(m3[0][0], 8.0f);
  EXPECT_EQ(m3[0][2], 9.3f);
  EXPECT_EQ(m3[2][0], 29.1f);
  EXPECT_EQ(m3[2][2], 28.3f);

  m3 = m3 - 1.0f;
  EXPECT_EQ(m3[0][0], 7.0f);
  EXPECT_EQ(m3[0][2], 8.3f);
  EXPECT_EQ(m3[2][0], 28.1f);
  EXPECT_EQ(m3[2][2], 27.3f);

  m3 = 1.0f - m3;
  EXPECT_EQ(m3[0][0], -6.0f);
  EXPECT_EQ(m3[0][2], -7.3f);
  EXPECT_EQ(m3[2][0], -27.1f);
  EXPECT_EQ(m3[2][2], -26.3f);
}

TEST(math_matrix_types, MultiplyOperator)
{
  float3x3 m3(float3(1.0f), float3(2.0f), float3(2.0f));

  m3 = m3 * 2;
  EXPECT_EQ(m3[0][0], 2.0f);
  EXPECT_EQ(m3[2][2], 4.0f);

  m3 = 2 * m3;
  EXPECT_EQ(m3[0][0], 4.0f);
  EXPECT_EQ(m3[2][2], 8.0f);

  m3 *= 2;
  EXPECT_EQ(m3[0][0], 8.0f);
  EXPECT_EQ(m3[2][2], 16.0f);
}

TEST(math_matrix_types, MatrixMultiplyOperator)
{
  float2x2 a(float2(1, 2), float2(3, 4));
  float2x2 b(float2(5, 6), float2(7, 8));

  float2x2 result = a * b;
  EXPECT_EQ(result[0][0], b[0][0] * a[0][0] + b[0][1] * a[1][0]);
  EXPECT_EQ(result[0][1], b[0][0] * a[0][1] + b[0][1] * a[1][1]);
  EXPECT_EQ(result[1][0], b[1][0] * a[0][0] + b[1][1] * a[1][0]);
  EXPECT_EQ(result[1][1], b[1][0] * a[0][1] + b[1][1] * a[1][1]);

  result = a;
  result *= b;
  EXPECT_EQ(result[0][0], b[0][0] * a[0][0] + b[0][1] * a[1][0]);
  EXPECT_EQ(result[0][1], b[0][0] * a[0][1] + b[0][1] * a[1][1]);
  EXPECT_EQ(result[1][0], b[1][0] * a[0][0] + b[1][1] * a[1][0]);
  EXPECT_EQ(result[1][1], b[1][0] * a[0][1] + b[1][1] * a[1][1]);

  /* Test SSE2 implementation. */
  float4x4 result2 = float4x4(2) * float4x4(6);
  EXPECT_EQ(result2, float4x4(12));

  float3x3 result3 = float3x3(2) * float3x3(6);
  EXPECT_EQ(result3, float3x3(12));
}

TEST(math_matrix_types, VectorMultiplyOperator)
{
  float3x2 mat(float2(1, 2), float2(3, 4), float2(5, 6));

  float2 result = mat * float3(7, 8, 9);
  EXPECT_EQ(result[0], 1 * 7 + 3 * 8 + 5 * 9);
  EXPECT_EQ(result[1], 2 * 7 + 4 * 8 + 6 * 9);
}

}  // namespace blender::tests
