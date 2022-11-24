/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "BLI_math_mat_types.hh"

namespace blender::tests {

using float3x3 = mat_base<float, 3, 3>;
using float2x2 = mat_base<float, 2, 2>;
using float3x2 = mat_base<float, 3, 2>;

using namespace blender::math;

TEST(math_mat_types, ScalarConstructor)
{
  float2x2 m(5.0f);
  EXPECT_EQ(m[0][0], 5.0f);
  EXPECT_EQ(m[1][1], 5.0f);
  EXPECT_EQ(m[0][1], 0.0f);
  EXPECT_EQ(m[1][0], 0.0f);
}

TEST(math_mat_types, VectorConstructor)
{
  float3x2 m({1.0f, 2.0f}, {3.0f, 4.0f}, {5.0f, 6.0f});
  EXPECT_EQ(m[0][0], 1.0f);
  EXPECT_EQ(m[0][1], 2.0f);
  EXPECT_EQ(m[1][0], 3.0f);
  EXPECT_EQ(m[1][1], 4.0f);
  EXPECT_EQ(m[2][0], 5.0f);
  EXPECT_EQ(m[2][1], 6.0f);
}

TEST(math_mat_types, ComponentMasking)
{
  float3x3 m3({1.1f, 1.2f, 1.3f}, {2.1f, 2.2f, 2.3f}, {3.1f, 3.2f, 3.3f});
  float2x2 m2(m3);
  EXPECT_EQ(m2[0][0], 1.1f);
  EXPECT_EQ(m2[1][1], 1.2f);
  EXPECT_EQ(m2[0][1], 2.1f);
  EXPECT_EQ(m2[1][0], 2.2f);
}

TEST(math_mat_types, PointerConversion)
{
  float array[4] = {1.0f, 2.0f, 3.0f, 4.0f};
  float2x2 m2(array);
  EXPECT_EQ(m2[0][0], 1.0f);
  EXPECT_EQ(m2[0][1], 2.0f);
  EXPECT_EQ(m2[1][0], 3.0f);
  EXPECT_EQ(m2[1][1], 4.0f);
}

TEST(math_mat_types, PointerArrayConversion)
{
  float array[2][2] = {{1.0f, 2.0f}, {3.0f, 4.0f}};
  float(*ptr)[2] = array;
  float2x2 m2(ptr);
  EXPECT_EQ(m2[0][0], 1.0f);
  EXPECT_EQ(m2[0][1], 2.0f);
  EXPECT_EQ(m2[1][0], 3.0f);
  EXPECT_EQ(m2[1][1], 4.0f);
}

}  // namespace blender::tests
