/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup bli
 */

#include "BLI_math_mat_types.hh"

/* Eigen gives annoying huge amount of warnings here, silence them! */
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wlogical-op"
#endif

#include <Eigen/Core>
#include <Eigen/Dense>

namespace blender::math {

using float2x2 = mat_base<float, 2, 2>;
using float3x3 = mat_base<float, 3, 3>;
using float4x4 = mat_base<float, 4, 4>;
using double2x2 = mat_base<double, 2, 2>;
using double3x3 = mat_base<double, 3, 3>;
using double4x4 = mat_base<double, 4, 4>;

using Eigen::Map;
using Eigen::Matrix;

template<> float determinant(const float2x2 &mat)
{
  return Map<const Matrix<float, 2, 2>>((const float *)mat).determinant();
}

template<> float determinant(const float3x3 &mat)
{
  return Map<const Matrix<float, 3, 3>>((const float *)mat).determinant();
}

template<> float determinant(const float4x4 &mat)
{
  return Map<const Matrix<float, 4, 4>>((const float *)mat).determinant();
}

template<> double determinant(const double2x2 &mat)
{
  return Map<const Matrix<double, 2, 2>>((const double *)mat).determinant();
}

template<> double determinant(const double3x3 &mat)
{
  return Map<const Matrix<double, 3, 3>>((const double *)mat).determinant();
}

template<> double determinant(const double4x4 &mat)
{
  return Map<const Matrix<double, 4, 4>>((const double *)mat).determinant();
}

template<typename T, int Size>
inline mat_base<T, Size, Size> inverse_impl(mat_base<T, Size, Size> &mat)
{
  mat_base<T, Size, Size> result;
  Map<Matrix<T, Size, Size>> M((T *)mat);
  Map<Matrix<T, Size, Size>> R((T *)result);
  bool is_invertible = true;
  M.computeInverseWithCheck(R, is_invertible, 0.0f);
  if (!is_invertible) {
    R = R.Zero();
  }
  return result;
}

template<> float2x2 inverse(float2x2 &mat)
{
  return inverse_impl(mat);
}

template<> float3x3 inverse(float3x3 &mat)
{
  return inverse_impl(mat);
}

template<> float4x4 inverse(float4x4 &mat)
{
  return inverse_impl(mat);
}

template<> double2x2 inverse(double2x2 &mat)
{
  return inverse_impl(mat);
}

template<> double3x3 inverse(double3x3 &mat)
{
  return inverse_impl(mat);
}

template<> double4x4 inverse(double4x4 &mat)
{
  return inverse_impl(mat);
}

}  // namespace blender::math
