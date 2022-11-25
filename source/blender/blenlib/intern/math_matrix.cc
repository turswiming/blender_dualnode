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

using float2x2 = blender::mat_base<float, 2, 2>;
using float3x3 = blender::mat_base<float, 3, 3>;
using float4x4 = blender::mat_base<float, 4, 4>;
using double2x2 = blender::mat_base<double, 2, 2>;
using double3x3 = blender::mat_base<double, 3, 3>;
using double4x4 = blender::mat_base<double, 4, 4>;

namespace blender::math {

using Eigen::Map;
using Eigen::Matrix;

/* -------------------------------------------------------------------- */
/** \name Determinant
 * \{ */

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

/** \} */

/* -------------------------------------------------------------------- */
/** \name Inverse
 * \{ */

template<typename T, int Size>
inline mat_base<T, Size, Size> inverse_impl(const mat_base<T, Size, Size> &mat)
{
  mat_base<T, Size, Size> result;
  Map<const Matrix<T, Size, Size>> M((const T *)mat);
  Map<Matrix<T, Size, Size>> R((T *)result);
  bool is_invertible = true;
  M.computeInverseWithCheck(R, is_invertible, 0.0f);
  if (!is_invertible) {
    R = R.Zero();
  }
  return result;
}

template<> float2x2 inverse(const float2x2 &mat)
{
  return inverse_impl(mat);
}

template<> float3x3 inverse(const float3x3 &mat)
{
  return inverse_impl(mat);
}

template<> float4x4 inverse(const float4x4 &mat)
{
  return inverse_impl(mat);
}

template<> double2x2 inverse(const double2x2 &mat)
{
  return inverse_impl(mat);
}

template<> double3x3 inverse(const double3x3 &mat)
{
  return inverse_impl(mat);
}

template<> double4x4 inverse(const double4x4 &mat)
{
  return inverse_impl(mat);
}

/** \} */

}  // namespace blender::math

/* -------------------------------------------------------------------- */
/** \name Matrix multiplication
 *
 * Manual optimization using SSE2.
 *
 * \{ */

namespace blender {

#ifdef BLI_HAVE_SSE2
template<> float3x3 operator*(const float3x3 &a, const float3x3 &b)
{
  union {
    float3x3 result;
    /* Protect from buffer overrun. */
    mat_base<float, 3, 4> _pad0;
  };

  __m128 A0 = _mm_loadu_ps(a[0]);
  __m128 A1 = _mm_loadu_ps(a[1]);
  /** \note This will load one extra float after the end of \arg a. But result is discarded. */
  __m128 A2 = _mm_loadu_ps(a[2]);

  for (int i = 0; i < 3; i++) {
    __m128 B0 = _mm_set1_ps(b[i][0]);
    __m128 B1 = _mm_set1_ps(b[i][1]);
    /** \note This will load one extra float after the end of \arg a. But result is discarded. */
    __m128 B2 = _mm_set1_ps(b[i][2]);

    __m128 sum = _mm_add_ps(_mm_add_ps(_mm_mul_ps(B0, A0), _mm_mul_ps(B1, A1)),
                            _mm_mul_ps(B2, A2));

    /** \note This will load one extra float after the end of `result`. */
    _mm_storeu_ps(result[i], sum);
  }
  return result;
}

template<> float4x4 operator*(const float4x4 &a, const float4x4 &b)
{
  float4x4 result;

  __m128 A0 = _mm_loadu_ps(a[0]);
  __m128 A1 = _mm_loadu_ps(a[1]);
  __m128 A2 = _mm_loadu_ps(a[2]);
  __m128 A3 = _mm_loadu_ps(a[3]);

  for (int i = 0; i < 4; i++) {
    __m128 B0 = _mm_set1_ps(b[i][0]);
    __m128 B1 = _mm_set1_ps(b[i][1]);
    __m128 B2 = _mm_set1_ps(b[i][2]);
    __m128 B3 = _mm_set1_ps(b[i][3]);

    __m128 sum = _mm_add_ps(_mm_add_ps(_mm_mul_ps(B0, A0), _mm_mul_ps(B1, A1)),
                            _mm_add_ps(_mm_mul_ps(B2, A2), _mm_mul_ps(B3, A3)));

    _mm_storeu_ps(result[i], sum);
  }
  return result;
}
#endif

}  // namespace blender

/** \} */
