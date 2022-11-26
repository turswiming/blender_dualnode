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
#include <Eigen/Eigenvalues>

using float2x2 = blender::mat_base<float, 2, 2>;
using float3x3 = blender::mat_base<float, 3, 3>;
using float4x4 = blender::mat_base<float, 4, 4>;
using double2x2 = blender::mat_base<double, 2, 2>;
using double3x3 = blender::mat_base<double, 3, 3>;
using double4x4 = blender::mat_base<double, 4, 4>;

namespace blender::math {

using Eigen::Map;
using Eigen::Matrix;
using Eigen::Stride;

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

template<> bool is_negative(const float4x4 &mat)
{
  /* Don't use determinant(float4x4) as only the 3x3 components are needed
   * when the matrix is used as a transformation to represent location/scale/rotation. */
  return Map<const Matrix<float, 3, 3>, 0, Stride<4, 1>>((const float *)mat).determinant() <
         float(0);
}

template<> bool is_negative(const double4x4 &mat)
{
  /* Don't use determinant(float4x4) as only the 3x3 components are needed
   * when the matrix is used as a transformation to represent location/scale/rotation. */
  return Map<const Matrix<double, 3, 3>, 0, Stride<4, 1>>((const double *)mat).determinant() <
         double(0);
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

/* -------------------------------------------------------------------- */
/** \name Polar Decomposition
 * \{ */

/**
 * Right polar decomposition:
 *     M = UP
 *
 * U is the 'rotation'-like component, the closest orthogonal matrix to M.
 * P is the 'scaling'-like component, defined in U space.
 *
 * See https://en.wikipedia.org/wiki/Polar_decomposition for more.
 */
template<typename T>
static void polar_decompose(const mat_base<T, 3, 3> &mat3,
                            mat_base<T, 3, 3> &r_U,
                            mat_base<T, 3, 3> &r_P)
{
  /* From svd decomposition (M = WSV*), we have:
   *     U = WV*
   *     P = VSV*
   */
  mat_3x3<T> W, V;
  vec_base<T, 3> S_val;

  {
    using namespace Eigen;
    using MatrixT = Matrix<T, 3, 3>;
    using VectorT = Matrix<T, 3, 1>;
    /* Blender and Eigen matrices are both column-major.
     * Since our matrix is squared, we can use thinU/V. */
    /** WORKAROUND:
     * (ComputeThinU | ComputeThinV) must be set as runtime parameters in Eigen < 3.4.0.
     * But this requires the matrix type to be dynamic to avoid an assert.
     */
    using MatrixDynamicT = Matrix<T, Dynamic, Dynamic>;
    JacobiSVD<MatrixDynamicT, NoQRPreconditioner> svd(
        Map<const MatrixDynamicT>((const T *)mat3, 3, 3), ComputeThinU | ComputeThinV);

    (Map<MatrixT>(W)) = svd.matrixU();
    (Map<VectorT>(S_val)) = svd.singularValues();
    (Map<MatrixT>(V)) = svd.matrixV();
  }

  mat_3x3<T> S = mat_3x3<T>::from_scale(S_val);
  mat_3x3<T> Vt = transpose(V);

  r_U = W * Vt;
  r_P = (V * S) * Vt;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Interpolate
 * \{ */

template<typename T>
mat_base<T, 3, 3> interpolate_impl(const mat_base<T, 3, 3> &A, const mat_base<T, 3, 3> &B, T t)
{
  /* 'Rotation' component ('U' part of polar decomposition,
   * the closest orthogonal matrix to M3 rot/scale
   * transformation matrix), spherically interpolated. */
  mat_3x3<T> U_A, U_B;
  /* 'Scaling' component ('P' part of polar decomposition, i.e. scaling in U-defined space),
   * linearly interpolated. */
  mat_3x3<T> P_A, P_B;

  math::polar_decompose(A, U_A, P_A);
  math::polar_decompose(B, U_B, P_B);

  /* Quaternions cannot represent an axis flip. If such a singularity is detected, choose a
   * different decomposition of the matrix that still satisfies A = U_A * P_A but which has a
   * positive determinant and thus no axis flips. This resolves T77154.
   *
   * Note that a flip of two axes is just a rotation of 180 degrees around the third axis, and
   * three flipped axes are just an 180 degree rotation + a single axis flip. It is thus sufficient
   * to solve this problem for single axis flips. */
  if (is_negative(U_A)) {
    U_A = -U_A;
    P_A = -P_A;
  }
  if (is_negative(U_B)) {
    U_B = -U_B;
    P_B = -P_B;
  }

  rotation::Quaternion<T> quat_A = U_A.to_quaternion();
  rotation::Quaternion<T> quat_B = U_B.to_quaternion();
  rotation::Quaternion<T> quat = math::interpolate(quat_A, quat_B, t);
  mat_3x3<T> U = mat_3x3<T>::from_rotation(quat);

  mat_3x3<T> P = math::interpolate_linear(P_A, P_B, t);
  /* And we reconstruct rot/scale matrix from interpolated polar components */
  return U * P;
}

float3x3 interpolate(const float3x3 &A, const float3x3 &B, float t)
{
  return interpolate_impl(A, B, t);
}

double3x3 interpolate(const double3x3 &A, const double3x3 &B, double t)
{
  return interpolate_impl(A, B, t);
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

/** \} */

}  // namespace blender

/** \} */
