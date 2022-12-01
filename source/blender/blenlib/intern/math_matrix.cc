/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup bli
 */

#include "BLI_math_matrix.hh"
#include "BLI_math_rotation_new.hh"

/* Eigen gives annoying huge amount of warnings here, silence them! */
#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic ignored "-Wlogical-op"
#endif

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using Eigen::Map;
using Eigen::Matrix;
using Eigen::Stride;

/* -------------------------------------------------------------------- */
/** \name Matrix multiplication
 * \{ */

namespace blender {

template<typename T, int NumCol, int NumRow>
mat_base<T, NumCol, NumRow> operator*(const mat_base<T, NumCol, NumRow> &a,
                                      const mat_base<T, NumCol, NumRow> &b)
{
  /* This is the reference implementation.
   * Subclass are free to overload it with vectorized / optimized code. */
  /** TODO(fclem): Only tested for square matrices. Might still contain bugs. */
  mat_base<T, NumCol, NumRow> result(0);
  unroll<NumCol>([&](auto c) {
    unroll<NumRow>([&](auto r) {
      /* This is vector multiplication. */
      result[c] += b[c][r] * a[r];
    });
  });
  return result;
}

template float2x2 operator*(const float2x2 &a, const float2x2 &b);
template float3x3 operator*(const float3x3 &a, const float3x3 &b);
template double2x2 operator*(const double2x2 &a, const double2x2 &b);
template double3x3 operator*(const double3x3 &a, const double3x3 &b);
template double4x4 operator*(const double4x4 &a, const double4x4 &b);

#ifdef BLI_HAVE_SSE2
template<> float4x4 operator*(const float4x4 &a, const float4x4 &b)
{
  float4x4 result;

  __m128 A0 = _mm_load_ps(a[0]);
  __m128 A1 = _mm_load_ps(a[1]);
  __m128 A2 = _mm_load_ps(a[2]);
  __m128 A3 = _mm_load_ps(a[3]);

  for (int i = 0; i < 4; i++) {
    __m128 B0 = _mm_set1_ps(b[i][0]);
    __m128 B1 = _mm_set1_ps(b[i][1]);
    __m128 B2 = _mm_set1_ps(b[i][2]);
    __m128 B3 = _mm_set1_ps(b[i][3]);

    __m128 sum = _mm_add_ps(_mm_add_ps(_mm_mul_ps(B0, A0), _mm_mul_ps(B1, A1)),
                            _mm_add_ps(_mm_mul_ps(B2, A2), _mm_mul_ps(B3, A3)));

    _mm_store_ps(result[i], sum);
  }
  return result;
}
#else
template<> float4x4 operator*(const float4x4 &a, const float4x4 &b)
{
  using namespace math;
  float4x4 T = transpose(b);

  float4x4 mat;
  mat.x = float4(dot(a.x, T.x), dot(a.x, T.y), dot(a.x, T.z), dot(a.x, T.w));
  mat.y = float4(dot(a.y, T.x), dot(a.y, T.y), dot(a.y, T.z), dot(a.y, T.w));
  mat.z = float4(dot(a.z, T.x), dot(a.z, T.y), dot(a.z, T.z), dot(a.z, T.w));
  mat.w = float4(dot(a.w, T.x), dot(a.w, T.y), dot(a.w, T.z), dot(a.w, T.w));
  return mat;
}
#endif

}  // namespace blender

/** \} */

namespace blender::math {

/* -------------------------------------------------------------------- */
/** \name Determinant
 * \{ */

template<typename T, int Size> T determinant(const mat_base<T, Size, Size> &mat)
{
  return Map<const Matrix<T, Size, Size>>((const T *)mat.ptr()).determinant();
}

template float determinant(const float2x2 &mat);
template float determinant(const float3x3 &mat);
template float determinant(const float4x4 &mat);
template double determinant(const double2x2 &mat);
template double determinant(const double3x3 &mat);
template double determinant(const double4x4 &mat);

template<typename T> bool is_negative(const mat_base<T, 4, 4> &mat)
{
  return Map<const Matrix<T, 3, 3>, 0, Stride<4, 1>>((const T *)mat.ptr()).determinant() < T(0);
}

template bool is_negative(const float4x4 &mat);
template bool is_negative(const double4x4 &mat);

/** \} */

/* -------------------------------------------------------------------- */
/** \name Inverse
 * \{ */

template<typename T, int Size> mat_base<T, Size, Size> invert(const mat_base<T, Size, Size> &mat)
{
  mat_base<T, Size, Size> result;
  Map<const Matrix<T, Size, Size>> M((const T *)mat.ptr());
  Map<Matrix<T, Size, Size>> R((T *)result.ptr());
  bool is_invertible = true;
  M.computeInverseWithCheck(R, is_invertible, 0.0f);
  if (!is_invertible) {
    R = R.Zero();
  }
  return result;
}

template float2x2 invert(const float2x2 &mat);
template float3x3 invert(const float3x3 &mat);
template float4x4 invert(const float4x4 &mat);
template double2x2 invert(const double2x2 &mat);
template double3x3 invert(const double3x3 &mat);
template double4x4 invert(const double4x4 &mat);

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
  mat_base<T, 3, 3> W, V;
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
        Map<const MatrixDynamicT>((const T *)mat3.ptr(), 3, 3), ComputeThinU | ComputeThinV);

    (Map<MatrixT>((T *)W.ptr())) = svd.matrixU();
    (Map<VectorT>(S_val)) = svd.singularValues();
    (Map<MatrixT>((T *)V.ptr())) = svd.matrixV();
  }

  mat_base<T, 3, 3> S = from_scale<mat_base<T, 3, 3>>(S_val);
  mat_base<T, 3, 3> Vt = transpose(V);

  r_U = W * Vt;
  r_P = (V * S) * Vt;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Interpolate
 * \{ */

template<typename T>
mat_base<T, 3, 3> interpolate(const mat_base<T, 3, 3> &A, const mat_base<T, 3, 3> &B, T t)
{
  using Mat3T = mat_base<T, 3, 3>;
  /* 'Rotation' component ('U' part of polar decomposition,
   * the closest orthogonal matrix to M3 rot/scale
   * transformation matrix), spherically interpolated. */
  Mat3T U_A, U_B;
  /* 'Scaling' component ('P' part of polar decomposition, i.e. scaling in U-defined space),
   * linearly interpolated. */
  Mat3T P_A, P_B;

  polar_decompose(A, U_A, P_A);
  polar_decompose(B, U_B, P_B);

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

  detail::Quaternion<T> quat_A = math::to_quaternion(U_A);
  detail::Quaternion<T> quat_B = math::to_quaternion(U_B);
  detail::Quaternion<T> quat = math::interpolate(quat_A, quat_B, t);
  Mat3T U = from_rotation<Mat3T>(quat);

  Mat3T P = interpolate_linear(P_A, P_B, t);
  /* And we reconstruct rot/scale matrix from interpolated polar components */
  return U * P;
}

template float3x3 interpolate(const float3x3 &A, const float3x3 &B, float t);
template double3x3 interpolate(const double3x3 &A, const double3x3 &B, double t);

template<typename T>
mat_base<T, 4, 4> interpolate(const mat_base<T, 4, 4> &A, const mat_base<T, 4, 4> &B, T t)
{
  mat_base<T, 4, 4> result = mat_base<T, 4, 4>(
      interpolate(mat_base<T, 3, 3>(A), mat_base<T, 3, 3>(B), t));

  /* Location component, linearly interpolated. */
  const auto &loc_a = static_cast<const mat_base<T, 4, 4> &>(A).location();
  const auto &loc_b = static_cast<const mat_base<T, 4, 4> &>(B).location();
  result.location() = interpolate(loc_a, loc_b, t);

  return result;
}

template float4x4 interpolate(const float4x4 &A, const float4x4 &B, float t);
template double4x4 interpolate(const double4x4 &A, const double4x4 &B, double t);

/** \} */

/* -------------------------------------------------------------------- */
/** \name Matrix to rotation
 * \{ */

namespace detail {

template<typename T>
void normalized_to_eul2(const mat_base<T, 3, 3> &mat,
                        detail::EulerXYZ<T> &eul1,
                        detail::EulerXYZ<T> &eul2)
{
  BLI_assert(math::is_unit_scale(mat));

  const float cy = std::hypot(mat[0][0], mat[0][1]);
  if (cy > 16.0f * FLT_EPSILON) {
    eul1.x = std::atan2(mat[1][2], mat[2][2]);
    eul1.y = std::atan2(-mat[0][2], cy);
    eul1.z = std::atan2(mat[0][1], mat[0][0]);

    eul2.x = std::atan2(-mat[1][2], -mat[2][2]);
    eul2.y = std::atan2(-mat[0][2], -cy);
    eul2.z = std::atan2(-mat[0][1], -mat[0][0]);
  }
  else {
    eul1.x = std::atan2(-mat[2][1], mat[1][1]);
    eul1.y = std::atan2(-mat[0][2], cy);
    eul1.z = 0.0f;

    eul2 = eul1;
  }
}

template void normalized_to_eul2(const float3x3 &mat,
                                 detail::EulerXYZ<float> &eul1,
                                 detail::EulerXYZ<float> &eul2);
template void normalized_to_eul2(const double3x3 &mat,
                                 detail::EulerXYZ<double> &eul1,
                                 detail::EulerXYZ<double> &eul2);

template<typename T> detail::Quaternion<T> normalized_to_quat_fast(const mat_base<T, 3, 3> &mat)
{
  BLI_assert(math::is_unit_scale(mat));
  /* Caller must ensure matrices aren't negative for valid results, see: T24291, T94231. */
  BLI_assert(!math::is_negative(mat));

  detail::Quaternion<T> q;

  /* Method outlined by Mike Day, ref: https://math.stackexchange.com/a/3183435/220949
   * with an additional `sqrtf(..)` for higher precision result.
   * Removing the `sqrt` causes tests to fail unless the precision is set to 1e-6 or larger. */

  if (mat[2][2] < 0.0f) {
    if (mat[0][0] > mat[1][1]) {
      const T trace = 1.0f + mat[0][0] - mat[1][1] - mat[2][2];
      T s = 2.0f * std::sqrt(trace);
      if (mat[1][2] < mat[2][1]) {
        /* Ensure W is non-negative for a canonical result. */
        s = -s;
      }
      q.y = 0.25f * s;
      s = 1.0f / s;
      q.x = (mat[1][2] - mat[2][1]) * s;
      q.z = (mat[0][1] + mat[1][0]) * s;
      q.w = (mat[2][0] + mat[0][2]) * s;
      if (UNLIKELY((trace == 1.0f) && (q.x == 0.0f && q.z == 0.0f && q.w == 0.0f))) {
        /* Avoids the need to normalize the degenerate case. */
        q.y = 1.0f;
      }
    }
    else {
      const T trace = 1.0f - mat[0][0] + mat[1][1] - mat[2][2];
      T s = 2.0f * std::sqrt(trace);
      if (mat[2][0] < mat[0][2]) {
        /* Ensure W is non-negative for a canonical result. */
        s = -s;
      }
      q.z = 0.25f * s;
      s = 1.0f / s;
      q.x = (mat[2][0] - mat[0][2]) * s;
      q.y = (mat[0][1] + mat[1][0]) * s;
      q.w = (mat[1][2] + mat[2][1]) * s;
      if (UNLIKELY((trace == 1.0f) && (q.x == 0.0f && q.y == 0.0f && q.w == 0.0f))) {
        /* Avoids the need to normalize the degenerate case. */
        q.z = 1.0f;
      }
    }
  }
  else {
    if (mat[0][0] < -mat[1][1]) {
      const T trace = 1.0f - mat[0][0] - mat[1][1] + mat[2][2];
      T s = 2.0f * std::sqrt(trace);
      if (mat[0][1] < mat[1][0]) {
        /* Ensure W is non-negative for a canonical result. */
        s = -s;
      }
      q.w = 0.25f * s;
      s = 1.0f / s;
      q.x = (mat[0][1] - mat[1][0]) * s;
      q.y = (mat[2][0] + mat[0][2]) * s;
      q.z = (mat[1][2] + mat[2][1]) * s;
      if (UNLIKELY((trace == 1.0f) && (q.x == 0.0f && q.y == 0.0f && q.z == 0.0f))) {
        /* Avoids the need to normalize the degenerate case. */
        q.w = 1.0f;
      }
    }
    else {
      /* NOTE(@campbellbarton): A zero matrix will fall through to this block,
       * needed so a zero scaled matrices to return a quaternion without rotation, see: T101848.
       */
      const T trace = 1.0f + mat[0][0] + mat[1][1] + mat[2][2];
      T s = 2.0f * std::sqrt(trace);
      q.x = 0.25f * s;
      s = 1.0f / s;
      q.y = (mat[1][2] - mat[2][1]) * s;
      q.z = (mat[2][0] - mat[0][2]) * s;
      q.w = (mat[0][1] - mat[1][0]) * s;
      if (UNLIKELY((trace == 1.0f) && (q.y == 0.0f && q.z == 0.0f && q.w == 0.0f))) {
        /* Avoids the need to normalize the degenerate case. */
        q.x = 1.0f;
      }
    }
  }

  BLI_assert(!(q.x < 0.0f));
  BLI_assert(math::is_unit_scale(vec_base<T, 4>(q)));
  return q;
}

template<typename T>
detail::Quaternion<T> normalized_to_quat_with_checks(const mat_base<T, 3, 3> &mat)
{
  const T det = math::determinant(mat);
  if (UNLIKELY(!isfinite(det))) {
    return detail::Quaternion<T>::identity();
  }
  else if (UNLIKELY(det < T(0))) {
    return normalized_to_quat_fast(-mat);
  }
  return normalized_to_quat_fast(mat);
}

template Quaternion<float> normalized_to_quat_with_checks(const float3x3 &mat);
template Quaternion<double> normalized_to_quat_with_checks(const double3x3 &mat);

}  // namespace detail

/** \} */

/* -------------------------------------------------------------------- */
/** \name Init Helpers.
 * \{ */

namespace detail {

template<typename T, int NumCol, int NumRow>
mat_base<T, NumCol, NumRow> from_rotation(const EulerXYZ<T> &rotation)
{
  using MatT = mat_base<T, NumCol, NumRow>;
  using IntermediateType = typename TypeTraits<T>::IntermediateType;
  IntermediateType ci = math::cos(rotation.x);
  IntermediateType cj = math::cos(rotation.y);
  IntermediateType ch = math::cos(rotation.z);
  IntermediateType si = math::sin(rotation.x);
  IntermediateType sj = math::sin(rotation.y);
  IntermediateType sh = math::sin(rotation.z);
  IntermediateType cc = ci * ch;
  IntermediateType cs = ci * sh;
  IntermediateType sc = si * ch;
  IntermediateType ss = si * sh;

  MatT mat;
  mat[0][0] = T(cj * ch);
  mat[1][0] = T(sj * sc - cs);
  mat[2][0] = T(sj * cc + ss);

  mat[0][1] = T(cj * sh);
  mat[1][1] = T(sj * ss + cc);
  mat[2][1] = T(sj * cs - sc);

  mat[0][2] = T(-sj);
  mat[1][2] = T(cj * si);
  mat[2][2] = T(cj * ci);
  return mat;
}

template<typename T, int NumCol, int NumRow>
mat_base<T, NumCol, NumRow> from_rotation(const Quaternion<T> &rotation)
{
  using MatT = mat_base<T, NumCol, NumRow>;
  using IntermediateType = typename TypeTraits<T>::IntermediateType;
  IntermediateType q0 = M_SQRT2 * IntermediateType(rotation.x);
  IntermediateType q1 = M_SQRT2 * IntermediateType(rotation.y);
  IntermediateType q2 = M_SQRT2 * IntermediateType(rotation.z);
  IntermediateType q3 = M_SQRT2 * IntermediateType(rotation.w);

  IntermediateType qda = q0 * q1;
  IntermediateType qdb = q0 * q2;
  IntermediateType qdc = q0 * q3;
  IntermediateType qaa = q1 * q1;
  IntermediateType qab = q1 * q2;
  IntermediateType qac = q1 * q3;
  IntermediateType qbb = q2 * q2;
  IntermediateType qbc = q2 * q3;
  IntermediateType qcc = q3 * q3;

  MatT mat;
  mat[0][0] = T(1.0 - qbb - qcc);
  mat[0][1] = T(qdc + qab);
  mat[0][2] = T(-qdb + qac);

  mat[1][0] = T(-qdc + qab);
  mat[1][1] = T(1.0 - qaa - qcc);
  mat[1][2] = T(qda + qbc);

  mat[2][0] = T(qdb + qac);
  mat[2][1] = T(-qda + qbc);
  mat[2][2] = T(1.0 - qaa - qbb);
  return mat;
}

template<typename T, int NumCol, int NumRow>
mat_base<T, NumCol, NumRow> from_rotation(const AxisAngle<T> &rotation)
{
  using MatT = mat_base<T, NumCol, NumRow>;
  using Vec3T = typename MatT::vec3_type;
  const T angle_sin = rotation.angle_sin();
  const T angle_cos = rotation.angle_cos();
  const Vec3T &axis = rotation.axis();

  BLI_assert(is_unit_scale(axis));

  T ico = (T(1) - angle_cos);
  Vec3T nsi = axis * angle_sin;

  Vec3T n012 = (axis * axis) * ico;
  T n_01 = (axis[0] * axis[1]) * ico;
  T n_02 = (axis[0] * axis[2]) * ico;
  T n_12 = (axis[1] * axis[2]) * ico;

  MatT mat = from_scale<MatT>(n012 + angle_cos);
  mat[0][1] = n_01 + nsi[2];
  mat[0][2] = n_02 - nsi[1];
  mat[1][0] = n_01 - nsi[2];
  mat[1][2] = n_12 + nsi[0];
  mat[2][0] = n_02 + nsi[1];
  mat[2][1] = n_12 - nsi[0];
  return mat;
}

template mat_base<float, 3, 3> from_rotation(const EulerXYZ<float> &rotation);
template mat_base<float, 4, 4> from_rotation(const EulerXYZ<float> &rotation);
template mat_base<float, 3, 3> from_rotation(const Quaternion<float> &rotation);
template mat_base<float, 4, 4> from_rotation(const Quaternion<float> &rotation);
template mat_base<float, 3, 3> from_rotation(const AxisAngle<float> &rotation);
template mat_base<float, 4, 4> from_rotation(const AxisAngle<float> &rotation);

}  // namespace detail

/** \} */

/* -------------------------------------------------------------------- */
/** \name Transform function.
 * \{ */

template<typename T>
vec_base<T, 3> transform_point(const mat_base<T, 3, 3> &mat, const vec_base<T, 3> &point)
{
  return mat * point;
}

template<typename T>
vec_base<T, 3> transform_point(const mat_base<T, 4, 4> &mat, const vec_base<T, 3> &point)
{
  /* TODO(fclem): mat3 view. */
  return vec_base<T, 3>(mat * vec_base<T, 4>(point, T(1)));
}

template<typename T>
vec_base<T, 3> transform_direction(const mat_base<T, 3, 3> &mat, const vec_base<T, 3> &direction)
{
  return mat * direction;
}

template<typename T>
vec_base<T, 3> transform_direction(const mat_base<T, 4, 4> &mat, const vec_base<T, 3> &direction)
{
  /* TODO(fclem): mat3 view. */
  return vec_base<T, 3>(mat * vec_base<T, 4>(direction, T(0)));
}

template<typename T>
vec_base<T, 3> project_point(const mat_base<T, 4, 4> &mat, const vec_base<T, 3> &point)
{
  vec_base<T, 4> tmp(point, T(0));
  tmp = mat * tmp;
  tmp /= tmp.w;
  return vec_base<T, 3>(tmp);
}

template float3 transform_point(const float3x3 &mat, const float3 &point);
template float3 transform_point(const float4x4 &mat, const float3 &point);
template float3 transform_direction(const float3x3 &mat, const float3 &direction);
template float3 transform_direction(const float4x4 &mat, const float3 &direction);
template float3 project_point(const float4x4 &mat, const float3 &direction);

/** \} */

/* -------------------------------------------------------------------- */
/** \name Projection Matrices.
 * \{ */

namespace projection {

template<typename T>
mat_base<T, 4, 4> orthographic(T left, T right, T bottom, T top, T near_clip, T far_clip)
{
  const T x_delta = right - left;
  const T y_delta = top - bottom;
  const T z_delta = far_clip - near_clip;

  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>::identity();
  if (x_delta != 0 && y_delta != 0 && z_delta != 0) {
    mat[0][0] = T(2.0) / x_delta;
    mat[3][0] = -(right + left) / x_delta;
    mat[1][1] = T(2.0) / y_delta;
    mat[3][1] = -(top + bottom) / y_delta;
    mat[2][2] = -T(2.0) / z_delta; /* NOTE: negate Z. */
    mat[3][2] = -(far_clip + near_clip) / z_delta;
  }
  return mat;
}

template<typename T>
mat_base<T, 4, 4> perspective(T left, T right, T bottom, T top, T near_clip, T far_clip)
{
  const T x_delta = right - left;
  const T y_delta = top - bottom;
  const T z_delta = far_clip - near_clip;

  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>::identity();
  if (x_delta != 0 && y_delta != 0 && z_delta != 0) {
    mat[0][0] = near_clip * T(2.0) / x_delta;
    mat[1][1] = near_clip * T(2.0) / y_delta;
    mat[2][0] = (right + left) / x_delta; /* NOTE: negate Z. */
    mat[2][1] = (top + bottom) / y_delta;
    mat[2][2] = -(far_clip + near_clip) / z_delta;
    mat[2][3] = -1.0f;
    mat[3][2] = (-2.0f * near_clip * far_clip) / z_delta;
  }
  return mat;
}

template float4x4 orthographic(
    float left, float right, float bottom, float top, float near_clip, float far_clip);
template float4x4 perspective(
    float left, float right, float bottom, float top, float near_clip, float far_clip);

}  // namespace projection

/** \} */

}  // namespace blender::math
