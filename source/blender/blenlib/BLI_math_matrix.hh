/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup bli
 */

#include "BLI_math_matrix_types.hh"
#include "BLI_math_rotation.hh"
#include "BLI_math_vector.hh"

namespace blender::math {

/* Returns true if each individual columns are unit scaled. Mainly for assert usage. */
template<typename T, int NumCol, int NumRow>
inline bool is_unit_scale(const mat_base<T, NumCol, NumRow> &v)
{
  for (int i = 0; i < NumCol; i++) {
    if (!is_unit_scale(v[i])) {
      return false;
    }
  }
  return true;
}

/* -------------------------------------------------------------------- */
/** \name Matrix Operations
 * \{ */

/**
 * Returns the inverse of a square matrix.
 */
template<typename T, int Size> mat_base<T, Size, Size> inverse(const mat_base<T, Size, Size> &mat);

/**
 * Flip the matrix around its diagonal.
 */
template<typename T, int NumCol, int NumRow>
inline mat_base<T, NumCol, NumRow> transpose(const mat_base<T, NumRow, NumCol> &mat)
{
  mat_base<T, NumCol, NumRow> result;
  for (int i = 0; i < NumCol; i++) {
    for (int j = 0; j < NumRow; j++) {
      result[i][j] = mat[j][i];
    }
  }
  return result;
}

template<typename T, int Size> T determinant(const mat_base<T, Size, Size> &mat);

/** \note linear interpolation for each component. */
template<typename T, int NumCol, int NumRow>
inline mat_base<T, NumCol, NumRow> interpolate_linear(const mat_base<T, NumCol, NumRow> &a,
                                                      const mat_base<T, NumCol, NumRow> &b,
                                                      T t)
{
  mat_base<T, NumCol, NumRow> result;
  unroll<NumCol>([&](auto c) { result[c] = interpolate(a[c], b[c], t); });
  return result;
}

/**
 * A polar-decomposition-based interpolation between matrix A and matrix B.
 *
 * \note This code is about five times slower as the 'naive' interpolation done by #blend_m3_m3m3
 * (it typically remains below 2 usec on an average i74700,
 * while #blend_m3_m3m3 remains below 0.4 usec).
 * However, it gives expected results even with non-uniformly scaled matrices,
 * see T46418 for an example.
 *
 * Based on "Matrix Animation and Polar Decomposition", by Ken Shoemake & Tom Duff
 *
 * \param R: Resulting interpolated matrix.
 * \param A: Input matrix which is totally effective with `t = 0.0`.
 * \param B: Input matrix which is totally effective with `t = 1.0`.
 * \param t: Interpolation factor.
 */
template<typename T>
mat_base<T, 3, 3> interpolate(const mat_base<T, 3, 3> &A, const mat_base<T, 3, 3> &B, T t);

/**
 * Complete transform matrix interpolation,
 * based on polar-decomposition-based interpolation from #interp_m3_m3m3.
 *
 * \param A: Input matrix which is totally effective with `t = 0.0`.
 * \param B: Input matrix which is totally effective with `t = 1.0`.
 * \param t: Interpolation factor.
 */
template<typename T>
mat_base<T, 4, 4> interpolate(const mat_base<T, 4, 4> &A, const mat_base<T, 4, 4> &B, T t);

/**
 * Normalize individually each column of the matrix.
 */
template<typename T> inline T normalize(const T &a)
{
  T result;
  unroll<T::col_len>([&](auto i) { result[i] = math::normalize(a[i]); });
  return result;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Compare / Test
 * \{ */

/**
 * Returns true if matrix has inverted handedness.
 */
template<typename T, int Size> bool is_negative(const mat_base<T, Size, Size> &mat)
{
  return determinant(mat) < T(0);
}
template<> bool is_negative(const mat_base<float, 4, 4> &mat);
template<> bool is_negative(const mat_base<double, 4, 4> &mat);

/**
 * Returns true if matrices are equal within the given limit.
 */
template<typename T, int NumCol, int NumRow>
inline bool compare(const mat_base<T, NumCol, NumRow> &a,
                    const mat_base<T, NumCol, NumRow> &b,
                    const T limit)
{
  for (int i = 0; i < NumCol; i++) {
    for (int j = 0; j < NumRow; j++) {
      if (std::abs(a[i][j] - b[i][j]) > limit) {
        return false;
      }
    }
  }
  return true;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name 3x3 init helpers.
 * \{ */

namespace mat3x3 {

template<typename T> mat_base<T, 3, 3> from_rotation(const rotation::EulerXYZ<T> &rotation)
{
  double ci = std::cos(rotation[0]);
  double cj = std::cos(rotation[1]);
  double ch = std::cos(rotation[2]);
  double si = std::sin(rotation[0]);
  double sj = std::sin(rotation[1]);
  double sh = std::sin(rotation[2]);
  double cc = ci * ch;
  double cs = ci * sh;
  double sc = si * ch;
  double ss = si * sh;

  mat_base<T, 3, 3> mat;
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

template<typename T> mat_base<T, 3, 3> from_rotation(const rotation::Quaternion<T> &rotation)
{
  double q0 = M_SQRT2 * double(rotation[0]);
  double q1 = M_SQRT2 * double(rotation[1]);
  double q2 = M_SQRT2 * double(rotation[2]);
  double q3 = M_SQRT2 * double(rotation[3]);

  double qda = q0 * q1;
  double qdb = q0 * q2;
  double qdc = q0 * q3;
  double qaa = q1 * q1;
  double qab = q1 * q2;
  double qac = q1 * q3;
  double qbb = q2 * q2;
  double qbc = q2 * q3;
  double qcc = q3 * q3;

  mat_base<T, 3, 3> mat;
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

template<typename T> mat_base<T, 3, 3> from_scale(const vec_base<T, 3> &scale)
{
  return mat_base<T, 3, 3>::from_diagonal(scale);
}

template<typename T>
mat_base<T, 3, 3> from_rot_scale(const rotation::EulerXYZ<T> &rotation,
                                 const vec_base<T, 3> &scale)
{
  return mat3x3::from_rotation(rotation) * mat3x3::from_scale(scale);
}

template<typename T>
mat_base<T, 3, 3> from_normalized_axis_data(const vec_base<T, 3> forward, const vec_base<T, 3> up)
{
  BLI_ASSERT_UNIT_V3(forward);
  BLI_ASSERT_UNIT_V3(up);

  mat_base<T, 3, 3> matrix;
  matrix.forward() = forward;
  /* Negate the cross product so that the resulting matrix has determinant 1 (instead of -1).
   * Without the negation, the result would be a so called improper rotation. That means it
   * contains a reflection. Such an improper rotation matrix could not be converted to another
   * representation of a rotation such as euler angles. */
  matrix.right() = -math::cross(forward, up);
  matrix.up() = up;
  return matrix;
}

}  // namespace mat3x3

/** \} */

/* -------------------------------------------------------------------- */
/** \name 4x4 init helpers.
 * \{ */

namespace mat4x4 {

template<typename T> mat_base<T, 4, 4> from_location(const vec_base<T, 3> location)
{
  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>::identity();
  mat.location() = location;
  return mat;
}

template<typename T> mat_base<T, 4, 4> from_rotation(const rotation::EulerXYZ<T> rotation)
{
  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>(mat3x3::from_rotation(rotation));
  return mat;
}

template<typename T> mat_base<T, 4, 4> from_scale(const vec_base<T, 3> scale)
{
  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>(mat3x3::from_scale(scale));
  return mat;
}

template<typename T>
mat_base<T, 4, 4> from_loc_rot(const vec_base<T, 3> location, const rotation::EulerXYZ<T> rotation)
{
  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>(mat3x3::from_rotation(rotation));
  mat.location() = location;
  return mat;
}

template<typename T>
mat_base<T, 4, 4> from_loc_rot_scale(const vec_base<T, 3> location,
                                     const rotation::EulerXYZ<T> rotation,
                                     const vec_base<T, 3> scale)
{
  mat_base<T, 4, 4> mat = mat_base<T, 4, 4>(mat3x3::from_rot_scale(rotation, scale));
  mat.location() = location;
  return mat;
}

template<typename T>
mat_base<T, 4, 4> from_normalized_axis_data(const vec_base<T, 3> location,
                                            const vec_base<T, 3> forward,
                                            const vec_base<T, 3> up)
{
  mat_base<T, 4, 4> matrix = mat_base<T, 4, 4>(mat3x3::from_normalized_axis_data(forward, up));
  matrix.location() = location;
  return matrix;
}

}  // namespace mat4x4

/** \} */

/* -------------------------------------------------------------------- */
/** \name Conversion function.
 * \{ */

/* Implementation details. */
namespace mat3x3 {

template<typename T>
void normalized_to_eul2(const mat_base<T, 3, 3> &mat,
                        rotation::EulerXYZ<T> &eul1,
                        rotation::EulerXYZ<T> &eul2);

template<typename T>
rotation::Quaternion<T> normalized_to_quat_with_checks(const mat_base<T, 3, 3> &mat);

}  // namespace mat3x3

template<typename T, bool Normalized = false>
inline rotation::EulerXYZ<T> to_euler(const mat_base<T, 3, 3> &mat)
{
  rotation::EulerXYZ<T> eul1, eul2;
  if constexpr (Normalized) {
    mat3x3::normalized_to_eul2(mat, eul1, eul2);
  }
  else {
    mat3x3::normalized_to_eul2(normalize(mat), eul1, eul2);
  }
  /* Return best, which is just the one with lowest values it in. */
  return (length_manhattan(vec_base<T, 3>(eul1)) > length_manhattan(vec_base<T, 3>(eul2))) ? eul2 :
                                                                                             eul1;
}

template<typename T, bool Normalized = false>
inline rotation::EulerXYZ<T> to_euler(const mat_base<T, 4, 4> &mat)
{
  /* TODO(fclem): Avoid the copy with 3x3 ref. */
  return to_euler<T, Normalized>(mat_base<T, 3, 3>(mat));
}

template<typename T, bool Normalized = false>
inline rotation::Quaternion<T> to_quaternion(const mat_base<T, 3, 3> &mat)
{
  using namespace math;
  if constexpr (Normalized) {
    return mat3x3::normalized_to_quat_with_checks(mat);
  }
  else {
    return mat3x3::normalized_to_quat_with_checks(normalize(mat));
  }
}

template<typename T, bool Normalized = false>
inline rotation::Quaternion<T> to_quaternion(const mat_base<T, 4, 4> &mat)
{
  /* TODO(fclem): Avoid the copy with 3x3 ref. */
  return to_quaternion<T, Normalized>(mat_base<T, 3, 3>(mat));
}

template<typename T, int NumCol, int NumRow>
inline vec_base<T, 3> to_scale(const mat_base<T, NumCol, NumRow> &mat)
{
  return {length(mat.forward()), length(mat.right()), length(mat.up())};
}

/** \} */

}  // namespace blender::math
