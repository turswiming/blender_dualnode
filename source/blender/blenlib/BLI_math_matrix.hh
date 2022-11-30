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
[[nodiscard]] inline bool is_unit_scale(const mat_base<T, NumCol, NumRow> &m)
{
  for (int i = 0; i < NumCol; i++) {
    if (!is_unit_scale(m[i])) {
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
template<typename T, int Size>
[[nodiscard]] mat_base<T, Size, Size> invert(const mat_base<T, Size, Size> &mat);

/**
 * Flip the matrix around its diagonal. Also flips dimensions for non square matrices.
 */
template<typename T, int NumCol, int NumRow>
[[nodiscard]] inline mat_base<T, NumCol, NumRow> transpose(const mat_base<T, NumRow, NumCol> &mat)
{
  mat_base<T, NumCol, NumRow> result;
  for (int i = 0; i < NumCol; i++) {
    for (int j = 0; j < NumRow; j++) {
      result[i][j] = mat[j][i];
    }
  }
  return result;
}

template<typename T, int Size> [[nodiscard]] T determinant(const mat_base<T, Size, Size> &mat);

/** Interpolate each component linearly. */
template<typename T, int NumCol, int NumRow>
[[nodiscard]] inline mat_base<T, NumCol, NumRow> interpolate_linear(
    const mat_base<T, NumCol, NumRow> &a, const mat_base<T, NumCol, NumRow> &b, T t)
{
  mat_base<T, NumCol, NumRow> result;
  unroll<NumCol>([&](auto c) { result[c] = interpolate(a[c], b[c], t); });
  return result;
}

/**
 * A polar-decomposition-based interpolation between matrix A and matrix B.
 *
 * \note This code is about five times slower than the 'naive' interpolation
 * (it typically remains below 2 usec on an average i74700,
 * while naive implementation remains below 0.4 usec).
 * However, it gives expected results even with non-uniformly scaled matrices,
 * see T46418 for an example.
 *
 * Based on "Matrix Animation and Polar Decomposition", by Ken Shoemake & Tom Duff
 *
 * \param A: Input matrix which is totally effective with `t = 0.0`.
 * \param B: Input matrix which is totally effective with `t = 1.0`.
 * \param t: Interpolation factor.
 */
template<typename T>
[[nodiscard]] mat_base<T, 3, 3> interpolate(const mat_base<T, 3, 3> &A,
                                            const mat_base<T, 3, 3> &B,
                                            T t);

/**
 * Complete transform matrix interpolation,
 * based on polar-decomposition-based interpolation from #interpolate<T, 3, 3>.
 *
 * \param A: Input matrix which is totally effective with `t = 0.0`.
 * \param B: Input matrix which is totally effective with `t = 1.0`.
 * \param t: Interpolation factor.
 */
template<typename T>
[[nodiscard]] mat_base<T, 4, 4> interpolate(const mat_base<T, 4, 4> &A,
                                            const mat_base<T, 4, 4> &B,
                                            T t);

/**
 * Normalize individually each column of the matrix.
 */
template<typename T> [[nodiscard]] inline T normalize(const T &a)
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
 *
 * \note It doesn't use determinant(mat4x4) as only the 3x3 components are needed
 * when the matrix is used as a transformation to represent location/scale/rotation.
 */
template<typename T, int Size> [[nodiscard]] bool is_negative(const mat_base<T, Size, Size> &mat)
{
  return determinant(mat) < T(0);
}
template<typename T> [[nodiscard]] bool is_negative(const mat_base<T, 4, 4> &mat);

/**
 * Returns true if matrices are equal within the given epsilon.
 */
template<typename T, int NumCol, int NumRow>
[[nodiscard]] inline bool compare(const mat_base<T, NumCol, NumRow> &a,
                                  const mat_base<T, NumCol, NumRow> &b,
                                  const T epsilon)
{
  for (int i = 0; i < NumCol; i++) {
    for (int j = 0; j < NumRow; j++) {
      if (math::abs(a[i][j] - b[i][j]) > epsilon) {
        return false;
      }
    }
  }
  return true;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Init helpers.
 * \{ */

/* Implementation details. */
namespace detail {

template<typename T, int NumCol, int NumRow>
[[nodiscard]] mat_base<T, NumCol, NumRow> from_rotation(const EulerXYZ<T> &rotation);

template<typename T, int NumCol, int NumRow>
[[nodiscard]] mat_base<T, NumCol, NumRow> from_rotation(const Quaternion<T> &rotation);

template<typename T, int NumCol, int NumRow>
[[nodiscard]] mat_base<T, NumCol, NumRow> from_rotation(const AxisAngle<T> &rotation);

}  // namespace detail

template<typename MatT>
[[nodiscard]] MatT from_location(const vec_base<typename MatT::base_type, 3> &location)
{
  MatT mat = MatT::identity();
  mat.location() = location;
  return mat;
}

/**
 * Create a matrix whose diagonal is defined by the given scale vector.
 * If vector dimension is lower than matrix diagonal, the missing terms are filled with ones.
 */
template<typename MatT, int ScaleDim>
[[nodiscard]] MatT from_scale(const vec_base<typename MatT::base_type, ScaleDim> &scale)
{
  BLI_STATIC_ASSERT(ScaleDim <= MatT::min_dim,
                    "Scale dimension should fit the matrix diagonal length.");
  MatT result(0);
  unroll<MatT::min_dim>(
      [&](auto i) { result[i][i] = (i < ScaleDim) ? scale[i] : typename MatT::base_type(1); });
  return result;
}

template<typename MatT, typename RotationT>
[[nodiscard]] MatT from_rotation(const RotationT &rotation)
{
  return detail::from_rotation<typename MatT::base_type, MatT::col_len, MatT::row_len>(rotation);
}

template<typename MatT, typename RotationT, typename VectorT>
[[nodiscard]] MatT from_rot_scale(const RotationT &rotation, const VectorT &scale)
{
  return from_rotation<MatT>(rotation) * from_scale<MatT>(scale);
}

template<typename MatT, typename RotationT, int ScaleDim>
[[nodiscard]] MatT from_loc_rot_scale(const vec_base<typename MatT::base_type, 3> &location,
                                      const RotationT &rotation,
                                      const vec_base<typename MatT::base_type, ScaleDim> &scale)
{
  using Mat3x3 = mat_base<typename MatT::base_type, 3, 3>;
  MatT mat = MatT(from_rot_scale<Mat3x3>(rotation, scale));
  mat.location() = location;
  return mat;
}

template<typename MatT, typename RotationT>
[[nodiscard]] MatT from_loc_rot(const vec_base<typename MatT::base_type, 3> &location,
                                const RotationT &rotation)
{
  using Mat3x3 = mat_base<typename MatT::base_type, 3, 3>;
  MatT mat = MatT(from_rotation<Mat3x3>(rotation));
  mat.location() = location;
  return mat;
}

template<typename MatT, typename VectorT>
[[nodiscard]] MatT from_normalized_axis_data(const VectorT forward, const VectorT up)
{
  BLI_ASSERT_UNIT_V3(forward);
  BLI_ASSERT_UNIT_V3(up);

  MatT matrix;
  matrix.forward() = forward;
  /* Negate the cross product so that the resulting matrix has determinant 1 (instead of -1).
   * Without the negation, the result would be a so called improper rotation. That means it
   * contains a reflection. Such an improper rotation matrix could not be converted to another
   * representation of a rotation such as euler angles. */
  matrix.right() = -math::cross(forward, up);
  matrix.up() = up;
  return matrix;
}

template<typename MatT, typename VectorT>
[[nodiscard]] MatT from_normalized_axis_data(const VectorT location,
                                             const VectorT forward,
                                             const VectorT up)
{
  using Mat3x3 = mat_base<typename MatT::base_type, 3, 3>;
  MatT matrix = MatT(from_normalized_axis_data<Mat3x3>(forward, up));
  matrix.location() = location;
  return matrix;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Conversion function.
 * \{ */

/* Implementation details. */
namespace detail {

template<typename T>
void normalized_to_eul2(const mat_base<T, 3, 3> &mat, EulerXYZ<T> &eul1, EulerXYZ<T> &eul2);

template<typename T>
[[nodiscard]] Quaternion<T> normalized_to_quat_with_checks(const mat_base<T, 3, 3> &mat);

}  // namespace detail

template<typename T, bool Normalized = false>
[[nodiscard]] inline EulerXYZ<T> to_euler(const mat_base<T, 3, 3> &mat)
{
  EulerXYZ<T> eul1, eul2;
  if constexpr (Normalized) {
    detail::normalized_to_eul2(mat, eul1, eul2);
  }
  else {
    detail::normalized_to_eul2(normalize(mat), eul1, eul2);
  }
  /* Return best, which is just the one with lowest values it in. */
  return (length_manhattan(vec_base<T, 3>(eul1)) > length_manhattan(vec_base<T, 3>(eul2))) ? eul2 :
                                                                                             eul1;
}

template<typename T, bool Normalized = false>
[[nodiscard]] inline rotation::EulerXYZ<T> to_euler(const mat_base<T, 4, 4> &mat)
{
  /* TODO(fclem): Avoid the copy with 3x3 ref. */
  return to_euler<T, Normalized>(mat_base<T, 3, 3>(mat));
}

template<typename T, bool Normalized = false>
[[nodiscard]] inline rotation::Quaternion<T> to_quaternion(const mat_base<T, 3, 3> &mat)
{
  using namespace math;
  if constexpr (Normalized) {
    return detail::normalized_to_quat_with_checks(mat);
  }
  else {
    return detail::normalized_to_quat_with_checks(normalize(mat));
  }
}

template<typename T, bool Normalized = false>
[[nodiscard]] inline rotation::Quaternion<T> to_quaternion(const mat_base<T, 4, 4> &mat)
{
  /* TODO(fclem): Avoid the copy with 3x3 ref. */
  return to_quaternion<T, Normalized>(mat_base<T, 3, 3>(mat));
}

template<typename T, int NumCol, int NumRow>
[[nodiscard]] inline vec_base<T, 3> to_scale(const mat_base<T, NumCol, NumRow> &mat)
{
  return {length(mat.x_axis()), length(mat.y_axis()), length(mat.z_axis())};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Transform function.
 * \{ */

template<typename T>
[[nodiscard]] vec_base<T, 3> transform_point(const mat_base<T, 3, 3> &mat,
                                             const vec_base<T, 3> &point);

template<typename T>
[[nodiscard]] vec_base<T, 3> transform_point(const mat_base<T, 4, 4> &mat,
                                             const vec_base<T, 3> &point);

template<typename T>
[[nodiscard]] vec_base<T, 3> transform_direction(const mat_base<T, 3, 3> &mat,
                                                 const vec_base<T, 3> &direction);

template<typename T>
[[nodiscard]] vec_base<T, 3> transform_direction(const mat_base<T, 4, 4> &mat,
                                                 const vec_base<T, 3> &direction);

template<typename T>
[[nodiscard]] vec_base<T, 3> project_point(const mat_base<T, 4, 4> &mat,
                                           const vec_base<T, 3> &point);

/** \} */

/* -------------------------------------------------------------------- */
/** \name Projection Matrices.
 * \{ */

namespace projection {

/**
 * \brief Create an orthographic projection matrix using OpenGL coordinate convention:
 * Maps each axis range to [-1..1] range for all axes.
 * The resulting matrix can be used with either #project_point or #transform_point.
 */
template<typename T>
[[nodiscard]] mat_base<T, 4, 4> orthographic(
    T left, T right, T bottom, T top, T near_clip, T far_clip);

/**
 * \brief Create a perspective projection matrix using OpenGL coordinate convention:
 * Maps each axis range to [-1..1] range for all axes.
 * `left`, `right`, `bottom`, `top` are frustum side distances at `z=near_clip`.
 * The resulting matrix can be used with #project_point.
 */
template<typename T>
[[nodiscard]] mat_base<T, 4, 4> perspective(
    T left, T right, T bottom, T top, T near_clip, T far_clip);

/**
 * \brief Create a perspective projection matrix using OpenGL coordinate convention:
 * Maps each axis range to [-1..1] range for all axes.
 * Uses field of view angles instead of plane distances.
 * The resulting matrix can be used with #project_point.
 */
template<typename T>
[[nodiscard]] mat_base<T, 4, 4> perspective_fov(
    T angle_left, T angle_right, T angle_bottom, T angle_top, T near_clip, T far_clip)
{
  mat_base<T, 4, 4> mat = perspective(math::tan(angle_left),
                                      math::tan(angle_right),
                                      math::tan(angle_bottom),
                                      math::tan(angle_top),
                                      near_clip,
                                      far_clip);
  mat[0][0] /= near_clip;
  mat[1][1] /= near_clip;
  return mat;
}

}  // namespace projection

/** \} */

}  // namespace blender::math
