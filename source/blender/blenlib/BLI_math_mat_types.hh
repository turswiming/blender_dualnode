/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup bli
 *
 * Template for matrix types.
 * This only works for tightly packed `T` without alignment padding.
 */

#include <array>
#include <cmath>
#include <iostream>
#include <type_traits>

#include "BLI_math_vec_types.hh"
#include "BLI_math_vector.hh"
#include "BLI_utildefines.h"

namespace blender {

#if 0
template<typename T> struct mat_view {
  T *data;
  IndexRange col_range;
  IndexRange row_range;

  mat_view(T *src, IndexRange column_range, IndexRange row_range)
  {
    data = src;
    col_range = column_range;
    row_range = row_range;
  }

  friend operator=()
  {
  }
};
#endif

template<typename T, int NumCol, int NumRow>
struct mat_base : public vec_struct_base<vec_base<T, NumRow>, NumCol> {

  using base_type = T;
  using col_type = vec_base<T, NumRow>;
  using row_type = vec_base<T, NumCol>;
  static constexpr int min_dim = (NumRow < NumCol) ? NumRow : NumCol;

  mat_base() = default;

  /** Initialize the diagonal of the matrix to this value and the rest with zero. Matches GLSL. */
  explicit mat_base(T value)
  {
    unroll<NumCol>([&](auto i) {
      (*this)[i] = col_type(0);
      (*this)[i][i] = value;
    });
  }

  template<typename U, BLI_ENABLE_IF((std::is_convertible_v<U, T>))>
  explicit mat_base(U value) : mat_base(T(value))
  {
  }

/* Workaround issue with template BLI_ENABLE_IF((Size == 2)) not working. */
#define BLI_ENABLE_IF_MAT(_size, _test) int S = _size, BLI_ENABLE_IF((S _test))

  template<BLI_ENABLE_IF_MAT(NumCol, == 2)> mat_base(col_type _x, col_type _y)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
  }

  template<BLI_ENABLE_IF_MAT(NumCol, == 3)> mat_base(col_type _x, col_type _y, col_type _z)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
    (*this)[2] = _z;
  }

  template<BLI_ENABLE_IF_MAT(NumCol, == 4)>
  mat_base(col_type _x, col_type _y, col_type _z, col_type _w)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
    (*this)[2] = _z;
    (*this)[3] = _w;
  }

  /** Masking. */

  template<typename U, int OtherNumCol, int OtherNumRow>
  explicit mat_base(const mat_base<U, OtherNumCol, OtherNumRow> &other)
  {
    if constexpr ((OtherNumRow >= NumRow) && (OtherNumCol >= NumCol)) {
      unroll<NumCol>([&](auto i) { (*this)[i] = col_type(other[i]); });
    }
    else {
      /* Allow enlarging following GLSL standard (i.e: mat4x4(mat3x3())). */
      unroll<NumCol>([&](auto i) {
        unroll<NumRow>([&](auto j) {
          if (i < OtherNumCol && j < OtherNumRow) {
            (*this)[i][j] = other[i][j];
          }
          else if (i == j) {
            (*this)[i][j] = T(1);
          }
          else {
            (*this)[i][j] = T(0);
          }
        });
      });
    }
  }

#undef BLI_ENABLE_IF_MAT

  /** Conversion from pointers (from C-style vectors). */

  mat_base(const T *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = reinterpret_cast<const col_type *>(ptr)[i]; });
  }

  template<typename U, BLI_ENABLE_IF((std::is_convertible_v<U, T>))>
  explicit mat_base(const U *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = ptr[i]; });
  }

  mat_base(const T (*ptr)[NumCol]) : mat_base(static_cast<const T *>(ptr[0]))
  {
  }

  /** Conversion from other matrix types. */

  template<typename U> explicit mat_base(const mat_base<U, NumRow, NumCol> &vec)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = col_type(vec[i]); });
  }

  /** C-style pointer dereference. */
  /** TODO(fclem): Ideally, theses implicit conversion should be removed once all the functions are
   * ported to the C++ API. */

  operator const T *() const
  {
    return reinterpret_cast<const T *>(this);
  }

  operator T *()
  {
    return reinterpret_cast<T *>(this);
  }

  using c_style_mat = T[NumCol][NumRow];

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  const c_style_mat &ptr() const
  {
    return *reinterpret_cast<c_style_mat *>(this);
  }

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  c_style_mat &ptr()
  {
    return *reinterpret_cast<c_style_mat *>(this);
  }

  /** Array access. */

  const col_type &operator[](int index) const
  {
    BLI_assert(index >= 0);
    BLI_assert(index < NumCol);
    return reinterpret_cast<const col_type *>(this)[index];
  }

  col_type &operator[](int index)
  {
    BLI_assert(index >= 0);
    BLI_assert(index < NumCol);
    return reinterpret_cast<col_type *>(this)[index];
  }

  /** Matrix operators. */

  friend mat_base operator+(const mat_base &a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b[i]; });
    return result;
  }

  friend mat_base operator+(const mat_base &a, T b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b; });
    return result;
  }

  friend mat_base operator+(T a, const mat_base &b)
  {
    return b + a;
  }

  mat_base &operator+=(const mat_base &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b[i]; });
    return *this;
  }

  mat_base &operator+=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b; });
    return *this;
  }

  friend mat_base operator-(const mat_base &a)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = -a[i]; });
    return result;
  }

  friend mat_base operator-(const mat_base &a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b[i]; });
    return result;
  }

  friend mat_base operator-(const mat_base &a, T b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b; });
    return result;
  }

  friend mat_base operator-(T a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a - b[i]; });
    return result;
  }

  mat_base &operator-=(const mat_base &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b[i]; });
    return *this;
  }

  mat_base &operator-=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b; });
    return *this;
  }

  /** IMPORTANT: This is matrix multiplication. Not per component. */
  friend mat_base operator*(const mat_base &a, const mat_base &b)
  {
    /** \note this is the reference implementation.
     * Subclass are free to overload it with vectorized / optimized code. */
    /** \note Only tested for square matrices. Might still contain bugs. */
    mat_base result = mat_base(0);
    unroll<NumCol>([&](auto c) {
      unroll<NumRow>([&](auto r) {
        /** \note this is vector multiplication. */
        result[c] += b[c][r] * a[r];
      });
    });
    return result;
  }

  /** IMPORTANT: This is per component multiplication. */
  template<typename FactorT> friend mat_base operator*(const mat_base &a, FactorT b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b; });
    return result;
  }

  /** IMPORTANT: This is per component multiplication. */
  template<typename FactorT> friend mat_base operator*(FactorT a, const mat_base &b)
  {
    return b * a;
  }

  /** IMPORTANT: This is matrix multiplication. Not per component. */
  mat_base &operator*=(const mat_base &b)
  {
    const mat_base &a = *this;
    *this = a * b;
    return *this;
  }

  /** IMPORTANT: This is per component multiplication. */
  template<typename FactorT> mat_base &operator*=(FactorT b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] *= b; });
    return *this;
  }

  /** Vector operators. */

  friend col_type operator*(const mat_base &a, const row_type &b)
  {
    /** \note this is the reference implementation.
     * Subclass are free to overload it with vectorized / optimized code. */
    col_type result(0);
    unroll<NumCol>([&](auto c) { result += b[c] * a[c]; });
    return result;
  }

  /** Compare. */

  friend bool operator==(const mat_base &a, const mat_base &b)
  {
    for (int i = 0; i < NumCol; i++) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const mat_base &a, const mat_base &b)
  {
    return !(a == b);
  }

  /** Misc */

  static mat_base identity()
  {
    return mat_base(1);
  }

  /** Initialize the diagonal of the matrix with the given vector. */
  static mat_base from_diagonal(vec_base<T, min_dim> vec)
  {
    mat_base result(0);
    unroll<min_dim>([&](auto i) { result[i][i] = vec[i]; });
    return result;
  }

  uint64_t hash() const
  {
    uint64_t h = 435109;
    unroll<NumCol * NumRow>([&](auto i) {
      T value = (static_cast<const T *>(this))[i];
      h = h * 33 + *reinterpret_cast<const as_uint_type<T> *>(&value);
    });
    return h;
  }

  friend std::ostream &operator<<(std::ostream &stream, const mat_base &mat)
  {
    stream << "(\n";
    for (int i = 0; i < NumCol; i++) {
      stream << "(";
      for (int j = 0; j < NumRow; j++) {
        /** NOTE: j and i are swapped to follow mathematical convention. */
        stream << mat[j][i];
        if (j < NumRow - 1) {
          stream << ", ";
        }
      }
      stream << ")";
      if (i < NumCol - 1) {
        stream << ",";
      }
      stream << "\n";
    }
    stream << ")\n";
    return stream;
  }
};

namespace math {

template<typename T, int Size> mat_base<T, Size, Size> inverse(mat_base<T, Size, Size> &mat);
template<> mat_base<float, 2, 2> inverse(mat_base<float, 2, 2> &mat);
template<> mat_base<float, 3, 3> inverse(mat_base<float, 3, 3> &mat);
template<> mat_base<float, 4, 4> inverse(mat_base<float, 4, 4> &mat);
template<> mat_base<double, 2, 2> inverse(mat_base<double, 2, 2> &mat);
template<> mat_base<double, 3, 3> inverse(mat_base<double, 3, 3> &mat);
template<> mat_base<double, 4, 4> inverse(mat_base<double, 4, 4> &mat);

/**
 * Matrix inversion can be implemented more efficiently for affine matrices.
 */
template<typename T> inline mat_base<T, 4, 4> inverse_affine(mat_base<T, 4, 4> &mat)
{
  BLI_assert(mat[0][3] == 0.0f && mat[1][3] == 0.0f && mat[2][3] == 0.0f && mat[3][3] == 1.0f);
  return inverse(mat);
}

template<typename T, int NumCol, int NumRow>
inline mat_base<T, NumCol, NumRow> transpose(mat_base<T, NumRow, NumCol> &mat)
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
template<> float determinant(const mat_base<float, 3, 3> &fmat);
template<> float determinant(const mat_base<float, 4, 4> &fmat);
template<> double determinant(const mat_base<double, 3, 3> &dmat);
template<> double determinant(const mat_base<double, 4, 4> &dmat);

template<typename T>
inline mat_base<T, 3, 3> interpolate(const mat_base<T, 3, 3> &a, const mat_base<T, 3, 3> &b, T t)
{
#if 0
  /* 'Rotation' component ('U' part of polar decomposition,
   * the closest orthogonal matrix to M3 rot/scale
   * transformation matrix), spherically interpolated. */
  float U_A[3][3], U_B[3][3], U[3][3];
  float quat_A[4], quat_B[4], quat[4];
  /* 'Scaling' component ('P' part of polar decomposition, i.e. scaling in U-defined space),
   * linearly interpolated. */
  float P_A[3][3], P_B[3][3], P[3][3];

  int i;

  mat3_polar_decompose(A, U_A, P_A);
  mat3_polar_decompose(B, U_B, P_B);

  /* Quaternions cannot represent an axis flip. If such a singularity is detected, choose a
   * different decomposition of the matrix that still satisfies A = U_A * P_A but which has a
   * positive determinant and thus no axis flips. This resolves T77154.
   *
   * Note that a flip of two axes is just a rotation of 180 degrees around the third axis, and
   * three flipped axes are just an 180 degree rotation + a single axis flip. It is thus sufficient
   * to solve this problem for single axis flips. */
  if (is_negative_m3(U_A)) {
    mul_m3_fl(U_A, -1.0f);
    mul_m3_fl(P_A, -1.0f);
  }
  if (is_negative_m3(U_B)) {
    mul_m3_fl(U_B, -1.0f);
    mul_m3_fl(P_B, -1.0f);
  }

  mat3_to_quat(quat_A, U_A);
  mat3_to_quat(quat_B, U_B);
  interp_qt_qtqt(quat, quat_A, quat_B, t);
  quat_to_mat3(U, quat);

  for (i = 0; i < 3; i++) {
    interp_v3_v3v3(P[i], P_A[i], P_B[i], t);
  }

  /* And we reconstruct rot/scale matrix from interpolated polar components */
  mul_m3_m3m3(R, U, P);
#endif
  /* TODO */
  return mat_base<T, 3, 3>(0);
}

/**
 * Complete transform matrix interpolation,
 * based on polar-decomposition-based interpolation from #interp_m3_m3m3.
 *
 * \param A: Input matrix which is totally effective with `t = 0.0`.
 * \param B: Input matrix which is totally effective with `t = 1.0`.
 * \param t: Interpolation factor.
 */
template<typename T>
inline mat_base<T, 4, 4> interpolate(const mat_base<T, 4, 4> &a, const mat_base<T, 4, 4> &b, T t)
{
  using mat_4x4 = mat_base<T, 4, 4>;
  using mat_3x3 = mat_base<T, 3, 3>;
  using vec_3 = vec_base<T, 3>;

  mat_4x4 result = mat_4x4(interpolate(mat_3x3(a), mat_3x3(b), t));

  /* Location component, linearly interpolated. */
  vec_3 location = *reinterpret_cast<vec_3 *>(&result[3][0]);
  vec_3 loc_a = *reinterpret_cast<vec_3 *>(&a[3][0]);
  vec_3 loc_b = *reinterpret_cast<vec_3 *>(&b[3][0]);
  location = mix(loc_a, loc_b, t);

  return result;
}

template<typename T, int NumCol, int NumRow>
inline mat_base<T, NumCol, NumRow> normalize(const mat_base<T, NumCol, NumRow> &a)
{
  mat_base<T, NumCol, NumRow> result;
  unroll<NumCol>([&](auto i) { result[i] = math::normalize(a[i]); });
  return result;
}

}  // namespace math

namespace experiment {

template<typename T> struct RotationEuler : public vec_base<T, 3> {
  /** TODO(fclem): Maybe add rotation order enum here? or in type? */
};

template<typename T> struct mat_3x3 : public mat_base<T, 3, 3> {
  using vec3_type = vec_base<T, 3>;

  /** Init Helpers. */

  static mat_3x3 from_rotation(const RotationEuler<T> rotation)
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

    mat_3x3 mat;
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

  static mat_3x3 from_scale(const vec3_type scale)
  {
    return mat_3x3::from_diagonal(scale);
  }

  static mat_3x3 from_scale(T scale)
  {
    return mat_3x3(scale);
  }

  static mat_3x3 from_rot_scale(const RotationEuler<T> rotation, const vec3_type scale)
  {
    return mat_3x3::from_rotation(rotation) * mat_3x3::from_scale(scale);
  }

  static mat_3x3 from_normalized_axis_data(const vec3_type forward, const vec3_type up)
  {
    BLI_ASSERT_UNIT_V3(forward);
    BLI_ASSERT_UNIT_V3(up);

    mat_3x3 matrix;
    matrix.forward() = forward;
    /* Negate the cross product so that the resulting matrix has determinant 1 (instead of -1).
     * Without the negation, the result would be a so called improper rotation. That means it
     * contains a reflection. Such an improper rotation matrix could not be converted to another
     * representation of a rotation such as euler angles. */
    matrix.left() = -math::cross(forward, up);
    matrix.up() = up;
    return matrix;
  }

  /** Access helpers. Using Blender coordinate system. */

  vec3_type &forward()
  {
    return (*this)[0];
  }

  vec3_type &left()
  {
    return (*this)[1];
  }

  vec3_type &up()
  {
    return (*this)[2];
  }

  /** Methods. */

  /* Assumes XYZ rotation order. */
  vec3_type to_euler() const
  {
    vec3_type eul1, eul2;
    normalize(*this).normalized_to_eul2(eul1, eul2);
    /* Return best, which is just the one with lowest values it in. */
    return (length_manhattan(eul1) > length_manhattan(eul2)) ? eul2 : eul1;
  }

  vec3_type to_scale() const
  {
    return {length(forward()), length(left()), length(up())};
  }

  bool is_negative() const
  {
    return determinant(*this) < 0.0f;
  }

 private:
  void normalized_to_eul2(RotationEuler<T> &eul1, RotationEuler<T> &eul2)
  {
    BLI_ASSERT_UNIT_M3(this->ptr());

    const float cy = hypotf((*this)[0][0], (*this)[0][1]);
    if (cy > 16.0f * FLT_EPSILON) {
      eul1[0] = std::atan2((*this)[1][2], (*this)[2][2]);
      eul1[1] = std::atan2(-(*this)[0][2], cy);
      eul1[2] = std::atan2((*this)[0][1], (*this)[0][0]);

      eul2[0] = std::atan2(-(*this)[1][2], -(*this)[2][2]);
      eul2[1] = std::atan2(-(*this)[0][2], -cy);
      eul2[2] = std::atan2(-(*this)[0][1], -(*this)[0][0]);
    }
    else {
      eul1[0] = std::atan2(-(*this)[2][1], (*this)[1][1]);
      eul1[1] = std::atan2(-(*this)[0][2], cy);
      eul1[2] = 0.0f;

      eul2 = eul1;
    }
  }
};

template<typename T> struct mat_4x4 : public mat_base<T, 4, 4> {
  using mat_3x3 = mat_base<T, 3, 3>;
  using vec3_type = vec_base<T, 3>;

  /** Init Helpers. */

  static mat_4x4 from_loc_rot_scale(const vec3_type location,
                                    const RotationEuler<T> rotation,
                                    const vec3_type scale)
  {
    mat_4x4 mat = mat_4x4(0);
    mat.location() = location;
    mat.as_3x3() = mat_3x3::from_rot_scale();
    return mat;
  }

  static mat_4x4 from_location(const vec3_type location)
  {
    mat_4x4 mat(1.0f);
    mat.location() = location;
    return mat;
  }

  static mat_4x4 from_normalized_axis_data(const vec3_type location,
                                           const vec3_type forward,
                                           const vec3_type up)
  {
    mat_4x4 matrix = mat_4x4(mat_3x3::from_normalized_axis_data(forward, up));
    matrix.location() = location;
    return matrix;
  }

  /** Access helpers. Using Blender coordinate system. */

  vec3_type &forward()
  {
    return *reinterpret_cast<vec3_type *>((*this)[0]);
  }

  vec3_type &left()
  {
    return *reinterpret_cast<vec3_type *>((*this)[1]);
  }

  vec3_type &up()
  {
    return *reinterpret_cast<vec3_type *>((*this)[2]);
  }

  vec3_type &location()
  {
    return *reinterpret_cast<vec3_type *>((*this)[3]);
  }

  /** Methods. */

  /* Assumes XYZ rotation order. */
  vec3_type to_euler() const
  {
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return mat_3x3(*this).to_euler();
  }

  vec3_type to_scale() const
  {
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return mat_3x3(*this).to_scale();
  }

  bool is_negative() const
  {
    /* Don't use determinant(float4x4) as only the 3x3 components are needed
     * when the matrix is used as a transformation to represent location/scale/rotation. */
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return determinant(mat_3x3(*this)) < 0.0f;
  }
};

}  // namespace experiment
}  // namespace blender
