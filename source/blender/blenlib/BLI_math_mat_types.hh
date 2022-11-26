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
#include "BLI_math_vector.h" /* For len_squared_v3 in BLI_ASSERT_UNIT_V3. */
#include "BLI_math_vector.hh"
#include "BLI_simd.h"
#include "BLI_utildefines.h"

namespace blender {

template<typename T, int NumCol, int NumRow>
struct mat_base : public vec_struct_base<vec_base<T, NumRow>, NumCol> {

  using base_type = T;
  using col_type = vec_base<T, NumRow>;
  using row_type = vec_base<T, NumCol>;
  static constexpr int min_dim = (NumRow < NumCol) ? NumRow : NumCol;
  static constexpr int col_len = NumCol;
  static constexpr int row_len = NumRow;

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

  explicit mat_base(const T *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = reinterpret_cast<const col_type *>(ptr)[i]; });
  }

  template<typename U, BLI_ENABLE_IF((std::is_convertible_v<U, T>))>
  explicit mat_base(const U *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = ptr[i]; });
  }

  explicit mat_base(const T (*ptr)[NumCol]) : mat_base(static_cast<const T *>(ptr[0]))
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
  /** \note defined outside of class. */
  // friend mat_base operator*(const mat_base &a, const mat_base &b);

  /** IMPORTANT: This is per component multiplication. */
  friend mat_base operator*(const mat_base &a, T b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b; });
    return result;
  }

  /** IMPORTANT: This is per component multiplication. */
  friend mat_base operator*(T a, const mat_base &b)
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
  mat_base &operator*=(T b)
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

  static mat_base zero()
  {
    return mat_base(0);
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

template<typename T, int NumCol, int NumRow>
mat_base<T, NumCol, NumRow> operator*(const mat_base<T, NumCol, NumRow> &a,
                                      const mat_base<T, NumCol, NumRow> &b);

#ifdef BLI_HAVE_SSE2
template<>
mat_base<float, 3, 3> operator*(const mat_base<float, 3, 3> &a, const mat_base<float, 3, 3> &b);
template<>
mat_base<float, 4, 4> operator*(const mat_base<float, 4, 4> &a, const mat_base<float, 4, 4> &b);
#endif

namespace math {

template<typename T, int Size> mat_base<T, Size, Size> inverse(const mat_base<T, Size, Size> &mat);

/**
 * Matrix inversion can be implemented more efficiently for affine matrices.
 */
template<typename T> inline mat_base<T, 4, 4> inverse_affine(const mat_base<T, 4, 4> &mat)
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

template<typename T, int Size> bool is_negative(const mat_base<T, Size, Size> &mat)
{
  return determinant(mat) < T(0);
}
template<> bool is_negative(const mat_base<float, 4, 4> &mat);
template<> bool is_negative(const mat_base<double, 4, 4> &mat);

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

template<typename T> inline T normalize(const T &a)
{
  T result;
  unroll<T::col_len>([&](auto i) { result[i] = math::normalize(a[i]); });
  return result;
}

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

/* TODO(fclem): Move to vector math. */
template<typename T, int Size>
inline bool compare(const vec_base<T, Size> &a, const vec_base<T, Size> &b, const T limit)
{
  for (int i = 0; i < Size; i++) {
    if (std::abs(a[i] - b[i]) > limit) {
      return false;
    }
  }
  return true;
}

template<typename T> struct AssertUnitEpsilon {
  static constexpr T value = T(BLI_ASSERT_UNIT_EPSILON_DB);
};

/* TODO(fclem): Move to vector math. */
template<typename T, int Size> inline bool is_unit_scale(const vec_base<T, Size> &v)
{
  /**
   * \note Checks are flipped so NAN doesn't assert.
   * This is done because we're making sure the value was normalized and in the case we
   * don't want NAN to be raising asserts since there is nothing to be done in that case.
   */
  const T test_unit = math::length_squared(v);
  return (!(std::abs(test_unit - T(1)) >= AssertUnitEpsilon<T>::value) ||
          !(std::abs(test_unit) >= AssertUnitEpsilon<T>::value));
}

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

}  // namespace math

namespace rotation {

template<typename T> struct EulerXYZ : public vec_base<T, 3> {
  /** TODO(fclem): Maybe add rotation order enum here? or in typename? */

  /** Inherit constructors. */
  using vec_base<T, 3>::vec_base;

  friend std::ostream &operator<<(std::ostream &stream, const EulerXYZ &mat)
  {
    return stream << "EulerXYZ" << static_cast<vec_base<T, 3>>(mat);
  }
};

template<typename T> struct Quaternion : public vec_base<T, 4> {
  /** Inherit constructors. */
  using vec_base<T, 4>::vec_base;

  static Quaternion<T> identity()
  {
    return {1, 0, 0, 0};
  }

  friend std::ostream &operator<<(std::ostream &stream, const Quaternion &mat)
  {
    return stream << "Quaternion" << static_cast<vec_base<T, 4>>(mat);
  }
};

}  // namespace rotation

namespace math {

template<typename U> struct AssertUnitEpsilon<rotation::Quaternion<U>> {
  static constexpr U value = AssertUnitEpsilon<U>::value * 10;
};

template<typename T> inline vec_base<T, 2> interp_dot_slerp(const T t, const T cosom)
{
  const T eps = 1e-4f;

  BLI_assert(IN_RANGE_INCL(cosom, T(-1.0001), T(1.0001)));

  vec_base<T, 2> w;
  /* within [-1..1] range, avoid aligned axis */
  if (LIKELY(std::abs(cosom) < (T(1.0) - eps))) {
    float omega, sinom;

    omega = std::acos(cosom);
    sinom = std::sin(omega);

    w[0] = std::sin((T(1.0) - t) * omega) / sinom;
    w[1] = std::sin(t * omega) / sinom;
  }
  else {
    /* fallback to lerp */
    w[0] = 1.0f - t;
    w[1] = t;
  }
  return w;
}

template<typename T>
inline rotation::Quaternion<T> interpolate(const rotation::Quaternion<T> &a,
                                           const rotation::Quaternion<T> &b,
                                           T t)
{
  BLI_assert(is_unit_scale(a));
  BLI_assert(is_unit_scale(b));

  vec_base<T, 4> quat = a;
  T cosom = dot(a, b);
  /* Rotate around shortest angle. */
  if (cosom < T(0)) {
    cosom = -cosom;
    quat = -quat;
  }

  vec_base<T, 2> w = interp_dot_slerp(t, cosom);

  return rotation::Quaternion<T>(w[0] * quat + w[1] * b);
}

}  // namespace math

template<typename T> struct mat_3x3 : public mat_base<T, 3, 3> {
  using vec3_type = vec_base<T, 3>;

  /** Inherit constructors. */
  using mat_base<T, 3, 3>::mat_base;

  mat_3x3(const mat_base<T, 3, 3> &base) : mat_base<T, 3, 3>(base)
  {
  }

  /** Init Helpers. */

  static mat_3x3 from_rotation(const rotation::EulerXYZ<T> rotation)
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

  static mat_3x3 from_rotation(const rotation::Quaternion<T> rotation)
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

    mat_3x3 mat;
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

  static mat_3x3 from_scale(const vec3_type scale)
  {
    return mat_3x3::from_diagonal(scale);
  }

  static mat_3x3 from_rot_scale(const rotation::EulerXYZ<T> rotation, const vec3_type scale)
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
    matrix.right() = -math::cross(forward, up);
    matrix.up() = up;
    return matrix;
  }

  /** Access helpers. Using Blender coordinate system. */

  vec3_type &forward()
  {
    return (*this)[0];
  }

  vec3_type &right()
  {
    return (*this)[1];
  }

  vec3_type &up()
  {
    return (*this)[2];
  }

  const vec3_type &forward() const
  {
    return (*this)[0];
  }

  const vec3_type &right() const
  {
    return (*this)[1];
  }

  const vec3_type &up() const
  {
    return (*this)[2];
  }

  /** Methods. */

  template<bool Normalized = false> rotation::EulerXYZ<T> to_euler() const
  {
    using namespace math;
    rotation::EulerXYZ<T> eul1, eul2;
    if constexpr (Normalized) {
      this->normalized_to_eul2(eul1, eul2);
    }
    else {
      normalize(*this).normalized_to_eul2(eul1, eul2);
    }
    /* Return best, which is just the one with lowest values it in. */
    return (length_manhattan(vec3_type(eul1)) > length_manhattan(vec3_type(eul2))) ? eul2 : eul1;
  }

  template<bool Normalized = false> rotation::Quaternion<T> to_quaternion() const
  {
    using namespace math;
    if constexpr (Normalized) {
      return this->normalized_to_quat_with_checks();
    }
    else {
      return normalize(*this).normalized_to_quat_with_checks();
    }
  }

  vec3_type to_scale() const
  {
    return {math::length(forward()), math::length(right()), math::length(up())};
  }

  friend std::ostream &operator<<(std::ostream &stream, const mat_3x3 &mat)
  {
    return stream << static_cast<mat_base<T, 3, 3>>(mat);
  }

 private:
  void normalized_to_eul2(rotation::EulerXYZ<T> &eul1, rotation::EulerXYZ<T> &eul2) const
  {
    BLI_assert(math::is_unit_scale(*this));

    const float cy = std::hypot((*this)[0][0], (*this)[0][1]);
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

  rotation::Quaternion<T> normalized_to_quat_fast() const
  {
    BLI_assert(math::is_unit_scale(*this));
    /* Caller must ensure matrices aren't negative for valid results, see: T24291, T94231. */
    BLI_assert(!math::is_negative(*this));

    const mat_3x3 &mat = *this;
    rotation::Quaternion<T> q;

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
        q[1] = 0.25f * s;
        s = 1.0f / s;
        q[0] = (mat[1][2] - mat[2][1]) * s;
        q[2] = (mat[0][1] + mat[1][0]) * s;
        q[3] = (mat[2][0] + mat[0][2]) * s;
        if (UNLIKELY((trace == 1.0f) && (q[0] == 0.0f && q[2] == 0.0f && q[3] == 0.0f))) {
          /* Avoids the need to normalize the degenerate case. */
          q[1] = 1.0f;
        }
      }
      else {
        const T trace = 1.0f - mat[0][0] + mat[1][1] - mat[2][2];
        T s = 2.0f * std::sqrt(trace);
        if (mat[2][0] < mat[0][2]) {
          /* Ensure W is non-negative for a canonical result. */
          s = -s;
        }
        q[2] = 0.25f * s;
        s = 1.0f / s;
        q[0] = (mat[2][0] - mat[0][2]) * s;
        q[1] = (mat[0][1] + mat[1][0]) * s;
        q[3] = (mat[1][2] + mat[2][1]) * s;
        if (UNLIKELY((trace == 1.0f) && (q[0] == 0.0f && q[1] == 0.0f && q[3] == 0.0f))) {
          /* Avoids the need to normalize the degenerate case. */
          q[2] = 1.0f;
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
        q[3] = 0.25f * s;
        s = 1.0f / s;
        q[0] = (mat[0][1] - mat[1][0]) * s;
        q[1] = (mat[2][0] + mat[0][2]) * s;
        q[2] = (mat[1][2] + mat[2][1]) * s;
        if (UNLIKELY((trace == 1.0f) && (q[0] == 0.0f && q[1] == 0.0f && q[2] == 0.0f))) {
          /* Avoids the need to normalize the degenerate case. */
          q[3] = 1.0f;
        }
      }
      else {
        /* NOTE(@campbellbarton): A zero matrix will fall through to this block,
         * needed so a zero scaled matrices to return a quaternion without rotation, see: T101848.
         */
        const T trace = 1.0f + mat[0][0] + mat[1][1] + mat[2][2];
        T s = 2.0f * std::sqrt(trace);
        q[0] = 0.25f * s;
        s = 1.0f / s;
        q[1] = (mat[1][2] - mat[2][1]) * s;
        q[2] = (mat[2][0] - mat[0][2]) * s;
        q[3] = (mat[0][1] - mat[1][0]) * s;
        if (UNLIKELY((trace == 1.0f) && (q[1] == 0.0f && q[2] == 0.0f && q[3] == 0.0f))) {
          /* Avoids the need to normalize the degenerate case. */
          q[0] = 1.0f;
        }
      }
    }

    BLI_assert(!(q[0] < 0.0f));
    BLI_assert(math::is_unit_scale(q));
    return q;
  }

  rotation::Quaternion<T> normalized_to_quat_with_checks() const
  {
    const T det = math::determinant(*this);
    if (UNLIKELY(!isfinite(det))) {
      return rotation::Quaternion<T>::identity();
    }
    else if (UNLIKELY(det < T(0))) {
      return mat_3x3(-(*this)).normalized_to_quat_fast();
    }
    return this->normalized_to_quat_fast();
  }
};

template<typename T> struct mat_4x4 : public mat_base<T, 4, 4> {
  using mat_3x3 = mat_3x3<T>;
  using vec3_type = vec_base<T, 3>;

  /** Inherit constructors. */
  using mat_base<T, 4, 4>::mat_base;

  mat_4x4(const mat_base<T, 4, 4> &base) : mat_base<T, 4, 4>(base)
  {
  }

  /** Init Helpers. */

  static mat_4x4 from_location(const vec3_type location)
  {
    mat_4x4 mat = mat_4x4::identity();
    mat.location() = location;
    return mat;
  }

  static mat_4x4 from_rotation(const rotation::EulerXYZ<T> rotation)
  {
    mat_4x4 mat = mat_4x4(mat_3x3::from_rotation(rotation));
    return mat;
  }

  static mat_4x4 from_scale(const vec3_type scale)
  {
    mat_4x4 mat = mat_4x4(mat_3x3::from_scale(scale));
    return mat;
  }

  static mat_4x4 from_loc_rot(const vec3_type location, const rotation::EulerXYZ<T> rotation)
  {
    mat_4x4 mat = mat_4x4(mat_3x3::from_rotation(rotation));
    mat.location() = location;
    return mat;
  }

  static mat_4x4 from_loc_rot_scale(const vec3_type location,
                                    const rotation::EulerXYZ<T> rotation,
                                    const vec3_type scale)
  {
    mat_4x4 mat = mat_4x4(mat_3x3::from_rot_scale(rotation, scale));
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
    return *reinterpret_cast<vec3_type *>(&(*this)[0]);
  }

  vec3_type &right()
  {
    return *reinterpret_cast<vec3_type *>(&(*this)[1]);
  }

  vec3_type &up()
  {
    return *reinterpret_cast<vec3_type *>(&(*this)[2]);
  }

  vec3_type &location()
  {
    return *reinterpret_cast<vec3_type *>(&(*this)[3]);
  }

  const vec3_type &forward() const
  {
    return *reinterpret_cast<const vec3_type *>(&(*this)[0]);
  }

  const vec3_type &right() const
  {
    return *reinterpret_cast<const vec3_type *>(&(*this)[1]);
  }

  const vec3_type &up() const
  {
    return *reinterpret_cast<const vec3_type *>(&(*this)[2]);
  }

  const vec3_type &location() const
  {
    return *reinterpret_cast<const vec3_type *>(&(*this)[3]);
  }

  /** Methods. */

  template<bool Normalized = false> rotation::EulerXYZ<T> to_euler() const
  {
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return mat_3x3(*this).template to_euler<Normalized>();
  }

  template<bool Normalized = false> rotation::Quaternion<T> to_quaternion() const
  {
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return mat_3x3(*this).template to_quaternion<Normalized>();
  }

  vec3_type to_scale() const
  {
    /* TODO(fclem): Avoid the copy with 3x3 ref. */
    return mat_3x3(*this).to_scale();
  }

  friend std::ostream &operator<<(std::ostream &stream, const mat_4x4 &mat)
  {
    return stream << static_cast<mat_base<T, 4, 4>>(mat);
  }
};

}  // namespace blender
