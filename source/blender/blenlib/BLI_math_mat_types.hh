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

template<typename T, int NumCol, int NumRow>
struct mat_base : public vec_struct_base<vec_base<T, NumRow>, NumCol> {

  using base_type = T;
  using col_type = vec_base<T, NumRow>;

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

  template<typename U,
           int OtherNumRow,
           int OtherNumCol,
           BLI_ENABLE_IF((OtherNumRow > NumRow) && (OtherNumCol > NumCol))>
  explicit mat_base(const mat_base<U, OtherNumRow, OtherNumCol> &other)
  {
    /* TODO(fclem): Allow enlarging following GLSL standard (i.e: mat4(mat3())). */
    unroll<NumCol>([&](auto i) { (*this)[i] = col_type(other[i]); });
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
    unroll<NumCol>([&](auto i) { (*this)[i] = T(vec[i]); });
  }

  /** C-style pointer dereference. */

  operator const T *() const
  {
    return reinterpret_cast<const T *>(this);
  }

  operator T *()
  {
    return reinterpret_cast<T *>(this);
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
#if 0

  /** Matrix operators. */

  friend mat_base operator+(const mat_base &a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b[i]; });
    return result;
  }

  mat_base &operator+=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b; });
  }

  mat_base &operator+=(mat_base b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b[i]; });
    return *this;
  }

  friend mat_base operator-(const mat_base &a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b[i]; });
    return result;
  }

  mat_base &operator-=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b; });
  }

  mat_base &operator-=(mat_base b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b[i]; });
    return *this;
  }

  /** IMPORTANT: This is matrix multiplication. Not per component. */
  friend mat_base operator*(const mat_base &a, const mat_base &b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b[i]; });
    return result;
  }

  /** IMPORTANT: This is matrix multiplication. Not per component. */
  mat_base &operator*=(mat_base b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] *= b[i]; });
    return *this;
  }

  /** IMPORTANT: This is per component multiplication. */
  mat_base &operator*=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] *= b; });
  }

  /** IMPORTANT: This is per component multiplication. */
  friend mat_base operator*(const mat_base &a, T b)
  {
    mat_base result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b; });
    return result;
  }

  /** Vector operators. */

  friend col_type operator*(const mat_base &a, const col_type &b)
  {
    /* TODO */
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

  template<int OtherNumRow, int OtherNumCol> struct mat_sub_ref {
    using mat_sub = mat_base<T, OtherNumRow, OtherNumCol>;

    mat_base &data;

    mat_sub_ref &operator+=(T b)
    {
      data.as_eigen().block(0, 0, OtherNumRow, OtherNumCol) += b;
    }
  };

  template<int OtherNumRow,
           int OtherNumCol,
           BLI_ENABLE_IF((OtherNumRow <= NumRow) && (OtherNumCol <= NumCol))>
  mat_sub_ref<OtherNumRow, OtherNumCol> &as()
  {
    return {*this};
  }

  /** Misc */

  static mat_base identity()
  {
    return mat_base(1);
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
    char fchar[16];
    stream << "(\n";
    for (int i = 0; i < NumCol; i++) {
      stream << "(";
      for (int j = 0; j < NumRow; j++) {
        /* Print with the same precision so we have aligned numbers. */
        /** NOTE: j and i are swapped to follow mathematical convention. */
        /** TODO(fclem): Use non float  */
        snprintf(fchar, sizeof(fchar), "%11.6f", mat[j][i]);
        stream << fchar;
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
#endif
};
#if 0
namespace experiment {

template<typename T> struct mat_4x4 : public mat_base<T, 4, 4> {
  using vec3_type = vec_base<T, 3>;
  /** Init Helpers. */

  /* Assumes an XYZ euler order. */
  static mat_4x4 from_loc_eul_scale(const vec3_type location,
                                    const vec3_type rotation,
                                    const vec3_type scale)
  {
    mat_4x4 mat;
    loc_eul_size_to_mat4(mat.values, location, rotation, scale);
    return mat;
  }

  static mat_4x4 from_location(const vec3_type location)
  {
    mat_4x4 mat(1.0f);
    mat.position() = location;
    return mat;
  }

  static mat_4x4 from_normalized_axis_data(const vec3_type location,
                                           const vec3_type forward,
                                           const vec3_type up)
  {
    BLI_ASSERT_UNIT_V3(forward);
    BLI_ASSERT_UNIT_V3(up);

    /* Negate the cross product so that the resulting matrix has determinant 1 (instead of -1).
     * Without the negation, the result would be a so called improper rotation. That means it
     * contains a reflection. Such an improper rotation matrix could not be converted to another
     * representation of a rotation such as euler angles. */
    const vec3_type cross = -math::cross(forward, up);

    mat_4x4 matrix;
    matrix.forward() = forward;
    matrix.left() = cross;
    matrix.up() = up;
    matrix.position() = location;

    matrix[0][3] = 0.0f;
    matrix[1][3] = 0.0f;
    matrix[2][3] = 0.0f;
    matrix[3][3] = 1.0f;

    return matrix;
  }

  /** C-style casting. */

  using c_style_float4x4 = float[4][4];
  c_style_float4x4 &ptr()
  {
    return *reinterpret_cast<c_style_float4x4 *>(this);
  }

  const c_style_float4x4 &ptr() const
  {
    return *reinterpret_cast<c_style_float4x4 *>(this);
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

  vec3_type &position()
  {
    return *reinterpret_cast<vec3_type *>((*this)[3]);
  }

  /** Methods. */

  /* Assumes XYZ rotation order. */
  vec3_type to_euler() const
  {
    vec3_type euler;
    mat4_to_eul(euler, values);
    return euler;
  }

  vec3_type to_scale() const
  {
    vec3_type scale;
    mat4_to_size(scale, values);
    return scale;
  }

  bool is_negative() const
  {
    /* Don't use determinant(float4x4) as only the 3x3 components are needed
     * when the matrix is used as a transformation to represent location/scale/rotation. */
    return determinant(this->as<3, 3>()) < 0.0f;
  }
};

namespace blender::math {

template<typename T, int Size> inline mat_base<T, Size, Size> inverse(mat_base<T, Size, Size> &mat)
{
  /* TODO: Eigen? */
}

/**
 * Matrix inversion can be implemented more efficiently for affine matrices.
 */
template<typename T> inline mat_base<T, 4, 4> inverse_affine(mat_base<T, 4, 4> &mat)
{
  BLI_assert(mat[0][3] == 0.0f && mat[1][3] == 0.0f && mat[2][3] == 0.0f && mat[3][3] == 1.0f);
  /* TODO: Eigen? */
}

template<typename T, int NumRow, int NumCol>
inline mat_base<T, NumCol, NumRow> transpose(mat_base<T, NumRow, NumCol> &mat)
{
  /* TODO: Eigen? */
}

template<typename T, int Size> inline T determinant(mat_base<T, Size, Size> &mat)
{
  /* TODO: Eigen? */
}

/* TODO(fclem): Should we keep interpolate as name? this isn't a real mix. */
template<typename T> inline mat_base<T, 4, 4> mix(mat_base<T, 4, 4> &a, mat_base<T, 4, 4> &b, T t)
{
  /* TODO: Port */
  // interp_m4_m4m4(result, a.values, b.values, t);
}

}  // namespace blender::math

}  // namespace experiment
#endif
}  // namespace blender
