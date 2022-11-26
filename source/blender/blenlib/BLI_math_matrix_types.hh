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
  using vec3_type = vec_base<T, 3>;
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

  /** Access helpers. Using Blender coordinate system. */

  vec3_type &forward()
  {
    BLI_STATIC_ASSERT(NumCol >= 1, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[0]);
  }

  vec3_type &right()
  {
    BLI_STATIC_ASSERT(NumCol >= 2, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[1]);
  }

  vec3_type &up()
  {
    BLI_STATIC_ASSERT(NumCol >= 3, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[2]);
  }

  vec3_type &location()
  {
    BLI_STATIC_ASSERT(NumCol >= 4, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[3]);
  }

  const vec3_type &forward() const
  {
    BLI_STATIC_ASSERT(NumCol >= 1, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[0]);
  }

  const vec3_type &right() const
  {
    BLI_STATIC_ASSERT(NumCol >= 2, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[1]);
  }

  const vec3_type &up() const
  {
    BLI_STATIC_ASSERT(NumCol >= 3, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[2]);
  }

  const vec3_type &location() const
  {
    BLI_STATIC_ASSERT(NumCol >= 4, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[3]);
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

using float2x2 = blender::mat_base<float, 2, 2>;
using float3x3 = blender::mat_base<float, 3, 3>;
using float4x4 = blender::mat_base<float, 4, 4>;
using double2x2 = blender::mat_base<double, 2, 2>;
using double3x3 = blender::mat_base<double, 3, 3>;
using double4x4 = blender::mat_base<double, 4, 4>;

}  // namespace blender
