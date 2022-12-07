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
#include "BLI_utildefines.h"
#include "BLI_utility_mixins.hh"

namespace blender {

template<typename T,
         int NumCol,
         int NumRow,
         int SrcNumCol,
         int SrcNumRow,
         int SrcStartCol,
         int SrcStartRow>
struct MatView;

template<typename T, int NumCol, int NumRow>
struct alignas(4 * sizeof(T)) MatBase : public vec_struct_base<vec_base<T, NumRow>, NumCol> {

  using base_type = T;
  using vec3_type = vec_base<T, 3>;
  using col_type = vec_base<T, NumRow>;
  using row_type = vec_base<T, NumCol>;
  static constexpr int min_dim = (NumRow < NumCol) ? NumRow : NumCol;
  static constexpr int col_len = NumCol;
  static constexpr int row_len = NumRow;

  MatBase() = default;

  /** Initialize the diagonal of the matrix to this value and the rest with zero. Matches GLSL. */
  explicit MatBase(T value)
  {
    unroll<NumCol>([&](auto i) {
      (*this)[i] = col_type(0);
      (*this)[i][i] = value;
    });
  }

  template<typename U, BLI_ENABLE_IF((std::is_convertible_v<U, T>))>
  explicit MatBase(U value) : MatBase(T(value))
  {
  }

/* Workaround issue with template BLI_ENABLE_IF((Size == 2)) not working. */
#define BLI_ENABLE_IF_MAT(_size, _test) int S = _size, BLI_ENABLE_IF((S _test))

  template<BLI_ENABLE_IF_MAT(NumCol, == 2)> MatBase(col_type _x, col_type _y)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
  }

  template<BLI_ENABLE_IF_MAT(NumCol, == 3)> MatBase(col_type _x, col_type _y, col_type _z)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
    (*this)[2] = _z;
  }

  template<BLI_ENABLE_IF_MAT(NumCol, == 4)>
  MatBase(col_type _x, col_type _y, col_type _z, col_type _w)
  {
    (*this)[0] = _x;
    (*this)[1] = _y;
    (*this)[2] = _z;
    (*this)[3] = _w;
  }

  /** Masking. */

  template<typename U, int OtherNumCol, int OtherNumRow>
  explicit MatBase(const MatBase<U, OtherNumCol, OtherNumRow> &other)
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

  explicit MatBase(const T *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = reinterpret_cast<const col_type *>(ptr)[i]; });
  }

  template<typename U, BLI_ENABLE_IF((std::is_convertible_v<U, T>))> explicit MatBase(const U *ptr)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = ptr[i]; });
  }

  explicit MatBase(const T (*ptr)[NumCol]) : MatBase(static_cast<const T *>(ptr[0]))
  {
  }

  /** Conversion from other matrix types. */

  template<typename U> explicit MatBase(const MatBase<U, NumRow, NumCol> &vec)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = col_type(vec[i]); });
  }

  /** C-style pointer dereference. */

  using c_style_mat = T[NumCol][NumRow];

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  const c_style_mat &ptr() const
  {
    return *reinterpret_cast<const c_style_mat *>(this);
  }

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  c_style_mat &ptr()
  {
    return *reinterpret_cast<c_style_mat *>(this);
  }

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  const T *base_ptr() const
  {
    return reinterpret_cast<const T *>(this);
  }

  /** \note Prevent implicit cast to types that could fit other pointer constructor. */
  T *base_ptr()
  {
    return reinterpret_cast<T *>(this);
  }

  /** View creation. */

  template<int ViewNumCol = NumCol,
           int ViewNumRow = NumRow,
           int SrcStartCol = 0,
           int SrcStartRow = 0>
  const MatView<T, ViewNumCol, ViewNumRow, NumCol, NumRow, SrcStartCol, SrcStartRow> view() const
  {
    return MatView<T, ViewNumCol, ViewNumRow, NumCol, NumRow, SrcStartCol, SrcStartRow>(
        const_cast<MatBase &>(*this));
  }

  template<int ViewNumCol = NumCol,
           int ViewNumRow = NumRow,
           int SrcStartCol = 0,
           int SrcStartRow = 0>
  MatView<T, ViewNumCol, ViewNumRow, NumCol, NumRow, SrcStartCol, SrcStartRow> view()
  {
    return MatView<T, ViewNumCol, ViewNumRow, NumCol, NumRow, SrcStartCol, SrcStartRow>(*this);
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

  vec3_type &x_axis()
  {
    BLI_STATIC_ASSERT(NumCol >= 1, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[0]);
  }

  vec3_type &y_axis()
  {
    BLI_STATIC_ASSERT(NumCol >= 2, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<vec3_type *>(&(*this)[1]);
  }

  vec3_type &z_axis()
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

  const vec3_type &x_axis() const
  {
    BLI_STATIC_ASSERT(NumCol >= 1, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[0]);
  }

  const vec3_type &y_axis() const
  {
    BLI_STATIC_ASSERT(NumCol >= 2, "Wrong Matrix dimension");
    BLI_STATIC_ASSERT(NumRow >= 3, "Wrong Matrix dimension");
    return *reinterpret_cast<const vec3_type *>(&(*this)[1]);
  }

  const vec3_type &z_axis() const
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

  friend MatBase operator+(const MatBase &a, const MatBase &b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b[i]; });
    return result;
  }

  friend MatBase operator+(const MatBase &a, T b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b; });
    return result;
  }

  friend MatBase operator+(T a, const MatBase &b)
  {
    return b + a;
  }

  MatBase &operator+=(const MatBase &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b[i]; });
    return *this;
  }

  MatBase &operator+=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b; });
    return *this;
  }

  friend MatBase operator-(const MatBase &a)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = -a[i]; });
    return result;
  }

  friend MatBase operator-(const MatBase &a, const MatBase &b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b[i]; });
    return result;
  }

  friend MatBase operator-(const MatBase &a, T b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b; });
    return result;
  }

  friend MatBase operator-(T a, const MatBase &b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a - b[i]; });
    return result;
  }

  MatBase &operator-=(const MatBase &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b[i]; });
    return *this;
  }

  MatBase &operator-=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b; });
    return *this;
  }

  /** Multiply two matrices using matrix multiplication. */
  MatBase operator*(const MatBase &b) const
  {
    const MatBase &a = *this;
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    /** TODO(fclem): Only tested for square matrices. Might still contain bugs. */
    MatBase<T, NumCol, NumRow> result(0);
    unroll<NumCol>([&](auto c) {
      unroll<NumRow>([&](auto r) {
        /* This is vector multiplication. */
        result[c] += b[c][r] * a[r];
      });
    });
    return result;
  }

  /** Multiply each component by a scalar. */
  friend MatBase operator*(const MatBase &a, T b)
  {
    MatBase result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b; });
    return result;
  }

  /** Multiply each component by a scalar. */
  friend MatBase operator*(T a, const MatBase &b)
  {
    return b * a;
  }

  /** Multiply two matrices using matrix multiplication. */
  MatBase &operator*=(const MatBase &b)
  {
    const MatBase &a = *this;
    *this = a * b;
    return *this;
  }

  /** Multiply each component by a scalar. */
  MatBase &operator*=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] *= b; });
    return *this;
  }

  /** Vector operators. */

  friend col_type operator*(const MatBase &a, const row_type &b)
  {
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    col_type result(0);
    unroll<NumCol>([&](auto c) { result += b[c] * a[c]; });
    return result;
  }

  /** Multiply by the transposed. */
  friend row_type operator*(const col_type &a, const MatBase &b)
  {
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    row_type result(0);
    unroll<NumCol>([&](auto c) { unroll<NumRow>([&](auto r) { result[c] += b[c][r] * a[r]; }); });
    return result;
  }

  /** Compare. */

  friend bool operator==(const MatBase &a, const MatBase &b)
  {
    for (int i = 0; i < NumCol; i++) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const MatBase &a, const MatBase &b)
  {
    return !(a == b);
  }

  /** Miscellaneous. */

  static MatBase identity()
  {
    return MatBase(1);
  }

  static MatBase zero()
  {
    return MatBase(0);
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

  friend std::ostream &operator<<(std::ostream &stream, const MatBase &mat)
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

template<typename T,
         /** The view dimensions. */
         int NumCol,
         int NumRow,
         /** The source matrix dimensions. */
         int SrcNumCol,
         int SrcNumRow,
         /** The base offset inside the source matrix. */
         int SrcStartCol,
         int SrcStartRow>
struct MatView : NonCopyable, NonMovable {
  using MatT = MatBase<T, NumCol, NumRow>;
  using SrcMatT = MatBase<T, SrcNumCol, SrcNumRow>;
  using col_type = vec_base<T, NumRow>;
  using row_type = vec_base<T, NumCol>;

 private:
  SrcMatT &src_;

 public:
  MatView() = delete;

  MatView(SrcMatT &src) : src_(src)
  {
    BLI_STATIC_ASSERT(SrcStartCol < SrcNumCol, "View does not fit source matrix dimensions");
    BLI_STATIC_ASSERT(SrcStartRow < SrcNumRow, "View does not fit source matrix dimensions");
  }

  /** Array access. */

  const col_type &operator[](int index) const
  {
    BLI_assert(index >= 0);
    BLI_assert(index < NumCol);
    return *reinterpret_cast<const col_type *>(&src_[index + SrcStartCol][SrcStartRow]);
  }

  col_type &operator[](int index)
  {
    BLI_assert(index >= 0);
    BLI_assert(index < NumCol);
    return *reinterpret_cast<col_type *>(&src_[index + SrcStartCol][SrcStartRow]);
  }

  /** Conversion back to matrix. */

  operator MatT() const
  {
    MatT mat;
    unroll<NumCol>([&](auto c) { mat[c] = (*this)[c]; });
    return mat;
  }

  /** Copy Assignment. */

  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  MatView &operator=(const MatView<T,
                                   NumCol,
                                   NumRow,
                                   OtherSrcNumCol,
                                   OtherSrcNumRow,
                                   OtherSrcStartCol,
                                   OtherSrcStartRow> &other)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] = other[i]; });
    return *this;
  }

  MatView &operator=(const MatT &other)
  {
    *this = other.template view();
    return *this;
  }

  /** Matrix operators. */

  friend MatT operator+(const MatView &a, T b)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] + b; });
    return result;
  }

  friend MatT operator+(T a, const MatView &b)
  {
    return b + a;
  }

  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  MatView &operator+=(const MatView<T,
                                    NumCol,
                                    NumRow,
                                    OtherSrcNumCol,
                                    OtherSrcNumRow,
                                    OtherSrcStartCol,
                                    OtherSrcStartRow> &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b[i]; });
    return *this;
  }

  MatView &operator+=(const MatT &b)
  {
    return *this += b.template view();
  }

  MatView &operator+=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] += b; });
    return *this;
  }

  friend MatT operator-(const MatView &a)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = -a[i]; });
    return result;
  }

  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  friend MatT operator-(const MatView &a,
                        const MatView<T,
                                      NumCol,
                                      NumRow,
                                      OtherSrcNumCol,
                                      OtherSrcNumRow,
                                      OtherSrcStartCol,
                                      OtherSrcStartRow> &b)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b[i]; });
    return result;
  }

  friend MatT operator-(const MatView &a, const MatT &b)
  {
    return a - b.template view();
  }

  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  friend MatT operator-(const MatView<T,
                                      NumCol,
                                      NumRow,
                                      OtherSrcNumCol,
                                      OtherSrcNumRow,
                                      OtherSrcStartCol,
                                      OtherSrcStartRow> &a,
                        const MatView &b)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b[i]; });
    return result;
  }

  friend MatT operator-(const MatT &a, const MatView &b)
  {
    return a.template view() - b;
  }

  friend MatT operator-(const MatView &a, T b)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] - b; });
    return result;
  }

  friend MatView operator-(T a, const MatView &b)
  {
    MatView result;
    unroll<NumCol>([&](auto i) { result[i] = a - b[i]; });
    return result;
  }

  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  MatView &operator-=(const MatView<T,
                                    NumCol,
                                    NumRow,
                                    OtherSrcNumCol,
                                    OtherSrcNumRow,
                                    OtherSrcStartCol,
                                    OtherSrcStartRow> &b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b[i]; });
    return *this;
  }

  MatView &operator-=(const MatT &b)
  {
    return *this -= b.template view();
  }

  MatView &operator-=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] -= b; });
    return *this;
  }

  /** Multiply two matrices using matrix multiplication. */
  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  MatT operator*(const MatView<T,
                               NumCol,
                               NumRow,
                               OtherSrcNumCol,
                               OtherSrcNumRow,
                               OtherSrcStartCol,
                               OtherSrcStartRow> &b) const
  {
    const MatView &a = *this;
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    /** TODO(fclem): Only tested for square matrices. Might still contain bugs. */
    MatT result(0);
    unroll<NumCol>([&](auto c) {
      unroll<NumRow>([&](auto r) {
        /* This is vector multiplication. */
        result[c] += b[c][r] * a[r];
      });
    });
    return result;
  }

  MatT operator*(const MatT &b) const
  {
    return *this * b.template view();
  }

  /** Multiply each component by a scalar. */
  friend MatT operator*(const MatView &a, T b)
  {
    MatT result;
    unroll<NumCol>([&](auto i) { result[i] = a[i] * b; });
    return result;
  }

  /** Multiply each component by a scalar. */
  friend MatT operator*(T a, const MatView &b)
  {
    return b * a;
  }

  /** Multiply two matrices using matrix multiplication. */
  template<int OtherSrcNumCol, int OtherSrcNumRow, int OtherSrcStartCol, int OtherSrcStartRow>
  MatView &operator*=(const MatView<T,
                                    NumCol,
                                    NumRow,
                                    OtherSrcNumCol,
                                    OtherSrcNumRow,
                                    OtherSrcStartCol,
                                    OtherSrcStartRow> &b)
  {
    const MatView &a = *this;
    *this = a * b;
    return *this;
  }

  MatView &operator*=(const MatT &b)
  {
    return *this *= b.template view();
  }

  /** Multiply each component by a scalar. */
  MatView &operator*=(T b)
  {
    unroll<NumCol>([&](auto i) { (*this)[i] *= b; });
    return *this;
  }

  /** Vector operators. */

  friend col_type operator*(const MatView &a, const row_type &b)
  {
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    col_type result(0);
    unroll<NumCol>([&](auto c) { result += b[c] * a[c]; });
    return result;
  }

  /** Multiply by the transposed. */
  friend row_type operator*(const col_type &a, const MatView &b)
  {
    /* This is the reference implementation.
     * Might be overloaded with vectorized / optimized code. */
    row_type result(0);
    unroll<NumCol>([&](auto c) { unroll<NumRow>([&](auto r) { result[c] += b[c][r] * a[r]; }); });
    return result;
  }

  /** Compare. */

  friend bool operator==(const MatView &a, const MatView &b)
  {
    for (int i = 0; i < NumCol; i++) {
      if (a[i] != b[i]) {
        return false;
      }
    }
    return true;
  }

  friend bool operator!=(const MatView &a, const MatView &b)
  {
    return !(a == b);
  }

  /** Miscellaneous. */

  friend std::ostream &operator<<(std::ostream &stream, const MatView &mat)
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

using float2x2 = MatBase<float, 2, 2>;
using float2x3 = MatBase<float, 2, 3>;
using float2x4 = MatBase<float, 2, 4>;
using float3x2 = MatBase<float, 3, 2>;
using float3x3 = MatBase<float, 3, 3>;
using float3x4 = MatBase<float, 3, 4>;
using float4x2 = MatBase<float, 4, 2>;
using float4x3 = MatBase<float, 4, 3>;
using float4x4 = MatBase<float, 4, 4>;

using double2x2 = MatBase<double, 2, 2>;
using double2x3 = MatBase<double, 2, 3>;
using double2x4 = MatBase<double, 2, 4>;
using double3x2 = MatBase<double, 3, 2>;
using double3x3 = MatBase<double, 3, 3>;
using double3x4 = MatBase<double, 3, 4>;
using double4x2 = MatBase<double, 4, 2>;
using double4x3 = MatBase<double, 4, 3>;
using double4x4 = MatBase<double, 4, 4>;

/* Specialization for SSE optimization. */
template<> float4x4 float4x4::operator*(const float4x4 &b) const;

extern template float2x2 float2x2::operator*(const float2x2 &b) const;
extern template float3x3 float3x3::operator*(const float3x3 &b) const;
extern template double2x2 double2x2::operator*(const double2x2 &b) const;
extern template double3x3 double3x3::operator*(const double3x3 &b) const;
extern template double4x4 double4x4::operator*(const double4x4 &b) const;

}  // namespace blender
