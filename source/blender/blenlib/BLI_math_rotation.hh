/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bli
 */

#include "BLI_math_vec_types.hh"
#include "BLI_math_vector.hh"

namespace blender::math::detail {
/**
 * Rotation Types
 *
 * These are wrappers around the vector types.
 * It gives more semantic informations allowing for overloaded functions depending on the rotation
 * type. It also prevent implicit cast from rotation to vector types.
 */

/* Forward declaration. */
template<typename T> struct AxisAngle;
template<typename T> struct Quaternion;

template<typename T> struct EulerXYZ {
  T x, y, z;

  EulerXYZ() = default;

  EulerXYZ(const T &x, const T &y, const T &z)
  {
    this->x = x;
    this->y = y;
    this->z = z;
  }

  EulerXYZ(const vec_base<T, 3> &vec) : EulerXYZ(UNPACK3(vec)){};

  static EulerXYZ identity()
  {
    return {0, 0, 0};
  }

  explicit operator vec_base<T, 3>() const
  {
    return {this->x, this->y, this->z};
  }

  explicit operator AxisAngle<T>() const;

  explicit operator Quaternion<T>() const;

  friend std::ostream &operator<<(std::ostream &stream, const EulerXYZ &rot)
  {
    return stream << "EulerXYZ" << static_cast<vec_base<T, 3>>(rot);
  }
};

template<typename T> struct Quaternion {
  T x, y, z, w;

  Quaternion() = default;

  Quaternion(const T &x, const T &y, const T &z, const T &w)
  {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
  }

  Quaternion(const vec_base<T, 4> &vec) : Quaternion(UNPACK4(vec)){};

  static Quaternion identity()
  {
    return {1, 0, 0, 0};
  }

  explicit operator vec_base<T, 4>() const
  {
    return {this->x, this->y, this->z, this->w};
  }

  explicit operator EulerXYZ<T>() const;

  explicit operator AxisAngle<T>() const;

  friend std::ostream &operator<<(std::ostream &stream, const Quaternion &rot)
  {
    return stream << "Quaternion" << static_cast<vec_base<T, 4>>(rot);
  }
};

template<typename T> struct AxisAngle {
  vec_base<T, 3> axis;
  T angle;

  explicit AxisAngle(){};

  AxisAngle(const vec_base<T, 3> &axis, T angle)
  {
    T length;
    const vec_base<T, 3> normalized_axis = math::normalize_and_get_length(axis, length);
    if (length > 0.0f) {
      this->axis = normalized_axis;
      this->angle = angle;
    }
    else {
      *this = identity();
    }
  }

  static AxisAngle<T> identity()
  {
    return {{0, 1, 0}, 0};
  }

  explicit operator Quaternion<T>() const;

  explicit operator EulerXYZ<T>() const;

  friend bool operator==(const AxisAngle &a, const AxisAngle &b)
  {
    return (a.axis == b.axis) && (a.angle == b.angle);
  }

  friend bool operator!=(const AxisAngle &a, const AxisAngle &b)
  {
    return (a != b);
  }

  friend std::ostream &operator<<(std::ostream &stream, const AxisAngle &rot)
  {
    return stream << "AxisAngle(axis=" << rot.axis << ", angle=" << rot.angle << ")";
  }
};

/**
 * A version of AxisAngle that expects axis to be already normalized.
 * Implicitly cast back to AxisAngle.
 */
template<typename T> struct AxisAngleNormalized : public AxisAngle<T> {
  AxisAngleNormalized(const vec_base<T, 3> &axis, T angle) : AxisAngle<T>()
  {
    BLI_assert(is_unit_scale(axis));
    this->axis = axis;
    this->angle = angle;
  }

  operator AxisAngle<T>() const
  {
    return *this;
  }
};

/**
 * Intermediate Types.
 *
 * Some functions need to have higher precision than standard floats.
 */
template<typename T> struct TypeTraits {
  using IntermediateType = T;
};
template<> struct TypeTraits<float> {
  using IntermediateType = double;
};

};  // namespace blender::math::detail

namespace blender {
/* Most common used types. */
using EulerXYZ = math::detail::EulerXYZ<float>;
using Quaternion = math::detail::Quaternion<float>;
using AxisAngle = math::detail::AxisAngle<float>;
using AxisAngleNormalized = math::detail::AxisAngleNormalized<float>;
};  // namespace blender

namespace blender::math {

template<typename U> struct AssertUnitEpsilon<detail::Quaternion<U>> {
  static constexpr U value = AssertUnitEpsilon<U>::value * 10;
};

/**
 * Generic function for implementing slerp
 * (quaternions and spherical vector coords).
 *
 * \param t: factor in [0..1]
 * \param cosom: dot product from normalized vectors/quats.
 * \param r_w: calculated weights.
 */
template<typename T> inline vec_base<T, 2> interpolate_dot_slerp(const T t, const T cosom)
{
  const T eps = T(1e-4);

  BLI_assert(IN_RANGE_INCL(cosom, T(-1.0001), T(1.0001)));

  vec_base<T, 2> w;
  /* Within [-1..1] range, avoid aligned axis. */
  if (LIKELY(math::abs(cosom) < (T(1) - eps))) {
    const T omega = math::acos(cosom);
    const T sinom = math::sin(omega);

    w[0] = math::sin((T(1) - t) * omega) / sinom;
    w[1] = math::sin(t * omega) / sinom;
  }
  else {
    /* Fallback to lerp */
    w[0] = T(1) - t;
    w[1] = t;
  }
  return w;
}

template<typename T>
inline detail::Quaternion<T> interpolate(const detail::Quaternion<T> &a,
                                         const detail::Quaternion<T> &b,
                                         T t)
{
  using Vec4T = vec_base<T, 4>;
  BLI_assert(is_unit_scale(Vec4T(a)));
  BLI_assert(is_unit_scale(Vec4T(b)));

  Vec4T quat = Vec4T(a);
  T cosom = dot(Vec4T(a), Vec4T(b));
  /* Rotate around shortest angle. */
  if (cosom < T(0)) {
    cosom = -cosom;
    quat = -quat;
  }

  vec_base<T, 2> w = interpolate_dot_slerp(t, cosom);

  return detail::Quaternion<T>(w[0] * quat + w[1] * Vec4T(b));
}

/**
 * Rotate the unit-length \a direction around the unit-length \a axis by the \a angle.
 */
float3 rotate_direction_around_axis(const float3 &direction, const float3 &axis, float angle);

/**
 * Rotate any arbitrary \a vector around the \a center position, with a unit-length \a axis
 * and the specified \a angle.
 */
float3 rotate_around_axis(const float3 &vector,
                          const float3 &center,
                          const float3 &axis,
                          float angle);

}  // namespace blender::math
