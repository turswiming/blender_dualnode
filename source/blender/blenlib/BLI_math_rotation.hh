/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

/** \file
 * \ingroup bli
 */

#include "BLI_math_vec_types.hh"

namespace blender::rotation {
/**
 * Rotation Types
 *
 * These are wrappers around the vector types.
 * It gives more semantic informations allowing for overloaded functions depending on the rotation
 * type. It also prevent implicit cast from rotation to vector types.
 */

template<typename T> struct EulerXYZ : public vec_base<T, 3> {
  using vec_base<T, 3>::vec_base;

  static EulerXYZ<T> identity()
  {
    return {0, 0, 0};
  }

  friend std::ostream &operator<<(std::ostream &stream, const EulerXYZ &mat)
  {
    return stream << "EulerXYZ" << static_cast<vec_base<T, 3>>(mat);
  }
};

template<typename T> struct Quaternion : public vec_base<T, 4> {
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

};  // namespace blender::rotation

namespace blender {
/* Most common used types. */
using EulerXYZ = rotation::EulerXYZ<float>;
using Quaternion = rotation::Quaternion<float>;
};  // namespace blender

namespace blender::math {

template<typename U> struct AssertUnitEpsilon<rotation::Quaternion<U>> {
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

  vec_base<T, 2> w = interpolate_dot_slerp(t, cosom);

  return rotation::Quaternion<T>(w[0] * quat + w[1] * b);
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
