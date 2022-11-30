/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bli
 */

#include "BLI_math_base.h"
#include "BLI_math_matrix.hh"
#include "BLI_math_rotation.hh"
#include "BLI_math_vector.h"
#include "BLI_math_vector.hh"

namespace blender::rotation {

#ifdef DEBUG
#  define BLI_ASSERT_UNIT_QUATERNION(_q) \
    { \
      auto &rot_vec = static_cast<const vec_base<T, 4>>(_q); \
      T quat_length = math::length_squared(rot_vec); \
      if (!(quat_length == 0 || (math::abs(quat_length - 1) < 0.0001))) { \
        std::cout << "Warning! " << __func__ << " called with non-normalized quaternion: size " \
                  << quat_length << " *** report a bug ***\n"; \
      } \
    }
#else
#  define BLI_ASSERT_UNIT_QUATERNION(_q)
#endif

template<typename T> Quaternion<T>::operator EulerXYZ<T>()
{
  using Mat3T = mat_base<T, 3, 3>;
  const Quaternion<T> &quat = *this;
  BLI_ASSERT_UNIT_QUATERNION(quat)
  Mat3T unit_mat = math::from_rotation<Mat3T>(quat);
  return math::to_euler<T, true>(unit_mat);
}

template<typename T> EulerXYZ<T>::operator Quaternion<T>()
{
  const EulerXYZ<T> &eul = *this;
  const T ti = eul[0] * T(0.5);
  const T tj = eul[1] * T(0.5);
  const T th = eul[2] * T(0.5);
  const T ci = math::cos(ti);
  const T cj = math::cos(tj);
  const T ch = math::cos(th);
  const T si = math::sin(ti);
  const T sj = math::sin(tj);
  const T sh = math::sin(th);
  const T cc = ci * ch;
  const T cs = ci * sh;
  const T sc = si * ch;
  const T ss = si * sh;

  Quaternion<T> quat;
  quat[0] = cj * cc + sj * ss;
  quat[1] = cj * sc - sj * cs;
  quat[2] = cj * ss + sj * cc;
  quat[3] = cj * cs - sj * sc;
  return quat;
}

template<typename T> Quaternion<T>::operator AxisAngle<T>()
{
  const Quaternion<T> &quat = *this;
  BLI_ASSERT_UNIT_QUATERNION(quat)

  /* Calculate angle/2, and sin(angle/2). */
  float ha = math::acos(quat[0]);
  float si = math::sin(ha);

  AxisAngle<T> rot;
  /* From half-angle to angle. */
  rot.angle = ha * 2;
  /* Prevent division by zero for axis conversion. */
  if (math::abs(si) < 0.0005) {
    si = 1.0f;
  }

  rot.axis = vec_base<T, 3>(quat[1], quat[2], quat[3]) / si;
  if (math::is_zero(rot.axis)) {
    rot.axis[1] = 1.0f;
  }
  return rot;
}

template<typename T> AxisAngle<T>::operator Quaternion<T>()
{
  const AxisAngle<T> &rot = *this;
  BLI_assert(math::is_unit_scale(rot.axis));
  const T phi = 0.5f * rot.angle;
  const T cosine = math::cos(phi);
  const T sine = math::sin(phi);

  Quaternion<T> quat;
  quat[0] = cosine;
  quat[1] = axis[0] * sine;
  quat[2] = axis[1] * sine;
  quat[3] = axis[2] * sine;
  return quat;
}

template<typename T> EulerXYZ<T>::operator AxisAngle<T>()
{
  /* Use quaternions as intermediate representation for now... */
  return AxisAngle<T>(Quaternion<T>(*this));
}

template<typename T> AxisAngle<T>::operator EulerXYZ<T>()
{
  /* Use quaternions as intermediate representation for now... */
  return EulerXYZ<T>(Quaternion<T>(*this));
}

template Quaternion<float>::operator EulerXYZ<float>();
template EulerXYZ<float>::operator Quaternion<float>();
template Quaternion<float>::operator AxisAngle<float>();
template AxisAngle<float>::operator Quaternion<float>();
template EulerXYZ<float>::operator AxisAngle<float>();
template AxisAngle<float>::operator EulerXYZ<float>();

template Quaternion<double>::operator EulerXYZ<double>();
template EulerXYZ<double>::operator Quaternion<double>();
template Quaternion<double>::operator AxisAngle<double>();
template AxisAngle<double>::operator Quaternion<double>();
template EulerXYZ<double>::operator AxisAngle<double>();
template AxisAngle<double>::operator EulerXYZ<double>();

}  // namespace blender::rotation

namespace blender::math {

float3 rotate_direction_around_axis(const float3 &direction, const float3 &axis, const float angle)
{
  BLI_ASSERT_UNIT_V3(direction);
  BLI_ASSERT_UNIT_V3(axis);

  const float3 axis_scaled = axis * math::dot(direction, axis);
  const float3 diff = direction - axis_scaled;
  const float3 cross = math::cross(axis, diff);

  return axis_scaled + diff * std::cos(angle) + cross * std::sin(angle);
}

float3 rotate_around_axis(const float3 &vector,
                          const float3 &center,
                          const float3 &axis,
                          const float angle)

{
  float3 result = vector - center;
  float mat[3][3];
  axis_angle_normalized_to_mat3(mat, axis, angle);
  mul_m3_v3(mat, result);
  return result + center;
}

}  // namespace blender::math
