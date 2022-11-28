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

template<typename T> explicit EulerXYZ<T> &operator EulerXYZ<T>(const Quaternion<T> &rot)
{
  using namespace math::mat3x3;
  BLI_ASSERT_UNIT_QUATERNION(rot)
  mat_base<T, 3, 3> unit_mat = from_rotation(rot);
  return to_euler<true>(unit_mat);
}

template<typename T> explicit Quaternion<T> &operator Quaternion<T>(const EulerXYZ<T> &rot)
{
  T ti = eul[0] * T(0.5);
  T tj = eul[1] * T(0.5);
  T th = eul[2] * T(0.5);
  T ci = math::cos(ti);
  T cj = math::cos(tj);
  T ch = math::cos(th);
  T si = math::sin(ti);
  T sj = math::sin(tj);
  T sh = math::sin(th);
  T cc = ci * ch;
  T cs = ci * sh;
  T sc = si * ch;
  T ss = si * sh;

  Quaternion<T> quat;
  quat[0] = cj * cc + sj * ss;
  quat[1] = cj * sc - sj * cs;
  quat[2] = cj * ss + sj * cc;
  quat[3] = cj * cs - sj * sc;
  return quat;
}

template<typename T> explicit AxisAngle<T> &operator AxisAngle<T>(const Quaternion<T> &rot)
{
  BLI_ASSERT_UNIT_QUATERNION(rot)

  /* Calculate angle/2, and sin(angle/2). */
  float ha = math::acos(q[0]);
  float si = math::sin(ha);
  /* From half-angle to angle. */
  *angle = ha * 2;
  /* Prevent division by zero for axis conversion. */
  if (fabsf(si) < 0.0005f) {
    si = 1.0f;
  }

  axis[0] = q[1] / si;
  axis[1] = q[2] / si;
  axis[2] = q[3] / si;
  if (is_zero_v3(axis)) {
    axis[1] = 1.0f;
  }
  return;
}

template<typename T> explicit Quaternion<T> &operator Quaternion<T>(const AxisAngle<T> &rot)
{
  BLI_assert(math::is_unit_scale(rot.axis()));
  const T phi = 0.5f * rot.angle();
  const T cosine = math::cos(phi);
  const T sine = math::sin(phi);
  Quaternion<T> quat;
  quat[0] = cosine;
  for (int i = 1; i < 4; i++) {
    quat[i] = axis[i] * sine;
  }
  return;
}

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
