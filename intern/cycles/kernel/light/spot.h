/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/light/common.h"

CCL_NAMESPACE_BEGIN

ccl_device float spot_light_attenuation(float3 dir, float spot_angle, float spot_smooth, float3 N)
{
  float attenuation = dot(dir, N);

  if (attenuation <= spot_angle) {
    attenuation = 0.0f;
  }
  else {
    float t = attenuation - spot_angle;

    if (t < spot_smooth && spot_smooth != 0.0f)
      attenuation *= smoothstepf(t / spot_smooth);
  }

  return attenuation;
}

template<bool in_volume_segment>
ccl_device_inline bool spot_light_sample(const ccl_global KernelLight *klight,
                                         const float randu,
                                         const float randv,
                                         const float3 P,
                                         ccl_private LightSample *ls)
{
  ls->P = make_float3(klight->co[0], klight->co[1], klight->co[2]);

  const float3 center = make_float3(klight->co[0], klight->co[1], klight->co[2]);
  const float radius = klight->spot.radius;
  const float3 dir = make_float3(klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);
  /* disk oriented normal */
  const float3 lightN = normalize(P - center);
  ls->P = center;

  if (radius > 0.0f) {
    /* disk light */
    ls->P += disk_light_sample(lightN, randu, randv) * radius;
  }

  const float invarea = klight->spot.invarea;
  ls->pdf = invarea;

  ls->D = normalize_len(ls->P - P, &ls->t);
  /* we set the light normal to the outgoing direction to support texturing */
  ls->Ng = -ls->D;

  ls->eval_fac = (0.25f * M_1_PI_F) * invarea;

  /* spot light attenuation */
  ls->eval_fac *= spot_light_attenuation(
      dir, klight->spot.spot_angle, klight->spot.spot_smooth, -ls->D);
  if (!in_volume_segment && ls->eval_fac == 0.0f) {
    return false;
  }

  float2 uv = map_to_sphere(ls->Ng);
  ls->u = uv.x;
  ls->v = uv.y;

  ls->pdf *= lamp_light_pdf(lightN, -ls->D, ls->t);
  return true;
}

ccl_device_forceinline void spot_light_update_position(const ccl_global KernelLight *klight,
                                                       ccl_private LightSample *ls,
                                                       const float3 P)
{
  ls->D = normalize_len(ls->P - P, &ls->t);
  ls->Ng = -ls->D;

  float2 uv = map_to_sphere(ls->Ng);
  ls->u = uv.x;
  ls->v = uv.y;

  float invarea = klight->spot.invarea;
  ls->eval_fac = (0.25f * M_1_PI_F) * invarea;
  ls->pdf = invarea;

  /* spot light attenuation */
  float3 dir = make_float3(klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);
  ls->eval_fac *= spot_light_attenuation(
      dir, klight->spot.spot_angle, klight->spot.spot_smooth, ls->Ng);
}

ccl_device_inline bool spot_light_intersect(const ccl_global KernelLight *klight,
                                            const ccl_private Ray *ccl_restrict ray,
                                            ccl_private float *t)
{
  /* Spot/Disk light. */
  const float3 lightP = make_float3(klight->co[0], klight->co[1], klight->co[2]);
  const float radius = klight->spot.radius;
  if (radius == 0.0f) {
    return false;
  }
  /* disk oriented normal */
  const float3 lightN = normalize(ray->P - lightP);
  /* One sided. */
  if (dot(ray->D, lightN) >= 0.0f) {
    return false;
  }

  float3 P;
  return ray_disk_intersect(ray->P, ray->D, ray->tmin, ray->tmax, lightP, lightN, radius, &P, t);
}

ccl_device_inline bool spot_light_sample_from_intersection(
    const ccl_global KernelLight *klight,
    ccl_private const Intersection *ccl_restrict isect,
    const float3 ray_P,
    const float3 ray_D,
    ccl_private LightSample *ccl_restrict ls)
{
  const float3 center = make_float3(klight->co[0], klight->co[1], klight->co[2]);
  const float3 dir = make_float3(klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);
  /* the normal of the oriented disk */
  const float3 lightN = normalize(ray_P - center);
  /* We set the light normal to the outgoing direction to support texturing. */
  ls->Ng = -ls->D;

  float invarea = klight->spot.invarea;
  ls->eval_fac = (0.25f * M_1_PI_F) * invarea;
  ls->pdf = invarea;

  /* spot light attenuation */
  ls->eval_fac *= spot_light_attenuation(
      dir, klight->spot.spot_angle, klight->spot.spot_smooth, -ls->D);

  if (ls->eval_fac == 0.0f) {
    return false;
  }

  float2 uv = map_to_sphere(ls->Ng);
  ls->u = uv.x;
  ls->v = uv.y;

  /* compute pdf */
  if (ls->t != FLT_MAX) {
    ls->pdf *= lamp_light_pdf(lightN, -ls->D, ls->t);
  }
  else {
    ls->pdf = 0.f;
  }

  return true;
}

ccl_device_inline float spot_light_tree_weight(const ccl_global KernelLight *klight,
                                               const float3 P,
                                               const float3 N)
{
  const float3 light_P = make_float3(klight->co[0], klight->co[1], klight->co[2]);

  const float radius = klight->spot.radius;
  const float cos_theta = klight->spot.spot_angle;
  const float theta = fast_acosf(cos_theta);
  const float3 light_dir = make_float3(
      klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);

  const float h1 = radius * fast_sinf(theta);
  const float d1 = radius * cos_theta;
  const float h2 = d1 / fast_tanf(theta);

  const float3 apex = light_P - (h1 + h2) * light_dir;
  const float3 apex_to_point = normalize(P - apex);

  return (dot(apex_to_point, light_dir) < cos_theta) ? 0.0f : 1.0f;
}

CCL_NAMESPACE_END
