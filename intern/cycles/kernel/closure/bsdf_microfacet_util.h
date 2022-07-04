/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/closure/bsdf_util.h"

CCL_NAMESPACE_BEGIN

/* GGX microfacet with Smith shadow-masking from:
 *
 * Microfacet Models for Refraction through Rough Surfaces
 * B. Walter, S. R. Marschner, H. Li, K. E. Torrance, EGSR 2007
 *
 * VNDF sampling as well as D and lambda terms from:
 *
 * Sampling the GGX Distribution of Visible Normals.
 * E. Heitz and E. d'Eon, JCGT Vol. 7, No. 4, 2018.
 * https://jcgt.org/published/0007/04/01/
 *
 * Also see for more details on marking-shadowing:
 * Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs.
 * E. Heitz, JCGT Vol. 3, No. 2, 2014.
 * https://jcgt.org/published/0003/02/03/ */

ccl_device_forceinline float microfacet_ggx_lambda(const float cosTheta, const float alpha2)
{
  float tanTheta2 = (1 - sqr(cosTheta)) / sqr(cosTheta);
  return 0.5f * (safe_sqrtf(1 + alpha2 * tanTheta2) - 1);
}

ccl_device_forceinline float microfacet_ggx_lambda_aniso(const float3 w,
                                                         const float alpha_x,
                                                         const float alpha_y)
{
  return 0.5f * (safe_sqrtf(1 + (sqr(w.x * alpha_x) + sqr(w.y * alpha_y)) / sqr(w.z)) - 1);
}

ccl_device_forceinline float microfacet_ggx_D(const float cosThetaM, const float alpha2)
{
  float cosThetaM2 = sqr(cosThetaM);
  float tanThetaM2 = (1 - cosThetaM2) / cosThetaM2;
  return alpha2 / (M_PI_F * sqr(cosThetaM2 * (alpha2 + tanThetaM2)));
}

ccl_device_forceinline float microfacet_ggx_D_aniso(const float3 m,
                                                    const float alpha_x,
                                                    const float alpha_y)
{
  return 1 /
         (M_PI_F * alpha_x * alpha_y * sqr(sqr(m.x / alpha_x) + sqr(m.y / alpha_y) + sqr(m.z)));
}

ccl_device_forceinline float microfacet_GTR1_D(const float cosThetaM, const float alpha2)
{
  if (alpha2 >= 1.0f)
    return M_1_PI_F;
  float t = 1.0f + (alpha2 - 1.0f) * sqr(cosThetaM);
  return (alpha2 - 1.0f) / (M_PI_F * logf(alpha2) * t);
}

ccl_device_forceinline float3 microfacet_ggx_sample_vndf(
    const float3 V, const float alpha_x, const float alpha_y, const float U1, const float U2)
{
  /* Section 3.2: Transforming the view direction to the hemisphere configuration. */
  float3 Vh = normalize(make_float3(alpha_x * V.x, alpha_y * V.y, V.z));
  /* Section 4.1: Orthonormal basis (with special case if cross product is zero). */
  float lensq = sqr(Vh.x) + sqr(Vh.y);
  float3 T1 = lensq > 1e-7f ? make_float3(-Vh.y, Vh.x, 0.0f) / sqrtf(lensq) :
                              make_float3(1.0f, 0.0f, 0.0f);
  float3 T2 = cross(Vh, T1);
  /* Section 4.2: Parameterization of the projected area. */
  float2 t = concentric_sample_disk(U1, U2);
  t.y = mix(safe_sqrtf(1.0f - sqr(t.x)), t.y, 0.5f * (1.0f + Vh.z));
  /* Section 4.3: Reprojection onto hemisphere. */
  float3 Mh = t.x * T1 + t.y * T2 + safe_sqrtf(1.0f - len_squared(t)) * Vh;
  /* Section 3.4: Transforming the normal back to the ellipsoid configuration. */
  return normalize(make_float3(alpha_x * Mh.x, alpha_y * Mh.y, max(0.0f, Mh.z)));
}

/* Albedo correction.
 * TODO: Use proper lookup table code for this data. */

ccl_device_forceinline float microfacet_ggx_glass_E(float mu, float rough, float ior)
{
  bool inv_table = (ior < 1.0f);
  if (inv_table) {
    ior = 1.0f / ior;
  }

  rough = saturatef(1 - rough - 1.0f/32.0f) * 16.0f;
  mu = saturatef(mu - 1.0f/32.0f) * 16.0f;
  ior = saturatef(sqrtf(0.5f * (ior - 1.0f)) - 1.0f/32.0f) * 16.0f;

  int rough_i = min(15, (int)rough);
  int rough_i1 = min(15, rough_i + 1);
  int mu_i = min(15, (int)mu);
  int mu_i1 = min(15, mu_i + 1);
  int ior_i = min(15, (int)ior);
  int ior_i1 = min(15, ior_i + 1);

  rough -= rough_i;
  mu -= mu_i;
  ior -= ior_i;

  auto &table = inv_table ? table_ggx_glass_inv_E : table_ggx_glass_E;
  float a = lerp(table[ior_i][rough_i][mu_i], table[ior_i][rough_i][mu_i1], mu);
  float b = lerp(table[ior_i][rough_i1][mu_i], table[ior_i][rough_i1][mu_i1], mu);
  float c = lerp(table[ior_i1][rough_i][mu_i], table[ior_i1][rough_i][mu_i1], mu);
  float d = lerp(table[ior_i1][rough_i1][mu_i], table[ior_i1][rough_i1][mu_i1], mu);

  return saturatef(lerp(lerp(a, b, rough), lerp(c, d, rough), ior));
}

ccl_device_forceinline float microfacet_ggx_dielectric_E(float mu, float rough, float ior)
{
  mu = lerp(mu, inverse_lerp(1.0f, fresnel_dielectric_cos(1.0f, ior), fresnel_dielectric_cos(mu, ior)), 0.5f);

  bool inv_table = (ior < 1.0f);
  if (inv_table) {
    ior = 1.0f / ior;
  }

  rough = saturatef(1 - rough - 1.0f/32.0f) * 16.0f;
  mu = saturatef(mu - 1.0f/32.0f) * 16.0f;
  ior = saturatef(sqrtf(0.5f * (ior - 1.0f)) - 1.0f/32.0f) * 16.0f;

  int rough_i = min(15, (int)rough);
  int rough_i1 = min(15, rough_i + 1);
  int mu_i = min(15, (int)mu);
  int mu_i1 = min(15, mu_i + 1);
  int ior_i = min(15, (int)ior);
  int ior_i1 = min(15, ior_i + 1);

  rough -= rough_i;
  mu -= mu_i;
  ior -= ior_i;

  auto &table = inv_table ? table_ggx_dielectric_inv_E : table_ggx_dielectric_E;
  float a = lerp(table[ior_i][rough_i][mu_i], table[ior_i][rough_i][mu_i1], mu);
  float b = lerp(table[ior_i][rough_i1][mu_i], table[ior_i][rough_i1][mu_i1], mu);
  float c = lerp(table[ior_i1][rough_i][mu_i], table[ior_i1][rough_i][mu_i1], mu);
  float d = lerp(table[ior_i1][rough_i1][mu_i], table[ior_i1][rough_i1][mu_i1], mu);

  return saturatef(lerp(lerp(a, b, rough), lerp(c, d, rough), ior));
}

ccl_device_forceinline float microfacet_ggx_E(float mu, float rough)
{
  rough = saturatef(1 - rough - 1.0f/64.0f) * 32.0f;
  mu = saturatef(mu - 1.0f/64.0f) * 32.0f;

  int rough_i = min(31, (int)rough);
  int rough_i1 = min(31, rough_i + 1);
  int mu_i = min(31, (int)mu);
  int mu_i1 = min(31, mu_i + 1);

  rough -= rough_i;
  mu -= mu_i;

  float a = lerp(table_ggx_E[rough_i][mu_i], table_ggx_E[rough_i][mu_i1], mu);
  float b = lerp(table_ggx_E[rough_i1][mu_i], table_ggx_E[rough_i1][mu_i1], mu);
  return saturatef(lerp(a, b, rough));
}

ccl_device_forceinline float microfacet_ggx_E_avg(float rough)
{
  rough = saturatef(1 - rough - 1.0f/64.0f) * 32.0f;
  int rough_i = min(31, (int)rough);
  int rough_i1 = min(31, rough_i + 1);
  rough -= rough_i;
  return saturatef(lerp(table_ggx_E_avg[rough_i], table_ggx_E_avg[rough_i1], rough));
}

ccl_device_inline float3 metallic_Fss(float3 F0, float3 F90, float falloff)
{
  /* Fss for lerp(F0, F90, (1-cosNI)^falloff) */
  return lerp(F0, F90, 2.0f / (sqr(falloff) + 3*falloff + 2));
}

ccl_device_inline float3 schlick_fresnel_Fss(float3 F0)
{
  /* TODO validate using multiGGX code */
  return (one_float3() + 20.0f * saturate(F0)) * (1.0f / 21.0f);
}

/* TODO Imageworks source */
ccl_device_inline float dielectric_fresnel_Fss(float eta)
{
  /* TODO validate using multiGGX code */
  float f;
  if (eta < 1.0f) {
    f = 0.997118f + eta * (0.1014f - eta * (0.965241f + eta * 0.130607f));
  }
  else {
    f = (eta - 1.0f) / (4.08567f + 1.00071f * eta);
  }
  return f;
}

CCL_NAMESPACE_END
