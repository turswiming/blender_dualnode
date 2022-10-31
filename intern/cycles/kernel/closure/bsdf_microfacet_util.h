/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/closure/bsdf_util.h"

#include "kernel/util/lookup_table.h"

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

/* Albedo correction. */

ccl_device_forceinline float microfacet_ggx_glass_E(KernelGlobals kg,
                                                    float mu,
                                                    float rough,
                                                    float ior)
{
  bool inv_table = (ior < 1.0f);
  int offset = inv_table ? kernel_data.tables.ggx_glass_inv_E_offset :
                           kernel_data.tables.ggx_glass_E_offset;

  float x = mu, y = 1 - rough;
  float z = sqrtf(0.5f * ((inv_table ? 1.0f / ior : ior) - 1.0f));
  return lookup_table_read_3D(kg, x, y, z, offset, 16, 16, 16);
}

ccl_device_forceinline float microfacet_ggx_dielectric_E(KernelGlobals kg,
                                                         float mu,
                                                         float rough,
                                                         float ior)
{
  bool inv_table = (ior < 1.0f);
  int offset = inv_table ? kernel_data.tables.ggx_dielectric_inv_E_offset :
                           kernel_data.tables.ggx_dielectric_E_offset;

  float macro_fresnel = fresnel_dielectric_cos(mu, ior);
  float F0 = fresnel_dielectric_F0(ior);
  float x = mix(mu, inverse_lerp(1.0f, F0, macro_fresnel), 0.5f);
  float y = 1 - rough;
  float z = sqrtf(0.5f * ((inv_table ? 1.0f / ior : ior) - 1.0f));

  return lookup_table_read_3D(kg, x, y, z, offset, 16, 16, 16);
}

ccl_device_forceinline float microfacet_ggx_E(KernelGlobals kg, float mu, float rough)
{
  return lookup_table_read_2D(kg, mu, 1 - rough, kernel_data.tables.ggx_E_offset, 32, 32);
}

ccl_device_forceinline float microfacet_ggx_E_avg(KernelGlobals kg, float rough)
{
  return lookup_table_read(kg, 1 - rough, kernel_data.tables.ggx_E_avg_offset, 32);
}

ccl_device_forceinline float clearcoat_E(KernelGlobals kg, float mu, float rough)
{
  float x = mu, y = 1 - rough;
  float table = lookup_table_read_2D(kg, x, y, kernel_data.tables.ggx_clearcoat_E_offset, 16, 16);
  return table * fresnel_dielectric_cos(mu, 1.5f);
}

ccl_device_inline Spectrum fresnel_metallic_Fss(Spectrum F0, Spectrum B)
{
  return saturate(mix(F0, one_spectrum(), 1.0f / 21.0f) - B * (1.0f / 126.0f));
}

ccl_device_inline Spectrum schlick_fresnel_Fss(Spectrum F0)
{
  return saturate(mix(F0, one_spectrum(), 1.0f / 21.0f));
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

/* Source for thin film iridescence model:
 * "A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence"
 * by Laurent Belcour and Pascal Barla (SIGGRAPH 2017)
 * https://belcour.github.io/blog/research/publication/2017/05/01/brdf-thin-film.html
 */

ccl_device_inline void fresnel_dielectric_complex(
    float cosTheta1, float eta1, float eta2, ccl_private float2 *R, ccl_private float2 *phi)
{

  float sinTheta1Sqr = 1.0f - sqr(cosTheta1);

  float etaSqr = sqr(eta1 / eta2);

  /* TODO: Simplify. */
  if (etaSqr * sinTheta1Sqr > 1) {
    phi->x = 2.0f * atanf(-etaSqr * sqrtf(sinTheta1Sqr - 1.0f / etaSqr) / cosTheta1);
    phi->y = 2.0f * atanf(-sqrtf(sinTheta1Sqr - 1.0f / etaSqr) / cosTheta1);
    *R = one_float2();
  }
  else {
    float cosTheta2 = sqrtf(1.0f - etaSqr * sinTheta1Sqr);
    float rx = (eta2 * cosTheta1 - eta1 * cosTheta2) / (eta2 * cosTheta1 + eta1 * cosTheta2);
    float ry = (eta1 * cosTheta1 - eta2 * cosTheta2) / (eta1 * cosTheta1 + eta2 * cosTheta2);
    phi->x = (rx < 0.0f) ? M_PI_F : 0.0f;
    phi->y = (ry < 0.0f) ? M_PI_F : 0.0f;
    *R = sqr(make_float2(rx, ry));
  }
}

ccl_device_inline float thin_film_sensitivity_component(
    float phase, float shift, float val, float pos, float var)
{
  return val * cosf(pos * phase + shift) * expf(-var * sqr(phase));
}

ccl_device_inline float3 thin_film_sensitivity(float opd, float shift)
{
  float phase = M_2PI_F * opd * 1.0e-9f;
  return make_float3(
      thin_film_sensitivity_component(phase, shift, 0.8465900f, 1.6810e+06f, 4.3278e+09f) +
          thin_film_sensitivity_component(phase, shift, 0.1538683f, 2.2399e+06f, 4.5282e+09f),
      thin_film_sensitivity_component(phase, shift, 1.0002218f, 1.7953e+06f, 9.3046e+09f),
      thin_film_sensitivity_component(phase, shift, 1.0011225f, 2.2084e+06f, 6.6121e+09f));
}

ccl_device_inline float3 fresnel_dielectric_thin_film(float cosThetaOuter,
                                                      float iorInner,
                                                      float iorOuter,
                                                      float thickness)
{
  float cosThetaInner = sqrtf(1.0f - sqr(1.0f / iorOuter) * (1.0f - sqr(cosThetaOuter)));

  float2 R12, phi12;
  fresnel_dielectric_complex(cosThetaOuter, 1.0f, iorOuter, &R12, &phi12);

  float2 R23, phi23;
  fresnel_dielectric_complex(cosThetaInner, iorOuter, iorInner, &R23, &phi23);

  float2 T121 = one_float2() - R12;
  float2 R123 = R12 * R23;
  float2 Rs = sqr(T121) * R23 / (one_float2() - R123);
  float2 Cm = Rs - T121;

  float opd = thickness * cosThetaInner;
  float2 phi2 = make_float2(M_PI_F, M_PI_F) - phi12 + phi23;
  float2 r123 = sqrt(R123);

  float3 I = average(R12 + Rs) * thin_film_sensitivity(0.0f, 0.0f);
  for (int m = 1; m < 4; m++) {
    Cm *= r123;
    float3 SmSh = thin_film_sensitivity(m * opd, m * phi2.x);
    float3 SmPh = thin_film_sensitivity(m * opd, m * phi2.y);
    I += Cm.x * SmSh + Cm.y * SmPh;
  }

  kernel_assert(isfinite_safe(I));

  /* CIE XYZ -> sRGB incl. whitepoint adaption E -> D65
   * TODO: Somehow properly handle color management here?
   * Does OCIO have anything for handling reflectances? */
  return saturate(make_float3(3.17482107f, -0.98993255f, 0.06047909f) * I.x +
                  make_float3(-1.69745555f, 1.94998695f, -0.20891802f) * I.y +
                  make_float3(-0.47752224f, 0.04004475f, 1.14851336f) * I.z);
}

CCL_NAMESPACE_END
