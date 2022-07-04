/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/sample/lcg.h"

CCL_NAMESPACE_BEGIN

ccl_device_inline float3
microfacet_ggx_glass_albedo_scaling(KernelGlobals kg,
                                    ccl_private const ShaderData *sd,
                                    ccl_private const MicrofacetBsdf *bsdf,
                                    const float3 Fss)
{
  float mu = dot(sd->I, bsdf->N);
  float rough = sqrtf(sqrtf(bsdf->alpha_x * bsdf->alpha_y));
  float E = microfacet_ggx_glass_E(kg, mu, rough, bsdf->ior);

  /* Close enough for glass, coloring here is unphysical anyways and it's unclear how to
   * approximate it better. */
  float3 Fms = Fss;

  return one_float3() + Fms * ((1.0f - E) / E);
  /* TODO: Ensure that increase in weight does not mess up glossy color, albedo etc. passes */
}

/* Currently no non-albedo-scaled version is implemented, could easily be added
 * but would still break compatibility with the old glass due to the microfacet Fresnel. */

ccl_device int bsdf_microfacet_multi_ggx_glass_setup(KernelGlobals kg,
                                                     ccl_private MicrofacetBsdf *bsdf,
                                                     ccl_private const ShaderData *sd,
                                                     const float3 color)
{
  bsdf->extra = NULL;

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = bsdf->alpha_x;

  bsdf->type = CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_ID;

  bsdf->weight *= microfacet_ggx_glass_albedo_scaling(kg, sd, bsdf, saturate(color));

  return SD_BSDF | SD_BSDF_HAS_EVAL | SD_BSDF_NEEDS_LCG;
}

ccl_device int bsdf_microfacet_multi_ggx_glass_fresnel_setup(KernelGlobals kg,
                                                             ccl_private MicrofacetBsdf *bsdf,
                                                             ccl_private const ShaderData *sd)
{
  bsdf->extra->cspec0 = saturate(bsdf->extra->cspec0);

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = bsdf->alpha_x;

  bsdf->type = CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID;

  bsdf_microfacet_fresnel_color(sd, bsdf);

  float3 Fss = schlick_fresnel_Fss(bsdf->extra->cspec0);
  bsdf->weight *= microfacet_ggx_glass_albedo_scaling(kg, sd, bsdf, Fss);

  return SD_BSDF | SD_BSDF_HAS_EVAL | SD_BSDF_NEEDS_LCG;
}

ccl_device float3 bsdf_microfacet_ggx_glass_eval_reflect(ccl_private const ShaderClosure *sc,
                                                         const float3 I,
                                                         const float3 omega_in,
                                                         ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  kernel_assert(alpha_x * alpha_y > 1e-7f);
  float3 N = bsdf->N;

  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);
  if (cosNI <= 0 || cosNO <= 0) {
    *pdf = 0.0f;
    return zero_float3();
  }

  float3 m = normalize(omega_in + I);
  float alpha2 = alpha_x * alpha_y;
  float D = microfacet_ggx_D(dot(N, m), alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  float common = D * 0.25f / cosNO;

  float F = fresnel_dielectric_cos(dot(m, I), bsdf->ior);
  float out = F * common / (1 + lambdaO + lambdaI);
  *pdf = common / (1 + lambdaO);

  float3 eval = make_float3(out, out, out);
  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    eval *= reflection_color(bsdf, omega_in, m);
  }
  return eval;
}

ccl_device float3 bsdf_microfacet_ggx_glass_eval_transmit(ccl_private const ShaderClosure *sc,
                                                          const float3 I,
                                                          const float3 omega_in,
                                                          ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  kernel_assert(alpha_x * alpha_y > 1e-7f);
  float eta = bsdf->ior;
  float3 N = bsdf->N;

  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);
  if (cosNO <= 0 || cosNI >= 0) {
    *pdf = 0.0f;
    return zero_float3();
  }

  float3 ht = -(eta * omega_in + I);
  float3 m = normalize(ht);
  float cosMO = dot(m, I);
  float cosMI = dot(m, omega_in);

  float F = fresnel_dielectric_cos(cosMO, eta);
  if (F == 1.0f) {
    /* TIR */
    *pdf = 0.0f;
    return zero_float3();
  }

  float alpha2 = alpha_x * alpha_y;
  float D = microfacet_ggx_D(dot(N, m), alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  float Ht2 = dot(ht, ht);

  float common = fabsf(cosMI * cosMO) * D * sqr(eta) / (cosNO * Ht2);
  float out = (1.0f - F) * common / (1 + lambdaO + lambdaI);
  *pdf = common / (1 + lambdaO);

  float3 eval = make_float3(out, out, out);
  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    eval *= bsdf->extra->color;
  }
  return eval;
}

ccl_device int bsdf_microfacet_ggx_glass_sample(ccl_private const ShaderClosure *sc,
                                                float3 Ng,
                                                float3 I,
                                                float3 dIdx,
                                                float3 dIdy,
                                                float randu,
                                                float randv,
                                                ccl_private float3 *eval,
                                                ccl_private float3 *omega_in,
                                                ccl_private float3 *domega_in_dx,
                                                ccl_private float3 *domega_in_dy,
                                                ccl_private float *pdf,
                                                ccl_private uint *lcg_state)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  kernel_assert(alpha_x * alpha_y > 1e-7f);
  float eta = bsdf->ior;
  float3 N = bsdf->N;
  int label;

  float cosNO = dot(N, I);
  if (cosNO <= 0) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 X, Y, Z = N;
  make_orthonormals(Z, &X, &Y);

  /* importance sampling with distribution of visible normals. vectors are
   * transformed to local space before and after */
  float3 local_O = make_float3(dot(X, I), dot(Y, I), cosNO);
  float3 local_m = microfacet_ggx_sample_vndf(local_O, alpha_x, alpha_y, randu, randv);

  float3 m = X * local_m.x + Y * local_m.y + Z * local_m.z;
  float cosThetaM = local_m.z;

  float cosMO = dot(m, I);
  if (cosMO <= 0.0f) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  float3 R, T;
#ifdef __RAY_DIFFERENTIALS__
  float3 dRdx, dRdy, dTdx, dTdy;
#endif
  bool inside; /* Will never be inside, we already checked cosMO */
  float fresnel = fresnel_dielectric(eta,
                                     m,
                                     I,
                                     &R,
                                     &T,
#ifdef __RAY_DIFFERENTIALS__
                                     dIdx,
                                     dIdy,
                                     &dRdx,
                                     &dRdy,
                                     &dTdx,
                                     &dTdy,
#endif
                                     &inside);

  float randw = lcg_step_float(lcg_state);
  bool do_reflect = randw < fresnel;

  /* Common microfacet model terms. */
  float alpha2 = alpha_x * alpha_y;
  float D = microfacet_ggx_D(cosThetaM, alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);

  float cosNI, common;
  if (do_reflect) {
    cosNI = dot(N, R);
    if (cosNI <= 0.0f || dot(Ng, R) <= 0.0f) {
      *pdf = 0.0f;
      return LABEL_NONE;
    }

    label = LABEL_REFLECT | LABEL_GLOSSY;
    *omega_in = R;

#ifdef __RAY_DIFFERENTIALS__
    *domega_in_dx = dRdx;
    *domega_in_dy = dRdy;
#endif

    common = fresnel * D * 0.25f / cosNO;
  }
  else {
    cosNI = dot(N, T);
    if (cosNI >= 0.0f || dot(Ng, T) >= 0.0f) {
      *pdf = 0.0f;
      return LABEL_NONE;
    }

    label = LABEL_TRANSMIT | LABEL_GLOSSY;
    *omega_in = T;
#ifdef __RAY_DIFFERENTIALS__
    *domega_in_dx = dTdx;
    *domega_in_dy = dTdy;
#endif

    float cosMI = dot(m, *omega_in);
    float Ht2 = sqr(eta * cosMI + cosMO);

    common = (1.0f - fresnel) * D * fabsf(cosMI * cosMO) * sqr(eta) / (cosNO * Ht2);
  }

  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);
  float out = common / (1 + lambdaO + lambdaI);
  *pdf = common / (1 + lambdaO);
  *eval = make_float3(out, out, out);

  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    *eval *= (label & LABEL_REFLECT) ? reflection_color(bsdf, *omega_in, m) : bsdf->extra->color;
  }
  return label;
}

CCL_NAMESPACE_END
