/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

/* DISNEY PRINCIPLED SHEEN BRDF
 *
 * Shading model by Brent Burley (Disney): "Physically Based Shading at Disney" (2012)
 */

#include "kernel/closure/bsdf_util.h"

CCL_NAMESPACE_BEGIN

typedef struct PrincipledSheenBsdf {
  SHADER_CLOSURE_BASE;
  float avg_value;
  float roughness;
} PrincipledSheenBsdf;

static_assert(sizeof(ShaderClosure) >= sizeof(PrincipledSheenBsdf),
              "PrincipledSheenBsdf is too large!");

ccl_device_inline float calculate_avg_principled_sheen_brdf(float3 N, float3 I)
{
  /* To compute the average, we set the half-vector to the normal, resulting in
   * NdotI = NdotL = NdotV = LdotH */
  float NdotI = dot(N, I);
  if (NdotI < 0.0f) {
    return 0.0f;
  }

  return schlick_fresnel(NdotI) * NdotI;
}

ccl_device float3
calculate_principled_sheen_brdf(float3 N, float3 V, float3 L, float3 H, ccl_private float *pdf)
{
  float NdotL = dot(N, L);
  float NdotV = dot(N, V);

  if (NdotL < 0 || NdotV < 0) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  float LdotH = dot(L, H);

  float value = schlick_fresnel(LdotH) * NdotL;

  return make_float3(value, value, value);
}

/* Based on
 * https://dassaultsystemes-technology.github.io/EnterprisePBRShadingModel/spec-2022x.md.html#components/sheen.
 */
ccl_device_inline float sheen_v2_lambda(float mu, float w)
{
  float a = lerp(11.9095f, 13.7000f, w);
  float b = lerp(4.68753f, 2.92754f, w);
  float c = lerp(0.33467f, 0.28670f, w);
  float d = lerp(-2.22664f, -0.81757f, w);
  float e = lerp(-1.76591f, -1.22466f, w);

  float exponent;
  if (mu < 0.5f) {
    exponent = a / (1 + b * powf(mu, c)) + d * mu + e;
  }
  else {
    exponent = 2 * a / (1 + b * exp2(-c)) - a / (1 + b * powf(1 - mu, c)) + d * mu + e;
  }
  return expf(exponent);
}

ccl_device_inline float3 sheen_v2_eval(float3 N, float3 V, float3 L, float3 H, float r, float *pdf)
{
  float cosNH = dot(N, H), cosNV = dot(N, V), cosNL = dot(N, L);

  if (cosNH < 0 || cosNV < 0 || cosNL < 0) {
    return zero_float3();
  }

  /* Evaluate microfacet distribution. */
  float sinTheta2 = 1 - sqr(cosNH);
  float invR = 1 / r;
  float D = M_1_2PI_F * (2 + invR) * powf(sinTheta2, 0.5f * invR);

  /* Evaluate shadowing-masking terms. */
  float w = -1.59612f / (1 + 0.20375f * fast_safe_powf(r, -0.55825f)) + 1.32805f;
  float lambdaV = sheen_v2_lambda(cosNV, w);
  float lambdaL = sheen_v2_lambda(cosNL, w);

  /* Soften shadow terminator. */
  lambdaL = fast_safe_powf(lambdaL, 1 + 2 * sqr(sqr(sqr(1 - cosNL))));

  /* Combined microfacet BSDF.
   * Usual form is F*D*G/(4*cosNL*cosNV), but here we have no Fresnel, we skip dividing by cosNL
   * since Cycles convention is returning BSDF*cosNL, and we use the combined shadowing-masking
   * term G=1/(1+lambdaV+lambdaL).
   */
  float val = D / (4 * cosNV * (1 + lambdaV + lambdaL));
  return make_float3(val, val, val);
}

ccl_device_forceinline float sheen_v2_E(float mu, float rough)
{
  /* TODO: Cleanup. */
  rough = saturatef(1 - sqrtf(rough) - 1.0f / 64.0f) * 32.0f;
  mu = saturatef(mu - 1.0f / 64.0f) * 32.0f;

  int rough_i = min(31, (int)rough);
  int rough_i1 = min(31, rough_i + 1);
  int mu_i = min(31, (int)mu);
  int mu_i1 = min(31, mu_i + 1);

  rough -= rough_i;
  mu -= mu_i;

  float a = lerp(table_sheen_E[rough_i][mu_i], table_sheen_E[rough_i][mu_i1], mu);
  float b = lerp(table_sheen_E[rough_i1][mu_i], table_sheen_E[rough_i1][mu_i1], mu);
  return lerp(a, b, rough);
}

ccl_device int bsdf_principled_sheen_setup(ccl_private const ShaderData *sd,
                                           ccl_private PrincipledSheenBsdf *bsdf)
{
  bsdf->type = CLOSURE_BSDF_PRINCIPLED_SHEEN_ID;
  bsdf->avg_value = calculate_avg_principled_sheen_brdf(bsdf->N, sd->I);
  bsdf->sample_weight *= bsdf->avg_value;
  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device int bsdf_principled_sheen_v2_setup(ccl_private const ShaderData *sd,
                                              ccl_private PrincipledSheenBsdf *bsdf)
{
  // TODO: Also expose as separate node. Add enum to Velvet BSDF maybe?
  bsdf->type = CLOSURE_BSDF_PRINCIPLED_SHEEN_V2_ID;

  bsdf->roughness = clamp(bsdf->roughness, 1e-3f, 1.0f);

  bsdf->avg_value = sheen_v2_E(dot(bsdf->N, sd->I), bsdf->roughness);
  bsdf->sample_weight *= bsdf->avg_value;

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device float3 bsdf_principled_sheen_eval_reflect(ccl_private const ShaderClosure *sc,
                                                     const float3 I,
                                                     const float3 omega_in,
                                                     ccl_private float *pdf)
{
  ccl_private const PrincipledSheenBsdf *bsdf = (ccl_private const PrincipledSheenBsdf *)sc;

  float3 N = bsdf->N;
  float3 V = I;         // outgoing
  float3 L = omega_in;  // incoming
  float3 H = normalize(L + V);

  if (dot(N, omega_in) > 0.0f) {
    *pdf = fmaxf(dot(N, omega_in), 0.0f) * M_1_PI_F;
    if (bsdf->type == CLOSURE_BSDF_PRINCIPLED_SHEEN_V2_ID) {
      return sheen_v2_eval(N, V, L, H, bsdf->roughness, pdf);
    }
    else {
      return calculate_principled_sheen_brdf(N, V, L, H, pdf);
    }
  }
  else {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }
}

ccl_device float3 bsdf_principled_sheen_eval_transmit(ccl_private const ShaderClosure *sc,
                                                      const float3 I,
                                                      const float3 omega_in,
                                                      ccl_private float *pdf)
{
  *pdf = 0.0f;
  return make_float3(0.0f, 0.0f, 0.0f);
}

ccl_device int bsdf_principled_sheen_sample(ccl_private const ShaderClosure *sc,
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
                                            ccl_private float *pdf)
{
  ccl_private const PrincipledSheenBsdf *bsdf = (ccl_private const PrincipledSheenBsdf *)sc;

  float3 N = bsdf->N;

  sample_cos_hemisphere(N, randu, randv, omega_in, pdf);

  if (dot(Ng, *omega_in) > 0) {
    float3 H = normalize(I + *omega_in);

    if (bsdf->type == CLOSURE_BSDF_PRINCIPLED_SHEEN_V2_ID) {
      *eval = sheen_v2_eval(N, I, *omega_in, H, bsdf->roughness, pdf);
    }
    else {
      *eval = calculate_principled_sheen_brdf(N, I, *omega_in, H, pdf);
    }

#ifdef __RAY_DIFFERENTIALS__
    // TODO: find a better approximation for the diffuse bounce
    *domega_in_dx = -((2 * dot(N, dIdx)) * N - dIdx);
    *domega_in_dy = -((2 * dot(N, dIdy)) * N - dIdy);
#endif
  }
  else {
    *eval = make_float3(0.0f, 0.0f, 0.0f);
    *pdf = 0.0f;
  }
  return LABEL_REFLECT | LABEL_DIFFUSE;
}

CCL_NAMESPACE_END
