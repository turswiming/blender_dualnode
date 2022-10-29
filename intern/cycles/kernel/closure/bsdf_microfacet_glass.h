/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/sample/lcg.h"

#include "util/hash.h"

CCL_NAMESPACE_BEGIN

ccl_device_inline Spectrum
microfacet_ggx_glass_albedo_scaling(KernelGlobals kg,
                                    ccl_private const ShaderData *sd,
                                    ccl_private const MicrofacetBsdf *bsdf,
                                    const Spectrum Fss)
{
  float mu = dot(sd->I, bsdf->N);
  float rough = sqrtf(sqrtf(bsdf->alpha_x * bsdf->alpha_y));
  float E = microfacet_ggx_glass_E(kg, mu, rough, bsdf->ior);

  /* Close enough for glass, coloring here is unphysical anyways and it's unclear how to
   * approximate it better. */
  Spectrum Fms = Fss;

  return one_spectrum() + Fms * ((1.0f - E) / E);
  /* TODO: Ensure that increase in weight does not mess up glossy color, albedo etc. passes */
}

/* Currently no non-albedo-scaled version is implemented, could easily be added
 * but would still break compatibility with the old glass due to the microfacet Fresnel. */

ccl_device int bsdf_microfacet_multi_ggx_glass_setup(KernelGlobals kg,
                                                     ccl_private MicrofacetBsdf *bsdf,
                                                     ccl_private const ShaderData *sd,
                                                     const Spectrum color)
{
  bsdf->extra = NULL;

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = bsdf->alpha_x;

  bsdf->type = CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_ID;

  bsdf->weight *= microfacet_ggx_glass_albedo_scaling(kg, sd, bsdf, saturate(color));

  return SD_BSDF | SD_BSDF_HAS_EVAL;
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

  Spectrum Fss = schlick_fresnel_Fss(bsdf->extra->cspec0);
  bsdf->weight *= microfacet_ggx_glass_albedo_scaling(kg, sd, bsdf, Fss);

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device_inline Spectrum glass_fresnel(ccl_private const MicrofacetBsdf *bsdf,
                                         float cosTheta,
                                         float *reflect_pdf)
{
  if (bsdf->T.x > 0.0f) {
    Spectrum F = fresnel_dielectric_thin_film(cosTheta, bsdf->ior, bsdf->T.y, bsdf->T.x);
    *reflect_pdf = average(F);
    return F;
  }
  else {
    *reflect_pdf = fresnel_dielectric_cos(cosTheta, bsdf->ior);
    return make_spectrum(*reflect_pdf);
  }
}

ccl_device Spectrum bsdf_microfacet_ggx_glass_eval_reflect(ccl_private const MicrofacetBsdf *bsdf,
                                                           const float3 N,
                                                           const float3 I,
                                                           const float3 omega_in,
                                                           ccl_private float *pdf,
                                                           const float alpha_x,
                                                           const float alpha_y,
                                                           const float cosNO,
                                                           const float cosNI)
{
  if (cosNI <= 0 || cosNO <= 0) {
    *pdf = 0.0f;
    return zero_spectrum();
  }

  float alpha2 = alpha_x * alpha_y;
  float3 m = normalize(omega_in + I);
  float D = microfacet_ggx_D(dot(N, m), alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  float common = D * 0.25f / cosNO;

  float reflect_pdf;
  Spectrum F = glass_fresnel(bsdf, dot(m, I), &reflect_pdf);

  Spectrum eval = F * common / (1 + lambdaO + lambdaI);
  *pdf = reflect_pdf * common / (1 + lambdaO);

  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    eval *= reflection_color(bsdf, omega_in, m);
  }
  return eval;
}

ccl_device Spectrum bsdf_microfacet_ggx_glass_eval_transmit(ccl_private const MicrofacetBsdf *bsdf,
                                                            const float3 N,
                                                            const float3 I,
                                                            const float3 omega_in,
                                                            ccl_private float *pdf,
                                                            const float alpha_x,
                                                            const float alpha_y,
                                                            const float cosNO,
                                                            const float cosNI)
{

  if (cosNO <= 0 || cosNI >= 0) {
    *pdf = 0.0f;
    return zero_spectrum();
  }

  float eta = bsdf->ior;
  float3 ht = -(eta * omega_in + I);
  float3 m = normalize(ht);
  float cosMO = dot(m, I);
  float cosMI = dot(m, omega_in);

  float reflect_pdf;
  Spectrum F = glass_fresnel(bsdf, cosMO, &reflect_pdf);

  if (F == one_spectrum()) {
    /* TIR */
    *pdf = 0.0f;
    return zero_spectrum();
  }

  float alpha2 = alpha_x * alpha_y;
  float D = microfacet_ggx_D(dot(N, m), alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  float Ht2 = dot(ht, ht);

  float common = fabsf(cosMI * cosMO) * D * sqr(eta) / (cosNO * Ht2);
  Spectrum eval = (one_spectrum() - F) * common / (1 + lambdaO + lambdaI);
  *pdf = (1.0f - reflect_pdf) / (1 + lambdaO);

  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    eval *= bsdf->extra->color;
  }
  return eval;
}

ccl_device Spectrum bsdf_microfacet_ggx_glass_eval(ccl_private const ShaderClosure *sc,
                                                   const float3 I,
                                                   const float3 omega_in,
                                                   ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  const float alpha_x = bsdf->alpha_x;
  const float alpha_y = bsdf->alpha_y;
  const float3 N = bsdf->N;
  const float cosNO = dot(N, I);
  const float cosNI = dot(N, omega_in);

  if (alpha_x * alpha_y <= 1e-7f) {
    *pdf = 0.0f;
    return zero_spectrum();
  }

  return (cosNI < 0.0f) ? bsdf_microfacet_ggx_glass_eval_transmit(
                              bsdf, N, I, omega_in, pdf, alpha_x, alpha_y, cosNO, cosNI) :
                          bsdf_microfacet_ggx_glass_eval_reflect(
                              bsdf, N, I, omega_in, pdf, alpha_x, alpha_y, cosNO, cosNI);
}

ccl_device int bsdf_microfacet_ggx_glass_sample(ccl_private const ShaderClosure *sc,
                                                float3 Ng,
                                                float3 I,
                                                float randu,
                                                float randv,
                                                ccl_private Spectrum *eval,
                                                ccl_private float3 *omega_in,
                                                ccl_private float *pdf,
                                                ccl_private float2 *sampled_roughness,
                                                ccl_private float *eta)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  float3 N = bsdf->N;
  int label;

  *sampled_roughness = make_float2(alpha_x, alpha_y);
  *eta = bsdf->ior;  // TODO: Do we need to invert in case of refraction?

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
  bool inside; /* Will never be inside, we already checked cosMO */
  fresnel_dielectric(bsdf->ior, m, I, &R, &T, &inside);

  /* TODO: Avoid recomputing Fresnel here if there's no thin-film?
   * fresnel_dielectric already does it otherwise. */
  float reflect_pdf;
  Spectrum F = glass_fresnel(bsdf, cosMO, &reflect_pdf);

  // TODO: Somehow get a properly stratified value here, this causes considerable noise
  float randw = hash_float2_to_float(make_float2(randu, randv));
  bool do_reflect = randw < reflect_pdf;

  float alpha2 = alpha_x * alpha_y;
  if (alpha2 <= 1e-7f) {
    /* Specular case, just return some high number for MIS */
    if (do_reflect) {
      *pdf = reflect_pdf * 1e6f;
      *eval = F * 1e6f;
      *omega_in = R;
    }
    else {
      *pdf = (1.0f - reflect_pdf) * 1e6f;
      *eval = (one_spectrum() - F) * 1e6f;
      *omega_in = T;
    }

    if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
      *eval *= do_reflect ? reflection_color(bsdf, *omega_in, m) : bsdf->extra->color;
    }

    return LABEL_SINGULAR | (do_reflect ? LABEL_REFLECT : LABEL_TRANSMIT);
  }

  /* Common microfacet model terms. */
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
    *pdf = reflect_pdf;
    *eval = F;

    /* Compute microfacet common term. */
    common = D * 0.25f / cosNO;
  }
  else {
    cosNI = dot(N, T);
    if (cosNI >= 0.0f || dot(Ng, T) >= 0.0f) {
      *pdf = 0.0f;
      return LABEL_NONE;
    }

    label = LABEL_TRANSMIT | LABEL_GLOSSY;
    *omega_in = T;
    *pdf = 1.0f - reflect_pdf;
    *eval = one_spectrum() - F;

    /* Compute microfacet common term. */
    float cosMI = dot(m, *omega_in);
    float Ht2 = sqr(bsdf->ior * cosMI + cosMO);
    common = D * fabsf(cosMI * cosMO) * sqr(bsdf->ior) / (cosNO * Ht2);
  }

  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);
  *eval *= common / (1 + lambdaO + lambdaI);
  *pdf *= common / (1 + lambdaO);

  if (bsdf->type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_FRESNEL_ID) {
    *eval *= do_reflect ? reflection_color(bsdf, *omega_in, m) : bsdf->extra->color;
  }
  return label;
}

CCL_NAMESPACE_END
