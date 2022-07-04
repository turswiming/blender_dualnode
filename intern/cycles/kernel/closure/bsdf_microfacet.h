/* SPDX-License-Identifier: BSD-3-Clause
 *
 * Adapted from Open Shading Language
 * Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
 * All Rights Reserved.
 *
 * Modifications Copyright 2011-2022 Blender Foundation. */

#pragma once

#include "kernel/closure/bsdf_microfacet_util.h"
#include "kernel/closure/bsdf_util.h"

CCL_NAMESPACE_BEGIN

typedef struct MicrofacetExtra {
  float3 color, cspec0;
  float3 fresnel_color;
} MicrofacetExtra;

typedef struct MicrofacetBsdf {
  SHADER_CLOSURE_BASE;

  float alpha_x, alpha_y, ior;
  ccl_private MicrofacetExtra *extra;
  float3 T;
} MicrofacetBsdf;

static_assert(sizeof(ShaderClosure) >= sizeof(MicrofacetBsdf), "MicrofacetBsdf is too large!");

/* Calculate the reflection color
 *
 * If fresnel is used, the color is an interpolation of the F0 color and white
 * with respect to the fresnel
 *
 * Else it is simply white
 */
ccl_device_forceinline float3 reflection_color(ccl_private const MicrofacetBsdf *bsdf,
                                               float3 L,
                                               float3 H)
{
  if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_FRESNEL_ID) {
    return interpolate_fresnel_color(L, H, bsdf->ior, bsdf->extra->cspec0);
  }
  else if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
    return interpolate_fresnel_color(L, H, bsdf->ior, make_float3(0.04f, 0.04f, 0.04f));
  }
  else {
    return one_float3();
  }
}

ccl_device_forceinline void bsdf_microfacet_fresnel_color(ccl_private const ShaderData *sd,
                                                          ccl_private MicrofacetBsdf *bsdf)
{
  float3 average_fresnel = reflection_color(bsdf, sd->I, bsdf->N);
  bsdf->sample_weight *= average(average_fresnel);

  if (bsdf->extra) {
    bsdf->extra->fresnel_color = average_fresnel;
  }
}

ccl_device_inline float3 microfacet_ggx_albedo_scaling(ccl_private const MicrofacetBsdf *bsdf,
                                                       ccl_private const ShaderData *sd,
                                                       const float3 Fss)
{
  float mu = dot(sd->I, bsdf->N);
  float rough = sqrtf(sqrtf(bsdf->alpha_x * bsdf->alpha_y));
  float E = microfacet_ggx_E(mu, rough);

  float E_avg = microfacet_ggx_E_avg(rough);
  /* Fms here is based on the appendix of
   * https://blog.selfshadow.com/publications/s2017-shading-course/imageworks/s2017_pbs_imageworks_slides_v2.pdf,
   * with one Fss cancelled out since this is just a multiplier on top of
   * the single-scattering BSDF, which already contains one bounce of Fresnel. */
  float3 Fms = Fss * E_avg / (one_float3() - Fss * (1.0f - E_avg));

  return one_float3() + Fms * ((1.0f - E) / E);
  /* TODO: Ensure that increase in weight does not mess up glossy color, albedo etc. passes */
}

ccl_device int bsdf_microfacet_ggx_setup(ccl_private MicrofacetBsdf *bsdf)
{
  bsdf->extra = NULL;

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = saturatef(bsdf->alpha_y);

  bsdf->type = CLOSURE_BSDF_MICROFACET_GGX_ID;

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

/* Required to maintain OSL interface. */
ccl_device int bsdf_microfacet_ggx_isotropic_setup(ccl_private MicrofacetBsdf *bsdf)
{
  bsdf->alpha_y = bsdf->alpha_x;

  return bsdf_microfacet_ggx_setup(bsdf);
}

ccl_device int bsdf_microfacet_multi_ggx_setup(ccl_private MicrofacetBsdf *bsdf,
                                               ccl_private const ShaderData *sd,
                                               const float3 color)
{
  bsdf->weight *= microfacet_ggx_albedo_scaling(bsdf, sd, saturate(color));
  return bsdf_microfacet_ggx_setup(bsdf);
}

ccl_device int bsdf_microfacet_ggx_fresnel_setup(ccl_private MicrofacetBsdf *bsdf,
                                                 ccl_private const ShaderData *sd)
{
  bsdf->extra->cspec0 = saturate(bsdf->extra->cspec0);

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = saturatef(bsdf->alpha_y);

  bsdf->type = CLOSURE_BSDF_MICROFACET_GGX_FRESNEL_ID;

  bsdf_microfacet_fresnel_color(sd, bsdf);

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device int bsdf_microfacet_multi_ggx_fresnel_setup(ccl_private MicrofacetBsdf *bsdf,
                                                       ccl_private const ShaderData *sd)
{
  float3 Fss = schlick_fresnel_Fss(bsdf->extra->cspec0);
  bsdf->weight *= microfacet_ggx_albedo_scaling(bsdf, sd, Fss);
  return bsdf_microfacet_ggx_fresnel_setup(bsdf, sd);
}

ccl_device int bsdf_microfacet_ggx_clearcoat_setup(ccl_private MicrofacetBsdf *bsdf,
                                                   ccl_private const ShaderData *sd)
{
  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = bsdf->alpha_x;

  bsdf->type = CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID;

  bsdf_microfacet_fresnel_color(sd, bsdf);

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device int bsdf_microfacet_ggx_refraction_setup(ccl_private MicrofacetBsdf *bsdf)
{
  bsdf->extra = NULL;

  bsdf->alpha_x = saturatef(bsdf->alpha_x);
  bsdf->alpha_y = bsdf->alpha_x;

  bsdf->type = CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;

  return SD_BSDF | SD_BSDF_HAS_EVAL;
}

ccl_device void bsdf_microfacet_ggx_blur(ccl_private ShaderClosure *sc, float roughness)
{
  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)sc;

  bsdf->alpha_x = fmaxf(roughness, bsdf->alpha_x);
  bsdf->alpha_y = fmaxf(roughness, bsdf->alpha_y);
}

ccl_device float3 bsdf_microfacet_ggx_eval_reflect(ccl_private const ShaderClosure *sc,
                                                   const float3 I,
                                                   const float3 omega_in,
                                                   ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  float alpha2 = alpha_x * alpha_y;
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;

  if (m_refractive || alpha2 <= 1e-7f) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  /* Warning: Cycles' naming is misleading here!
   * I is the incoming direction in a ray-tracing sense, but in the shading context,
   * it is actually the outgoing direction since it points towards the camera.
   * Therefore, in the BSDF code, I is referred to as O and omega_in is referred to as I
   * in order to be consistent with papers.
   */

  /* Ensure that both direction are in the upper hemisphere */
  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);
  if (cosNI <= 0 || cosNO <= 0) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  /* Compute half vector */
  float3 m = normalize(omega_in + I);

  float D, lambdaO, lambdaI;
  if (alpha_x == alpha_y) {
    /* Isotropic case */

    if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
      /* use GTR1 for clearcoat */
      D = microfacet_GTR1_D(dot(N, m), alpha2);

      /* the alpha value for clearcoat is a fixed 0.25 => alpha2 = 0.25 * 0.25 */
      alpha2 = 0.0625f;
    }
    else {
      /* use GTR2 otherwise */
      D = microfacet_ggx_D(dot(N, m), alpha2);
    }

    lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
    lambdaI = microfacet_ggx_lambda(cosNI, alpha2);
  }
  else {
    /* Anisotropic case */
    float3 X, Y, Z = N;
    make_orthonormals_tangent(Z, bsdf->T, &X, &Y);

    /* Transform vectors into local coordinate space */
    float3 local_m = make_float3(dot(X, m), dot(Y, m), dot(Z, m));
    float3 local_O = make_float3(dot(X, I), dot(Y, I), cosNO);
    float3 local_I = make_float3(dot(X, omega_in), dot(Y, omega_in), cosNI);

    D = microfacet_ggx_D_aniso(local_m, alpha_x, alpha_y);
    lambdaO = microfacet_ggx_lambda_aniso(local_O, alpha_x, alpha_y);
    lambdaI = microfacet_ggx_lambda_aniso(local_I, alpha_x, alpha_y);
  }

  /* The full BSDF is (see e.g. eq. 20 in Walter et al. 2007):
   * f(i, o) = F(i, m) * G(i, o) * D(m) / (4*cosNI*cosNO).
   *
   * Here, F is the fresnel reflection term, G is the masking-shadowing term,
   * D is the microfacet distribution and cosNI/cosNO are cosines of angles.
   *
   * For G, this implementation uses the non-separable form of the Smith
   * masking-shadowing term, so G is defined in terms of a function Lambda:
   * G(i, o) = 1 / (1 + Lambda(i) + Lambda(o)).
   *
   * In Cycles, BSDF evaluation actually returns f(i, o)*cosNI, so one term
   * in the BSDFs denominator cancels out.
   *
   * The PDF of VNDF sampling is D(m) * G1(o) / (4*cosNO), where G1(o) is
   * 1 / (1 + Lambda(o)).
   */

  /* Evaluate BSDF */
  float common = D * 0.25f / cosNO;
  float3 F = reflection_color(bsdf, omega_in, m);
  float3 out = F * common / (1 + lambdaO + lambdaI);
  *pdf = common / (1 + lambdaO);

  return out;
}

ccl_device float3 bsdf_microfacet_ggx_eval_transmit(ccl_private const ShaderClosure *sc,
                                                    const float3 I,
                                                    const float3 omega_in,
                                                    ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  float alpha2 = alpha_x * alpha_y;
  float m_eta = bsdf->ior;
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;

  if (!m_refractive || alpha2 <= 1e-7f) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  /* Ensure that both directions are in the expected hemispheres. */
  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);
  if (cosNO <= 0 || cosNI >= 0) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  /* Compute half vector */
  float3 ht = -(m_eta * omega_in + I);
  float3 m = normalize(ht);
  float cosMO = dot(m, I);
  float cosMI = dot(m, omega_in);

  /* Evaluate microfacet model */
  float D = microfacet_ggx_D(dot(N, m), alpha2);
  float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  /* Evaluate BSDF */
  float Ht2 = dot(ht, ht);
  float common = fabsf(cosMI * cosMO) * D * sqr(m_eta) / (cosNO * Ht2);
  float out = common / (1 + lambdaO + lambdaI);
  *pdf = common / (1 + lambdaO);

  return make_float3(out, out, out);
}

ccl_device int bsdf_microfacet_ggx_sample(ccl_private const ShaderClosure *sc,
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
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  float alpha2 = alpha_x * alpha_y;
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;

  /* Ensure that the view direction is in the upper hemisphere. */
  float cosNO = dot(N, I);
  if (cosNO <= 0) {
    *pdf = 0.0f;
    return LABEL_NONE;
  }

  /* Form local coordinate frame */
  float3 X, Y, Z = N;
  if (alpha_x == alpha_y)
    make_orthonormals(Z, &X, &Y);
  else
    make_orthonormals_tangent(Z, bsdf->T, &X, &Y);

  /* Sample distribution of visible normals to find the microfacet normal.
   * Sampling happens in the local frame. */
  float3 local_O = make_float3(dot(X, I), dot(Y, I), cosNO);
  float3 local_m = microfacet_ggx_sample_vndf(local_O, alpha_x, alpha_y, randu, randv);
  float3 m = X * local_m.x + Y * local_m.y + Z * local_m.z;
  float cosThetaM = local_m.z;

  /* Reflection or refraction? */
  if (!m_refractive) {
    /* Compute reflected direction and ensure that it is in the upper hemisphere.
     * Also check if the microfacet is masked (in that case, we'd hit it from the backside). */
    float cosMO = dot(m, I);
    *omega_in = 2 * cosMO * m - I;
    if (cosMO <= 0 || dot(Ng, *omega_in) <= 0) {
      *pdf = 0.0f;
      return LABEL_NONE;
    }

    float3 F = reflection_color(bsdf, *omega_in, m);
    if (alpha2 <= 1e-7f) {
      /* Specular case, just return some high number for MIS */
      *pdf = 1e6f;
      *eval = make_float3(1e6f, 1e6f, 1e6f) * F;
      return LABEL_REFLECT | LABEL_SINGULAR;
    }

    /* Evaluate microfacet model. */
    float D, lambdaO, lambdaI;
    if (alpha_x == alpha_y) {
      /* Isotropic case */

      if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
        /* use GTR1 for clearcoat */
        D = microfacet_GTR1_D(cosThetaM, alpha2);

        /* the alpha value for clearcoat is a fixed 0.25 => alpha2 = 0.25 * 0.25 */
        alpha2 = 0.0625f;

        /* recalculate lambdaO */
        lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
      }
      else {
        /* use GTR2 otherwise */
        D = microfacet_ggx_D(cosThetaM, alpha2);
      }

      float cosNI = dot(N, *omega_in);
      lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
      lambdaI = microfacet_ggx_lambda(cosNI, alpha2);
    }
    else {
      /* Anisotropic case */
      D = microfacet_ggx_D_aniso(local_m, alpha_x, alpha_y);

      float3 local_I = make_float3(dot(X, *omega_in), dot(Y, *omega_in), dot(N, *omega_in));
      lambdaO = microfacet_ggx_lambda_aniso(local_O, alpha_x, alpha_y);
      lambdaI = microfacet_ggx_lambda_aniso(local_I, alpha_x, alpha_y);
    }

    /* See bsdf_microfacet_ggx_eval_reflect for derivation. */
    float common = D * 0.25f / cosNO;
    *pdf = common / (1 + lambdaO);
    *eval = common * F / (1 + lambdaO + lambdaI);

#ifdef __RAY_DIFFERENTIALS__
    *domega_in_dx = (2 * dot(m, dIdx)) * m - dIdx;
    *domega_in_dy = (2 * dot(m, dIdy)) * m - dIdy;
#endif

    return LABEL_REFLECT | LABEL_GLOSSY;
  }
  else {
    /* Compute refracted direction */
    float3 R, T;
#ifdef __RAY_DIFFERENTIALS__
    float3 dRdx, dRdy, dTdx, dTdy;
#endif
    float m_eta = bsdf->ior, fresnel;
    bool inside;

    fresnel = fresnel_dielectric(m_eta,
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

    /* Ensure that the microfacet is nor masked and that we don't encounter TIR */
    if (inside || fresnel == 1.0f) {
      *pdf = 0.0f;
      return LABEL_NONE;
    }

    *omega_in = T;
#ifdef __RAY_DIFFERENTIALS__
    *domega_in_dx = dTdx;
    *domega_in_dy = dTdy;
#endif

    if (alpha2 <= 1e-7f || fabsf(m_eta - 1.0f) < 1e-4f) {
      /* some high number for MIS */
      *pdf = 1e6f;
      *eval = make_float3(1e6f, 1e6f, 1e6f);
      return LABEL_TRANSMIT | LABEL_SINGULAR;
    }

    /* Evaluate microfacet model */
    float D = microfacet_ggx_D(cosThetaM, alpha2);
    float cosNI = dot(N, *omega_in);
    float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
    float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

    /* Evaluate BSDF */
    float cosMI = dot(m, *omega_in);
    float cosMO = dot(m, I);
    float Ht2 = sqr(m_eta * cosMI + cosMO);
    float common = fabsf(cosMI * cosMO) * D * sqr(m_eta) / (cosNO * Ht2);
    float out = common / (1 + lambdaO + lambdaI);
    *pdf = common / (1 + lambdaO);

    *eval = make_float3(out, out, out);

    return LABEL_TRANSMIT | LABEL_GLOSSY;
  }
}

CCL_NAMESPACE_END
