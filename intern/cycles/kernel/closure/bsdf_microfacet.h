/* SPDX-License-Identifier: BSD-3-Clause
 *
 * Adapted from Open Shading Language
 * Copyright (c) 2009-2010 Sony Pictures Imageworks Inc., et al.
 * All Rights Reserved.
 *
 * Modifications Copyright 2011-2022 Blender Foundation. */

#pragma once

#include "kernel/closure/bsdf_util.h"

CCL_NAMESPACE_BEGIN

typedef struct MicrofacetExtra {
  float3 color, cspec0;
  float3 fresnel_color;
  float clearcoat;
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
  float3 F = make_float3(1.0f, 1.0f, 1.0f);
  bool use_fresnel = (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_FRESNEL_ID ||
                      bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID);
  if (use_fresnel) {
    float F0 = fresnel_dielectric_cos(1.0f, bsdf->ior);

    F = interpolate_fresnel_color(L, H, bsdf->ior, F0, bsdf->extra->cspec0);
  }

  return F;
}

ccl_device_forceinline void bsdf_microfacet_fresnel_color(ccl_private const ShaderData *sd,
                                                          ccl_private MicrofacetBsdf *bsdf)
{
  kernel_assert(CLOSURE_IS_BSDF_MICROFACET_FRESNEL(bsdf->type));

  float F0 = fresnel_dielectric_cos(1.0f, bsdf->ior);
  bsdf->extra->fresnel_color = interpolate_fresnel_color(
      sd->I, bsdf->N, bsdf->ior, F0, bsdf->extra->cspec0);

  if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
    bsdf->extra->fresnel_color *= 0.25f * bsdf->extra->clearcoat;
  }

  bsdf->sample_weight *= average(bsdf->extra->fresnel_color);
}

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
 * https://jcgt.org/published/0003/02/03/
 *
 * Anisotropy is only supported for reflection currently, but adding it for
 * transmission is just a matter of copying code from reflection if needed. */

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

ccl_device int bsdf_microfacet_ggx_clearcoat_setup(ccl_private MicrofacetBsdf *bsdf,
                                                   ccl_private const ShaderData *sd)
{
  bsdf->extra->cspec0 = saturate(bsdf->extra->cspec0);

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
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;

  if (m_refractive || alpha_x * alpha_y <= 1e-7f) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  /* Warning: Cycles' naming is misleading here!
   * I is the incoming direction in a ray-tracing sense, but in the shading context,
   * it is actually the outgoing direction since it points towards the camera.
   * Therefore, in the BSDF code, I is referred to as O and omega_in is referred to as I
   * in order to be consistent with papers.
   */
  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);

  if (cosNI > 0 && cosNO > 0) {
    /* get half vector */
    float3 m = normalize(omega_in + I);
    float alpha2 = alpha_x * alpha_y;
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

    /* The full BSDF is (see e.g. eq. 20 in Walter et al. 2017):
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

    float common = D * 0.25f / cosNO;

    float3 F = reflection_color(bsdf, omega_in, m);
    if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
      F *= 0.25f * bsdf->extra->clearcoat;
    }

    float3 out = F * common / (1 + lambdaO + lambdaI);
    *pdf = common / (1 + lambdaO);

    return out;
  }

  return make_float3(0.0f, 0.0f, 0.0f);
}

ccl_device float3 bsdf_microfacet_ggx_eval_transmit(ccl_private const ShaderClosure *sc,
                                                    const float3 I,
                                                    const float3 omega_in,
                                                    ccl_private float *pdf)
{
  ccl_private const MicrofacetBsdf *bsdf = (ccl_private const MicrofacetBsdf *)sc;
  float alpha_x = bsdf->alpha_x;
  float alpha_y = bsdf->alpha_y;
  float m_eta = bsdf->ior;
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;

  if (!m_refractive || alpha_x * alpha_y <= 1e-7f) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f);
  }

  float cosNO = dot(N, I);
  float cosNI = dot(N, omega_in);

  if (cosNO <= 0 || cosNI >= 0) {
    *pdf = 0.0f;
    return make_float3(0.0f, 0.0f, 0.0f); /* vectors on same side -- not possible */
  }
  /* compute half-vector of the refraction (eq. 16) */
  float3 ht = -(m_eta * omega_in + I);
  float3 Ht = normalize(ht);
  float cosHO = dot(Ht, I);
  float cosHI = dot(Ht, omega_in);

  float D, lambdaO, lambdaI;

  /* eq. 33: first we calculate D(m) with m=Ht: */
  float alpha2 = alpha_x * alpha_y;
  D = microfacet_ggx_D(dot(N, Ht), alpha2);

  lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
  lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

  /* probability */
  float Ht2 = dot(ht, ht);

  /* eq. 2 in distribution of visible normals sampling
   * pm = Dw = G1o * dot(m, I) * D / dot(N, I); */

  /* out = fabsf(cosHI * cosHO) * (m_eta * m_eta) * G * D / (cosNO * Ht2)
   * pdf = pm * (m_eta * m_eta) * fabsf(cosHI) / Ht2 */
  float common = fabsf(cosHI * cosHO) * D * sqr(m_eta) / (cosNO * Ht2);
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
  bool m_refractive = bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID;
  float3 N = bsdf->N;
  int label;

  float cosNO = dot(N, I);
  if (cosNO > 0) {
    float3 X, Y, Z = N;

    if (alpha_x == alpha_y)
      make_orthonormals(Z, &X, &Y);
    else
      make_orthonormals_tangent(Z, bsdf->T, &X, &Y);

    /* importance sampling with distribution of visible normals. vectors are
     * transformed to local space before and after */
    float3 local_O = make_float3(dot(X, I), dot(Y, I), cosNO);
    float3 local_m = microfacet_ggx_sample_vndf(local_O, alpha_x, alpha_y, randu, randv);

    float3 m = X * local_m.x + Y * local_m.y + Z * local_m.z;
    float cosThetaM = local_m.z;

    /* reflection or refraction? */
    if (!m_refractive) {
      float cosMO = dot(m, I);
      label = LABEL_REFLECT | LABEL_GLOSSY;

      if (cosMO > 0) {
        /* eq. 39 - compute actual reflected direction */
        *omega_in = 2 * cosMO * m - I;

        if (dot(Ng, *omega_in) > 0) {
          if (alpha_x * alpha_y <= 1e-7f) {
            /* Specular case, just return some high number for MIS */
            *pdf = 1e6f;
            *eval = make_float3(1e6f, 1e6f, 1e6f);

            bool use_fresnel = (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_FRESNEL_ID ||
                                bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID);

            /* if fresnel is used, calculate the color with reflection_color(...) */
            if (use_fresnel) {
              *eval *= reflection_color(bsdf, *omega_in, m);
            }

            label = LABEL_REFLECT | LABEL_SINGULAR;
          }
          else {
            /* Evaluate microfacet model. */
            float alpha2 = alpha_x * alpha_y;
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

              float3 local_I = make_float3(
                  dot(X, *omega_in), dot(Y, *omega_in), dot(N, *omega_in));
              lambdaO = microfacet_ggx_lambda_aniso(local_O, alpha_x, alpha_y);
              lambdaI = microfacet_ggx_lambda_aniso(local_I, alpha_x, alpha_y);
            }

            /* See bsdf_microfacet_ggx_eval_reflect for derivation. */
            float common = D * 0.25f / cosNO;
            *pdf = common / (1 + lambdaO);

            float3 F = reflection_color(bsdf, *omega_in, m);

            *eval = common * F / (1 + lambdaO + lambdaI);
          }

          if (bsdf->type == CLOSURE_BSDF_MICROFACET_GGX_CLEARCOAT_ID) {
            *eval *= 0.25f * bsdf->extra->clearcoat;
          }

#ifdef __RAY_DIFFERENTIALS__
          *domega_in_dx = (2 * dot(m, dIdx)) * m - dIdx;
          *domega_in_dy = (2 * dot(m, dIdy)) * m - dIdy;
#endif
        }
        else {
          *eval = make_float3(0.0f, 0.0f, 0.0f);
          *pdf = 0.0f;
        }
      }
    }
    else {
      label = LABEL_TRANSMIT | LABEL_GLOSSY;

      /* CAUTION: the i and o variables are inverted relative to the paper
       * eq. 39 - compute actual refractive direction */
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

      if (!inside && fresnel != 1.0f) {

        *omega_in = T;
#ifdef __RAY_DIFFERENTIALS__
        *domega_in_dx = dTdx;
        *domega_in_dy = dTdy;
#endif

        if (alpha_x * alpha_y <= 1e-7f || fabsf(m_eta - 1.0f) < 1e-4f) {
          /* some high number for MIS */
          *pdf = 1e6f;
          *eval = make_float3(1e6f, 1e6f, 1e6f);
          label = LABEL_TRANSMIT | LABEL_SINGULAR;
        }
        else {
          /* eq. 33 */
          float alpha2 = alpha_x * alpha_y;
          float D = microfacet_ggx_D(cosThetaM, alpha2);

          /* eval BRDF*cosNI */
          float cosNI = dot(N, *omega_in);

          /* eq. 34: now calculate G1(i,m) */
          float lambdaO = microfacet_ggx_lambda(cosNO, alpha2);
          float lambdaI = microfacet_ggx_lambda(cosNI, alpha2);

          /* eq. 21 */
          float cosHI = dot(m, *omega_in);
          float cosHO = dot(m, I);
          float Ht2 = m_eta * cosHI + cosHO;
          Ht2 *= Ht2;

          /* see eval function for derivation */
          float common = D * (m_eta * m_eta) / (cosNO * Ht2);
          float out = fabsf(cosHI * cosHO) * common / (1 + lambdaO + lambdaI);
          *pdf = cosHO * fabsf(cosHI) * common / (1 + lambdaO);

          *eval = make_float3(out, out, out);
        }
      }
      else {
        *eval = make_float3(0.0f, 0.0f, 0.0f);
        *pdf = 0.0f;
      }
    }
  }
  else {
    label = (m_refractive) ? LABEL_TRANSMIT | LABEL_GLOSSY : LABEL_REFLECT | LABEL_GLOSSY;
  }
  return label;
}

CCL_NAMESPACE_END
