/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

/* Functions to evaluate shaders. */

#pragma once

#include "kernel/closure/alloc.h"
#include "kernel/closure/bsdf.h"
#include "kernel/closure/bsdf_util.h"
#include "kernel/closure/emissive.h"

#include "kernel/integrator/guiding.h"

#ifdef __SVM__
#  include "kernel/svm/svm.h"
#endif
#ifdef __OSL__
#  include "kernel/osl/osl.h"
#endif

CCL_NAMESPACE_BEGIN

#define RIS_COSINE
#if OPENPGL_VERSION_MINOR >= 5
#  define RIS_INCOMMING_RADIANCE
#endif
/* Guiding */

#ifdef __PATH_GUIDING__

ccl_device float surface_shader_average_sample_weight_squared_roughness(
    ccl_private const ShaderData *sd)
{
  float avg_roughness = 0.0f;
  float sum_sample_weight = 0.0f;
  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (!CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
      continue;
    }
    avg_roughness += sc->sample_weight * bsdf_get_specular_roughness_squared(sc);
    sum_sample_weight += sc->sample_weight;
  }

  avg_roughness = avg_roughness > 0.f ? avg_roughness / sum_sample_weight : 0.f;
  return avg_roughness;
}

ccl_device_inline void surface_shader_prepare_guiding(KernelGlobals kg,
                                                      IntegratorState state,
                                                      ccl_private ShaderData *sd,
                                                      ccl_private const RNGState *rng_state)
{
  /* Have any BSDF to guide? */
  if (!(kernel_data.integrator.use_surface_guiding && (sd->flag & SD_BSDF_HAS_EVAL))) {
    state->guiding.use_surface_guiding = false;
    return;
  }

  const float surface_guiding_probability = kernel_data.integrator.surface_guiding_probability;
  const int guiding_directional_sampling_type =
      kernel_data.integrator.guiding_directional_sampling_type;
  const float guiding_roughness_threshold = kernel_data.integrator.guiding_roughness_threshold;
  float rand_bsdf_guiding = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);

  /* Compute proportion of diffuse BSDF and BSSRDFs .*/
  float diffuse_sampling_fraction = 0.0f;
  float bssrdf_sampling_fraction = 0.0f;
  float bsdf_bssrdf_sampling_sum = 0.0f;

  bool fully_opaque = true;

  for (int i = 0; i < sd->num_closure; i++) {
    ShaderClosure *sc = &sd->closure[i];
    if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
      const float sweight = sc->sample_weight;
      kernel_assert(sweight >= 0.0f);

      bsdf_bssrdf_sampling_sum += sweight;
      if (CLOSURE_IS_BSDF_DIFFUSE(sc->type) && sc->type < CLOSURE_BSDF_TRANSLUCENT_ID) {
        diffuse_sampling_fraction += sweight;
      }
      if (CLOSURE_IS_BSSRDF(sc->type)) {
        bssrdf_sampling_fraction += sweight;
      }

      if (CLOSURE_IS_BSDF_TRANSPARENT(sc->type) || CLOSURE_IS_BSDF_TRANSMISSION(sc->type)) {
        fully_opaque = false;
      }
    }
  }

  if (bsdf_bssrdf_sampling_sum > 0.0f) {
    diffuse_sampling_fraction /= bsdf_bssrdf_sampling_sum;
    bssrdf_sampling_fraction /= bsdf_bssrdf_sampling_sum;
  }

  float avg_roughness = surface_shader_average_sample_weight_squared_roughness(sd);
  avg_roughness = safe_sqrtf(avg_roughness);
  /* Init guiding (diffuse BSDFs only for now). */
  if (!fully_opaque || avg_roughness < guiding_roughness_threshold * guiding_roughness_threshold ||
      ((guiding_directional_sampling_type == GUIDING_DIRECTIONAL_SAMPLING_TYPE_PRODUCT) &&
       (diffuse_sampling_fraction <= 0.f)) ||
      !guiding_bsdf_init(kg, state, sd->P, sd->N, rand_bsdf_guiding)) {
    state->guiding.use_surface_guiding = false;
    state->guiding.surface_guiding_sampling_prob = 0.0f;
    return;
  }

  state->guiding.use_surface_guiding = true;
  if (kernel_data.integrator.guiding_directional_sampling_type ==
      GUIDING_DIRECTIONAL_SAMPLING_TYPE_PRODUCT) {
    state->guiding.surface_guiding_sampling_prob = surface_guiding_probability *
                                                   diffuse_sampling_fraction;
  }
  else if (kernel_data.integrator.guiding_directional_sampling_type ==
           GUIDING_DIRECTIONAL_SAMPLING_TYPE_RIS) {
    state->guiding.surface_guiding_sampling_prob = surface_guiding_probability;
  }
  else {  // GUIDING_DIRECTIONAL_SAMPLING_TYPE_ROUGHNESS
    state->guiding.surface_guiding_sampling_prob = surface_guiding_probability * avg_roughness;
  }
  state->guiding.bssrdf_sampling_prob = bssrdf_sampling_fraction;
  state->guiding.sample_surface_guiding_rand = rand_bsdf_guiding;

  kernel_assert(state->guiding.surface_guiding_sampling_prob > 0.0f &&
                state->guiding.surface_guiding_sampling_prob <= 1.0f);
}
#endif

ccl_device_inline void surface_shader_prepare_closures(KernelGlobals kg,
                                                       ConstIntegratorState state,
                                                       ccl_private ShaderData *sd,
                                                       const uint32_t path_flag)
{
  /* Filter out closures. */
  if (kernel_data.integrator.filter_closures) {
    if (kernel_data.integrator.filter_closures & FILTER_CLOSURE_EMISSION) {
      sd->closure_emission_background = zero_spectrum();
    }

    if (kernel_data.integrator.filter_closures & FILTER_CLOSURE_DIRECT_LIGHT) {
      sd->flag &= ~SD_BSDF_HAS_EVAL;
    }

    if (path_flag & PATH_RAY_CAMERA) {
      for (int i = 0; i < sd->num_closure; i++) {
        ccl_private ShaderClosure *sc = &sd->closure[i];

        if ((CLOSURE_IS_BSDF_DIFFUSE(sc->type) &&
             (kernel_data.integrator.filter_closures & FILTER_CLOSURE_DIFFUSE)) ||
            (CLOSURE_IS_BSDF_GLOSSY(sc->type) &&
             (kernel_data.integrator.filter_closures & FILTER_CLOSURE_GLOSSY)) ||
            (CLOSURE_IS_BSDF_TRANSMISSION(sc->type) &&
             (kernel_data.integrator.filter_closures & FILTER_CLOSURE_TRANSMISSION))) {
          sc->type = CLOSURE_NONE_ID;
          sc->sample_weight = 0.0f;
        }
        else if ((CLOSURE_IS_BSDF_TRANSPARENT(sc->type) &&
                  (kernel_data.integrator.filter_closures & FILTER_CLOSURE_TRANSPARENT))) {
          sc->type = CLOSURE_HOLDOUT_ID;
          sc->sample_weight = 0.0f;
          sd->flag |= SD_HOLDOUT;
        }
      }
    }
  }

  /* Defensive sampling.
   *
   * We can likely also do defensive sampling at deeper bounces, particularly
   * for cases like a perfect mirror but possibly also others. This will need
   * a good heuristic. */
  if (INTEGRATOR_STATE(state, path, bounce) + INTEGRATOR_STATE(state, path, transparent_bounce) ==
          0 &&
      sd->num_closure > 1) {
    float sum = 0.0f;

    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private ShaderClosure *sc = &sd->closure[i];
      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
        sum += sc->sample_weight;
      }
    }

    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private ShaderClosure *sc = &sd->closure[i];
      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
        sc->sample_weight = max(sc->sample_weight, 0.125f * sum);
      }
    }
  }

  /* Filter glossy.
   *
   * Blurring of bsdf after bounces, for rays that have a small likelihood
   * of following this particular path (diffuse, rough glossy) */
  if (kernel_data.integrator.filter_glossy != FLT_MAX
#ifdef __MNEE__
      && !(INTEGRATOR_STATE(state, path, mnee) & PATH_MNEE_VALID)
#endif
  ) {
    float blur_pdf = kernel_data.integrator.filter_glossy *
                     INTEGRATOR_STATE(state, path, min_ray_pdf);

    if (blur_pdf < 1.0f) {
      float blur_roughness = sqrtf(1.0f - blur_pdf) * 0.5f;

      for (int i = 0; i < sd->num_closure; i++) {
        ccl_private ShaderClosure *sc = &sd->closure[i];
        if (CLOSURE_IS_BSDF(sc->type)) {
          bsdf_blur(kg, sc, blur_roughness);
        }
      }
    }
  }
}

/* BSDF */
#if 0
ccl_device_inline void surface_shader_validate_bsdf_sample(const KernelGlobals kg,
                                                           const ShaderClosure *sc,
                                                           const float3 wo,
                                                           const int org_label,
                                                           const float2 org_roughness,
                                                           const float org_eta)
{
  /* Validate the the bsdf_label and bsdf_roughness_eta functions
   * by estimating the values after a bsdf sample. */
  const int comp_label = bsdf_label(kg, sc, wo);
  kernel_assert(org_label == comp_label);

  float2 comp_roughness;
  float comp_eta;
  bsdf_roughness_eta(kg, sc, &comp_roughness, &comp_eta);
  kernel_assert(org_eta == comp_eta);
  kernel_assert(org_roughness.x == comp_roughness.x);
  kernel_assert(org_roughness.y == comp_roughness.y);
}
#endif

ccl_device_forceinline bool _surface_shader_exclude(ClosureType type, uint light_shader_flags)
{
  if (!(light_shader_flags & SHADER_EXCLUDE_ANY)) {
    return false;
  }
  if (light_shader_flags & SHADER_EXCLUDE_DIFFUSE) {
    if (CLOSURE_IS_BSDF_DIFFUSE(type)) {
      return true;
    }
  }
  if (light_shader_flags & SHADER_EXCLUDE_GLOSSY) {
    if (CLOSURE_IS_BSDF_GLOSSY(type)) {
      return true;
    }
  }
  if (light_shader_flags & SHADER_EXCLUDE_TRANSMIT) {
    if (CLOSURE_IS_BSDF_TRANSMISSION(type)) {
      return true;
    }
  }
  return false;
}

ccl_device_inline float _surface_shader_bsdf_eval_mis(KernelGlobals kg,
                                                      ccl_private ShaderData *sd,
                                                      const float3 wo,
                                                      ccl_private const ShaderClosure *skip_sc,
                                                      ccl_private BsdfEval *result_eval,
                                                      float sum_pdf,
                                                      float sum_sample_weight,
                                                      const uint light_shader_flags)
{
  /* This is the veach one-sample model with balance heuristic,
   * some PDF factors drop out when using balance heuristic weighting. */
  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (sc == skip_sc) {
      continue;
    }

    if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
      if (CLOSURE_IS_BSDF(sc->type)) {
        float bsdf_pdf = 0.0f;
        Spectrum eval = bsdf_eval(kg, sd, sc, wo, &bsdf_pdf);

        if (bsdf_pdf != 0.0f) {
          if (!_surface_shader_exclude(sc->type, light_shader_flags)) {
            bsdf_eval_accum(result_eval, sc->type, eval * sc->weight);
          }
          sum_pdf += bsdf_pdf * sc->sample_weight;
        }
      }

      sum_sample_weight += sc->sample_weight;
    }
  }

  return (sum_sample_weight > 0.0f) ? sum_pdf / sum_sample_weight : 0.0f;
}

ccl_device_inline float surface_shader_bsdf_eval_pdfs(const KernelGlobals kg,
                                                      ccl_private ShaderData *sd,
                                                      const float3 wo,
                                                      ccl_private BsdfEval *result_eval,
                                                      ccl_private float *pdfs,
                                                      const uint light_shader_flags)
{
  /* This is the veach one-sample model with balance heuristic, some pdf
   * factors drop out when using balance heuristic weighting. */
  float sum_pdf = 0.0f;
  float sum_sample_weight = 0.0f;
  bsdf_eval_init(result_eval, CLOSURE_NONE_ID, zero_spectrum());
  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
      if (CLOSURE_IS_BSDF(sc->type) && !_surface_shader_exclude(sc->type, light_shader_flags)) {
        float bsdf_pdf = 0.0f;
        Spectrum eval = bsdf_eval(kg, sd, sc, wo, &bsdf_pdf);
        kernel_assert(bsdf_pdf >= 0.0f);
        if (bsdf_pdf != 0.0f) {
          bsdf_eval_accum(result_eval, sc->type, eval * sc->weight);
          sum_pdf += bsdf_pdf * sc->sample_weight;
          kernel_assert(bsdf_pdf * sc->sample_weight >= 0.0f);
          pdfs[i] = bsdf_pdf * sc->sample_weight;
        }
        else {
          pdfs[i] = 0.0f;
        }
      }
      else {
        pdfs[i] = 0.0f;
      }

      sum_sample_weight += sc->sample_weight;
    }
    else {
      pdfs[i] = 0.0f;
    }
  }
  if (sum_pdf > 0.0f) {
    for (int i = 0; i < sd->num_closure; i++) {
      pdfs[i] /= sum_pdf;
    }
  }

  return (sum_sample_weight > 0.0f) ? sum_pdf / sum_sample_weight : 0.0f;
}

#ifndef __KERNEL_CUDA__
ccl_device
#else
ccl_device_inline
#endif
    float
    surface_shader_bsdf_eval(KernelGlobals kg,
                             IntegratorState state,
                             ccl_private ShaderData *sd,
                             const float3 wo,
                             ccl_private BsdfEval *bsdf_eval,
                             const uint light_shader_flags)
{
  bsdf_eval_init(bsdf_eval, CLOSURE_NONE_ID, zero_spectrum());

  float pdf = _surface_shader_bsdf_eval_mis(
      kg, sd, wo, NULL, bsdf_eval, 0.0f, 0.0f, light_shader_flags);

#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
  if (pdf > 0.f && state->guiding.use_surface_guiding) {
    const float guiding_sampling_prob = state->guiding.surface_guiding_sampling_prob;
    const float bssrdf_sampling_prob = state->guiding.bssrdf_sampling_prob;
    const float guide_pdf = guiding_bsdf_pdf(kg, state, wo);

    if (kernel_data.integrator.guiding_directional_sampling_type ==
        GUIDING_DIRECTIONAL_SAMPLING_TYPE_RIS) {
      pdf = (0.5f * guide_pdf * (1.0f - bssrdf_sampling_prob)) + 0.5f * pdf;
    }
    else {
      pdf = (guiding_sampling_prob * guide_pdf * (1.0f - bssrdf_sampling_prob)) +
            (1.0f - guiding_sampling_prob) * pdf;
    }
  }
#endif

  return pdf;
}

/* Randomly sample a BSSRDF or BSDF proportional to ShaderClosure.sample_weight. */
ccl_device_inline ccl_private const ShaderClosure *surface_shader_bsdf_bssrdf_pick(
    ccl_private const ShaderData *ccl_restrict sd, ccl_private float2 *rand_bsdf)
{
  int sampled = 0;

  if (sd->num_closure > 1) {
    /* Pick a BSDF or based on sample weights. */
    float sum = 0.0f;

    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private const ShaderClosure *sc = &sd->closure[i];

      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
        sum += sc->sample_weight;
      }
    }

    float r = (*rand_bsdf).x * sum;
    float partial_sum = 0.0f;

    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private const ShaderClosure *sc = &sd->closure[i];

      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
        float next_sum = partial_sum + sc->sample_weight;

        if (r < next_sum) {
          sampled = i;

          /* Rescale to reuse for direction sample, to better preserve stratification. */
          (*rand_bsdf).x = (r - partial_sum) / sc->sample_weight;
          break;
        }

        partial_sum = next_sum;
      }
    }
  }

  return &sd->closure[sampled];
}

/* Return weight for picked BSSRDF. */
ccl_device_inline Spectrum
surface_shader_bssrdf_sample_weight(ccl_private const ShaderData *ccl_restrict sd,
                                    ccl_private const ShaderClosure *ccl_restrict bssrdf_sc)
{
  Spectrum weight = bssrdf_sc->weight;

  if (sd->num_closure > 1) {
    float sum = 0.0f;
    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private const ShaderClosure *sc = &sd->closure[i];

      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
        sum += sc->sample_weight;
      }
    }
    weight *= sum / bssrdf_sc->sample_weight;
  }

  return weight;
}

#ifdef __PATH_GUIDING__
/* Sample direction for picked BSDF, and return evaluation and pdf for all
 * BSDFs combined using MIS. */

ccl_device int surface_shader_bsdf_guided_sample_closure_mis(KernelGlobals kg,
                                                             IntegratorState state,
                                                             ccl_private ShaderData *sd,
                                                             ccl_private const ShaderClosure *sc,
                                                             const float2 rand_bsdf,
                                                             ccl_private BsdfEval *bsdf_eval,
                                                             ccl_private float3 *wo,
                                                             ccl_private float *bsdf_pdf,
                                                             ccl_private float *unguided_bsdf_pdf,
                                                             ccl_private float2 *sampled_rougness,
                                                             ccl_private float *eta)
{
  /* BSSRDF should already have been handled elsewhere. */
  kernel_assert(CLOSURE_IS_BSDF(sc->type));

  const bool use_surface_guiding = state->guiding.use_surface_guiding;
  const float guiding_sampling_prob = state->guiding.surface_guiding_sampling_prob;
  const float bssrdf_sampling_prob = state->guiding.bssrdf_sampling_prob;

  /* Decide between sampling guiding distribution and BSDF. */
  bool sample_guiding = false;
  float rand_bsdf_guiding = state->guiding.sample_surface_guiding_rand;

  if (use_surface_guiding && rand_bsdf_guiding < guiding_sampling_prob) {
    sample_guiding = true;
    rand_bsdf_guiding /= guiding_sampling_prob;
  }
  else {
    rand_bsdf_guiding -= guiding_sampling_prob;
    rand_bsdf_guiding /= (1.0f - guiding_sampling_prob);
  }

  /* Initialize to zero. */
  int label = LABEL_NONE;
  Spectrum eval = zero_spectrum();
  bsdf_eval_init(bsdf_eval, CLOSURE_NONE_ID, eval);

  *unguided_bsdf_pdf = 0.0f;
  float guide_pdf = 0.0f;

  if (sample_guiding) {
    /* Sample guiding distribution. */
    guide_pdf = guiding_bsdf_sample(kg, state, rand_bsdf, wo);
    *bsdf_pdf = 0.0f;

    if (guide_pdf != 0.0f) {
      float unguided_bsdf_pdfs[MAX_CLOSURE];

      *unguided_bsdf_pdf = surface_shader_bsdf_eval_pdfs(
          kg, sd, *wo, bsdf_eval, unguided_bsdf_pdfs, 0);
      *bsdf_pdf = (guiding_sampling_prob * guide_pdf * (1.0f - bssrdf_sampling_prob)) +
                  ((1.0f - guiding_sampling_prob) * (*unguided_bsdf_pdf));
      float sum_pdfs = 0.0f;

      if (*unguided_bsdf_pdf > 0.0f) {
        int idx = -1;
        for (int i = 0; i < sd->num_closure; i++) {
          sum_pdfs += unguided_bsdf_pdfs[i];
          if (rand_bsdf_guiding <= sum_pdfs) {
            idx = i;
            break;
          }
        }

        kernel_assert(idx >= 0);
        /* Set the default idx to the last in the list.
         * in case of numerical problems and rand_bsdf_guiding is just >=1.0f and
         * the sum of all unguided_bsdf_pdfs is just < 1.0f. */
        idx = (rand_bsdf_guiding > sum_pdfs) ? sd->num_closure - 1 : idx;

        label = bsdf_label(kg, &sd->closure[idx], *wo);
      }
      else {
        *bsdf_pdf = 0.0f;
        *unguided_bsdf_pdf = 0.0f;
      }
    }

    kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);

    *sampled_rougness = make_float2(1.0f, 1.0f);
    *eta = 1.0f;
  }
  else {
    /* Sample BSDF. */
    *bsdf_pdf = 0.0f;
    label = bsdf_sample(
        kg, sd, sc, rand_bsdf.x, rand_bsdf.y, &eval, wo, unguided_bsdf_pdf, sampled_rougness, eta);
#  if 0
    if (*unguided_bsdf_pdf > 0.0f) {
      surface_shader_validate_bsdf_sample(kg, sc, *wo, label, sampled_roughness, eta);
    }
#  endif

    if (*unguided_bsdf_pdf != 0.0f) {
      bsdf_eval_init(bsdf_eval, sc->type, eval * sc->weight);

      kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);

      if (sd->num_closure > 1) {
        float sweight = sc->sample_weight;
        *unguided_bsdf_pdf = _surface_shader_bsdf_eval_mis(
            kg, sd, *wo, sc, bsdf_eval, (*unguided_bsdf_pdf) * sweight, sweight, 0);
        kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);
      }
      *bsdf_pdf = *unguided_bsdf_pdf;

      if (use_surface_guiding) {
        guide_pdf = guiding_bsdf_pdf(kg, state, *wo);
        *bsdf_pdf *= 1.0f - guiding_sampling_prob;
        *bsdf_pdf += guiding_sampling_prob * guide_pdf * (1.0f - bssrdf_sampling_prob);
      }
    }

    kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);
  }

  return label;
}

ccl_device int surface_shader_bsdf_guided_sample_closure_ris(KernelGlobals kg,
                                                             IntegratorState state,
                                                             ccl_private ShaderData *sd,
                                                             ccl_private const ShaderClosure *sc,
                                                             const float2 rand_bsdf,
                                                             ccl_private const RNGState *rng_state,
                                                             ccl_private BsdfEval *bsdf_eval,
                                                             ccl_private float3 *omega_in,
                                                             ccl_private float *bsdf_pdf,
                                                             ccl_private float *mis_pdf,
                                                             ccl_private float *unguided_bsdf_pdf,
                                                             ccl_private float2 *sampled_rougness,
                                                             ccl_private float *eta)
{
  /* BSSRDF should already have been handled elsewhere. */
  kernel_assert(CLOSURE_IS_BSDF(sc->type));

  const bool use_surface_guiding = state->guiding.use_surface_guiding;
  const float guiding_sampling_prob = state->guiding.surface_guiding_sampling_prob;
  const float bssrdf_sampling_prob = state->guiding.bssrdf_sampling_prob;

  /* Decide between sampling guiding distribution and BSDF. */
  float rand_bsdf_guiding = state->guiding.sample_surface_guiding_rand;

  /* Initialize to zero. */
  int label = LABEL_NONE;
  Spectrum eval = zero_spectrum();
  bsdf_eval_init(bsdf_eval, CLOSURE_NONE_ID, eval);

  *unguided_bsdf_pdf = 0.0f;
  float guide_pdf = 0.0f;

  if (use_surface_guiding && guiding_sampling_prob > 0.f) {
    /* Sample guiding distribution. */
    float2 rand_guiding_bsdf_ris[2];
    rand_guiding_bsdf_ris[0] = path_state_rng_2D(kg, rng_state, PRNG_SURFACE_RIS_GUIDING_0);
    rand_guiding_bsdf_ris[1] = path_state_rng_2D(kg, rng_state, PRNG_SURFACE_RIS_GUIDING_1);

    float pi_factor = 2.0f;

    float2 sampled_rougness_ris[2];
    float eta_ris[2];
    int label_ris[2];

    float3 omega_in_ris[2];
    float ris_pdfs[2] = {0.f, 0.f};
    float guide_pdfs[2] = {0.f, 0.f};
#  ifdef RIS_INCOMMING_RADIANCE
    float incomming_radiance_pdfs[2] = {0.f, 0.f};
#  endif
    float bsdf_pdfs[2] = {0.f, 0.f};
    float cosines[2] = {0.f, 0.f};
    BsdfEval bsdf_evals[2];
    float avg_bsdf_evals[2] = {0.f, 0.f};
    Spectrum evals[2] = {zero_spectrum(), zero_spectrum()};

    int ris_idx = 0;

    // RIS0 - sample BSDF
    label_ris[0] = bsdf_sample(kg,
                               sd,
                               sc,
                               rand_guiding_bsdf_ris[0].x,
                               rand_guiding_bsdf_ris[0].y,
                               &evals[0],
                               &omega_in_ris[0],
                               &bsdf_pdfs[0],
                               &sampled_rougness_ris[0],
                               &eta_ris[0]);

    bsdf_eval_init(&bsdf_evals[0], sc->type, evals[0] * sc->weight);
    if (bsdf_pdfs[0] > 0.f) {
      cosines[0] = max(0.01f, fabsf(dot(sd->N, omega_in_ris[0])));
      if (sd->num_closure > 1) {
        float sweight = sc->sample_weight;
        bsdf_pdfs[0] = _surface_shader_bsdf_eval_mis(
            kg, sd, omega_in_ris[0], sc, &bsdf_evals[0], (bsdf_pdfs[0]) * sweight, sweight, 0);
        kernel_assert(reduce_min(bsdf_eval_sum(&bsdf_evals[0])) >= 0.0f);
      }
      avg_bsdf_evals[0] = (bsdf_evals[0].sum[0] + bsdf_evals[0].sum[1] + bsdf_evals[0].sum[2]) /
                          3.0f;
      guide_pdfs[0] = guiding_bsdf_pdf(kg, state, omega_in_ris[0]);
      guide_pdfs[0] *= (1.0f - bssrdf_sampling_prob);
#  ifdef RIS_INCOMMING_RADIANCE
      incomming_radiance_pdfs[0] = guiding_surface_incomming_radiance_pdf(
          kg, state, omega_in_ris[0]);
#  endif
      bsdf_pdfs[0] = max(0.f, bsdf_pdfs[0]);
    }

    // RIS1 - sample guiding
    float unguided_bsdf_pdfs[MAX_CLOSURE];
    bsdf_eval_init(&bsdf_evals[1], CLOSURE_NONE_ID, eval);
    guide_pdfs[1] = guiding_bsdf_sample(kg, state, rand_guiding_bsdf_ris[1], &omega_in_ris[1]);
    guide_pdfs[1] *= (1.0f - bssrdf_sampling_prob);
#  ifdef RIS_INCOMMING_RADIANCE
    incomming_radiance_pdfs[1] = guiding_surface_incomming_radiance_pdf(
        kg, state, omega_in_ris[1]);
#  endif
    cosines[1] = max(0.01f, fabsf(dot(sd->N, omega_in_ris[1])));
    bsdf_pdfs[1] = surface_shader_bsdf_eval_pdfs(
        kg, sd, omega_in_ris[1], &bsdf_evals[1], unguided_bsdf_pdfs, 0);
    label_ris[1] = label_ris[0];
    avg_bsdf_evals[1] = (bsdf_evals[1].sum[0] + bsdf_evals[1].sum[1] + bsdf_evals[1].sum[2]) /
                        3.0f;
    bsdf_pdfs[1] = max(0.f, bsdf_pdfs[1]);

    int num_samples = 0;
    float sum_ris_pdfs = 0.f;
    if (avg_bsdf_evals[0] > 0.f && bsdf_pdfs[0] > 1e-10f && guide_pdfs[0] > 0.f) {
#  ifdef RIS_COSINE
#    ifdef RIS_INCOMMING_RADIANCE
      ris_pdfs[0] = (avg_bsdf_evals[0] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * incomming_radiance_pdfs[0])))) /
                    (0.5f * (bsdf_pdfs[0] + guide_pdfs[0]));
#    else
      ris_pdfs[0] = (avg_bsdf_evals[0] / cosines[0] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * guide_pdfs[0])))) /
                    (0.5f * (bsdf_pdfs[0] + guide_pdfs[0]));
#    endif
#  else
      ris_pdfs[0] = (avg_bsdf_evals[0] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * guide_pdfs[0])))) /
                    (0.5f * (bsdf_pdfs[0] + guide_pdfs[0]));
#  endif
      sum_ris_pdfs += ris_pdfs[0];
      num_samples++;
    }
    else {
      ris_pdfs[0] = 0.f;
    }
    assert(sum_ris_pdfs >= 0.f);

    if (avg_bsdf_evals[1] > 0.f && bsdf_pdfs[1] > 1e-10f && guide_pdfs[1] > 0.f) {
#  ifdef RIS_COSINE
#    ifdef RIS_INCOMMING_RADIANCE
      ris_pdfs[1] = (avg_bsdf_evals[1] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * incomming_radiance_pdfs[1])))) /
                    (0.5f * (bsdf_pdfs[1] + guide_pdfs[1]));
#    else
      ris_pdfs[1] = (avg_bsdf_evals[1] / cosines[1] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * guide_pdfs[1])))) /
                    (0.5f * (bsdf_pdfs[1] + guide_pdfs[1]));
#    endif
#  else
      ris_pdfs[1] = (avg_bsdf_evals[1] *
                     ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                       (guiding_sampling_prob * guide_pdfs[1])))) /
                    (0.5f * (bsdf_pdfs[1] + guide_pdfs[1]));
#  endif
      sum_ris_pdfs += ris_pdfs[1];
      num_samples++;
    }
    else {
      ris_pdfs[1] = 0.f;
    }
    assert(sum_ris_pdfs >= 0.f);

    if (num_samples == 0 || !(sum_ris_pdfs > 1e-10f)) {
      *bsdf_pdf = 0.0f;
      *mis_pdf = 0.0f;
      return label;
    }

    float rand_ris_select = rand_bsdf_guiding * sum_ris_pdfs;

    float sum_ris = 0.0f;
    for (int i = 0; i < 2; i++) {
      sum_ris += ris_pdfs[i];
      if (rand_ris_select <= sum_ris) {
        ris_idx = i;
        break;
      }
    }

    assert(sum_ris_pdfs >= 0.f);
    assert(ris_idx < 2);

#  ifdef RIS_COSINE
#    ifdef RIS_INCOMMING_RADIANCE
    guide_pdf = (avg_bsdf_evals[ris_idx] *
                 ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                   (guiding_sampling_prob * incomming_radiance_pdfs[ris_idx])))) *
                (float(2) / sum_ris_pdfs);
#    else
    guide_pdf = (avg_bsdf_evals[ris_idx] / cosines[ris_idx] *
                 ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                   (guiding_sampling_prob * guide_pdfs[ris_idx])))) *
                (float(2) / sum_ris_pdfs);
#    endif
#  else
    guide_pdf = (avg_bsdf_evals[ris_idx] *
                 ((((1.0f - guiding_sampling_prob) * (1.0f / (pi_factor * float(M_PI)))) +
                   (guiding_sampling_prob * guide_pdfs[ris_idx])))) *
                (float(2) / sum_ris_pdfs);
#  endif
    *unguided_bsdf_pdf = bsdf_pdfs[ris_idx];
    *mis_pdf = 0.5f * (bsdf_pdfs[ris_idx] + guide_pdfs[ris_idx]);
    *bsdf_pdf = guide_pdf;

    *omega_in = omega_in_ris[ris_idx];
    label = label_ris[ris_idx];

    *sampled_rougness = sampled_rougness_ris[ris_idx];
    *eta = eta_ris[ris_idx];
    *bsdf_eval = bsdf_evals[ris_idx];

    assert(std::isfinite(guide_pdf));
    assert(std::isfinite(*bsdf_pdf));

    if (!(*bsdf_pdf > 1e-10f)) {
      *bsdf_pdf = 0.0f;
      *mis_pdf = 0.0f;
      return label;
    }

    assert(*bsdf_pdf > 0.f);
    assert(*bsdf_pdf >= 1e-20f);
    assert(guide_pdf >= 0.f);

    /// select label sampled_roughness and eta
    if (ris_idx == 1) {
      // if (!(label | LABEL_SINGULAR) && !(label | LABEL_TRANSPARENT)) {
      float sum_pdfs = 0.0f;
      float rnd = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_RIS_GUIDING_2);

      if (bsdf_pdfs[1] > 0.0f) {
        int idx = -1;
        for (int i = 0; i < sd->num_closure; i++) {
          sum_pdfs += unguided_bsdf_pdfs[i];
          if (rnd <= sum_pdfs) {
            idx = i;
            break;
          }
        }
        // kernel_assert(idx >= 0);
        /* Set the default idx to the last in the list.
         * in case of numerical problems and rand_bsdf_guiding is just >=1.0f and
         * the sum of all unguided_bsdf_pdfs is just < 1.0f. */
        idx = (rnd > sum_pdfs) ? sd->num_closure - 1 : idx;

        label = bsdf_label(kg, &sd->closure[idx], *omega_in);
        bsdf_roughness_eta(kg, &sd->closure[idx], sampled_rougness, eta);
      }
    }

    assert(std::isfinite(*bsdf_pdf));
    assert(*bsdf_pdf >= 0.f);
    kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);
  }
  else {
    /* Sample BSDF. */
    *bsdf_pdf = 0.0f;
    label = bsdf_sample(kg,
                        sd,
                        sc,
                        rand_bsdf.x,
                        rand_bsdf.y,
                        &eval,
                        omega_in,
                        unguided_bsdf_pdf,
                        sampled_rougness,
                        eta);
#  if 0
    if (*unguided_bsdf_pdf > 0.0f) {
      surface_shader_validate_bsdf_sample(kg, sc, *omega_in, label, sampled_roughness, eta);
    }
#  endif

    if (*unguided_bsdf_pdf != 0.0f) {
      bsdf_eval_init(bsdf_eval, sc->type, eval * sc->weight);

      kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);

      if (sd->num_closure > 1) {
        float sweight = sc->sample_weight;
        *unguided_bsdf_pdf = _surface_shader_bsdf_eval_mis(
            kg, sd, *omega_in, sc, bsdf_eval, (*unguided_bsdf_pdf) * sweight, sweight, 0);
        kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);
      }
      *bsdf_pdf = *unguided_bsdf_pdf;
      *mis_pdf = *bsdf_pdf;
    }

    kernel_assert(reduce_min(bsdf_eval_sum(bsdf_eval)) >= 0.0f);
  }

  return label;
}

#endif

/* Sample direction for picked BSDF, and return evaluation and pdf for all
 * BSDFs combined using MIS. */
ccl_device int surface_shader_bsdf_sample_closure(KernelGlobals kg,
                                                  ccl_private ShaderData *sd,
                                                  ccl_private const ShaderClosure *sc,
                                                  const float2 rand_bsdf,
                                                  ccl_private BsdfEval *bsdf_eval,
                                                  ccl_private float3 *wo,
                                                  ccl_private float *pdf,
                                                  ccl_private float2 *sampled_roughness,
                                                  ccl_private float *eta)
{
  /* BSSRDF should already have been handled elsewhere. */
  kernel_assert(CLOSURE_IS_BSDF(sc->type));

  int label;
  Spectrum eval = zero_spectrum();

  *pdf = 0.0f;
  label = bsdf_sample(
      kg, sd, sc, rand_bsdf.x, rand_bsdf.y, &eval, wo, pdf, sampled_roughness, eta);

  if (*pdf != 0.0f) {
    bsdf_eval_init(bsdf_eval, sc->type, eval * sc->weight);

    if (sd->num_closure > 1) {
      float sweight = sc->sample_weight;
      *pdf = _surface_shader_bsdf_eval_mis(kg, sd, *wo, sc, bsdf_eval, *pdf * sweight, sweight, 0);
    }
  }
  else {
    bsdf_eval_init(bsdf_eval, sc->type, zero_spectrum());
  }

  return label;
}

ccl_device float surface_shader_average_roughness(ccl_private const ShaderData *sd)
{
  float roughness = 0.0f;
  float sum_weight = 0.0f;

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF(sc->type)) {
      /* sqrt once to undo the squaring from multiplying roughness on the
       * two axes, and once for the squared roughness convention. */
      float weight = fabsf(average(sc->weight));
      roughness += weight * sqrtf(safe_sqrtf(bsdf_get_roughness_squared(sc)));
      sum_weight += weight;
    }
  }

  return (sum_weight > 0.0f) ? roughness / sum_weight : 0.0f;
}

ccl_device Spectrum surface_shader_transparency(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  if (sd->flag & SD_HAS_ONLY_VOLUME) {
    return one_spectrum();
  }
  else if (sd->flag & SD_TRANSPARENT) {
    return sd->closure_transparent_extinction;
  }
  else {
    return zero_spectrum();
  }
}

ccl_device void surface_shader_disable_transparency(KernelGlobals kg, ccl_private ShaderData *sd)
{
  if (sd->flag & SD_TRANSPARENT) {
    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private ShaderClosure *sc = &sd->closure[i];

      if (sc->type == CLOSURE_BSDF_TRANSPARENT_ID) {
        sc->sample_weight = 0.0f;
        sc->weight = zero_spectrum();
      }
    }

    sd->flag &= ~SD_TRANSPARENT;
  }
}

ccl_device Spectrum surface_shader_alpha(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  Spectrum alpha = one_spectrum() - surface_shader_transparency(kg, sd);

  alpha = saturate(alpha);

  return alpha;
}

ccl_device Spectrum surface_shader_diffuse(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  Spectrum eval = zero_spectrum();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_DIFFUSE(sc->type) || CLOSURE_IS_BSSRDF(sc->type))
      eval += sc->weight;
  }

  return eval;
}

ccl_device Spectrum surface_shader_glossy(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  Spectrum eval = zero_spectrum();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_GLOSSY(sc->type))
      eval += sc->weight;
  }

  return eval;
}

ccl_device Spectrum surface_shader_transmission(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  Spectrum eval = zero_spectrum();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_TRANSMISSION(sc->type))
      eval += sc->weight;
  }

  return eval;
}

ccl_device float3 surface_shader_average_normal(KernelGlobals kg, ccl_private const ShaderData *sd)
{
  float3 N = zero_float3();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];
    if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type))
      N += sc->N * fabsf(average(sc->weight));
  }

  return (is_zero(N)) ? sd->N : normalize(N);
}

ccl_device Spectrum surface_shader_ao(KernelGlobals kg,
                                      ccl_private const ShaderData *sd,
                                      const float ao_factor,
                                      ccl_private float3 *N_)
{
  Spectrum eval = zero_spectrum();
  float3 N = zero_float3();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_DIFFUSE(sc->type)) {
      ccl_private const DiffuseBsdf *bsdf = (ccl_private const DiffuseBsdf *)sc;
      eval += sc->weight * ao_factor;
      N += bsdf->N * fabsf(average(sc->weight));
    }
  }

  *N_ = (is_zero(N)) ? sd->N : normalize(N);
  return eval;
}

#ifdef __SUBSURFACE__
ccl_device float3 surface_shader_bssrdf_normal(ccl_private const ShaderData *sd)
{
  float3 N = zero_float3();

  for (int i = 0; i < sd->num_closure; i++) {
    ccl_private const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSSRDF(sc->type)) {
      ccl_private const Bssrdf *bssrdf = (ccl_private const Bssrdf *)sc;
      float avg_weight = fabsf(average(sc->weight));

      N += bssrdf->N * avg_weight;
    }
  }

  return (is_zero(N)) ? sd->N : normalize(N);
}
#endif /* __SUBSURFACE__ */

/* Constant emission optimization */

ccl_device bool surface_shader_constant_emission(KernelGlobals kg,
                                                 int shader,
                                                 ccl_private Spectrum *eval)
{
  int shader_index = shader & SHADER_MASK;
  int shader_flag = kernel_data_fetch(shaders, shader_index).flags;

  if (shader_flag & SD_HAS_CONSTANT_EMISSION) {
    const float3 emission_rgb = make_float3(
        kernel_data_fetch(shaders, shader_index).constant_emission[0],
        kernel_data_fetch(shaders, shader_index).constant_emission[1],
        kernel_data_fetch(shaders, shader_index).constant_emission[2]);
    *eval = rgb_to_spectrum(emission_rgb);

    return true;
  }

  return false;
}

/* Background */

ccl_device Spectrum surface_shader_background(ccl_private const ShaderData *sd)
{
  if (sd->flag & SD_EMISSION) {
    return sd->closure_emission_background;
  }
  else {
    return zero_spectrum();
  }
}

/* Emission */

ccl_device Spectrum surface_shader_emission(ccl_private const ShaderData *sd)
{
  if (sd->flag & SD_EMISSION) {
    return emissive_simple_eval(sd->Ng, sd->wi) * sd->closure_emission_background;
  }
  else {
    return zero_spectrum();
  }
}

/* Holdout */

ccl_device Spectrum surface_shader_apply_holdout(KernelGlobals kg, ccl_private ShaderData *sd)
{
  Spectrum weight = zero_spectrum();

  /* For objects marked as holdout, preserve transparency and remove all other
   * closures, replacing them with a holdout weight. */
  if (sd->object_flag & SD_OBJECT_HOLDOUT_MASK) {
    if ((sd->flag & SD_TRANSPARENT) && !(sd->flag & SD_HAS_ONLY_VOLUME)) {
      weight = one_spectrum() - sd->closure_transparent_extinction;

      for (int i = 0; i < sd->num_closure; i++) {
        ccl_private ShaderClosure *sc = &sd->closure[i];
        if (!CLOSURE_IS_BSDF_TRANSPARENT(sc->type)) {
          sc->type = NBUILTIN_CLOSURES;
        }
      }

      sd->flag &= ~(SD_CLOSURE_FLAGS - (SD_TRANSPARENT | SD_BSDF));
    }
    else {
      weight = one_spectrum();
    }
  }
  else {
    for (int i = 0; i < sd->num_closure; i++) {
      ccl_private const ShaderClosure *sc = &sd->closure[i];
      if (CLOSURE_IS_HOLDOUT(sc->type)) {
        weight += sc->weight;
      }
    }
  }

  return weight;
}

/* Surface Evaluation */

template<uint node_feature_mask, typename ConstIntegratorGenericState>
ccl_device void surface_shader_eval(KernelGlobals kg,
                                    ConstIntegratorGenericState state,
                                    ccl_private ShaderData *ccl_restrict sd,
                                    ccl_global float *ccl_restrict buffer,
                                    uint32_t path_flag,
                                    bool use_caustics_storage = false)
{
  /* If path is being terminated, we are tracing a shadow ray or evaluating
   * emission, then we don't need to store closures. The emission and shadow
   * shader data also do not have a closure array to save GPU memory. */
  int max_closures;
  if (path_flag & (PATH_RAY_TERMINATE | PATH_RAY_SHADOW | PATH_RAY_EMISSION)) {
    max_closures = 0;
  }
  else {
    max_closures = use_caustics_storage ? CAUSTICS_MAX_CLOSURE : kernel_data.max_closures;
  }

  sd->num_closure = 0;
  sd->num_closure_left = max_closures;

#ifdef __OSL__
  if (kernel_data.kernel_features & KERNEL_FEATURE_OSL) {
    osl_eval_nodes<SHADER_TYPE_SURFACE>(kg, state, sd, path_flag);
  }
  else
#endif
  {
#ifdef __SVM__
    svm_eval_nodes<node_feature_mask, SHADER_TYPE_SURFACE>(kg, state, sd, buffer, path_flag);
#else
    if (sd->object == OBJECT_NONE) {
      sd->closure_emission_background = make_spectrum(0.8f);
      sd->flag |= SD_EMISSION;
    }
    else {
      ccl_private DiffuseBsdf *bsdf = (ccl_private DiffuseBsdf *)bsdf_alloc(
          sd, sizeof(DiffuseBsdf), make_spectrum(0.8f));
      if (bsdf != NULL) {
        bsdf->N = sd->N;
        sd->flag |= bsdf_diffuse_setup(bsdf);
      }
    }
#endif
  }
}

CCL_NAMESPACE_END
