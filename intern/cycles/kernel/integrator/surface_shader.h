/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

/* Functions to evaluate shaders. */

#pragma once

#include "kernel/closure/alloc.h"
#include "kernel/closure/bsdf.h"
#include "kernel/closure/bsdf_util.h"
#include "kernel/closure/emissive.h"

#include "kernel/integrator/guiding.h"

#include "kernel/svm/svm.h"

#ifdef __OSL__
#  include "kernel/osl/shader.h"
#endif

CCL_NAMESPACE_BEGIN

/* Guiding */

#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
ccl_device_inline void surface_shader_prepare_guiding(KernelGlobals kg,
                                                      IntegratorState state,
                                                      ccl_private ShaderData *sd,
                                                      ccl_private const RNGState *rng_state)
{
  const bool guiding = kernel_data.integrator.use_guiding;
  const bool surface_guiding = kernel_data.integrator.use_surface_guiding;
  const float surface_guiding_probability = kernel_data.integrator.surface_guiding_probability;

  float diffuse_sampling_fraction = 0.f;
  float bssrdf_sampling_fraction = 0.f;
  float bsdf_bssrdf_sampling_sum = 0.f;
  float guiding_sampling_prob = 0.f;

  float grand = 0.f;

  bool useGuiding = false;

  if (guiding && surface_guiding) {

    grand = path_state_rng_1D(kg, rng_state, PRNG_SURFACE_BSDF_GUIDING);

    for (int i = 0; i < sd->num_closure; i++) {
      ShaderClosure *sc = &sd->closure[i];
      if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {

        const float sweight = sc->sample_weight;
        kernel_assert(sweight >= 0.f);
        bsdf_bssrdf_sampling_sum += sweight;
        if (CLOSURE_IS_BSDF_DIFFUSE(sc->type) && sc->type < CLOSURE_BSDF_TRANSLUCENT_ID) {
          diffuse_sampling_fraction += sweight;
        }
        if (CLOSURE_IS_BSSRDF(sc->type)) {
          bssrdf_sampling_fraction += sweight;
        }
      }
    }

    if (bsdf_bssrdf_sampling_sum > 0.f) {
      diffuse_sampling_fraction /= bsdf_bssrdf_sampling_sum;
      bssrdf_sampling_fraction /= bsdf_bssrdf_sampling_sum;
    }

    if (diffuse_sampling_fraction > 0.f) {
      useGuiding = state->guiding.surface_sampling_distribution->Init(
          kg->opgl_guiding_field, guiding_point3f(sd->P), grand, true);

      if (useGuiding) {
        state->guiding.surface_sampling_distribution->ApplyCosineProduct(guiding_point3f(sd->N));
        guiding_sampling_prob = surface_guiding_probability * diffuse_sampling_fraction;
      }
    }
  }

  kernel_assert(guiding_sampling_prob >= 0.f && guiding_sampling_prob <= 1.0f);
  state->guiding.surface_guiding_sampling_prob = guiding_sampling_prob;
  state->guiding.bssrdf_sampling_prob = bssrdf_sampling_fraction;
  state->guiding.use_surface_guiding = useGuiding;
  state->guiding.sample_surface_guiding_rand = grand;
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

ccl_device_inline bool surface_shader_is_transmission(ccl_private const ShaderData *sd,
                                                      const float3 omega_in)
{
  return dot(sd->N, omega_in) < 0.0f;
}

ccl_device_inline bool surface_shader_is_transmission3(const ShaderClosure *sc,
                                                       const float3 omega_in)
{
  return dot(sc->N, omega_in) < 0.0f;
}

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
                                                      const float3 omega_in,
                                                      const bool is_transmission,
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
      if (CLOSURE_IS_BSDF(sc->type) && !_surface_shader_exclude(sc->type, light_shader_flags)) {
        float bsdf_pdf = 0.0f;
        Spectrum eval = bsdf_eval(kg, sd, sc, omega_in, is_transmission, &bsdf_pdf);

        if (bsdf_pdf != 0.0f) {
          bsdf_eval_accum(result_eval, sc->type, eval * sc->weight);
          sum_pdf += bsdf_pdf * sc->sample_weight;
        }
      }

      sum_sample_weight += sc->sample_weight;
    }
  }

  return (sum_sample_weight > 0.0f) ? sum_pdf / sum_sample_weight : 0.0f;
}

ccl_device_inline float surface_shader_bsdf_eval_pdfs(const KernelGlobals kg,
                                                      ShaderData *sd,
                                                      const float3 omega_in,
                                                      const bool is_transmission,
                                                      BsdfEval *result_eval,
                                                      float *pdfs,
                                                      const uint light_shader_flags)
{
  float sum_pdf = 0.0f;
  float sum_sample_weight = 0.0f;
  bsdf_eval_init(result_eval, CLOSURE_NONE_ID, zero_spectrum());
  /* this is the veach one-sample model with balance heuristic, some pdf
   * factors drop out when using balance heuristic weighting */
  for (int i = 0; i < sd->num_closure; i++) {
    const ShaderClosure *sc = &sd->closure[i];

    if (CLOSURE_IS_BSDF_OR_BSSRDF(sc->type)) {
      if (CLOSURE_IS_BSDF(sc->type) && !_surface_shader_exclude(sc->type, light_shader_flags)) {
        float bsdf_pdf = 0.0f;
        Spectrum eval = bsdf_eval(kg, sd, sc, omega_in, is_transmission, &bsdf_pdf);
        kernel_assert(bsdf_pdf >= 0.f);
        if (bsdf_pdf != 0.0f) {
          bsdf_eval_accum(result_eval, sc->type, eval * sc->weight);
          sum_pdf += bsdf_pdf * sc->sample_weight;
          kernel_assert(bsdf_pdf * sc->sample_weight >= 0.f);
          pdfs[i] = bsdf_pdf * sc->sample_weight;
        }
        else {
          pdfs[i] = 0.f;
        }
      }
      else {
        pdfs[i] = 0.f;
      }

      sum_sample_weight += sc->sample_weight;
    }
    else {
      pdfs[i] = 0.f;
    }
  }
  if (sum_pdf > 0.f) {
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
                             ccl_private ShaderData *sd,
                             const float3 omega_in,
                             const bool is_transmission,
                             ccl_private BsdfEval *bsdf_eval,
                             const uint light_shader_flags)
{
  bsdf_eval_init(bsdf_eval, CLOSURE_NONE_ID, zero_spectrum());

  return _surface_shader_bsdf_eval_mis(
      kg, sd, omega_in, is_transmission, NULL, bsdf_eval, 0.0f, 0.0f, light_shader_flags);
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

#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
/* Sample direction for picked BSDF, and return evaluation and pdf for all
 * BSDFs combined using MIS. */

ccl_device int surface_shader_bsdf_guided_sample_closure(KernelGlobals kg,
                                                         IntegratorState state,
                                                         ShaderData *sd,
                                                         const ShaderClosure *sc,
                                                         const float2 rand_bsdf,
                                                         BsdfEval *bsdf_eval,
                                                         float3 *omega_in,
                                                         float *guided_bsdf_pdf,
                                                         float *bsdf_pdf,
                                                         float2 *sampled_rougness,
                                                         float *eta)
{
  /* BSSRDF should already have been handled elsewhere. */
  kernel_assert(CLOSURE_IS_BSDF(sc->type));
  bool useGuiding = kernel_data.integrator.use_guiding && state->guiding.use_surface_guiding;
  const float guiding_sampling_prob = state->guiding.surface_guiding_sampling_prob;
  float randw = state->guiding.sample_surface_guiding_rand;

  bool sampleGuding = false;
  if (useGuiding && randw < guiding_sampling_prob) {
    sampleGuding = true;
    randw /= guiding_sampling_prob;
  }
  else {
    randw -= guiding_sampling_prob;
    randw /= (1.0f - guiding_sampling_prob);
  }

  // used for debugging can be removed later
  Spectrum bsdf_weight = bsdf_eval_sum(bsdf_eval);

  int label = LABEL_NONE;
  Spectrum eval = zero_spectrum();
  bsdf_eval_init(bsdf_eval, CLOSURE_NONE_ID, eval);

  *bsdf_pdf = 0.f;
  float guide_pdf = 0.f;

  float bsdf_pdfs[MAX_CLOSURE];

  float bssrdf_sampling_prob = state->guiding.bssrdf_sampling_prob;

  if (sampleGuding) {
    // sample guiding distr.
    pgl_vec3f pglWo;
    pgl_point2f pglSample = openpgl::cpp::Point2(rand_bsdf.x, rand_bsdf.y);
    guide_pdf = state->guiding.surface_sampling_distribution->SamplePDF(pglSample, pglWo);
    *omega_in = make_float3(pglWo.x, pglWo.y, pglWo.z);
    *guided_bsdf_pdf = 0.0f;
    if (guide_pdf != 0.0f) {
      // TODO: update is_transmission when closure is picked
      const bool is_transmission = surface_shader_is_transmission(sd, *omega_in);
      *bsdf_pdf = surface_shader_bsdf_eval_pdfs(
          kg, sd, *omega_in, is_transmission, bsdf_eval, bsdf_pdfs, 0);
      *guided_bsdf_pdf = (guiding_sampling_prob * guide_pdf * (1.0f - bssrdf_sampling_prob)) +
                         ((1.f - guiding_sampling_prob) * (*bsdf_pdf));
      float sumPDFs = 0.0f;

      if (*bsdf_pdf > 0.f) {
        int idx = -1;
        for (int i = 0; i < sd->num_closure; i++) {
          sumPDFs += bsdf_pdfs[i];
          if (randw <= sumPDFs) {
            idx = i;
            break;
          }
        }

        kernel_assert(idx >= 0);
        // set the default idx to the last in the list
        // in case of numerical problems and randw is just >=1.0f and
        // the sum of all bsdf_pdfs is just < 1.0f
        idx = (randw > sumPDFs) ? sd->num_closure - 1 : idx;

        label = bsdf_label(kg, &sd->closure[idx], is_transmission);
      }
    }

    bsdf_weight = bsdf_eval_sum(bsdf_eval);
    kernel_assert(reduce_min(bsdf_weight) >= 0.0f);

    *sampled_rougness = make_float2(1.f, 1.f);
    *eta = 1.f;
  }
  else {
    *guided_bsdf_pdf = 0.0f;
    label = bsdf_sample(
        kg, sd, sc, rand_bsdf.x, rand_bsdf.y, &eval, omega_in, bsdf_pdf, sampled_rougness, eta);
#  ifdef WITH_CYCLES_DEBUG
    ///////
    // validation code to test the bsdf_label function
    ///////
    bool is_transmission3 = surface_shader_bsdf_is_transmission3(sc, *omega_in);
    int label2 = bsdf_label(kg, sc, is_transmission3);

    if (*bsdf_pdf > 0.f && label != label2) {
      std::cout << "LABEL ERROR: " << std::endl;
      std::cout << "LABEL:  reflect = " << (label & LABEL_REFLECT)
                << "\t transmit = " << (label & LABEL_TRANSMIT)
                << "\t diffuse = " << (label & LABEL_DIFFUSE)
                << "\t glossy = " << (label & LABEL_GLOSSY)
                << "\t singular = " << (label & LABEL_SINGULAR) << std::endl;
      std::cout << "LABEL2: reflect = " << (label2 & LABEL_REFLECT)
                << "\t transmit = " << (label2 & LABEL_TRANSMIT)
                << "\t diffuse = " << (label2 & LABEL_DIFFUSE)
                << "\t glossy = " << (label2 & LABEL_GLOSSY)
                << "\t singular = " << (label2 & LABEL_SINGULAR) << std::endl;
      // int label2 = bsdf_label(kg, sc, is_transmission3);
    }

    float2 sampled_rougness2;
    float eta2;
    bsdf_roughness_eta(kg, sc, &sampled_rougness2, &eta2);
    if (*bsdf_pdf > 0.f && (*eta != eta2 || sampled_rougness->x != sampled_rougness2.x ||
                            sampled_rougness->y != sampled_rougness2.y)) {
      std::cout << "ROUGHNESS ETA ERROR: " << std::endl;
      std::cout << "ETA:  eta = " << *eta << "\t eta2 = " << eta2
                << "\t rough.x = " << sampled_rougness->x << "\t rough.y = " << sampled_rougness->y
                << "\t rough2.x = " << sampled_rougness2.x
                << "\t rough2.y = " << sampled_rougness2.y << std::endl;
      bsdf_roughness_eta(kg, sc, &sampled_rougness2, &eta2);
    }

    ////////
#  endif

    if (*bsdf_pdf != 0.0f) {
      bsdf_eval_init(bsdf_eval, sc->type, eval * sc->weight);

      bsdf_weight = bsdf_eval_sum(bsdf_eval);
      kernel_assert(reduce_min(bsdf_weight) >= 0.0f);

      if (sd->num_closure > 1) {
        const bool is_transmission = surface_shader_is_transmission(sd, *omega_in);
        float sweight = sc->sample_weight;
        // BsdfEval bsdf_eval_old = *bsdf_eval;
        *bsdf_pdf = _surface_shader_bsdf_eval_mis(
            kg, sd, *omega_in, is_transmission, sc, bsdf_eval, (*bsdf_pdf) * sweight, sweight, 0);
        bsdf_weight = bsdf_eval_sum(bsdf_eval);
        kernel_assert(reduce_min(bsdf_weight) >= 0.0f);
      }
      *guided_bsdf_pdf = *bsdf_pdf;

      if (useGuiding) {
        const float3 omega = *omega_in;
        pgl_vec3f pglWo = openpgl::cpp::Vector3(omega[0], omega[1], omega[2]);
        guide_pdf = state->guiding.surface_sampling_distribution->PDF(pglWo);
        *guided_bsdf_pdf *= 1.0f - guiding_sampling_prob;
        *guided_bsdf_pdf += guiding_sampling_prob * guide_pdf * (1.0f - bssrdf_sampling_prob);
      }
    }
    else {
      bsdf_eval_init(bsdf_eval, sc->type, zero_spectrum());
    }

    bsdf_weight = bsdf_eval_sum(bsdf_eval);
    kernel_assert(reduce_min(bsdf_weight) >= 0.0f);
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
                                                  ccl_private float3 *omega_in,
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
      kg, sd, sc, rand_bsdf.x, rand_bsdf.y, &eval, omega_in, pdf, sampled_roughness, eta);

  if (*pdf != 0.0f) {
    bsdf_eval_init(bsdf_eval, sc->type, eval * sc->weight);

    if (sd->num_closure > 1) {
      const bool is_transmission = surface_shader_is_transmission(sd, *omega_in);
      float sweight = sc->sample_weight;
      *pdf = _surface_shader_bsdf_eval_mis(
          kg, sd, *omega_in, is_transmission, sc, bsdf_eval, *pdf * sweight, sweight, 0);
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
    return emissive_simple_eval(sd->Ng, sd->I) * sd->closure_emission_background;
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
  if (kg->osl) {
    if (sd->object == OBJECT_NONE && sd->lamp == LAMP_NONE) {
      OSLShader::eval_background(kg, state, sd, path_flag);
    }
    else {
      OSLShader::eval_surface(kg, state, sd, path_flag);
    }
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
