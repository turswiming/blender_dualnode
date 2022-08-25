/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#ifdef WITH_CYCLES_DEBUG
#  include "kernel/closure/alloc.h"
#  include "kernel/closure/bsdf.h"
#  include "kernel/film/write_passes.h"
#endif

CCL_NAMESPACE_BEGIN

/* Utilities. */

#if defined(__PATH_GUIDING__)
static pgl_vec3f guiding_vec3f(const float3 v)
{
  return openpgl::cpp::Vector3(v.x, v.y, v.z);
}

static pgl_point3f guiding_point3f(const float3 v)
{
  return openpgl::cpp::Point3(v.x, v.y, v.z);
}
#endif

/* Path recording for guiding. */

ccl_device_forceinline void guiding_record_light_surface_segment(
    KernelGlobals kg, IntegratorState state, ccl_private const Intersection *ccl_restrict isect)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f zero = guiding_vec3f(zero_float3());
  const pgl_vec3f one = guiding_vec3f(one_float3());
  const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
  const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
  const float3 P = ray_P + isect->t * ray_D;

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, guiding_point3f(P));
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, guiding_vec3f(-ray_D));
  openpgl::cpp::SetNormal(state->guiding.path_segment, guiding_vec3f(-ray_D));
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, guiding_vec3f(ray_D));
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, 1.0f);
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, one);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, one);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0f);
#endif
}

ccl_device_forceinline void guiding_record_surface_segment(KernelGlobals kg,
                                                           IntegratorState state,
                                                           ccl_private const ShaderData *sd)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f zero = guiding_vec3f(zero_float3());
  const pgl_vec3f one = guiding_vec3f(one_float3());

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, guiding_point3f(sd->P));
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, guiding_vec3f(sd->I));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, one);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_surface_bounce(KernelGlobals kg,
                                                          IntegratorState state,
                                                          ccl_private const ShaderData *sd,
                                                          const Spectrum weight,
                                                          const float pdf,
                                                          const float3 N,
                                                          const float3 omega_in,
                                                          const float2 roughness,
                                                          const float eta)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float min_roughness = safe_sqrtf(fminf(roughness.x, roughness.y));
  const bool is_delta = (min_roughness == 0.0f);
  const float3 weight_rgb = spectrum_to_rgb(weight);
  const float3 normal = clamp(N, -one_float3(), one_float3());

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, guiding_vec3f(one_float3()));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetNormal(state->guiding.path_segment, guiding_vec3f(normal));
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, guiding_vec3f(omega_in));
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, guiding_vec3f(weight_rgb));
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, is_delta);
  openpgl::cpp::SetEta(state->guiding.path_segment, eta);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, min_roughness);
#endif
}

ccl_device_forceinline void guiding_record_bssrdf_segment(KernelGlobals kg,
                                                          IntegratorState state,
                                                          const float3 P,
                                                          const float3 I)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f zero = guiding_vec3f(zero_float3());
  const pgl_vec3f one = guiding_vec3f(one_float3());

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, guiding_point3f(P));
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, guiding_vec3f(I));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, one);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_bssrdf_bounce(KernelGlobals kg,
                                                         IntegratorState state,
                                                         const Spectrum weight,
                                                         const float pdf,
                                                         const float3 N,
                                                         const float3 omega_in)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 weight_rgb = spectrum_to_rgb(weight);
  const float3 normal = clamp(N, -one_float3(), one_float3());

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, guiding_vec3f(one_float3()));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetNormal(state->guiding.path_segment, guiding_vec3f(normal));
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, guiding_vec3f(omega_in));
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, guiding_vec3f(weight_rgb));
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, false);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0f);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, 1.0f);
#endif
}

ccl_device_forceinline void guiding_record_volume_segment(KernelGlobals kg,
                                                          IntegratorState state,
                                                          const float3 P,
                                                          const float3 I)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f zero = guiding_vec3f(zero_float3());
  const pgl_vec3f one = guiding_vec3f(one_float3());

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();

  openpgl::cpp::SetPosition(state->guiding.path_segment, guiding_point3f(P));
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, guiding_vec3f(I));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, zero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, one);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_volume_bounce(KernelGlobals kg,
                                                         IntegratorState state,
                                                         ccl_private const ShaderData *sd,
                                                         const Spectrum weight,
                                                         const float pdf,
                                                         const float3 omega_in,
                                                         const float roughness)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 weight_rgb = spectrum_to_rgb(weight);
  const float3 normal = make_float3(0.0f, 0.0f, 1.0f);

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, guiding_vec3f(one_float3()));
  openpgl::cpp::SetNormal(state->guiding.path_segment, guiding_vec3f(normal));
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, guiding_vec3f(omega_in));
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, guiding_vec3f(weight_rgb));
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, false);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.f);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, roughness);
#endif
}

ccl_device_forceinline void guiding_record_background(KernelGlobals kg,
                                                      IntegratorState state,
                                                      const Spectrum L,
                                                      const float mis_weight)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 L_rgb = spectrum_to_rgb(L);
  const float3 ray_P = INTEGRATOR_STATE(state, ray, P);
  const float3 ray_D = INTEGRATOR_STATE(state, ray, D);
  const float3 P = ray_P + (1e6f) * ray_D;
  const float3 normal = make_float3(0.0f, 0.0f, 1.0f);

  openpgl::cpp::PathSegment background_segment;
  openpgl::cpp::SetPosition(&background_segment, guiding_vec3f(P));
  openpgl::cpp::SetNormal(&background_segment, guiding_vec3f(normal));
  openpgl::cpp::SetDirectionOut(&background_segment, guiding_vec3f(-ray_D));
  openpgl::cpp::SetDirectContribution(&background_segment, guiding_vec3f(L_rgb));
  openpgl::cpp::SetMiWeight(&background_segment, mis_weight);
  state->guiding.path_segment_storage->AddSegment(background_segment);
#endif
}

ccl_device_forceinline void guiding_record_surface_emission(KernelGlobals kg,
                                                            IntegratorState state,
                                                            const Spectrum Le,
                                                            const float mis_weight)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 Le_rgb = spectrum_to_rgb(Le);

  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, guiding_vec3f(Le_rgb));
  openpgl::cpp::SetMiWeight(state->guiding.path_segment, mis_weight);
#endif
}

ccl_device_forceinline void guiding_record_direct_light(KernelGlobals kg,
                                                        IntegratorShadowState state)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  if (state->shadow_path.path_segment) {
    const Spectrum Lo = INTEGRATOR_STATE(state, shadow_path, scattered_contribution);
    const float3 Lo_rgb = spectrum_to_rgb(Lo);
    openpgl::cpp::AddScatteredContribution(state->shadow_path.path_segment, guiding_vec3f(Lo_rgb));
  }
#endif
}

ccl_device_forceinline void guiding_record_continuation_probability(
    KernelGlobals kg, IntegratorState state, const float continuation_probability)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  if (state->guiding.path_segment) {
    openpgl::cpp::SetRussianRouletteProbability(state->guiding.path_segment,
                                                continuation_probability);
  }
#endif
}

/* Path guiding debug render passes. */

ccl_device_forceinline void guiding_write_debug_passes(KernelGlobals kg,
                                                       IntegratorStateCPU *state,
                                                       ccl_private const ShaderData *sd,
                                                       ccl_global float *ccl_restrict
                                                           render_buffer)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4 && defined(WITH_CYCLES_DEBUG)
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  if (INTEGRATOR_STATE(state, path, bounce) != 0) {
    return;
  }

  const uint32_t render_pixel_index = INTEGRATOR_STATE(state, path, render_pixel_index);
  const uint64_t render_buffer_offset = (uint64_t)render_pixel_index *
                                        kernel_data.film.pass_stride;
  ccl_global float *buffer = render_buffer + render_buffer_offset;

  if (kernel_data.film.pass_guiding_probability != PASS_UNUSED) {
    float guiding_prob = state->guiding.surface_guiding_sampling_prob;
    kernel_write_pass_float(buffer + kernel_data.film.pass_guiding_probability, guiding_prob);
  }

  if (kernel_data.film.pass_guiding_avg_roughness != PASS_UNUSED) {
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

    kernel_write_pass_float(buffer + kernel_data.film.pass_guiding_avg_roughness, avg_roughness);
  }
#endif
}

CCL_NAMESPACE_END
