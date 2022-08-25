/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#ifdef WITH_CYCLES_DEBUG
#  include "kernel/closure/alloc.h"
#  include "kernel/closure/bsdf.h"
#  include "kernel/film/write_passes.h"
#endif

CCL_NAMESPACE_BEGIN

/* Path recording for guiding. */

ccl_device_forceinline void guiding_record_light_surface_segment(
    KernelGlobals kg, IntegratorState state, ccl_private const Intersection *ccl_restrict isect)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f pglZero = openpgl::cpp::Vector3(0.f, 0.f, 0.f);
  const pgl_vec3f pglOne = openpgl::cpp::Vector3(1.f, 1.f, 1.f);
  float3 ray_P = INTEGRATOR_STATE(state, ray, P);
  float3 ray_D = INTEGRATOR_STATE(state, ray, D);
  float3 p = ray_P + isect->t * ray_D;
  pgl_point3f pglP = openpgl::cpp::Point3(p[0], p[1], p[2]);
  pgl_vec3f pglWi = openpgl::cpp::Vector3(-ray_D[0], -ray_D[1], -ray_D[2]);
  pgl_vec3f pglWo = openpgl::cpp::Vector3(ray_D[0], ray_D[1], ray_D[2]);

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, pglP);
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, pglWi);
  openpgl::cpp::SetNormal(state->guiding.path_segment, pglWi);
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, pglWo);
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, 1.0f);
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, pglOne);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, pglOne);
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

  const pgl_vec3f pglZero = openpgl::cpp::Vector3(0.f, 0.f, 0.f);
  const pgl_vec3f pglOne = openpgl::cpp::Vector3(1.f, 1.f, 1.f);
  pgl_point3f pglP = openpgl::cpp::Point3(sd->P.x, sd->P.y, sd->P.z);
  pgl_vec3f pglWi = openpgl::cpp::Vector3(sd->I.x, sd->I.y, sd->I.z);

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, pglP);
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, pglWi);
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, pglOne);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_surface_bounce(KernelGlobals kg,
                                                          IntegratorState state,
                                                          ccl_private const ShaderData *sd,
                                                          const Spectrum bsdf_weight,
                                                          const float bsdf_pdf,
                                                          const float3 bsdf_shading_normal,
                                                          const float3 bsdf_omega_in,
                                                          const float2 bsdf_roughness,
                                                          const float bsdf_eta)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  float roughness = fminf(bsdf_roughness.x, bsdf_roughness.y);
  if (roughness > 0.0f) {
    roughness = sqrt(roughness);
  }

  const bool bsdf_is_delta = roughness > 0.f ? false : true;
  const float3 bsdf_weight_rgb = spectrum_to_rgb(bsdf_weight);

  pgl_vec3f pglWo = openpgl::cpp::Vector3(bsdf_omega_in[0], bsdf_omega_in[1], bsdf_omega_in[2]);
  pgl_vec3f pglBSDFWeight = openpgl::cpp::Vector3(
      bsdf_weight_rgb[0], bsdf_weight_rgb[1], bsdf_weight_rgb[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(clamp(bsdf_shading_normal[0], -1.f, 1.f),
                                              clamp(bsdf_shading_normal[1], -1.f, 1.f),
                                              clamp(bsdf_shading_normal[2], -1.f, 1.f));

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment,
                                       openpgl::cpp::Vector3(1.0f, 1.0f, 1.0f));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetNormal(state->guiding.path_segment, pglNormal);
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, pglWo);
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, bsdf_pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, pglBSDFWeight);
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, bsdf_is_delta);
  openpgl::cpp::SetEta(state->guiding.path_segment, bsdf_eta);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, roughness);
#endif
}

ccl_device_forceinline void guiding_record_bssrdf_segment(KernelGlobals kg,
                                                          IntegratorState state,
                                                          const float3 &P,
                                                          const float3 &I)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f pglZero = openpgl::cpp::Vector3(0.f, 0.f, 0.f);
  const pgl_vec3f pglOne = openpgl::cpp::Vector3(1.f, 1.f, 1.f);
  pgl_point3f pglP = openpgl::cpp::Point3(P.x, P.y, P.z);
  pgl_vec3f pglWi = openpgl::cpp::Vector3(I.x, I.y, I.z);

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();
  openpgl::cpp::SetPosition(state->guiding.path_segment, pglP);
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, pglWi);
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, pglOne);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_bssrdf_bounce(KernelGlobals kg,
                                                         IntegratorState state,
                                                         const Spectrum bssrdf_weight,
                                                         const float bssrdf_pdf,
                                                         const float3 bssrdf_shading_normal,
                                                         const float3 bssrdf_omega_in)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 bssrdf_weight_rgb = spectrum_to_rgb(bssrdf_weight);

  pgl_vec3f pglWo = openpgl::cpp::Vector3(
      bssrdf_omega_in[0], bssrdf_omega_in[1], bssrdf_omega_in[2]);
  pgl_vec3f pglBSSRDFWeight = openpgl::cpp::Vector3(
      bssrdf_weight_rgb[0], bssrdf_weight_rgb[1], bssrdf_weight_rgb[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(clamp(bssrdf_shading_normal[0], -1.f, 1.f),
                                              clamp(bssrdf_shading_normal[1], -1.f, 1.f),
                                              clamp(bssrdf_shading_normal[2], -1.f, 1.f));

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment,
                                       openpgl::cpp::Vector3(1.0f, 1.0f, 1.0f));
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, false);
  openpgl::cpp::SetNormal(state->guiding.path_segment, pglNormal);
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, pglWo);
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, bssrdf_pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, pglBSSRDFWeight);
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, false);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0f);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, 1.0f);
#endif
}

ccl_device_forceinline void guiding_record_volume_segment(KernelGlobals kg,
                                                          IntegratorState state,
                                                          const float3 &P,
                                                          const float3 &I)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 1
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const pgl_vec3f pglZero = openpgl::cpp::Vector3(0.f, 0.f, 0.f);
  const pgl_vec3f pglOne = openpgl::cpp::Vector3(1.f, 1.f, 1.f);
  pgl_point3f pglP = openpgl::cpp::Point3(P.x, P.y, P.z);
  pgl_vec3f pglWi = openpgl::cpp::Vector3(I.x, I.y, I.z);

  state->guiding.path_segment = state->guiding.path_segment_storage->NextSegment();

  openpgl::cpp::SetPosition(state->guiding.path_segment, pglP);
  openpgl::cpp::SetDirectionOut(state->guiding.path_segment, pglWi);
  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetScatteredContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment, pglZero);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment, pglOne);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.0);
#endif
}

ccl_device_forceinline void guiding_record_volume_bounce(KernelGlobals kg,
                                                         IntegratorState state,
                                                         ccl_private const ShaderData *sd,
                                                         const Spectrum phase_weight,
                                                         const float phase_pdf,
                                                         const float3 phase_omega_in,
                                                         const float phase_roughness)
{
#if defined(__PATH_GUIDING__) && PATH_GUIDING_LEVEL >= 4
  if (!kernel_data.integrator.use_guiding) {
    return;
  }

  const float3 phase_weight_rgb = spectrum_to_rgb(phase_weight);

  pgl_vec3f pglWo = openpgl::cpp::Vector3(phase_omega_in[0], phase_omega_in[1], phase_omega_in[2]);
  pgl_vec3f pglPhaseWeight = openpgl::cpp::Vector3(
      phase_weight_rgb[0], phase_weight_rgb[1], phase_weight_rgb[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(0.f, 0.f, 1.f);

  kernel_assert(state->guiding.path_segment != nullptr);

  openpgl::cpp::SetVolumeScatter(state->guiding.path_segment, true);
  openpgl::cpp::SetTransmittanceWeight(state->guiding.path_segment,
                                       openpgl::cpp::Vector3(1.0f, 1.0f, 1.0f));
  openpgl::cpp::SetNormal(state->guiding.path_segment, pglNormal);
  openpgl::cpp::SetDirectionIn(state->guiding.path_segment, pglWo);
  openpgl::cpp::SetPDFDirectionIn(state->guiding.path_segment, phase_pdf);
  openpgl::cpp::SetScatteringWeight(state->guiding.path_segment, pglPhaseWeight);
  openpgl::cpp::SetIsDelta(state->guiding.path_segment, false);
  openpgl::cpp::SetEta(state->guiding.path_segment, 1.f);
  openpgl::cpp::SetRoughness(state->guiding.path_segment, phase_roughness);
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
  const float3 background_pos = state->ray.P + (1e6f) * state->ray.D;

  pgl_point3f pglEnvPos = openpgl::cpp::Vector3(
      background_pos[0], background_pos[1], background_pos[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(0.0f, 0.0f, 1.0f);
  pgl_vec3f pglDirectionOut = openpgl::cpp::Vector3(
      -state->ray.D.x, -state->ray.D.y, -state->ray.D.z);

  openpgl::cpp::PathSegment background_segment;
  openpgl::cpp::SetPosition(&background_segment, pglEnvPos);
  openpgl::cpp::SetNormal(&background_segment, pglNormal);
  openpgl::cpp::SetDirectionOut(&background_segment, pglDirectionOut);
  openpgl::cpp::SetDirectContribution(&background_segment,
                                      openpgl::cpp::Vector3(L_rgb[0], L_rgb[1], L_rgb[2]));
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

  openpgl::cpp::SetDirectContribution(state->guiding.path_segment,
                                      openpgl::cpp::Vector3(Le_rgb[0], Le_rgb[1], Le_rgb[2]));
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
    openpgl::cpp::AddScatteredContribution(state->shadow_path.path_segment,
                                           openpgl::cpp::Vector3(Lo_rgb[0], Lo_rgb[1], Lo_rgb[2]));
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