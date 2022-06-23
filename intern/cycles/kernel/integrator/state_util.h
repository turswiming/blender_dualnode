/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include <iostream>

#include "kernel/integrator/state.h"

#include "kernel/util/differential.h"

#if defined(PATH_GUIDING_DEBUG_PASS)
#  include "kernel/closure/alloc.h"
#  include "kernel/closure/bsdf.h"
#  include "kernel/film/write_passes.h"
#endif
CCL_NAMESPACE_BEGIN

/* Ray */

ccl_device_forceinline void integrator_state_write_ray(KernelGlobals kg,
                                                       IntegratorState state,
                                                       ccl_private const Ray *ccl_restrict ray)
{
  INTEGRATOR_STATE_WRITE(state, ray, P) = ray->P;
  INTEGRATOR_STATE_WRITE(state, ray, D) = ray->D;
  INTEGRATOR_STATE_WRITE(state, ray, t) = ray->t;
  INTEGRATOR_STATE_WRITE(state, ray, time) = ray->time;
  INTEGRATOR_STATE_WRITE(state, ray, dP) = ray->dP;
  INTEGRATOR_STATE_WRITE(state, ray, dD) = ray->dD;
}

ccl_device_forceinline void integrator_state_read_ray(KernelGlobals kg,
                                                      ConstIntegratorState state,
                                                      ccl_private Ray *ccl_restrict ray)
{
  ray->P = INTEGRATOR_STATE(state, ray, P);
  ray->D = INTEGRATOR_STATE(state, ray, D);
  ray->t = INTEGRATOR_STATE(state, ray, t);
  ray->time = INTEGRATOR_STATE(state, ray, time);
  ray->dP = INTEGRATOR_STATE(state, ray, dP);
  ray->dD = INTEGRATOR_STATE(state, ray, dD);
}

/* Shadow Ray */

ccl_device_forceinline void integrator_state_write_shadow_ray(
    KernelGlobals kg, IntegratorShadowState state, ccl_private const Ray *ccl_restrict ray)
{
  INTEGRATOR_STATE_WRITE(state, shadow_ray, P) = ray->P;
  INTEGRATOR_STATE_WRITE(state, shadow_ray, D) = ray->D;
  INTEGRATOR_STATE_WRITE(state, shadow_ray, t) = ray->t;
  INTEGRATOR_STATE_WRITE(state, shadow_ray, time) = ray->time;
  INTEGRATOR_STATE_WRITE(state, shadow_ray, dP) = ray->dP;
}

ccl_device_forceinline void integrator_state_read_shadow_ray(KernelGlobals kg,
                                                             ConstIntegratorShadowState state,
                                                             ccl_private Ray *ccl_restrict ray)
{
  ray->P = INTEGRATOR_STATE(state, shadow_ray, P);
  ray->D = INTEGRATOR_STATE(state, shadow_ray, D);
  ray->t = INTEGRATOR_STATE(state, shadow_ray, t);
  ray->time = INTEGRATOR_STATE(state, shadow_ray, time);
  ray->dP = INTEGRATOR_STATE(state, shadow_ray, dP);
  ray->dD = differential_zero_compact();
}

/* Intersection */

ccl_device_forceinline void integrator_state_write_isect(
    KernelGlobals kg, IntegratorState state, ccl_private const Intersection *ccl_restrict isect)
{
  INTEGRATOR_STATE_WRITE(state, isect, t) = isect->t;
  INTEGRATOR_STATE_WRITE(state, isect, u) = isect->u;
  INTEGRATOR_STATE_WRITE(state, isect, v) = isect->v;
  INTEGRATOR_STATE_WRITE(state, isect, object) = isect->object;
  INTEGRATOR_STATE_WRITE(state, isect, prim) = isect->prim;
  INTEGRATOR_STATE_WRITE(state, isect, type) = isect->type;
}

ccl_device_forceinline void integrator_state_read_isect(
    KernelGlobals kg, ConstIntegratorState state, ccl_private Intersection *ccl_restrict isect)
{
  isect->prim = INTEGRATOR_STATE(state, isect, prim);
  isect->object = INTEGRATOR_STATE(state, isect, object);
  isect->type = INTEGRATOR_STATE(state, isect, type);
  isect->u = INTEGRATOR_STATE(state, isect, u);
  isect->v = INTEGRATOR_STATE(state, isect, v);
  isect->t = INTEGRATOR_STATE(state, isect, t);
}

ccl_device_forceinline VolumeStack integrator_state_read_volume_stack(ConstIntegratorState state,
                                                                      int i)
{
  VolumeStack entry = {INTEGRATOR_STATE_ARRAY(state, volume_stack, i, object),
                       INTEGRATOR_STATE_ARRAY(state, volume_stack, i, shader)};
  return entry;
}

ccl_device_forceinline void integrator_state_write_volume_stack(IntegratorState state,
                                                                int i,
                                                                VolumeStack entry)
{
  INTEGRATOR_STATE_ARRAY_WRITE(state, volume_stack, i, object) = entry.object;
  INTEGRATOR_STATE_ARRAY_WRITE(state, volume_stack, i, shader) = entry.shader;
}

ccl_device_forceinline bool integrator_state_volume_stack_is_empty(KernelGlobals kg,
                                                                   ConstIntegratorState state)
{
  return (kernel_data.kernel_features & KERNEL_FEATURE_VOLUME) ?
             INTEGRATOR_STATE_ARRAY(state, volume_stack, 0, shader) == SHADER_NONE :
             true;
}

/* Shadow Intersection */

ccl_device_forceinline void integrator_state_write_shadow_isect(
    IntegratorShadowState state,
    ccl_private const Intersection *ccl_restrict isect,
    const int index)
{
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, t) = isect->t;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, u) = isect->u;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, v) = isect->v;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, object) = isect->object;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, prim) = isect->prim;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_isect, index, type) = isect->type;
}

ccl_device_forceinline void integrator_state_read_shadow_isect(
    ConstIntegratorShadowState state,
    ccl_private Intersection *ccl_restrict isect,
    const int index)
{
  isect->prim = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, prim);
  isect->object = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, object);
  isect->type = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, type);
  isect->u = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, u);
  isect->v = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, v);
  isect->t = INTEGRATOR_STATE_ARRAY(state, shadow_isect, index, t);
}

ccl_device_forceinline void integrator_state_copy_volume_stack_to_shadow(
    KernelGlobals kg, IntegratorShadowState shadow_state, ConstIntegratorState state)
{
  if (kernel_data.kernel_features & KERNEL_FEATURE_VOLUME) {
    int index = 0;
    int shader;
    do {
      shader = INTEGRATOR_STATE_ARRAY(state, volume_stack, index, shader);

      INTEGRATOR_STATE_ARRAY_WRITE(shadow_state, shadow_volume_stack, index, object) =
          INTEGRATOR_STATE_ARRAY(state, volume_stack, index, object);
      INTEGRATOR_STATE_ARRAY_WRITE(shadow_state, shadow_volume_stack, index, shader) = shader;

      ++index;
    } while (shader != OBJECT_NONE);
  }
}

ccl_device_forceinline void integrator_state_copy_volume_stack(KernelGlobals kg,
                                                               IntegratorState to_state,
                                                               ConstIntegratorState state)
{
  if (kernel_data.kernel_features & KERNEL_FEATURE_VOLUME) {
    int index = 0;
    int shader;
    do {
      shader = INTEGRATOR_STATE_ARRAY(state, volume_stack, index, shader);

      INTEGRATOR_STATE_ARRAY_WRITE(to_state, volume_stack, index, object) = INTEGRATOR_STATE_ARRAY(
          state, volume_stack, index, object);
      INTEGRATOR_STATE_ARRAY_WRITE(to_state, volume_stack, index, shader) = shader;

      ++index;
    } while (shader != OBJECT_NONE);
  }
}

ccl_device_forceinline VolumeStack
integrator_state_read_shadow_volume_stack(ConstIntegratorShadowState state, int i)
{
  VolumeStack entry = {INTEGRATOR_STATE_ARRAY(state, shadow_volume_stack, i, object),
                       INTEGRATOR_STATE_ARRAY(state, shadow_volume_stack, i, shader)};
  return entry;
}

ccl_device_forceinline bool integrator_state_shadow_volume_stack_is_empty(
    KernelGlobals kg, ConstIntegratorShadowState state)
{
  return (kernel_data.kernel_features & KERNEL_FEATURE_VOLUME) ?
             INTEGRATOR_STATE_ARRAY(state, shadow_volume_stack, 0, shader) == SHADER_NONE :
             true;
}

ccl_device_forceinline void integrator_state_write_shadow_volume_stack(IntegratorShadowState state,
                                                                       int i,
                                                                       VolumeStack entry)
{
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_volume_stack, i, object) = entry.object;
  INTEGRATOR_STATE_ARRAY_WRITE(state, shadow_volume_stack, i, shader) = entry.shader;
}

#if defined(__KERNEL_GPU__)
ccl_device_inline void integrator_state_copy_only(KernelGlobals kg,
                                                  ConstIntegratorState to_state,
                                                  ConstIntegratorState state)
{
  int index;

  /* Rely on the compiler to optimize out unused assignments and `while(false)`'s. */

#  define KERNEL_STRUCT_BEGIN(name) \
    index = 0; \
    do {

#  define KERNEL_STRUCT_MEMBER(parent_struct, type, name, feature) \
    if (kernel_integrator_state.parent_struct.name != nullptr) { \
      kernel_integrator_state.parent_struct.name[to_state] = \
          kernel_integrator_state.parent_struct.name[state]; \
    }

#  define KERNEL_STRUCT_ARRAY_MEMBER(parent_struct, type, name, feature) \
    if (kernel_integrator_state.parent_struct[index].name != nullptr) { \
      kernel_integrator_state.parent_struct[index].name[to_state] = \
          kernel_integrator_state.parent_struct[index].name[state]; \
    }

#  define KERNEL_STRUCT_END(name) \
    } \
    while (false) \
      ;

#  define KERNEL_STRUCT_END_ARRAY(name, cpu_array_size, gpu_array_size) \
    ++index; \
    } \
    while (index < gpu_array_size) \
      ;

#  define KERNEL_STRUCT_VOLUME_STACK_SIZE kernel_data.volume_stack_size

#  include "kernel/integrator/state_template.h"

#  undef KERNEL_STRUCT_BEGIN
#  undef KERNEL_STRUCT_MEMBER
#  undef KERNEL_STRUCT_ARRAY_MEMBER
#  undef KERNEL_STRUCT_END
#  undef KERNEL_STRUCT_END_ARRAY
#  undef KERNEL_STRUCT_VOLUME_STACK_SIZE
}

ccl_device_inline void integrator_state_move(KernelGlobals kg,
                                             ConstIntegratorState to_state,
                                             ConstIntegratorState state)
{
  integrator_state_copy_only(kg, to_state, state);

  INTEGRATOR_STATE_WRITE(state, path, queued_kernel) = 0;
}

ccl_device_inline void integrator_shadow_state_copy_only(KernelGlobals kg,
                                                         ConstIntegratorShadowState to_state,
                                                         ConstIntegratorShadowState state)
{
  int index;

  /* Rely on the compiler to optimize out unused assignments and `while(false)`'s. */

#  define KERNEL_STRUCT_BEGIN(name) \
    index = 0; \
    do {

#  define KERNEL_STRUCT_MEMBER(parent_struct, type, name, feature) \
    if (kernel_integrator_state.parent_struct.name != nullptr) { \
      kernel_integrator_state.parent_struct.name[to_state] = \
          kernel_integrator_state.parent_struct.name[state]; \
    }

#  define KERNEL_STRUCT_ARRAY_MEMBER(parent_struct, type, name, feature) \
    if (kernel_integrator_state.parent_struct[index].name != nullptr) { \
      kernel_integrator_state.parent_struct[index].name[to_state] = \
          kernel_integrator_state.parent_struct[index].name[state]; \
    }

#  define KERNEL_STRUCT_END(name) \
    } \
    while (false) \
      ;

#  define KERNEL_STRUCT_END_ARRAY(name, cpu_array_size, gpu_array_size) \
    ++index; \
    } \
    while (index < gpu_array_size) \
      ;

#  define KERNEL_STRUCT_VOLUME_STACK_SIZE kernel_data.volume_stack_size

#  include "kernel/integrator/shadow_state_template.h"

#  undef KERNEL_STRUCT_BEGIN
#  undef KERNEL_STRUCT_MEMBER
#  undef KERNEL_STRUCT_ARRAY_MEMBER
#  undef KERNEL_STRUCT_END
#  undef KERNEL_STRUCT_END_ARRAY
#  undef KERNEL_STRUCT_VOLUME_STACK_SIZE
}

ccl_device_inline void integrator_shadow_state_move(KernelGlobals kg,
                                                    ConstIntegratorState to_state,
                                                    ConstIntegratorState state)
{
  integrator_shadow_state_copy_only(kg, to_state, state);

  INTEGRATOR_STATE_WRITE(state, shadow_path, queued_kernel) = 0;
}

#endif

/* NOTE: Leaves kernel scheduling information untouched. Use INIT semantic for one of the paths
 * after this function. */
ccl_device_inline IntegratorState integrator_state_shadow_catcher_split(KernelGlobals kg,
                                                                        IntegratorState state)
{
#if defined(__KERNEL_GPU__)
  ConstIntegratorState to_state = atomic_fetch_and_add_uint32(
      &kernel_integrator_state.next_main_path_index[0], 1);

  integrator_state_copy_only(kg, to_state, state);
#else
  IntegratorStateCPU *ccl_restrict to_state = state + 1;

  /* Only copy the required subset for performance. */
  to_state->path = state->path;
  to_state->ray = state->ray;
  to_state->isect = state->isect;
  integrator_state_copy_volume_stack(kg, to_state, state);
#endif

  return to_state;
}

#ifdef __KERNEL_CPU__
ccl_device_inline int integrator_state_bounce(ConstIntegratorState state, const int)
{
  return INTEGRATOR_STATE(state, path, bounce);
}

ccl_device_inline int integrator_state_bounce(ConstIntegratorShadowState state, const int)
{
  return INTEGRATOR_STATE(state, shadow_path, bounce);
}

ccl_device_inline int integrator_state_diffuse_bounce(ConstIntegratorState state, const int)
{
  return INTEGRATOR_STATE(state, path, diffuse_bounce);
}

ccl_device_inline int integrator_state_diffuse_bounce(ConstIntegratorShadowState state, const int)
{
  return INTEGRATOR_STATE(state, shadow_path, diffuse_bounce);
}

ccl_device_inline int integrator_state_glossy_bounce(ConstIntegratorState state, const int)
{
  return INTEGRATOR_STATE(state, path, glossy_bounce);
}

ccl_device_inline int integrator_state_glossy_bounce(ConstIntegratorShadowState state, const int)
{
  return INTEGRATOR_STATE(state, shadow_path, glossy_bounce);
}

ccl_device_inline int integrator_state_transmission_bounce(ConstIntegratorState state, const int)
{
  return INTEGRATOR_STATE(state, path, transmission_bounce);
}

ccl_device_inline int integrator_state_transmission_bounce(ConstIntegratorShadowState state,
                                                           const int)
{
  return INTEGRATOR_STATE(state, shadow_path, transmission_bounce);
}

ccl_device_inline int integrator_state_transparent_bounce(ConstIntegratorState state, const int)
{
  return INTEGRATOR_STATE(state, path, transparent_bounce);
}

ccl_device_inline int integrator_state_transparent_bounce(ConstIntegratorShadowState state,
                                                          const int)
{
  return INTEGRATOR_STATE(state, shadow_path, transparent_bounce);
}
#else
ccl_device_inline int integrator_state_bounce(ConstIntegratorShadowState state,
                                              const uint32_t path_flag)
{
  return (path_flag & PATH_RAY_SHADOW) ? INTEGRATOR_STATE(state, shadow_path, bounce) :
                                         INTEGRATOR_STATE(state, path, bounce);
}

ccl_device_inline int integrator_state_diffuse_bounce(ConstIntegratorShadowState state,
                                                      const uint32_t path_flag)
{
  return (path_flag & PATH_RAY_SHADOW) ? INTEGRATOR_STATE(state, shadow_path, diffuse_bounce) :
                                         INTEGRATOR_STATE(state, path, diffuse_bounce);
}

ccl_device_inline int integrator_state_glossy_bounce(ConstIntegratorShadowState state,
                                                     const uint32_t path_flag)
{
  return (path_flag & PATH_RAY_SHADOW) ? INTEGRATOR_STATE(state, shadow_path, glossy_bounce) :
                                         INTEGRATOR_STATE(state, path, glossy_bounce);
}

ccl_device_inline int integrator_state_transmission_bounce(ConstIntegratorShadowState state,
                                                           const uint32_t path_flag)
{
  return (path_flag & PATH_RAY_SHADOW) ?
             INTEGRATOR_STATE(state, shadow_path, transmission_bounce) :
             INTEGRATOR_STATE(state, path, transmission_bounce);
}

ccl_device_inline int integrator_state_transparent_bounce(ConstIntegratorShadowState state,
                                                          const uint32_t path_flag)
{
  return (path_flag & PATH_RAY_SHADOW) ? INTEGRATOR_STATE(state, shadow_path, transparent_bounce) :
                                         INTEGRATOR_STATE(state, path, transparent_bounce);
}
#endif

#if defined(__PATH_GUIDING__)
#  if PATH_GUIDING_LEVEL >= 1
ccl_device_forceinline void guiding_new_virtual_light_segment(
    IntegratorState state, ccl_private const Intersection *ccl_restrict isect)
{
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
}

ccl_device_forceinline void guiding_new_surface_segment(IntegratorState state,
                                                        const ShaderData *sd)
{
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
}

ccl_device_forceinline void guiding_add_bsdf_data(IntegratorState state,
                                                  const ShaderData *sd,
                                                  const float3 bsdf_weight,
                                                  const float bsdf_pdf,
                                                  const float3 bsdf_shading_normal,
                                                  const float3 bsdf_omega_in,
                                                  const float2 bsdf_roughness,
                                                  const float bsdf_eta,
                                                  const bool bsdf_is_delta)
{
  float roughness = fminf(bsdf_roughness.x, bsdf_roughness.y);
  if (roughness > 0.0f)
    roughness = sqrt(roughness);

  pgl_vec3f pglWo = openpgl::cpp::Vector3(bsdf_omega_in[0], bsdf_omega_in[1], bsdf_omega_in[2]);
  pgl_vec3f pglBSDFWeight = openpgl::cpp::Vector3(bsdf_weight[0], bsdf_weight[1], bsdf_weight[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(clamp(bsdf_shading_normal[0], -1.f, 1.f),
                                              clamp(bsdf_shading_normal[1], -1.f, 1.f),
                                              clamp(bsdf_shading_normal[2], -1.f, 1.f));

  assert(state->guiding.path_segment != nullptr);

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
}

ccl_device_forceinline void guiding_new_bssrdf_segment(IntegratorState state,
                                                       const float3 &P,
                                                       const float3 &I)
{
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
}

ccl_device_forceinline void guiding_add_bssrdf_data(IntegratorState state,
                                                    const float3 bssrdf_weight,
                                                    const float bssrdf_pdf,
                                                    const float3 bssrdf_shading_normal,
                                                    const float3 bssrdf_omega_in)
{
  pgl_vec3f pglWo = openpgl::cpp::Vector3(
      bssrdf_omega_in[0], bssrdf_omega_in[1], bssrdf_omega_in[2]);
  pgl_vec3f pglBSSRDFWeight = openpgl::cpp::Vector3(
      bssrdf_weight[0], bssrdf_weight[1], bssrdf_weight[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(clamp(bssrdf_shading_normal[0], -1.f, 1.f),
                                              clamp(bssrdf_shading_normal[1], -1.f, 1.f),
                                              clamp(bssrdf_shading_normal[2], -1.f, 1.f));

  assert(state->guiding.path_segment != nullptr);

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
}

ccl_device_forceinline void guiding_new_volume_segment(IntegratorState state,
                                                       const float3 &P,
                                                       const float3 &I)
{
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
}

ccl_device_forceinline void guiding_add_phase_data(IntegratorState state,
                                                   const ShaderData *sd,
                                                   const float3 phase_weight,
                                                   const float phase_pdf,
                                                   const float3 phase_omega_in,
                                                   const float phase_roughness)
{
  pgl_vec3f pglWo = openpgl::cpp::Vector3(phase_omega_in[0], phase_omega_in[1], phase_omega_in[2]);
  pgl_vec3f pglPhaseWeight = openpgl::cpp::Vector3(
      phase_weight[0], phase_weight[1], phase_weight[2]);
  pgl_vec3f pglNormal = openpgl::cpp::Vector3(0.f, 0.f, 1.f);

  assert(state->guiding.path_segment != nullptr);

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
}

ccl_device_forceinline void guiding_add_background(IntegratorState state,
                                                   const float3 L,
                                                   const float mis_weight)
{
  float3 background_pos = state->ray.P + (1e6f) * state->ray.D;
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
                                      openpgl::cpp::Vector3(L[0], L[1], L[2]));
  openpgl::cpp::SetMiWeight(&background_segment, mis_weight);
  state->guiding.path_segment_storage->AddSegment(background_segment);
}

ccl_device_forceinline void guiding_add_direct_contribution(IntegratorState state,
                                                            const float3 Le,
                                                            const float mis_weight)
{
  openpgl::cpp::SetDirectContribution(state->guiding.path_segment,
                                      openpgl::cpp::Vector3(Le[0], Le[1], Le[2]));
  openpgl::cpp::SetMiWeight(state->guiding.path_segment, mis_weight);
}

ccl_device_forceinline void guiding_add_scattered_contribution(IntegratorShadowState state,
                                                               const float3 Lo)
{
  if (state->shadow_path.path_segment) {
    openpgl::cpp::AddScatteredContribution(state->shadow_path.path_segment,
                                           openpgl::cpp::Vector3(Lo[0], Lo[1], Lo[2]));
  }
}

ccl_device_forceinline void guiding_set_continuation_probability(
    IntegratorState state, const float continuation_probability)
{
  if (state->guiding.path_segment) {
    openpgl::cpp::SetRussianRouletteProbability(state->guiding.path_segment,
                                                continuation_probability);
  }
}

#  endif

#  if defined(PATH_GUIDING_DEBUG_PASS)
/* Functions for writing guiding related debug information into separate frame buffers*/

ccl_device_forceinline void guiding_write_guiding_prob_buffer(const KernelGlobalsCPU *kg,
                                                              IntegratorStateCPU *state,
                                                              ccl_global float *ccl_restrict
                                                                  render_buffer)
{
  if (INTEGRATOR_STATE(state, path, bounce) == 0) {
    const uint32_t render_pixel_index = INTEGRATOR_STATE(state, path, render_pixel_index);
    const uint64_t render_buffer_offset = (uint64_t)render_pixel_index *
                                          kernel_data.film.pass_stride;
    ccl_global float *buffer = render_buffer + render_buffer_offset;
    float guiding_prob = state->guiding.surface_guiding_sampling_prob;
    if (kernel_data.film.pass_opgl_guiding_prob != PASS_UNUSED) {
      kernel_write_pass_float(buffer + kernel_data.film.pass_opgl_guiding_prob, guiding_prob);
    }
  }
}

ccl_device_forceinline void guiding_write_avg_roughness_buffer(const KernelGlobalsCPU *kg,
                                                               IntegratorStateCPU *state,
                                                               ccl_private const ShaderData *sd,
                                                               ccl_global float *ccl_restrict
                                                                   render_buffer)
{
  if (INTEGRATOR_STATE(state, path, bounce) == 0) {
    const uint32_t render_pixel_index = INTEGRATOR_STATE(state, path, render_pixel_index);
    const uint64_t render_buffer_offset = (uint64_t)render_pixel_index *
                                          kernel_data.film.pass_stride;
    ccl_global float *buffer = render_buffer + render_buffer_offset;
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

    if (kernel_data.film.pass_opgl_avg_roughness != PASS_UNUSED) {
      kernel_write_pass_float(buffer + kernel_data.film.pass_opgl_avg_roughness, avg_roughness);
    }
  }
}
#  endif
#endif

CCL_NAMESPACE_END
