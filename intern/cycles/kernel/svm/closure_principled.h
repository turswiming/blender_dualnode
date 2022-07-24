/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

CCL_NAMESPACE_BEGIN

/* Principled v1 components */

ccl_device_inline void principled_v1_diffuse(ccl_private ShaderData *sd,
                                             float3 weight,
                                             float3 base_color,
                                             float diffuse_weight,
                                             float3 N,
                                             float roughness)
{
  if (diffuse_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  ccl_private PrincipledDiffuseBsdf *bsdf = (ccl_private PrincipledDiffuseBsdf *)bsdf_alloc(
      sd, sizeof(PrincipledDiffuseBsdf), diffuse_weight * base_color * weight);

  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->roughness = roughness;

  /* setup bsdf */
  sd->flag |= bsdf_principled_diffuse_setup(bsdf, PRINCIPLED_DIFFUSE_FULL);
}

ccl_device_inline void principled_v1_diffuse_sss(ccl_private ShaderData *sd,
                                                 ccl_private float *stack,
                                                 float3 weight,
                                                 int path_flag,
                                                 uint data_1,
                                                 uint data_2,
                                                 float3 base_color,
                                                 float diffuse_weight,
                                                 float3 N,
                                                 float roughness)
{
#ifdef __SUBSURFACE__
  uint method, subsurface_offset, aniso_offset, radius_offset;
  uint color_offset, ior_offset, dummy;
  svm_unpack_node_uchar4(data_1, &method, &subsurface_offset, &aniso_offset, &radius_offset);
  svm_unpack_node_uchar4(data_2, &color_offset, &ior_offset, &dummy, &dummy);

  float subsurface = stack_load_float(stack, subsurface_offset);
  float subsurface_anisotropy = stack_load_float(stack, aniso_offset);
  float subsurface_ior = stack_load_float(stack, ior_offset);
  float3 subsurface_color = stack_load_float3(stack, color_offset);
  float3 subsurface_radius = stack_load_float3(stack, radius_offset);

  float3 mixed_ss_base_color = mix(base_color, subsurface_color, subsurface);

  /* disable in case of diffuse ancestor, can't see it well then and
   * adds considerably noise due to probabilities of continuing path
   * getting lower and lower */
  if (path_flag & PATH_RAY_DIFFUSE_ANCESTOR) {
    subsurface = 0.0f;
  }

  /* diffuse */
  if (fabsf(average(mixed_ss_base_color)) > CLOSURE_WEIGHT_CUTOFF) {
    if (subsurface > CLOSURE_WEIGHT_CUTOFF) {
      float3 subsurf_weight = weight * mixed_ss_base_color * diffuse_weight;
      ccl_private Bssrdf *bssrdf = bssrdf_alloc(sd, subsurf_weight);

      if (bssrdf == NULL) {
        return;
      }

      bssrdf->radius = subsurface_radius * subsurface;
      bssrdf->albedo = mixed_ss_base_color;
      bssrdf->N = N;
      bssrdf->roughness = roughness;

      /* Clamps protecting against bad/extreme and non physical values. */
      subsurface_ior = clamp(subsurface_ior, 1.01f, 3.8f);
      bssrdf->anisotropy = clamp(subsurface_anisotropy, 0.0f, 0.9f);

      /* setup bsdf */
      sd->flag |= bssrdf_setup(sd, bssrdf, (ClosureType)method, subsurface_ior);
    }
    else {
      principled_v1_diffuse(sd, weight, mixed_ss_base_color, diffuse_weight, N, roughness);
    }
  }
#else
  /* diffuse */
  principled_v1_diffuse(sd, weight, base_color, diffuse_weight, N, roughness);
#endif
}

ccl_device_inline void principled_v1_specular(KernelGlobals kg,
                                              ccl_private ShaderData *sd,
                                              ccl_private float *stack,
                                              float3 weight,
                                              int path_flag,
                                              ClosureType distribution,
                                              uint data,
                                              float3 base_color,
                                              float3 N,
                                              float specular_weight,
                                              float metallic,
                                              float roughness,
                                              float specular_tint)
{
#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE)) {
    return;
  }
#endif

  uint specular_offset, aniso_offset, rotation_offset, tangent_offset;
  svm_unpack_node_uchar4(data, &specular_offset, &aniso_offset, &rotation_offset, &tangent_offset);

  float specular = stack_load_float(stack, specular_offset);

  if ((specular_weight <= CLOSURE_WEIGHT_CUTOFF) ||
      (specular + metallic <= CLOSURE_WEIGHT_CUTOFF)) {
    return;
  }

  float anisotropic = stack_load_float(stack, aniso_offset);
  float3 T = zero_float3();
  if (stack_valid(tangent_offset)) {
    T = stack_load_float3(stack, tangent_offset);
    T = rotate_around_axis(T, N, stack_load_float(stack, rotation_offset) * M_2PI_F);
  }

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), specular_weight * weight);
  if (bsdf == NULL) {
    return;
  }
  ccl_private MicrofacetExtra *extra = (ccl_private MicrofacetExtra *)closure_alloc_extra(
      sd, sizeof(MicrofacetExtra));
  if (extra == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->ior = (2.0f / (1.0f - safe_sqrtf(0.08f * specular))) - 1.0f;
  bsdf->T = T;
  bsdf->extra = extra;

  float aspect = safe_sqrtf(1.0f - anisotropic * 0.9f);

  bsdf->alpha_x = sqr(roughness) / aspect;
  bsdf->alpha_y = sqr(roughness) * aspect;

  // normalize lum. to isolate hue+sat
  float m_cdlum = linear_rgb_to_gray(kg, base_color);
  float3 m_ctint = m_cdlum > 0.0f ? base_color / m_cdlum : one_float3();
  float3 specular_color = mix(one_float3(), m_ctint, specular_tint);

  bsdf->extra->cspec0 = mix(specular * 0.08f * specular_color, base_color, metallic);
  bsdf->extra->color = base_color;

  /* setup bsdf */
  if (distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID ||
      roughness <= 0.075f) /* use single-scatter GGX */
    sd->flag |= bsdf_microfacet_ggx_fresnel_setup(bsdf, sd);
  else /* use multi-scatter GGX */
    sd->flag |= bsdf_microfacet_multi_ggx_fresnel_setup(kg, bsdf, sd);
}

ccl_device_inline void principled_v1_glass_refl(ccl_private ShaderData *sd,
                                                float3 weight,
                                                float3 base_color,
                                                float reflection_weight,
                                                float3 N,
                                                float roughness,
                                                float ior,
                                                float specular_tint)
{
  if (reflection_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), reflection_weight * weight);
  if (bsdf == NULL) {
    return;
  }

  ccl_private MicrofacetExtra *extra = (ccl_private MicrofacetExtra *)closure_alloc_extra(
      sd, sizeof(MicrofacetExtra));
  if (extra == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
  bsdf->extra = extra;

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);
  bsdf->ior = ior;

  bsdf->extra->color = base_color;
  bsdf->extra->cspec0 = mix(one_float3(), base_color, specular_tint);

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_ggx_fresnel_setup(bsdf, sd);
}

ccl_device_inline void principled_v1_glass_refr(ccl_private ShaderData *sd,
                                                float3 weight,
                                                float3 base_color,
                                                float refraction_weight,
                                                float3 N,
                                                float roughness,
                                                float ior)
{
  if (refraction_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), base_color * weight * refraction_weight);
  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
  bsdf->extra = NULL;

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);
  bsdf->ior = ior;

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_ggx_refraction_setup(bsdf);
}

ccl_device_inline void principled_v1_glass_single(KernelGlobals kg,
                                                  ccl_private ShaderData *sd,
                                                  ccl_private float *stack,
                                                  float3 weight,
                                                  int path_flag,
                                                  ClosureType distribution,
                                                  uint data,
                                                  float3 base_color,
                                                  float glass_weight,
                                                  float3 N,
                                                  float roughness,
                                                  float specular_tint)
{
  if (glass_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  uint transmission_roughness_offset, eta_offset, dummy;
  svm_unpack_node_uchar4(data, &eta_offset, &dummy, &dummy, &transmission_roughness_offset);
  float transmission_roughness = stack_load_float(stack, transmission_roughness_offset);
  float eta = fmaxf(stack_load_float(stack, eta_offset), 1e-5f);

  /* calculate ior */
  float ior = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

  // calculate fresnel for refraction
  float fresnel = fresnel_dielectric_cos(dot(N, sd->I), ior);

  /* reflection */
  float reflection_weight = glass_weight * fresnel;
#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE)) {
    reflection_weight = 0.0f;
  }
#endif
  principled_v1_glass_refl(
      sd, weight, base_color, reflection_weight, N, roughness, ior, specular_tint);

  /* refraction */
  /* TODO: MNEE ensured that this is always >0, is that correct?? */
  float refraction_weight = glass_weight * (1.0f - fresnel);
#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_refractive && (path_flag & PATH_RAY_DIFFUSE)) {
    refraction_weight = 0.0f;
  }
#endif
  if (distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID)
    transmission_roughness = 1.0f - (1.0f - roughness) * (1.0f - transmission_roughness);
  else
    transmission_roughness = roughness;
  principled_v1_glass_refr(
      sd, weight, base_color, refraction_weight, N, transmission_roughness, ior);
}

ccl_device_inline void principled_v1_glass_multi(KernelGlobals kg,
                                                 ccl_private ShaderData *sd,
                                                 ccl_private float *stack,
                                                 float3 weight,
                                                 int path_flag,
                                                 uint data,
                                                 float3 base_color,
                                                 float glass_weight,
                                                 float3 N,
                                                 float roughness,
                                                 float specular_tint)
{
#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_reflective && !kernel_data.integrator.caustics_refractive &&
      (path_flag & PATH_RAY_DIFFUSE)) {
    return;
  }
#endif

  if (glass_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  uint eta_offset, dummy;
  svm_unpack_node_uchar4(data, &eta_offset, &dummy, &dummy, &dummy);
  float eta = fmaxf(stack_load_float(stack, eta_offset), 1e-5f);

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), glass_weight * weight);
  if (bsdf == NULL) {
    return;
  }
  ccl_private MicrofacetExtra *extra = (ccl_private MicrofacetExtra *)closure_alloc_extra(
      sd, sizeof(MicrofacetExtra));
  if (extra == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->extra = extra;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);
  bsdf->ior = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

  bsdf->extra->color = base_color;
  bsdf->extra->cspec0 = mix(one_float3(), base_color, specular_tint);

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_multi_ggx_glass_fresnel_setup(kg, bsdf, sd);
}

ccl_device_inline void principled_v1_sheen(KernelGlobals kg,
                                           ccl_private ShaderData *sd,
                                           ccl_private float *stack,
                                           float3 weight,
                                           uint data,
                                           float3 base_color,
                                           float diffuse_weight,
                                           float3 N)
{
  uint sheen_offset, sheen_tint_offset, dummy;
  svm_unpack_node_uchar4(data, &dummy, &sheen_offset, &sheen_tint_offset, &dummy);
  float sheen = stack_load_float(stack, sheen_offset);

  float sheen_weight = diffuse_weight * sheen;
  if (sheen_weight <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  // normalize lum. to isolate hue+sat
  float m_cdlum = linear_rgb_to_gray(kg, base_color);
  float3 m_ctint = m_cdlum > 0.0f ? base_color / m_cdlum : one_float3();

  /* color of the sheen component */
  float sheen_tint = stack_load_float(stack, sheen_tint_offset);
  float3 sheen_color = mix(one_float3(), m_ctint, sheen_tint);

  ccl_private PrincipledSheenBsdf *bsdf = (ccl_private PrincipledSheenBsdf *)bsdf_alloc(
      sd, sizeof(PrincipledSheenBsdf), sheen_weight * sheen_color * weight);

  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;

  /* setup bsdf */
  sd->flag |= bsdf_principled_sheen_setup(sd, bsdf);
}

ccl_device_inline void principled_v1_clearcoat(KernelGlobals kg,
                                               ccl_private ShaderData *sd,
                                               ccl_private float *stack,
                                               float3 weight,
                                               int path_flag,
                                               uint data)
{
#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE)) {
    return;
  }
#endif

  uint clearcoat_offset, roughness_offset, normal_offset, dummy;
  svm_unpack_node_uchar4(data, &clearcoat_offset, &roughness_offset, &normal_offset, &dummy);
  float clearcoat = stack_load_float(stack, clearcoat_offset);

  if (clearcoat <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

  float roughness = stack_load_float(stack, roughness_offset);
  float3 N = stack_valid(normal_offset) ? stack_load_float3(stack, normal_offset) : sd->N;
  if (!(sd->type & PRIMITIVE_CURVE)) {
    N = ensure_valid_reflection(sd->Ng, sd->I, N);
  }

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), 0.25f * weight * clearcoat);

  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
  bsdf->ior = 1.5f;
  bsdf->extra = NULL;

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_ggx_clearcoat_setup(bsdf, sd);
}

/* Principled v2 components */

ccl_device_inline void principled_v2_diffuse_sss(ccl_private ShaderData *sd,
                                                 ccl_private float *stack,
                                                 float3 weight,
                                                 int path_flag,
                                                 uint data,
                                                 float3 base_color,
                                                 float ior,
                                                 float3 N)
{
  if (reduce_max(weight * base_color) <= CLOSURE_WEIGHT_CUTOFF) {
    return;
  }

#ifdef __SUBSURFACE__
  uint method, scale_offset, aniso_offset, radius_offset;
  svm_unpack_node_uchar4(data, &scale_offset, &aniso_offset, &radius_offset, &method);

  float aniso = stack_load_float(stack, aniso_offset);
  float3 radius = stack_load_float3(stack, radius_offset) * stack_load_float(stack, scale_offset);

  /* Fall back to diffuse if there has been a diffuse bounce before or the radius is too small. */
  if ((path_flag & PATH_RAY_DIFFUSE_ANCESTOR) == 0 && reduce_max(radius) > 1e-7f) {
    ccl_private Bssrdf *bssrdf = bssrdf_alloc(sd, base_color * weight);

    if (bssrdf == NULL) {
      return;
    }

    bssrdf->radius = radius;
    bssrdf->albedo = base_color;
    bssrdf->N = N;
    bssrdf->roughness = FLT_MAX;

    /* Clamps protecting against bad/extreme and non physical values. */
    bssrdf->anisotropy = clamp(aniso, 0.0f, 0.9f);

    /* setup bsdf */
    sd->flag |= bssrdf_setup(sd, bssrdf, (ClosureType)method, clamp(ior, 1.01f, 3.8f));

    return;
  }
#endif

  ccl_private DiffuseBsdf *bsdf = (ccl_private DiffuseBsdf *)bsdf_alloc(
      sd, sizeof(DiffuseBsdf), base_color * weight);

  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;

  /* setup bsdf */
  sd->flag |= bsdf_diffuse_setup(bsdf);
}

ccl_device_inline float3 principled_v2_clearcoat(KernelGlobals kg,
                                                 ccl_private ShaderData *sd,
                                                 ccl_private float *stack,
                                                 float3 weight,
                                                 int path_flag,
                                                 uint data)
{
  uint clearcoat_offset, roughness_offset, tint_offset, normal_offset;
  svm_unpack_node_uchar4(data, &clearcoat_offset, &roughness_offset, &tint_offset, &normal_offset);

  float3 N = stack_valid(normal_offset) ? stack_load_float3(stack, normal_offset) : sd->N;
  float3 tint = saturate(stack_load_float3(stack, tint_offset));
  if (tint != one_float3()) {
    /* Tint is normalized to perpendicular incidence.
     * Therefore, if we define the coating thickness as length 1, the length along the ray is
     * t = sqrt(1+tan^2(angle(N, I))) = sqrt(1+tan^2(acos(dotNI))) = 1 / dotNI.
     * From Beer's law, we have T = exp(-sigma_e * t).
     * Therefore, tint = exp(-sigma_e * 1) (per def.), so -sigma_e = log(tint).
     * From this, T = exp(log(tint) * t) = exp(log(tint)) ^ t = tint ^ t;
     *
     * Note that this is only an approximation - in particular, the exit path of the
     * light that bounces off the main layer is not accounted for.
     * Ideally, we should be attenuating by sqrt(tint) both ways. Currently, the model
     * is accurate for perfect mirrors since both entry and exit have the same cosNI there,
     * but if e.g. the base is diffuse and we look at it from a grazing angle, in reality
     * the cosNI of the exit bounce will be much higher on average, so the tint would be
     * less extreme.
     * TODO: Maybe account for this by setting
     * OD := 0.5*OD + 0.5*mix(1.59, OD, metallic * (1 - roughness)),
     * where 1.59 is the closest numerical fit for average optical depth on lambertian reflectors.
     * That way, mirrors preserve their look, but diffuse-ish objects have a more natural behavior.
     */
    float cosNI = dot(sd->I, N);
    /* Refract incoming direction into clearcoat material, which has a fixed IOR of 1.5.
     * TIR is no concern here since we're always coming from the outside. */
    float cosNT = sqrtf(1.0f - sqr(1.0f / 1.5f) * (1 - sqr(cosNI)));
    float optical_depth = 1.0f / cosNT;
    tint = pow(tint, optical_depth);
  }

  float clearcoat = saturatef(stack_load_float(stack, clearcoat_offset));
  if (clearcoat <= CLOSURE_WEIGHT_CUTOFF) {
    return tint;
  }

#ifdef __CAUSTICS_TRICKS__
  if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE)) {
    return tint;
  }
#endif

  float roughness = saturatef(stack_load_float(stack, roughness_offset));
  if (!(sd->type & PRIMITIVE_CURVE)) {
    N = ensure_valid_reflection(sd->Ng, sd->I, N);
  }

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), clearcoat * weight);

  if (bsdf == NULL) {
    return tint;
  }

  bsdf->N = N;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
  bsdf->ior = 1.5f;
  bsdf->extra = NULL;

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_ggx_clearcoat_v2_setup(kg, bsdf, sd);

  return tint * (1.0f - clearcoat_E(kg, dot(sd->I, N), roughness) * clearcoat);
}

ccl_device_inline float principled_v2_sheen(KernelGlobals kg,
                                            ccl_private ShaderData *sd,
                                            ccl_private float *stack,
                                            float3 weight,
                                            float3 N,
                                            uint data)
{
  uint sheen_offset, sheen_tint_offset, sheen_roughness_offset, dummy;
  svm_unpack_node_uchar4(data, &dummy, &sheen_offset, &sheen_tint_offset, &sheen_roughness_offset);

  float sheen = stack_load_float(stack, sheen_offset);
  if (sheen <= CLOSURE_WEIGHT_CUTOFF) {
    return 1.0f;
  }

  float roughness = stack_load_float(stack, sheen_roughness_offset);
  ccl_private PrincipledSheenBsdf *bsdf = (ccl_private PrincipledSheenBsdf *)bsdf_alloc(
      sd, sizeof(PrincipledSheenBsdf), sheen * weight);  // TODO include tint

  if (bsdf == NULL) {
    return 1.0f;
  }

  bsdf->N = N;
  bsdf->roughness = sqr(roughness);

  /* setup bsdf */
  sd->flag |= bsdf_principled_sheen_v2_setup(kg, sd, bsdf);

  return 1.0f - sheen * bsdf->avg_value;  // TODO include tint
}

ccl_device_inline float principled_v2_specular(KernelGlobals kg,
                                               ccl_private ShaderData *sd,
                                               ccl_private float *stack,
                                               float3 weight,
                                               float3 base_color,
                                               float roughness,
                                               float metallic,
                                               float ior,
                                               float transmission,
                                               float3 N,
                                               uint data1,
                                               uint data2)
{
  // TODO Handle caustics flag

  if (metallic + (1.0f - transmission) <= CLOSURE_WEIGHT_CUTOFF) {
    return 0.0f;
  }

  uint edge_offset, dummy;
  uint aniso_offset, rotation_offset, tangent_offset;
  svm_unpack_node_uchar4(data1, &dummy, &edge_offset, &dummy, &dummy);
  svm_unpack_node_uchar4(data2, &aniso_offset, &rotation_offset, &tangent_offset, &dummy);

  /* This function handles two specular components:
   * 1. Metallic: The overall energy is given by the metallic input
   * 2. Dielectric opaque: The overall energy is given by (1-metallic)*(1-transmission)
   * Both of these are handled by one closure, which adds up two Fresnel terms.
   * On top of that, there is also the transmissive component, that is handled by the glass code.
   */

  float anisotropic = stack_load_float(stack, aniso_offset);
  float aspect = safe_sqrtf(1.0f - anisotropic * 0.9f);
  float3 T = zero_float3();
  if (stack_valid(tangent_offset)) {
    T = stack_load_float3(stack, tangent_offset);
    T = rotate_around_axis(T, N, stack_load_float(stack, rotation_offset) * M_2PI_F);
  }
  float3 edge_color = stack_load_float3(stack, edge_offset);

  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), weight);
  if (bsdf == NULL) {
    return 0.0f;
  }
  ccl_private MicrofacetExtrav2 *extra = (ccl_private MicrofacetExtrav2 *)closure_alloc_extra(
      sd, sizeof(MicrofacetExtrav2));
  if (extra == NULL) {
    return 0.0f;
  }

  bsdf->N = N;
  bsdf->ior = ior;
  bsdf->T = T;
  bsdf->extra = (MicrofacetExtra *)extra;

  bsdf->alpha_x = sqr(roughness) / aspect;
  bsdf->alpha_y = sqr(roughness) * aspect;

  extra->metal_base = base_color;
  extra->metal_edge_factor = metallic_edge_factor(base_color, edge_color);

  float dielectric = (1.0f - metallic) * (1.0f - transmission);
  sd->flag |= bsdf_microfacet_ggx_fresnel_v2_setup(kg, bsdf, sd, metallic, dielectric);

  return microfacet_ggx_dielectric_E(kg, dot(sd->I, N), roughness, ior);
}

ccl_device_inline void principled_v2_glass(KernelGlobals kg,
                                           ccl_private ShaderData *sd,
                                           float3 weight,
                                           float transmission,
                                           float roughness,
                                           float ior,
                                           float3 N)
{
  ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
      sd, sizeof(MicrofacetBsdf), transmission * weight);
  if (bsdf == NULL) {
    return;
  }

  bsdf->N = N;
  bsdf->T = make_float3(0.0f, 0.0f, 0.0f);

  bsdf->alpha_x = bsdf->alpha_y = sqr(roughness);
  bsdf->ior = (sd->flag & SD_BACKFACING) ? 1.0f / ior : ior;

  /* setup bsdf */
  sd->flag |= bsdf_microfacet_multi_ggx_glass_setup(kg, bsdf, sd, one_float3());
}

ccl_device void svm_node_closure_principled_v2(KernelGlobals kg,
                                               ccl_private ShaderData *sd,
                                               ccl_private float *stack,
                                               uint4 node_1,
                                               uint4 node_2,
                                               float mix_weight,
                                               int path_flag,
                                               int *offset)
{
  float3 weight = sd->svm_closure_weight * mix_weight;

  /* Load shared parameter data. */
  uint base_color_offset, normal_offset, dummy;
  uint roughness_offset, metallic_offset, ior_offset, transmission_offset;
  svm_unpack_node_uchar4(node_1.y, &dummy, &base_color_offset, &normal_offset, &dummy);
  svm_unpack_node_uchar4(
      node_1.z, &roughness_offset, &metallic_offset, &ior_offset, &transmission_offset);

  float3 base_color = stack_load_float3(stack, base_color_offset);
  float3 N = stack_valid(normal_offset) ? stack_load_float3(stack, normal_offset) : sd->N;
  if (!(sd->type & PRIMITIVE_CURVE)) {
    N = ensure_valid_reflection(sd->Ng, sd->I, N);
  }
  float roughness = saturatef(stack_load_float(stack, roughness_offset));
  float metallic = saturatef(stack_load_float(stack, metallic_offset));
  float ior = fmaxf(stack_load_float(stack, ior_offset), 1e-5f);
  float transmission = saturatef(stack_load_float(stack, transmission_offset));

  weight *= principled_v2_clearcoat(kg, sd, stack, weight, path_flag, node_2.w);
  weight *= principled_v2_sheen(kg, sd, stack, weight, N, node_2.z);

  float dielectric_albedo = principled_v2_specular(kg,
                                                   sd,
                                                   stack,
                                                   weight,
                                                   base_color,
                                                   roughness,
                                                   metallic,
                                                   ior,
                                                   transmission,
                                                   N,
                                                   node_2.x,
                                                   node_2.y);
  weight *= 1.0f - metallic;

  principled_v2_glass(kg, sd, weight, transmission, roughness, ior, N);

  weight *= (1.0f - transmission) * (1.0f - dielectric_albedo);

  principled_v2_diffuse_sss(sd, stack, weight, path_flag, node_1.w, base_color, ior, N);
}

ccl_device void svm_node_closure_principled(KernelGlobals kg,
                                            ccl_private ShaderData *sd,
                                            ccl_private float *stack,
                                            uint4 node_1,
                                            uint4 node_2,
                                            float mix_weight,
                                            int path_flag,
                                            int *offset)
{
  /* Load distribution type. */
  uint packed_distribution, dummy;
  svm_unpack_node_uchar4(node_2.x, &dummy, &dummy, &dummy, &packed_distribution);
  ClosureType distribution = (ClosureType)packed_distribution;

  if (distribution == CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_ID) {
    svm_node_closure_principled_v2(kg, sd, stack, node_1, node_2, mix_weight, path_flag, offset);
    return;
  }

  /* Load shared parameter data. */
  uint base_color_offset, normal_offset;
  uint roughness_offset, metallic_offset, transmission_offset, specular_tint_offset;
  svm_unpack_node_uchar4(node_1.y, &dummy, &base_color_offset, &normal_offset, &dummy);
  svm_unpack_node_uchar4(
      node_1.z, &roughness_offset, &metallic_offset, &transmission_offset, &specular_tint_offset);

  float3 base_color = stack_load_float3(stack, base_color_offset);
  float3 N = stack_valid(normal_offset) ? stack_load_float3(stack, normal_offset) : sd->N;
  if (!(sd->type & PRIMITIVE_CURVE)) {
    N = ensure_valid_reflection(sd->Ng, sd->I, N);
  }
  float roughness = saturatef(stack_load_float(stack, roughness_offset));
  float metallic = saturatef(stack_load_float(stack, metallic_offset));
  float transmission = saturatef(stack_load_float(stack, transmission_offset));
  float specular_tint = saturatef(stack_load_float(stack, specular_tint_offset));

  /* Calculate closure mix weights.
   * The combined BSDF is mix(mix(Diffuse+SSS, Glass, transmission), Metal, metallic). */
  float diffuse_weight = (1.0f - metallic) * (1.0f - transmission);
  transmission *= 1.0f - metallic;
  /* NOTE: The mixing here is incorrect, the specular lobe should be metallic + (1 - transmission)
   * since it models both the metallic specular as well as the non-glass dielectric specular.
   * This only affects materials mixing diffuse, glass AND metal though. */
  float specular_weight = (1.0f - transmission);
  float3 weight = sd->svm_closure_weight * mix_weight;

  /* Diffuse and subsurface */
  principled_v1_diffuse_sss(
      sd, stack, weight, path_flag, node_1.w, node_2.x, base_color, diffuse_weight, N, roughness);

  /* sheen */
  principled_v1_sheen(kg, sd, stack, weight, node_2.z, base_color, diffuse_weight, N);

  /* specular reflection */
  principled_v1_specular(kg,
                         sd,
                         stack,
                         weight,
                         path_flag,
                         distribution,
                         node_2.y,
                         base_color,
                         N,
                         specular_weight,
                         metallic,
                         roughness,
                         specular_tint);

  /* glass */
  if (roughness <= 5e-2f || distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID) {
    principled_v1_glass_single(kg,
                               sd,
                               stack,
                               weight,
                               path_flag,
                               distribution,
                               node_2.z,
                               base_color,
                               transmission,
                               N,
                               roughness,
                               specular_tint);
  }
  else {
    principled_v1_glass_multi(kg,
                              sd,
                              stack,
                              weight,
                              path_flag,
                              node_2.z,
                              base_color,
                              transmission,
                              N,
                              roughness,
                              specular_tint);
  }

  /* clearcoat */
  principled_v1_clearcoat(kg, sd, stack, weight, path_flag, node_2.w);
}

CCL_NAMESPACE_END
