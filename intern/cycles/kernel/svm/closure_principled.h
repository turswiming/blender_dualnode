/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

CCL_NAMESPACE_BEGIN

ccl_device void svm_node_closure_principled(KernelGlobals kg,
                                            ccl_private ShaderData *sd,
                                            ccl_private float *stack,
                                            uint4 data_node,
                                            float param1,
                                            float param2,
                                            float3 N,
                                            float mix_weight,
                                            int path_flag,
                                            int *offset)
{
  uint specular_offset, roughness_offset, specular_tint_offset, anisotropic_offset, sheen_offset,
      sheen_tint_offset, clearcoat_offset, clearcoat_roughness_offset, eta_offset,
      transmission_offset, anisotropic_rotation_offset, transmission_roughness_offset;
  uint4 data_node2 = read_node(kg, offset);

  float3 T = stack_load_float3(stack, data_node.y);
  svm_unpack_node_uchar4(data_node.z,
                         &specular_offset,
                         &roughness_offset,
                         &specular_tint_offset,
                         &anisotropic_offset);
  svm_unpack_node_uchar4(data_node.w,
                         &sheen_offset,
                         &sheen_tint_offset,
                         &clearcoat_offset,
                         &clearcoat_roughness_offset);
  svm_unpack_node_uchar4(data_node2.x,
                         &eta_offset,
                         &transmission_offset,
                         &anisotropic_rotation_offset,
                         &transmission_roughness_offset);

  // get Disney principled parameters
  float metallic = param1;
  float subsurface = param2;
  float specular = stack_load_float(stack, specular_offset);
  float roughness = stack_load_float(stack, roughness_offset);
  float specular_tint = stack_load_float(stack, specular_tint_offset);
  float anisotropic = stack_load_float(stack, anisotropic_offset);
  float sheen = stack_load_float(stack, sheen_offset);
  float sheen_tint = stack_load_float(stack, sheen_tint_offset);
  float clearcoat = stack_load_float(stack, clearcoat_offset);
  float clearcoat_roughness = stack_load_float(stack, clearcoat_roughness_offset);
  float transmission = stack_load_float(stack, transmission_offset);
  float anisotropic_rotation = stack_load_float(stack, anisotropic_rotation_offset);
  float transmission_roughness = stack_load_float(stack, transmission_roughness_offset);
  float eta = fmaxf(stack_load_float(stack, eta_offset), 1e-5f);

  ClosureType distribution = (ClosureType)data_node2.y;
  ClosureType subsurface_method = (ClosureType)data_node2.z;

  /* rotate tangent */
  if (anisotropic_rotation != 0.0f)
    T = rotate_around_axis(T, N, anisotropic_rotation * M_2PI_F);

  /* calculate ior */
  float ior = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

  // calculate fresnel for refraction
  float cosNO = dot(N, sd->I);
  float fresnel = fresnel_dielectric_cos(cosNO, ior);

  // calculate weights of the diffuse and specular part
  float diffuse_weight = (1.0f - saturatef(metallic)) * (1.0f - saturatef(transmission));

  float final_transmission = saturatef(transmission) * (1.0f - saturatef(metallic));
  float specular_weight = (1.0f - final_transmission);

  // get the base color
  uint4 data_base_color = read_node(kg, offset);
  float3 base_color = stack_valid(data_base_color.x) ?
                          stack_load_float3(stack, data_base_color.x) :
                          make_float3(__uint_as_float(data_base_color.y),
                                      __uint_as_float(data_base_color.z),
                                      __uint_as_float(data_base_color.w));

  // get the additional clearcoat normal and subsurface scattering radius
  uint4 data_cn_ssr = read_node(kg, offset);
  float3 clearcoat_normal = stack_valid(data_cn_ssr.x) ? stack_load_float3(stack, data_cn_ssr.x) :
                                                         sd->N;
  if (!(sd->type & PRIMITIVE_CURVE)) {
    clearcoat_normal = ensure_valid_reflection(sd->Ng, sd->I, clearcoat_normal);
  }
  float3 subsurface_radius = stack_valid(data_cn_ssr.y) ? stack_load_float3(stack, data_cn_ssr.y) :
                                                          make_float3(1.0f, 1.0f, 1.0f);
  float subsurface_ior = stack_valid(data_cn_ssr.z) ? stack_load_float(stack, data_cn_ssr.z) :
                                                      1.4f;
  float subsurface_anisotropy = stack_valid(data_cn_ssr.w) ?
                                    stack_load_float(stack, data_cn_ssr.w) :
                                    0.0f;

  // get the subsurface color
  uint4 data_subsurface_color = read_node(kg, offset);
  float3 subsurface_color = stack_valid(data_subsurface_color.x) ?
                                stack_load_float3(stack, data_subsurface_color.x) :
                                make_float3(__uint_as_float(data_subsurface_color.y),
                                            __uint_as_float(data_subsurface_color.z),
                                            __uint_as_float(data_subsurface_color.w));

  float3 weight = sd->svm_closure_weight * mix_weight;

#ifdef __SUBSURFACE__
  float3 mixed_ss_base_color = subsurface_color * subsurface + base_color * (1.0f - subsurface);
  float3 subsurf_weight = weight * mixed_ss_base_color * diffuse_weight;

  /* disable in case of diffuse ancestor, can't see it well then and
   * adds considerably noise due to probabilities of continuing path
   * getting lower and lower */
  if (path_flag & PATH_RAY_DIFFUSE_ANCESTOR) {
    subsurface = 0.0f;

    /* need to set the base color in this case such that the
     * rays get the correctly mixed color after transmitting
     * the object */
    base_color = mixed_ss_base_color;
  }

  /* diffuse */
  if (fabsf(average(mixed_ss_base_color)) > CLOSURE_WEIGHT_CUTOFF) {
    if (subsurface <= CLOSURE_WEIGHT_CUTOFF && diffuse_weight > CLOSURE_WEIGHT_CUTOFF) {
      float3 diff_weight = weight * base_color * diffuse_weight;

      ccl_private PrincipledDiffuseBsdf *bsdf = (ccl_private PrincipledDiffuseBsdf *)bsdf_alloc(
          sd, sizeof(PrincipledDiffuseBsdf), diff_weight);

      if (bsdf) {
        bsdf->N = N;
        bsdf->roughness = roughness;

        /* setup bsdf */
        sd->flag |= bsdf_principled_diffuse_setup(bsdf, PRINCIPLED_DIFFUSE_FULL);
      }
    }
    else if (subsurface > CLOSURE_WEIGHT_CUTOFF) {
      ccl_private Bssrdf *bssrdf = bssrdf_alloc(sd, subsurf_weight);

      if (bssrdf) {
        bssrdf->radius = subsurface_radius * subsurface;
        bssrdf->albedo = mixed_ss_base_color;
        bssrdf->N = N;
        bssrdf->roughness = roughness;

        /* Clamps protecting against bad/extreme and non physical values. */
        subsurface_ior = clamp(subsurface_ior, 1.01f, 3.8f);
        bssrdf->anisotropy = clamp(subsurface_anisotropy, 0.0f, 0.9f);

        /* setup bsdf */
        sd->flag |= bssrdf_setup(sd, bssrdf, subsurface_method, subsurface_ior);
      }
    }
  }
#else
  /* diffuse */
  if (diffuse_weight > CLOSURE_WEIGHT_CUTOFF) {
    float3 diff_weight = weight * base_color * diffuse_weight;

    ccl_private PrincipledDiffuseBsdf *bsdf = (ccl_private PrincipledDiffuseBsdf *)bsdf_alloc(
        sd, sizeof(PrincipledDiffuseBsdf), diff_weight);

    if (bsdf) {
      bsdf->N = N;
      bsdf->roughness = roughness;

      /* setup bsdf */
      sd->flag |= bsdf_principled_diffuse_setup(bsdf, PRINCIPLED_DIFFUSE_FULL);
    }
  }
#endif

  /* sheen */
  if (diffuse_weight > CLOSURE_WEIGHT_CUTOFF && sheen > CLOSURE_WEIGHT_CUTOFF) {
    float m_cdlum = linear_rgb_to_gray(kg, base_color);
    float3 m_ctint = m_cdlum > 0.0f ?
                         base_color / m_cdlum :
                         make_float3(1.0f, 1.0f, 1.0f);  // normalize lum. to isolate hue+sat

    /* color of the sheen component */
    float3 sheen_color = make_float3(1.0f, 1.0f, 1.0f) * (1.0f - sheen_tint) +
                         m_ctint * sheen_tint;

    float3 sheen_weight = weight * sheen * sheen_color * diffuse_weight;

    ccl_private PrincipledSheenBsdf *bsdf = (ccl_private PrincipledSheenBsdf *)bsdf_alloc(
        sd, sizeof(PrincipledSheenBsdf), sheen_weight);

    if (bsdf) {
      bsdf->N = N;

      /* setup bsdf */
      sd->flag |= bsdf_principled_sheen_setup(sd, bsdf);
    }
  }

  /* specular reflection */
#ifdef __CAUSTICS_TRICKS__
  if (kernel_data.integrator.caustics_reflective || (path_flag & PATH_RAY_DIFFUSE) == 0) {
#endif
    if (specular_weight > CLOSURE_WEIGHT_CUTOFF &&
        (specular > CLOSURE_WEIGHT_CUTOFF || metallic > CLOSURE_WEIGHT_CUTOFF)) {
      float3 spec_weight = weight * specular_weight;

      ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
          sd, sizeof(MicrofacetBsdf), spec_weight);
      ccl_private MicrofacetExtra *extra = (bsdf != NULL) ?
                                               (ccl_private MicrofacetExtra *)closure_alloc_extra(
                                                   sd, sizeof(MicrofacetExtra)) :
                                               NULL;

      if (bsdf && extra) {
        bsdf->N = N;
        bsdf->ior = (2.0f / (1.0f - safe_sqrtf(0.08f * specular))) - 1.0f;
        bsdf->T = T;
        bsdf->extra = extra;

        float aspect = safe_sqrtf(1.0f - anisotropic * 0.9f);
        float r2 = roughness * roughness;

        bsdf->alpha_x = r2 / aspect;
        bsdf->alpha_y = r2 * aspect;

        float m_cdlum = 0.3f * base_color.x + 0.6f * base_color.y +
                        0.1f * base_color.z;  // luminance approx.
        float3 m_ctint = m_cdlum > 0.0f ?
                             base_color / m_cdlum :
                             make_float3(1.0f, 1.0f, 1.0f);  // normalize lum. to isolate hue+sat
        float3 tmp_col = make_float3(1.0f, 1.0f, 1.0f) * (1.0f - specular_tint) +
                         m_ctint * specular_tint;

        bsdf->extra->cspec0 = (specular * 0.08f * tmp_col) * (1.0f - metallic) +
                              base_color * metallic;
        bsdf->extra->color = base_color;
        bsdf->extra->clearcoat = 0.0f;

        /* setup bsdf */
        if (distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID ||
            roughness <= 0.075f) /* use single-scatter GGX */
          sd->flag |= bsdf_microfacet_ggx_fresnel_setup(bsdf, sd);
        else /* use multi-scatter GGX */
          sd->flag |= bsdf_microfacet_multi_ggx_fresnel_setup(bsdf, sd);
      }
    }
#ifdef __CAUSTICS_TRICKS__
  }
#endif

  /* BSDF */
#ifdef __CAUSTICS_TRICKS__
  if (kernel_data.integrator.caustics_reflective || kernel_data.integrator.caustics_refractive ||
      (path_flag & PATH_RAY_DIFFUSE) == 0) {
#endif
    if (final_transmission > CLOSURE_WEIGHT_CUTOFF) {
      float3 glass_weight = weight * final_transmission;
      float3 cspec0 = base_color * specular_tint +
                      make_float3(1.0f, 1.0f, 1.0f) * (1.0f - specular_tint);

      if (roughness <= 5e-2f ||
          distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID) { /* use single-scatter GGX */
        float refl_roughness = roughness;

        /* reflection */
#ifdef __CAUSTICS_TRICKS__
        if (kernel_data.integrator.caustics_reflective || (path_flag & PATH_RAY_DIFFUSE) == 0)
#endif
        {
          ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
              sd, sizeof(MicrofacetBsdf), glass_weight * fresnel);
          ccl_private MicrofacetExtra *extra =
              (bsdf != NULL) ?
                  (ccl_private MicrofacetExtra *)closure_alloc_extra(sd, sizeof(MicrofacetExtra)) :
                  NULL;

          if (bsdf && extra) {
            bsdf->N = N;
            bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
            bsdf->extra = extra;

            bsdf->alpha_x = refl_roughness * refl_roughness;
            bsdf->alpha_y = refl_roughness * refl_roughness;
            bsdf->ior = ior;

            bsdf->extra->color = base_color;
            bsdf->extra->cspec0 = cspec0;
            bsdf->extra->clearcoat = 0.0f;

            /* setup bsdf */
            sd->flag |= bsdf_microfacet_ggx_fresnel_setup(bsdf, sd);
          }
        }

        /* refraction */
#ifdef __CAUSTICS_TRICKS__
        if (kernel_data.integrator.caustics_refractive || (path_flag & PATH_RAY_DIFFUSE) == 0)
#endif
        {
          /* This is to prevent MNEE from receiving a null BSDF. */
          float refraction_fresnel = fmaxf(0.0001f, 1.0f - fresnel);
          ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
              sd, sizeof(MicrofacetBsdf), base_color * glass_weight * refraction_fresnel);
          if (bsdf) {
            bsdf->N = N;
            bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
            bsdf->extra = NULL;

            if (distribution == CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID)
              transmission_roughness = 1.0f -
                                       (1.0f - refl_roughness) * (1.0f - transmission_roughness);
            else
              transmission_roughness = refl_roughness;

            bsdf->alpha_x = transmission_roughness * transmission_roughness;
            bsdf->alpha_y = transmission_roughness * transmission_roughness;
            bsdf->ior = ior;

            /* setup bsdf */
            sd->flag |= bsdf_microfacet_ggx_refraction_setup(bsdf);
          }
        }
      }
      else { /* use multi-scatter GGX */
        ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
            sd, sizeof(MicrofacetBsdf), glass_weight);
        ccl_private MicrofacetExtra *extra =
            (bsdf != NULL) ?
                (ccl_private MicrofacetExtra *)closure_alloc_extra(sd, sizeof(MicrofacetExtra)) :
                NULL;

        if (bsdf && extra) {
          bsdf->N = N;
          bsdf->extra = extra;
          bsdf->T = make_float3(0.0f, 0.0f, 0.0f);

          bsdf->alpha_x = roughness * roughness;
          bsdf->alpha_y = roughness * roughness;
          bsdf->ior = ior;

          bsdf->extra->color = base_color;
          bsdf->extra->cspec0 = cspec0;
          bsdf->extra->clearcoat = 0.0f;

          /* setup bsdf */
          sd->flag |= bsdf_microfacet_multi_ggx_glass_fresnel_setup(bsdf, sd);
        }
      }
    }
#ifdef __CAUSTICS_TRICKS__
  }
#endif

  /* clearcoat */
#ifdef __CAUSTICS_TRICKS__
  if (kernel_data.integrator.caustics_reflective || (path_flag & PATH_RAY_DIFFUSE) == 0) {
#endif
    if (clearcoat > CLOSURE_WEIGHT_CUTOFF) {
      ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
          sd, sizeof(MicrofacetBsdf), 0.25f * clearcoat * weight);

      if (bsdf) {
        bsdf->N = clearcoat_normal;
        bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
        bsdf->ior = 1.5f;
        bsdf->extra = NULL;

        bsdf->alpha_x = bsdf->alpha_y = sqr(clearcoat_roughness);

        /* setup bsdf */
        sd->flag |= bsdf_microfacet_ggx_clearcoat_setup(bsdf, sd);
      }
    }
#ifdef __CAUSTICS_TRICKS__
  }
#endif
}

CCL_NAMESPACE_END
