/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#ifdef __PRINCIPLED__
#  include "kernel/svm/closure_principled.h"
#endif

CCL_NAMESPACE_BEGIN

/* Closure Nodes */

ccl_device void svm_node_glass_setup(ccl_private ShaderData *sd,
                                     ccl_private MicrofacetBsdf *bsdf,
                                     int type,
                                     float eta,
                                     float roughness,
                                     bool refract)
{
  if (type == CLOSURE_BSDF_SHARP_GLASS_ID) {
    if (refract) {
      bsdf->alpha_y = 0.0f;
      bsdf->alpha_x = 0.0f;
      bsdf->ior = eta;
      sd->flag |= bsdf_refraction_setup(bsdf);
    }
    else {
      bsdf->alpha_y = 0.0f;
      bsdf->alpha_x = 0.0f;
      bsdf->ior = 0.0f;
      sd->flag |= bsdf_reflection_setup(bsdf);
    }
  }
  else if (type == CLOSURE_BSDF_MICROFACET_BECKMANN_GLASS_ID) {
    bsdf->alpha_x = roughness;
    bsdf->alpha_y = roughness;
    bsdf->ior = eta;

    if (refract)
      sd->flag |= bsdf_microfacet_beckmann_refraction_setup(bsdf);
    else
      sd->flag |= bsdf_microfacet_beckmann_setup(bsdf);
  }
  else {
    bsdf->alpha_x = roughness;
    bsdf->alpha_y = roughness;
    bsdf->ior = eta;

    if (refract)
      sd->flag |= bsdf_microfacet_ggx_refraction_setup(bsdf);
    else
      sd->flag |= bsdf_microfacet_ggx_setup(bsdf);
  }
}

ccl_device_inline int svm_node_closure_bsdf_skip(KernelGlobals kg, int offset, uint type)
{
  if (type == CLOSURE_BSDF_PRINCIPLED_ID) {
    /* Read all principled BSDF extra data to get the right offset. */
    read_node(kg, &offset);
    read_node(kg, &offset);
    read_node(kg, &offset);
    read_node(kg, &offset);
  }

  return offset;
}

template<uint node_feature_mask, ShaderType shader_type>
ccl_device_noinline int svm_node_closure_bsdf(KernelGlobals kg,
                                              ccl_private ShaderData *sd,
                                              ccl_private float *stack,
                                              uint4 node,
                                              uint32_t path_flag,
                                              int offset)
{
  uint type, param1_offset, param2_offset;

  uint mix_weight_offset;
  svm_unpack_node_uchar4(node.y, &type, &param1_offset, &param2_offset, &mix_weight_offset);
  float mix_weight = (stack_valid(mix_weight_offset) ? stack_load_float(stack, mix_weight_offset) :
                                                       1.0f);

  /* note we read this extra node before weight check, so offset is added */
  uint4 data_node = read_node(kg, &offset);

  /* Only compute BSDF for surfaces, transparent variable is shared with volume extinction. */
  IF_KERNEL_NODES_FEATURE(BSDF)
  {
    if ((shader_type != SHADER_TYPE_SURFACE) || mix_weight == 0.0f) {
      return svm_node_closure_bsdf_skip(kg, offset, type);
    }
  }
  else
  {
    return svm_node_closure_bsdf_skip(kg, offset, type);
  }

  float3 N = stack_valid(data_node.x) ? stack_load_float3(stack, data_node.x) : sd->N;
  if (!(sd->type & PRIMITIVE_CURVE)) {
    N = ensure_valid_reflection(sd->Ng, sd->I, N);
  }

  float param1 = (stack_valid(param1_offset)) ? stack_load_float(stack, param1_offset) :
                                                __uint_as_float(node.z);
  float param2 = (stack_valid(param2_offset)) ? stack_load_float(stack, param2_offset) :
                                                __uint_as_float(node.w);

  switch (type) {
#ifdef __PRINCIPLED__
    case CLOSURE_BSDF_PRINCIPLED_ID: {
      svm_node_closure_principled(
          kg, sd, stack, data_node, param1, param2, N, mix_weight, path_flag, &offset);
      break;
    }
#endif /* __PRINCIPLED__ */
    case CLOSURE_BSDF_DIFFUSE_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private OrenNayarBsdf *bsdf = (ccl_private OrenNayarBsdf *)bsdf_alloc(
          sd, sizeof(OrenNayarBsdf), weight);

      if (bsdf) {
        bsdf->N = N;

        float roughness = param1;

        if (roughness == 0.0f) {
          sd->flag |= bsdf_diffuse_setup((ccl_private DiffuseBsdf *)bsdf);
        }
        else {
          bsdf->roughness = roughness;
          sd->flag |= bsdf_oren_nayar_setup(bsdf);
        }
      }
      break;
    }
    case CLOSURE_BSDF_TRANSLUCENT_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private DiffuseBsdf *bsdf = (ccl_private DiffuseBsdf *)bsdf_alloc(
          sd, sizeof(DiffuseBsdf), weight);

      if (bsdf) {
        bsdf->N = N;
        sd->flag |= bsdf_translucent_setup(bsdf);
      }
      break;
    }
    case CLOSURE_BSDF_TRANSPARENT_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      bsdf_transparent_setup(sd, weight, path_flag);
      break;
    }
    case CLOSURE_BSDF_REFLECTION_ID:
    case CLOSURE_BSDF_MICROFACET_GGX_ID:
    case CLOSURE_BSDF_MICROFACET_BECKMANN_ID:
    case CLOSURE_BSDF_ASHIKHMIN_SHIRLEY_ID:
    case CLOSURE_BSDF_MICROFACET_MULTI_GGX_ID: {
#ifdef __CAUSTICS_TRICKS__
      if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE))
        break;
#endif
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
          sd, sizeof(MicrofacetBsdf), weight);

      if (!bsdf) {
        break;
      }

      float roughness = sqr(param1);

      bsdf->N = N;
      bsdf->ior = 0.0f;
      bsdf->extra = NULL;

      if (data_node.y == SVM_STACK_INVALID) {
        bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
        bsdf->alpha_x = roughness;
        bsdf->alpha_y = roughness;
      }
      else {
        bsdf->T = stack_load_float3(stack, data_node.y);

        /* rotate tangent */
        float rotation = stack_load_float(stack, data_node.z);
        if (rotation != 0.0f)
          bsdf->T = rotate_around_axis(bsdf->T, bsdf->N, rotation * M_2PI_F);

        /* compute roughness */
        float anisotropy = clamp(param2, -0.99f, 0.99f);
        if (anisotropy < 0.0f) {
          bsdf->alpha_x = roughness / (1.0f + anisotropy);
          bsdf->alpha_y = roughness * (1.0f + anisotropy);
        }
        else {
          bsdf->alpha_x = roughness * (1.0f - anisotropy);
          bsdf->alpha_y = roughness / (1.0f - anisotropy);
        }
      }

      /* setup bsdf */
      if (type == CLOSURE_BSDF_REFLECTION_ID)
        sd->flag |= bsdf_reflection_setup(bsdf);
      else if (type == CLOSURE_BSDF_MICROFACET_BECKMANN_ID)
        sd->flag |= bsdf_microfacet_beckmann_setup(bsdf);
      else if (type == CLOSURE_BSDF_MICROFACET_GGX_ID)
        sd->flag |= bsdf_microfacet_ggx_setup(bsdf);
      else if (type == CLOSURE_BSDF_MICROFACET_MULTI_GGX_ID) {
        kernel_assert(stack_valid(data_node.w));
        float3 color = stack_load_float3(stack, data_node.w);
        sd->flag |= bsdf_microfacet_multi_ggx_setup(bsdf, sd, color);
      }
      else {
        sd->flag |= bsdf_ashikhmin_shirley_setup(bsdf);
      }

      break;
    }
    case CLOSURE_BSDF_REFRACTION_ID:
    case CLOSURE_BSDF_MICROFACET_GGX_REFRACTION_ID:
    case CLOSURE_BSDF_MICROFACET_BECKMANN_REFRACTION_ID: {
#ifdef __CAUSTICS_TRICKS__
      if (!kernel_data.integrator.caustics_refractive && (path_flag & PATH_RAY_DIFFUSE))
        break;
#endif
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
          sd, sizeof(MicrofacetBsdf), weight);

      if (bsdf) {
        bsdf->N = N;
        bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
        bsdf->extra = NULL;

        float eta = fmaxf(param2, 1e-5f);
        eta = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

        /* setup bsdf */
        if (type == CLOSURE_BSDF_REFRACTION_ID) {
          bsdf->alpha_x = 0.0f;
          bsdf->alpha_y = 0.0f;
          bsdf->ior = eta;

          sd->flag |= bsdf_refraction_setup(bsdf);
        }
        else {
          float roughness = sqr(param1);
          bsdf->alpha_x = roughness;
          bsdf->alpha_y = roughness;
          bsdf->ior = eta;

          if (type == CLOSURE_BSDF_MICROFACET_BECKMANN_REFRACTION_ID)
            sd->flag |= bsdf_microfacet_beckmann_refraction_setup(bsdf);
          else
            sd->flag |= bsdf_microfacet_ggx_refraction_setup(bsdf);
        }
      }

      break;
    }
    case CLOSURE_BSDF_SHARP_GLASS_ID:
    case CLOSURE_BSDF_MICROFACET_GGX_GLASS_ID:
    case CLOSURE_BSDF_MICROFACET_BECKMANN_GLASS_ID: {
#ifdef __CAUSTICS_TRICKS__
      if (!kernel_data.integrator.caustics_reflective &&
          !kernel_data.integrator.caustics_refractive && (path_flag & PATH_RAY_DIFFUSE)) {
        break;
      }
#endif
      float3 weight = sd->svm_closure_weight * mix_weight;

      /* index of refraction */
      float eta = fmaxf(param2, 1e-5f);
      eta = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

      /* fresnel */
      float cosNO = dot(N, sd->I);
      float fresnel = fresnel_dielectric_cos(cosNO, eta);
      float roughness = sqr(param1);

      /* reflection */
#ifdef __CAUSTICS_TRICKS__
      if (kernel_data.integrator.caustics_reflective || (path_flag & PATH_RAY_DIFFUSE) == 0)
#endif
      {
        ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
            sd, sizeof(MicrofacetBsdf), weight * fresnel);

        if (bsdf) {
          bsdf->N = N;
          bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
          bsdf->extra = NULL;
          svm_node_glass_setup(sd, bsdf, type, eta, roughness, false);
        }
      }

      /* refraction */
#ifdef __CAUSTICS_TRICKS__
      if (kernel_data.integrator.caustics_refractive || (path_flag & PATH_RAY_DIFFUSE) == 0)
#endif
      {
        /* This is to prevent MNEE from receiving a null BSDF.
         * TODO: Doesn't this always enable the closure? */
        float refraction_fresnel = fmaxf(0.0001f, 1.0f - fresnel);
        ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
            sd, sizeof(MicrofacetBsdf), weight * refraction_fresnel);

        if (bsdf) {
          bsdf->N = N;
          bsdf->T = make_float3(0.0f, 0.0f, 0.0f);
          bsdf->extra = NULL;
          svm_node_glass_setup(sd, bsdf, type, eta, roughness, true);
        }
      }

      break;
    }
    case CLOSURE_BSDF_MICROFACET_MULTI_GGX_GLASS_ID: {
#ifdef __CAUSTICS_TRICKS__
      if (!kernel_data.integrator.caustics_reflective &&
          !kernel_data.integrator.caustics_refractive && (path_flag & PATH_RAY_DIFFUSE))
        break;
#endif
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private MicrofacetBsdf *bsdf = (ccl_private MicrofacetBsdf *)bsdf_alloc(
          sd, sizeof(MicrofacetBsdf), weight);
      if (!bsdf) {
        break;
      }

      /* TODO: Detect sharp, fallback. */

      bsdf->N = N;
      bsdf->extra = NULL;
      bsdf->T = make_float3(0.0f, 0.0f, 0.0f);

      float roughness = sqr(param1);
      bsdf->alpha_x = roughness;
      bsdf->alpha_y = roughness;
      float eta = fmaxf(param2, 1e-5f);
      bsdf->ior = (sd->flag & SD_BACKFACING) ? 1.0f / eta : eta;

      kernel_assert(stack_valid(data_node.z));
      float3 color = stack_load_float3(stack, data_node.z);

      /* setup bsdf */
      sd->flag |= bsdf_microfacet_multi_ggx_glass_setup(bsdf, sd, color);
      break;
    }
    case CLOSURE_BSDF_ASHIKHMIN_VELVET_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private VelvetBsdf *bsdf = (ccl_private VelvetBsdf *)bsdf_alloc(
          sd, sizeof(VelvetBsdf), weight);

      if (bsdf) {
        bsdf->N = N;

        bsdf->sigma = saturatef(param1);
        sd->flag |= bsdf_ashikhmin_velvet_setup(bsdf);
      }
      break;
    }
    case CLOSURE_BSDF_GLOSSY_TOON_ID:
#ifdef __CAUSTICS_TRICKS__
      if (!kernel_data.integrator.caustics_reflective && (path_flag & PATH_RAY_DIFFUSE))
        break;
      ATTR_FALLTHROUGH;
#endif
    case CLOSURE_BSDF_DIFFUSE_TOON_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private ToonBsdf *bsdf = (ccl_private ToonBsdf *)bsdf_alloc(
          sd, sizeof(ToonBsdf), weight);

      if (bsdf) {
        bsdf->N = N;
        bsdf->size = param1;
        bsdf->smooth = param2;

        if (type == CLOSURE_BSDF_DIFFUSE_TOON_ID)
          sd->flag |= bsdf_diffuse_toon_setup(bsdf);
        else
          sd->flag |= bsdf_glossy_toon_setup(bsdf);
      }
      break;
    }
#ifdef __HAIR__
    case CLOSURE_BSDF_HAIR_PRINCIPLED_ID: {
      uint4 data_node2 = read_node(kg, &offset);
      uint4 data_node3 = read_node(kg, &offset);
      uint4 data_node4 = read_node(kg, &offset);

      float3 weight = sd->svm_closure_weight * mix_weight;

      uint offset_ofs, ior_ofs, color_ofs, parametrization;
      svm_unpack_node_uchar4(data_node.y, &offset_ofs, &ior_ofs, &color_ofs, &parametrization);
      float alpha = stack_load_float_default(stack, offset_ofs, data_node.z);
      float ior = stack_load_float_default(stack, ior_ofs, data_node.w);

      uint coat_ofs, melanin_ofs, melanin_redness_ofs, absorption_coefficient_ofs;
      svm_unpack_node_uchar4(data_node2.x,
                             &coat_ofs,
                             &melanin_ofs,
                             &melanin_redness_ofs,
                             &absorption_coefficient_ofs);

      uint tint_ofs, random_ofs, random_color_ofs, random_roughness_ofs;
      svm_unpack_node_uchar4(
          data_node3.x, &tint_ofs, &random_ofs, &random_color_ofs, &random_roughness_ofs);

      const AttributeDescriptor attr_descr_random = find_attribute(kg, sd, data_node4.y);
      float random = 0.0f;
      if (attr_descr_random.offset != ATTR_STD_NOT_FOUND) {
        random = primitive_surface_attribute_float(kg, sd, attr_descr_random, NULL, NULL);
      }
      else {
        random = stack_load_float_default(stack, random_ofs, data_node3.y);
      }

      ccl_private PrincipledHairBSDF *bsdf = (ccl_private PrincipledHairBSDF *)bsdf_alloc(
          sd, sizeof(PrincipledHairBSDF), weight);
      if (bsdf) {
        ccl_private PrincipledHairExtra *extra = (ccl_private PrincipledHairExtra *)
            closure_alloc_extra(sd, sizeof(PrincipledHairExtra));

        if (!extra)
          break;

        /* Random factors range: [-randomization/2, +randomization/2]. */
        float random_roughness = stack_load_float_default(
            stack, random_roughness_ofs, data_node3.w);
        float factor_random_roughness = 1.0f + 2.0f * (random - 0.5f) * random_roughness;
        float roughness = param1 * factor_random_roughness;
        float radial_roughness = param2 * factor_random_roughness;

        /* Remap Coat value to [0, 100]% of Roughness. */
        float coat = stack_load_float_default(stack, coat_ofs, data_node2.y);
        float m0_roughness = 1.0f - clamp(coat, 0.0f, 1.0f);

        bsdf->N = N;
        bsdf->v = roughness;
        bsdf->s = radial_roughness;
        bsdf->m0_roughness = m0_roughness;
        bsdf->alpha = alpha;
        bsdf->eta = ior;
        bsdf->extra = extra;

        switch (parametrization) {
          case NODE_PRINCIPLED_HAIR_DIRECT_ABSORPTION: {
            float3 absorption_coefficient = stack_load_float3(stack, absorption_coefficient_ofs);
            bsdf->sigma = absorption_coefficient;
            break;
          }
          case NODE_PRINCIPLED_HAIR_PIGMENT_CONCENTRATION: {
            float melanin = stack_load_float_default(stack, melanin_ofs, data_node2.z);
            float melanin_redness = stack_load_float_default(
                stack, melanin_redness_ofs, data_node2.w);

            /* Randomize melanin. */
            float random_color = stack_load_float_default(stack, random_color_ofs, data_node3.z);
            random_color = clamp(random_color, 0.0f, 1.0f);
            float factor_random_color = 1.0f + 2.0f * (random - 0.5f) * random_color;
            melanin *= factor_random_color;

            /* Map melanin 0..inf from more perceptually linear 0..1. */
            melanin = -logf(fmaxf(1.0f - melanin, 0.0001f));

            /* Benedikt Bitterli's melanin ratio remapping. */
            float eumelanin = melanin * (1.0f - melanin_redness);
            float pheomelanin = melanin * melanin_redness;
            float3 melanin_sigma = bsdf_principled_hair_sigma_from_concentration(eumelanin,
                                                                                 pheomelanin);

            /* Optional tint. */
            float3 tint = stack_load_float3(stack, tint_ofs);
            float3 tint_sigma = bsdf_principled_hair_sigma_from_reflectance(tint,
                                                                            radial_roughness);

            bsdf->sigma = melanin_sigma + tint_sigma;
            break;
          }
          case NODE_PRINCIPLED_HAIR_REFLECTANCE: {
            float3 color = stack_load_float3(stack, color_ofs);
            bsdf->sigma = bsdf_principled_hair_sigma_from_reflectance(color, radial_roughness);
            break;
          }
          default: {
            /* Fallback to brownish hair, same as defaults for melanin. */
            kernel_assert(!"Invalid Principled Hair parametrization!");
            bsdf->sigma = bsdf_principled_hair_sigma_from_concentration(0.0f, 0.8054375f);
            break;
          }
        }

        sd->flag |= bsdf_principled_hair_setup(sd, bsdf);
      }
      break;
    }
    case CLOSURE_BSDF_HAIR_REFLECTION_ID:
    case CLOSURE_BSDF_HAIR_TRANSMISSION_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;

      ccl_private HairBsdf *bsdf = (ccl_private HairBsdf *)bsdf_alloc(
          sd, sizeof(HairBsdf), weight);

      if (bsdf) {
        bsdf->N = N;
        bsdf->roughness1 = param1;
        bsdf->roughness2 = param2;
        bsdf->offset = -stack_load_float(stack, data_node.z);

        if (stack_valid(data_node.y)) {
          bsdf->T = normalize(stack_load_float3(stack, data_node.y));
        }
        else if (!(sd->type & PRIMITIVE_CURVE)) {
          bsdf->T = normalize(sd->dPdv);
          bsdf->offset = 0.0f;
        }
        else
          bsdf->T = normalize(sd->dPdu);

        if (type == CLOSURE_BSDF_HAIR_REFLECTION_ID) {
          sd->flag |= bsdf_hair_reflection_setup(bsdf);
        }
        else {
          sd->flag |= bsdf_hair_transmission_setup(bsdf);
        }
      }

      break;
    }
#endif /* __HAIR__ */

#ifdef __SUBSURFACE__
    case CLOSURE_BSSRDF_BURLEY_ID:
    case CLOSURE_BSSRDF_RANDOM_WALK_ID:
    case CLOSURE_BSSRDF_RANDOM_WALK_FIXED_RADIUS_ID: {
      float3 weight = sd->svm_closure_weight * mix_weight;
      ccl_private Bssrdf *bssrdf = bssrdf_alloc(sd, weight);

      if (bssrdf) {
        /* disable in case of diffuse ancestor, can't see it well then and
         * adds considerably noise due to probabilities of continuing path
         * getting lower and lower */
        if (path_flag & PATH_RAY_DIFFUSE_ANCESTOR)
          param1 = 0.0f;

        bssrdf->radius = stack_load_float3(stack, data_node.z) * param1;
        bssrdf->albedo = sd->svm_closure_weight;
        bssrdf->N = N;
        bssrdf->roughness = FLT_MAX;

        const float subsurface_ior = clamp(param2, 1.01f, 3.8f);
        const float subsurface_anisotropy = stack_load_float(stack, data_node.w);
        bssrdf->anisotropy = clamp(subsurface_anisotropy, 0.0f, 0.9f);

        sd->flag |= bssrdf_setup(sd, bssrdf, (ClosureType)type, subsurface_ior);
      }

      break;
    }
#endif
    default:
      break;
  }

  return offset;
}

template<ShaderType shader_type>
ccl_device_noinline void svm_node_closure_volume(KernelGlobals kg,
                                                 ccl_private ShaderData *sd,
                                                 ccl_private float *stack,
                                                 uint4 node)
{
#ifdef __VOLUME__
  /* Only sum extinction for volumes, variable is shared with surface transparency. */
  if (shader_type != SHADER_TYPE_VOLUME) {
    return;
  }

  uint type, density_offset, anisotropy_offset;

  uint mix_weight_offset;
  svm_unpack_node_uchar4(node.y, &type, &density_offset, &anisotropy_offset, &mix_weight_offset);
  float mix_weight = (stack_valid(mix_weight_offset) ? stack_load_float(stack, mix_weight_offset) :
                                                       1.0f);

  if (mix_weight == 0.0f) {
    return;
  }

  float density = (stack_valid(density_offset)) ? stack_load_float(stack, density_offset) :
                                                  __uint_as_float(node.z);
  density = mix_weight * fmaxf(density, 0.0f);

  /* Compute scattering coefficient. */
  float3 weight = sd->svm_closure_weight;

  if (type == CLOSURE_VOLUME_ABSORPTION_ID) {
    weight = make_float3(1.0f, 1.0f, 1.0f) - weight;
  }

  weight *= density;

  /* Add closure for volume scattering. */
  if (type == CLOSURE_VOLUME_HENYEY_GREENSTEIN_ID) {
    ccl_private HenyeyGreensteinVolume *volume = (ccl_private HenyeyGreensteinVolume *)bsdf_alloc(
        sd, sizeof(HenyeyGreensteinVolume), weight);

    if (volume) {
      float anisotropy = (stack_valid(anisotropy_offset)) ?
                             stack_load_float(stack, anisotropy_offset) :
                             __uint_as_float(node.w);
      volume->g = anisotropy; /* g */
      sd->flag |= volume_henyey_greenstein_setup(volume);
    }
  }

  /* Sum total extinction weight. */
  volume_extinction_setup(sd, weight);
#endif
}

template<ShaderType shader_type>
ccl_device_noinline int svm_node_principled_volume(KernelGlobals kg,
                                                   ccl_private ShaderData *sd,
                                                   ccl_private float *stack,
                                                   uint4 node,
                                                   uint32_t path_flag,
                                                   int offset)
{
#ifdef __VOLUME__
  uint4 value_node = read_node(kg, &offset);
  uint4 attr_node = read_node(kg, &offset);

  /* Only sum extinction for volumes, variable is shared with surface transparency. */
  if (shader_type != SHADER_TYPE_VOLUME) {
    return offset;
  }

  uint density_offset, anisotropy_offset, absorption_color_offset, mix_weight_offset;
  svm_unpack_node_uchar4(
      node.y, &density_offset, &anisotropy_offset, &absorption_color_offset, &mix_weight_offset);
  float mix_weight = (stack_valid(mix_weight_offset) ? stack_load_float(stack, mix_weight_offset) :
                                                       1.0f);

  if (mix_weight == 0.0f) {
    return offset;
  }

  /* Compute density. */
  float primitive_density = 1.0f;
  float density = (stack_valid(density_offset)) ? stack_load_float(stack, density_offset) :
                                                  __uint_as_float(value_node.x);
  density = mix_weight * fmaxf(density, 0.0f);

  if (density > CLOSURE_WEIGHT_CUTOFF) {
    /* Density and color attribute lookup if available. */
    const AttributeDescriptor attr_density = find_attribute(kg, sd, attr_node.x);
    if (attr_density.offset != ATTR_STD_NOT_FOUND) {
      primitive_density = primitive_volume_attribute_float(kg, sd, attr_density);
      density = fmaxf(density * primitive_density, 0.0f);
    }
  }

  if (density > CLOSURE_WEIGHT_CUTOFF) {
    /* Compute scattering color. */
    float3 color = sd->svm_closure_weight;

    const AttributeDescriptor attr_color = find_attribute(kg, sd, attr_node.y);
    if (attr_color.offset != ATTR_STD_NOT_FOUND) {
      color *= primitive_volume_attribute_float3(kg, sd, attr_color);
    }

    /* Add closure for volume scattering. */
    ccl_private HenyeyGreensteinVolume *volume = (ccl_private HenyeyGreensteinVolume *)bsdf_alloc(
        sd, sizeof(HenyeyGreensteinVolume), color * density);
    if (volume) {
      float anisotropy = (stack_valid(anisotropy_offset)) ?
                             stack_load_float(stack, anisotropy_offset) :
                             __uint_as_float(value_node.y);
      volume->g = anisotropy;
      sd->flag |= volume_henyey_greenstein_setup(volume);
    }

    /* Add extinction weight. */
    float3 zero = make_float3(0.0f, 0.0f, 0.0f);
    float3 one = make_float3(1.0f, 1.0f, 1.0f);
    float3 absorption_color = max(sqrt(stack_load_float3(stack, absorption_color_offset)), zero);
    float3 absorption = max(one - color, zero) * max(one - absorption_color, zero);
    volume_extinction_setup(sd, (color + absorption) * density);
  }

  /* Compute emission. */
  if (path_flag & PATH_RAY_SHADOW) {
    /* Don't need emission for shadows. */
    return offset;
  }

  uint emission_offset, emission_color_offset, blackbody_offset, temperature_offset;
  svm_unpack_node_uchar4(
      node.z, &emission_offset, &emission_color_offset, &blackbody_offset, &temperature_offset);
  float emission = (stack_valid(emission_offset)) ? stack_load_float(stack, emission_offset) :
                                                    __uint_as_float(value_node.z);
  float blackbody = (stack_valid(blackbody_offset)) ? stack_load_float(stack, blackbody_offset) :
                                                      __uint_as_float(value_node.w);

  if (emission > CLOSURE_WEIGHT_CUTOFF) {
    float3 emission_color = stack_load_float3(stack, emission_color_offset);
    emission_setup(sd, emission * emission_color);
  }

  if (blackbody > CLOSURE_WEIGHT_CUTOFF) {
    float T = stack_load_float(stack, temperature_offset);

    /* Add flame temperature from attribute if available. */
    const AttributeDescriptor attr_temperature = find_attribute(kg, sd, attr_node.z);
    if (attr_temperature.offset != ATTR_STD_NOT_FOUND) {
      float temperature = primitive_volume_attribute_float(kg, sd, attr_temperature);
      T *= fmaxf(temperature, 0.0f);
    }

    T = fmaxf(T, 0.0f);

    /* Stefan-Boltzmann law. */
    float T4 = sqr(sqr(T));
    float sigma = 5.670373e-8f * 1e-6f / M_PI_F;
    float intensity = sigma * mix(1.0f, T4, blackbody);

    if (intensity > CLOSURE_WEIGHT_CUTOFF) {
      float3 blackbody_tint = stack_load_float3(stack, node.w);
      float3 bb = blackbody_tint * intensity *
                  rec709_to_rgb(kg, svm_math_blackbody_color_rec709(T));
      emission_setup(sd, bb);
    }
  }
#endif
  return offset;
}

ccl_device_noinline void svm_node_closure_emission(ccl_private ShaderData *sd,
                                                   ccl_private float *stack,
                                                   uint4 node)
{
  uint mix_weight_offset = node.y;
  float3 weight = sd->svm_closure_weight;

  if (stack_valid(mix_weight_offset)) {
    float mix_weight = stack_load_float(stack, mix_weight_offset);

    if (mix_weight == 0.0f)
      return;

    weight *= mix_weight;
  }

  emission_setup(sd, weight);
}

ccl_device_noinline void svm_node_closure_background(ccl_private ShaderData *sd,
                                                     ccl_private float *stack,
                                                     uint4 node)
{
  uint mix_weight_offset = node.y;
  float3 weight = sd->svm_closure_weight;

  if (stack_valid(mix_weight_offset)) {
    float mix_weight = stack_load_float(stack, mix_weight_offset);

    if (mix_weight == 0.0f)
      return;

    weight *= mix_weight;
  }

  background_setup(sd, weight);
}

ccl_device_noinline void svm_node_closure_holdout(ccl_private ShaderData *sd,
                                                  ccl_private float *stack,
                                                  uint4 node)
{
  uint mix_weight_offset = node.y;

  if (stack_valid(mix_weight_offset)) {
    float mix_weight = stack_load_float(stack, mix_weight_offset);

    if (mix_weight == 0.0f)
      return;

    closure_alloc(
        sd, sizeof(ShaderClosure), CLOSURE_HOLDOUT_ID, sd->svm_closure_weight * mix_weight);
  }
  else
    closure_alloc(sd, sizeof(ShaderClosure), CLOSURE_HOLDOUT_ID, sd->svm_closure_weight);

  sd->flag |= SD_HOLDOUT;
}

/* Closure Nodes */

ccl_device_inline void svm_node_closure_store_weight(ccl_private ShaderData *sd, float3 weight)
{
  sd->svm_closure_weight = weight;
}

ccl_device void svm_node_closure_set_weight(ccl_private ShaderData *sd, uint r, uint g, uint b)
{
  float3 weight = make_float3(__uint_as_float(r), __uint_as_float(g), __uint_as_float(b));
  svm_node_closure_store_weight(sd, weight);
}

ccl_device void svm_node_closure_weight(ccl_private ShaderData *sd,
                                        ccl_private float *stack,
                                        uint weight_offset)
{
  float3 weight = stack_load_float3(stack, weight_offset);
  svm_node_closure_store_weight(sd, weight);
}

ccl_device_noinline void svm_node_emission_weight(KernelGlobals kg,
                                                  ccl_private ShaderData *sd,
                                                  ccl_private float *stack,
                                                  uint4 node)
{
  uint color_offset = node.y;
  uint strength_offset = node.z;

  float strength = stack_load_float(stack, strength_offset);
  float3 weight = stack_load_float3(stack, color_offset) * strength;

  svm_node_closure_store_weight(sd, weight);
}

ccl_device_noinline void svm_node_mix_closure(ccl_private ShaderData *sd,
                                              ccl_private float *stack,
                                              uint4 node)
{
  /* fetch weight from blend input, previous mix closures,
   * and write to stack to be used by closure nodes later */
  uint weight_offset, in_weight_offset, weight1_offset, weight2_offset;
  svm_unpack_node_uchar4(
      node.y, &weight_offset, &in_weight_offset, &weight1_offset, &weight2_offset);

  float weight = stack_load_float(stack, weight_offset);
  weight = saturatef(weight);

  float in_weight = (stack_valid(in_weight_offset)) ? stack_load_float(stack, in_weight_offset) :
                                                      1.0f;

  if (stack_valid(weight1_offset))
    stack_store_float(stack, weight1_offset, in_weight * (1.0f - weight));
  if (stack_valid(weight2_offset))
    stack_store_float(stack, weight2_offset, in_weight * weight);
}

/* (Bump) normal */

ccl_device void svm_node_set_normal(KernelGlobals kg,
                                    ccl_private ShaderData *sd,
                                    ccl_private float *stack,
                                    uint in_direction,
                                    uint out_normal)
{
  float3 normal = stack_load_float3(stack, in_direction);
  sd->N = normal;
  stack_store_float3(stack, out_normal, normal);
}

CCL_NAMESPACE_END
