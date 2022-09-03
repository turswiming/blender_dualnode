#pragma once

#include "kernel/light/light.h"

CCL_NAMESPACE_BEGIN

/* to-do: this seems like a relative expensive computation, and we can make it a lot cheaper
 * by using a bounding sphere instead of a bounding box. This will be more inaccurate, but it
 * might be fine when used along with the adaptive splitting. */
ccl_device float light_tree_bounding_box_angle(const float3 bbox_min,
                                               const float3 bbox_max,
                                               const float3 P,
                                               const float3 point_to_centroid)
{
  /* Iterate through all 8 possible points of the bounding box. */
  float theta_u = 0;
  float3 corners[8];
  corners[0] = bbox_min;
  corners[1] = make_float3(bbox_min.x, bbox_min.y, bbox_max.z);
  corners[2] = make_float3(bbox_min.x, bbox_max.y, bbox_min.z);
  corners[3] = make_float3(bbox_min.x, bbox_max.y, bbox_max.z);
  corners[4] = make_float3(bbox_max.x, bbox_min.y, bbox_min.z);
  corners[5] = make_float3(bbox_max.x, bbox_min.y, bbox_max.z);
  corners[6] = make_float3(bbox_max.x, bbox_max.y, bbox_min.z);
  corners[7] = bbox_max;
  for (int i = 0; i < 8; ++i) {
    float3 point_to_corner = normalize(corners[i] - P);
    const float cos_theta_u = dot(point_to_centroid, point_to_corner);
    theta_u = fmaxf(fast_acosf(cos_theta_u), theta_u);
  }
  return theta_u;
}

/* This is the general function for calculating the importance of either a cluster or an emitter.
 * Both of the specialized functions obtain the necessary data before calling this function. */
ccl_device float light_tree_node_importance(const float3 P,
                                            const float3 N,
                                            const float3 bbox_min,
                                            const float3 bbox_max,
                                            const float3 bcone_axis,
                                            const float theta_o,
                                            const float theta_e,
                                            const float energy)
{
  const float3 centroid = 0.5f * bbox_min + 0.5f * bbox_max;
  const float3 point_to_centroid = normalize(centroid - P);

  /* Since we're not using the splitting heuristic, we clamp
   * the distance to half the radius of the cluster. */
  const float distance_squared = fmaxf(0.25f * len_squared(centroid - bbox_max),
                                       len_squared(centroid - P));

  /* to-do: should there be a different importance calculations for different surfaces?
   * opaque surfaces could just return 0 importance in this case. */
  const float theta = fast_acosf(dot(bcone_axis, -point_to_centroid));
  float theta_i = fast_acosf(dot(point_to_centroid, N));
  if (theta_i > M_PI_2_F) {
    theta_i = M_PI_F - theta_i;
  }
  const float theta_u = light_tree_bounding_box_angle(bbox_min, bbox_max, P, point_to_centroid);

  /* Avoid using cosine until needed. */
  const float theta_prime = fmaxf(theta - theta_o - theta_u, 0.0f);
  if (theta_prime >= theta_e) {
    return 0.0f;
  }
  const float cos_theta_prime = fast_cosf(theta_prime);

  float cos_theta_i_prime = 1.0f;
  if (theta_i - theta_u > 0.0f) {
    cos_theta_i_prime = fabsf(fast_cosf(theta_i - theta_u));
  }

  /* to-do: find a good approximation for this value. */
  const float f_a = 1.0f;
  float importance = f_a * cos_theta_i_prime * energy / distance_squared * cos_theta_prime;
  return importance;
}

/* This is uniformly sampling the reservoir for now. */
ccl_device float light_tree_emitter_reservoir_weight(KernelGlobals kg,
                                                     const float3 P,
                                                     const float3 N,
                                                     int emitter_index)
{
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);
  const int prim = kemitter->prim_id;

  /* Triangles are handled normally for now. */
  if (prim < 0) {
    const int lamp = -prim - 1;

    const ccl_global KernelLight *klight = &kernel_data_fetch(lights, lamp);
    float3 light_P = make_float3(klight->co[0], klight->co[1], klight->co[2]);

    /* We use a special calculation to check if a light is
     * within the bounds of a spot or area light. */
    if (klight->type == LIGHT_SPOT) {
      const float radius = klight->spot.radius;
      const float cos_theta = klight->spot.spot_angle;
      const float theta = fast_acosf(cos_theta);
      const float3 light_P = make_float3(klight->co[0], klight->co[1], klight->co[2]);
      const float3 light_dir = make_float3(
          klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);

      const float h1 = radius * fast_sinf(theta);
      const float d1 = radius * cos_theta;
      const float h2 = d1 / fast_tanf(theta);

      const float3 apex = light_P - (h1 + h2) * light_dir;
      const float3 apex_to_point = normalize(P - apex);
      if (dot(apex_to_point, light_dir) < cos_theta) {
        return 0.0f;
      }
    }
    else if (klight->type == LIGHT_AREA) {
      float3 axisu = make_float3(
          klight->area.axisu[0], klight->area.axisu[1], klight->area.axisu[2]);
      float3 axisv = make_float3(
          klight->area.axisv[0], klight->area.axisv[1], klight->area.axisv[2]);
      float3 Ng = make_float3(klight->area.dir[0], klight->area.dir[1], klight->area.dir[2]);
      bool is_round = (klight->area.invarea < 0.0f);

      if (dot(light_P - P, Ng) > 0.0f) {
        return 0.0f;
      }

      if (!is_round) {
        if (klight->area.tan_spread > 0.0f) {
          if (!light_spread_clamp_area_light(
                  P, Ng, &light_P, &axisu, &axisv, klight->area.tan_spread)) {
            return 0.0f;
          }
        }
      }
    }
  }

  return 1.0f;
}

ccl_device float light_tree_emitter_importance(KernelGlobals kg,
                                               const float3 P,
                                               const float3 N,
                                               int emitter_index)
{
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);

  const float3 bbox_min = make_float3(kemitter->bounding_box_min[0],
                                      kemitter->bounding_box_min[1],
                                      kemitter->bounding_box_min[2]);
  const float3 bbox_max = make_float3(kemitter->bounding_box_max[0],
                                      kemitter->bounding_box_max[1],
                                      kemitter->bounding_box_max[2]);
  const float3 bcone_axis = make_float3(kemitter->bounding_cone_axis[0],
                                        kemitter->bounding_cone_axis[1],
                                        kemitter->bounding_cone_axis[2]);

  return light_tree_node_importance(
      P, N, bbox_min, bbox_max, bcone_axis, kemitter->theta_o, kemitter->theta_e, kemitter->energy);
}

ccl_device bool light_tree_should_split(KernelGlobals kg,
                                         const float3 P,
                                         const ccl_global KernelLightTreeNode *knode)
{
  const float splitting_threshold = kernel_data.integrator.splitting_threshold;
  if (splitting_threshold == 0.0f) {
    return false;
  }
  else if (splitting_threshold == 1.0f) {
    return true;
  }
  
  const float3 bbox_min = make_float3(
      knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
  const float3 bbox_max = make_float3(
      knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);
  const float3 centroid = 0.5f * bbox_min + 0.5f * bbox_max;

  const float radius = len(bbox_max - centroid);
  const float distance = len(P - centroid);

  if (distance < radius) {
    return true;
  }

  const float a = distance - radius;
  const float b = distance + radius;

  const float E_g = 1.0f / (a * b);
  const float E_e = knode->energy;

  /* This is a simplified version of the expression given in the paper. */
  const float V_g = (b - a) * (b - a) * E_g * E_g * E_g / 3.0f;
  const float V_e = knode->energy_variance;

  const float total_variance = V_e * V_g + V_e * E_g * E_g + E_e * E_e * V_g;
  const float normalized_variance = sqrt(sqrt(1.0f / (1.0f + sqrt(total_variance))));
  return (normalized_variance < splitting_threshold);
}

ccl_device float light_tree_cluster_importance(KernelGlobals kg,
                                               const float3 P,
                                               const float3 N,
                                               const ccl_global KernelLightTreeNode *knode)
{
  const float3 bbox_min = make_float3(
      knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
  const float3 bbox_max = make_float3(
      knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);
  const float3 bcone_axis = make_float3(
      knode->bounding_cone_axis[0], knode->bounding_cone_axis[1], knode->bounding_cone_axis[2]);

  return light_tree_node_importance(
      P, N, bbox_min, bbox_max, bcone_axis, knode->theta_o, knode->theta_e, knode->energy);
}

ccl_device int light_tree_cluster_select_emitter(KernelGlobals kg,
                                                 float *randu,
                                                 const float3 P,
                                                 const float3 N,
                                                 const ccl_global KernelLightTreeNode *knode,
                                                 float *pdf_factor)
{
  float total_emitter_importance = 0.0f;
  for (int i = 0; i < knode->num_prims; i++) {
    const int prim_index = -knode->child_index + i;
    total_emitter_importance += light_tree_emitter_importance(kg, P, N, prim_index);
  }

  /* to-do: need to handle a case when total importance is 0. */
  if (total_emitter_importance == 0.0f) {
    return -1;
  }

  /* Once we have the total importance, we can normalize the CDF and sample it. */
  const float inv_total_importance = 1.0f / total_emitter_importance;
  float emitter_cdf = 0.0f;
  for (int i = 0; i < knode->num_prims; i++) {
    const int prim_index = -knode->child_index + i;
    /* to-do: is there any way to cache these values, so that recalculation isn't needed? */
    const float emitter_pdf = light_tree_emitter_importance(kg, P, N, prim_index) *
                              inv_total_importance;
    if (*randu < emitter_cdf + emitter_pdf) {
      *randu = (*randu - emitter_cdf) / emitter_pdf;
      *pdf_factor *= emitter_pdf;
      return prim_index;
    }
    emitter_cdf += emitter_pdf;
  }

  /* This point should never be reached. */
  assert(false);
  return -1;
}

/* to-do: for now, we're not going to worry about being in a volume,
 * but this is how the other function determines whether we're in a volume or not. */
template<bool in_volume_segment>
ccl_device bool light_tree_sample(KernelGlobals kg,
                                  ccl_private const RNGState *rng_state,
                                  float *randu,
                                  const float randv,
                                  const float time,
                                  const float3 N,
                                  const float3 P,
                                  const int bounce,
                                  const uint32_t path_flag,
                                  ccl_private LightSample *ls,
                                  float *pdf_factor)
{
  /* We keep track of the currently selected primitive and its weight,
   * as well as the total weight as part of the weighted reservoir sampling. */
  int current_light = -1;
  float current_weight = -1.0f;
  float total_weight = 0.0f;
  float current_pdf = 1.0f;

  /* We need a stack to substitute for recursion. */
  const int stack_size = 32;
  int stack[stack_size];
  float pdfs[stack_size];
  int stack_index = 0;
  stack[0] = 0;
  pdfs[0] = 1.0f;

  /* For now, we arbitrarily limit splitting to 8 so that it doesn't continuously split. */
  int split_count = 0;

  /* First traverse the light tree until a leaf node is reached.
   * Also keep track of the probability of traversing to a given node,
   * so that we can scale our PDF accordingly later. */
  while (stack_index >= 0) {
    const float pdf = pdfs[stack_index];
    const int index = stack[stack_index];
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, index);

    /* If we're at a leaf node, we choose a primitive. Otherwise, we check if we should split
     * or traverse down the light tree. */
    if (knode->child_index <= 0) {
      float light_probability = 1.0f;
      const int selected_light = light_tree_cluster_select_emitter(kg, randu, P, N, knode, &light_probability);

      if (selected_light < 0) {
        stack_index--;
        continue;
      }

      const float light_weight = light_tree_emitter_reservoir_weight(
          kg, P, N, selected_light);
      if (light_weight == 0.0f) {
        stack_index--;
        continue;
      }
      total_weight += light_weight;

      /* We compute the probability of switching to the new candidate sample,
       * otherwise we stick with the old one. */
      const float selection_probability = light_weight / total_weight;
      if (*randu < selection_probability) {
        *randu = *randu / selection_probability;
        current_light = selected_light;
        current_weight = light_weight;
        current_pdf = pdf * light_probability;
      }
      else {
        *randu = (*randu - selection_probability) / (1.0f - selection_probability);
      }

      stack_index--;
      continue;
    }

    /* At an interior node, the left child is directly after the parent,
     * while the right child is stored as the child index.
     * We adaptively split if the variance is high enough. */
    const int left_index = index + 1;
    const int right_index = knode->child_index;
    if (light_tree_should_split(kg, P, knode) &&
        split_count < 8 && 
        stack_index < stack_size - 1) {
      stack[stack_index] = left_index;
      pdfs[stack_index] = pdf;
      stack[stack_index + 1] = right_index;
      pdfs[stack_index + 1] = pdf;
      stack_index++;
      split_count++;
      continue;
    }

    /* If we don't split, then we need to choose sampling between the left or right child. */
    const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, left_index);
    const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes, right_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, left);
    const float right_importance = light_tree_cluster_importance(kg, P, N, right);
    const float total_importance = left_importance + right_importance;

    if (total_importance == 0.0f) {
      stack_index--;
      continue;
    }
    float left_probability = left_importance / (left_importance + right_importance);

    if (*randu < left_probability) {
      stack[stack_index] = left_index;
      *randu = *randu / left_probability;
      pdfs[stack_index] = pdf * left_probability;
    }
    else {
      stack[stack_index] = right_index;
      *randu = (*randu - left_probability) / (1.0f - left_probability);
      pdfs[stack_index] = pdf * (1.0f - left_probability);
    }
  }

  if (total_weight == 0.0f) {
    return false;
  }

  *pdf_factor *= current_pdf * current_weight / total_weight;

  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         current_light);

  /* to-do: this is the same code as light_distribution_sample, except the index is determined
   * differently. Would it be better to refactor this into a separate function? */
  const int prim = kemitter->prim_id;
  if (prim >= 0) {
    /* Mesh light. */
    const int object = kemitter->mesh_light.object_id;

    /* Exclude synthetic meshes from shadow catcher pass. */
    if ((path_flag & PATH_RAY_SHADOW_CATCHER_PASS) &&
        !(kernel_data_fetch(object_flag, object) & SD_OBJECT_SHADOW_CATCHER)) {
      return false;
    }

    const int shader_flag = kemitter->mesh_light.shader_flag;
    triangle_light_sample<in_volume_segment>(kg, prim, object, *randu, randv, time, ls, P);
    ls->shader |= shader_flag;

    return (ls->pdf > 0.0f);
  }

  const int lamp = -prim - 1;

  if (UNLIKELY(light_select_reached_max_bounces(kg, lamp, bounce))) {
    return false;
  }

  return light_sample<in_volume_segment>(kg, lamp, *randu, randv, P, path_flag, ls);
}

ccl_device float light_tree_distant_light_importance(KernelGlobals kg,
                                                     const float3 P,
                                                     const float3 N,
                                                     const int index)
{
  ccl_global const KernelLightTreeDistantEmitter *kdistant = &kernel_data_fetch(
      light_tree_distant_group, index);

  if (kdistant->energy == 0.0f) {
    return 0.0f;
  }

  const float3 light_axis = make_float3(
      kdistant->direction[0], kdistant->direction[1], kdistant->direction[2]);
  const float theta = fast_acosf(dot(N, light_axis));
  const float theta_i_prime = theta - kdistant->bounding_radius;

  float cos_theta_i_prime = 1.0f;
  if (theta_i_prime > M_PI_2_F) {
    return 0.0f;
  } else if (theta - kdistant->bounding_radius > 0.0f) {
    cos_theta_i_prime = fast_cosf(theta - kdistant->bounding_radius);
  }

  /* to-do: find a good value for this. */
  const float f_a = 1.0f;
  float importance = f_a * cos_theta_i_prime * kdistant->energy;
  return importance;
}

template<bool in_volume_segment>
ccl_device bool light_tree_sample_distant_lights(KernelGlobals kg,
                                                ccl_private const RNGState *rng_state,
                                                float *randu,
                                                const float randv,
                                                const float time,
                                                const float3 N,
                                                const float3 P,
                                                const int bounce,
                                                const uint32_t path_flag,
                                                ccl_private LightSample *ls,
                                                float *pdf_factor)
{
  const int num_distant_lights = kernel_data.integrator.num_distant_lights;
  float total_importance = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    total_importance += light_tree_distant_light_importance(kg, P, N, i);
  }
  const float inv_total_importance = 1.0f / total_importance;

  float light_cdf = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    const float light_pdf = light_tree_distant_light_importance(kg, P, N, i) *
                            inv_total_importance;
    if (*randu < light_cdf + light_pdf) {
      *randu = (*randu - light_cdf) / light_pdf;
      *pdf_factor *= light_pdf;
      ccl_global const KernelLightTreeDistantEmitter *kdistant = &kernel_data_fetch(
          light_tree_distant_group, i);

      const int lamp = kdistant->prim_id;

      if (UNLIKELY(light_select_reached_max_bounces(kg, lamp, bounce))) {
        return false;
      }

      return light_sample<in_volume_segment>(kg, lamp, *randu, randv, P, path_flag, ls);
    }
    light_cdf += light_pdf;
  }

  /* We should never reach this point. */
  assert(false);
  return -1;
}

/* We need to be able to find the probability of selecting a given light for MIS. */
ccl_device float light_tree_pdf(KernelGlobals kg,
                                ConstIntegratorState state,
                                const float3 P,
                                const float3 N,
                                const int prim)
{
  float distant_light_importance = light_tree_distant_light_importance(
      kg, P, N, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, kroot);
  }
  const float total_group_importance = light_tree_importance + distant_light_importance;
  assert(total_group_importance != 0.0f);
  float pdf = light_tree_importance / total_group_importance;

  const int emitter = (prim >= 0) ? kernel_data_fetch(triangle_to_tree, prim) :
                                    kernel_data_fetch(light_to_tree, ~prim);
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter);
  ccl_global const KernelLightTreeNode *kleaf = &kernel_data_fetch(light_tree_nodes,
                                                                   kemitter->parent_index);

  /* We generate a random number to use for selecting a light. */
  RNGState rng_state;
  path_state_rng_load(state, &rng_state);
  /* to-do: is this the correct macro to use? */
  float randu = path_state_rng_1D(kg, &rng_state, PRNG_LIGHT); 

  /* We traverse to the leaf node and
   * find the probability of selecting the target light. */
  const int stack_size = 32;
  int stack[stack_size];
  float pdfs[stack_size];
  int stack_index = 0;
  stack[0] = 0;
  pdfs[0] = 1.0f;

  int split_count = 0;

  float light_tree_pdf = 0.0f;
  float light_leaf_pdf = 0.0f;
  float total_weight = 0.0f;
  float target_weight = 0.0f;

  uint bit_trail = kleaf->bit_trail;
  while (stack_index >= 0) {
    const float pdf = pdfs[stack_index];
    const int index = stack[stack_index];
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, index);

    if (knode->child_index <= 0) {
      int selected_light = -1;
      float light_weight = 0.0f;

      /* If we're at the leaf node containing the light we need, 
       * then we iterate through the lights to find the target emitter.
       * Otherwise, we randomly select one. */
      if (index == emitter) {
        light_tree_pdf = pdf;

        float target_emitter_importance = 0.0f;
        float total_emitter_importance = 0.0f;
        for (int i = 0; i < knode->num_prims; i++) {
          const int prim_index = -knode->child_index + i;
          float light_importance = light_tree_emitter_importance(kg, P, N, prim_index);
          ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(
              light_tree_emitters, prim_index);
          if (kemitter->prim_id == prim) {
            selected_light = prim_index;
            light_weight = light_tree_emitter_reservoir_weight(kg, P, N, selected_light);
            target_weight = light_weight;
            target_emitter_importance = light_importance;
          }
          total_emitter_importance += light_importance;
        }

        light_leaf_pdf = target_emitter_importance / total_emitter_importance;
      }
      else {
        float light_probability = 1.0f;
        const int selected_light = light_tree_cluster_select_emitter(
            kg, &randu, P, N, knode, &light_probability);
        light_weight = light_tree_emitter_reservoir_weight(kg, P, N, selected_light);
      }

      if (selected_light < 0) {
        stack_index--;
        continue;
      }

      if (light_weight == 0.0f) {
        stack_index--;
        continue;
      }
      total_weight += light_weight;

      stack_index--;
      continue;
    }

    /* At an interior node, the left child is directly after the parent,
     * while the right child is stored as the child index.
     * We adaptively split if the variance is high enough. */
    const int left_index = index + 1;
    const int right_index = knode->child_index;
    if (light_tree_should_split(kg, P, knode) &&
        split_count < 8 &&
        stack_index < stack_size - 1) {
      stack[stack_index] = left_index;
      pdfs[stack_index] = pdf;
      stack[stack_index + 1] = right_index;
      pdfs[stack_index + 1] = pdf;
      stack_index++;
      split_count++;
      continue;
    }

    /* If we don't split, the bit trail determines whether we go left or right. */
    const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, left_index);
    const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes,
                                                                     right_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, left);
    const float right_importance = light_tree_cluster_importance(kg, P, N, right);
    const float total_importance = left_importance + right_importance;

    if (total_importance == 0.0f) {
      stack_index--;
      continue;
    }
    float left_probability = left_importance / (left_importance + right_importance);

    if (bit_trail & 1) {
      stack[stack_index] = left_index;
      pdfs[stack_index] = pdf * left_probability;
    }
    else {
      stack[stack_index] = right_index;
      pdfs[stack_index] = pdf * (1.0f - left_probability);
    }
  }
 
  pdf *= light_leaf_pdf * light_tree_pdf * target_weight / total_weight;
  return pdf;
}

ccl_device float distant_lights_pdf(KernelGlobals kg,
                                    const float3 P,
                                    const float3 N,
                                    const int prim)
{
  float distant_light_importance = light_tree_distant_light_importance(
      kg, P, N, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, kroot);
  }
  const float total_group_importance = light_tree_importance + distant_light_importance;
  assert(total_group_importance != 0.0f);
  float pdf = distant_light_importance / total_group_importance;

  /* The light_to_tree array doubles as a lookup table for
   * both the light tree as well as the distant lights group.*/
  const int distant_light = kernel_data_fetch(light_to_tree, prim);
  const int num_distant_lights = kernel_data.integrator.num_distant_lights;

  float emitter_importance = 0.0f;
  float total_importance = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    float importance = light_tree_distant_light_importance(kg, P, N, i);
    if (i == distant_light) {
      emitter_importance = importance;
    }
    total_importance += importance;
  }

  pdf *= emitter_importance / total_importance;
  return pdf;
}

ccl_device bool light_tree_sample_from_position(KernelGlobals kg,
                                                ccl_private const RNGState *rng_state,
                                                float randu,
                                                const float randv,
                                                const float time,
                                                const float3 P,
                                                const float3 N,
                                                const int bounce,
                                                const uint32_t path_flag,
                                                ccl_private LightSample *ls)
{
  /* to-do: with weighted reservoir sampling, we can also try picking a sample from the distant light group
   * and compare it to the sample from the light tree. */
  float distant_light_importance = light_tree_distant_light_importance(
      kg, P, N, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, kroot);
  }
  const float total_importance = light_tree_importance + distant_light_importance;

  if (total_importance == 0.0f) {
    return false;
  }

  const float light_tree_probability = light_tree_importance / total_importance;

  float pdf_factor = 1.0f;
  bool ret;
  if (randu < light_tree_probability) {
    randu = randu / light_tree_probability;
    pdf_factor *= light_tree_probability;
    ret = light_tree_sample<false>(
        kg, rng_state, &randu, randv, time, N, P, bounce, path_flag, ls, &pdf_factor);
  }
  else {
    randu = (randu - light_tree_probability) / (1.0f - light_tree_probability);
    pdf_factor *= (1.0f - light_tree_probability);
    ret = light_tree_sample_distant_lights<false>(
        kg, rng_state, &randu, randv, time, N, P, bounce, path_flag, ls, &pdf_factor);
  }

  ls->pdf *= pdf_factor;
  return ret;
}

CCL_NAMESPACE_END
