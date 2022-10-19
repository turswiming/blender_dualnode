#pragma once

#include "kernel/light/light.h"

CCL_NAMESPACE_BEGIN

/* This is the general function for calculating the importance of either a cluster or an emitter.
 * Both of the specialized functions obtain the necessary data before calling this function. */
ccl_device float light_tree_node_importance(const float3 P,
                                            const float3 N,
                                            const bool has_transmission,
                                            const float3 bbox_min,
                                            const float3 bbox_max,
                                            const float3 bcone_axis,
                                            const float theta_o,
                                            const float theta_e,
                                            const float energy)
{
  const float3 centroid = 0.5f * bbox_min + 0.5f * bbox_max;
  const float3 point_to_centroid = normalize(centroid - P);

  /* TODO: we're using the splitting heuristic now, do we still need to clamp the distance to half
   * the radius of the cluster? */
  const float distance_squared = fmaxf(0.25f * len_squared(centroid - bbox_max),
                                       len_squared(centroid - P));

  const float cos_theta = dot(bcone_axis, -point_to_centroid);
  const float cos_theta_i = has_transmission ? fabsf(dot(point_to_centroid, N)) :
                                               dot(point_to_centroid, N);

  bool bbox_is_behind_surface = !has_transmission && (cos_theta_i < 0);

  /* TODO: this seems like a relative expensive computation, and we can make it a lot cheaper
   * by using a bounding sphere instead of a bounding box. This will be more inaccurate, but it
   * might be fine when used along with the adaptive splitting. */
  float cos_theta_u = 1.0f;

  /* Iterate through all 8 possible points of the bounding box. */
  for (int i = 0; i < 8; ++i) {
    const float3 corner = make_float3((i & 1) ? bbox_max.x : bbox_min.x,
                                      (i & 2) ? bbox_max.y : bbox_min.y,
                                      (i & 4) ? bbox_max.z : bbox_min.z);

    /* Caculate the bounding box angle. */
    float3 point_to_corner = normalize(corner - P);
    cos_theta_u = fminf(cos_theta_u, dot(point_to_centroid, point_to_corner));

    /* Figure out whether or not the bounding box is in front or behind the shading point. */
    if (bbox_is_behind_surface && dot(point_to_corner, N) >= 0) {
      bbox_is_behind_surface = false;
    }
  }

  /* If the node is guaranteed to be behind the surface we're sampling, and the surface is opaque,
   * then we can give the node an importance of 0 as it contributes nothing to the surface. */
  if (bbox_is_behind_surface) {
    return 0.0f;
  }

  const float sin_theta_u = safe_sqrtf(1.0f - sqr(cos_theta_u));

  /* cos(theta - theta_u) */
  const float sin_theta = safe_sqrtf(1.0f - sqr(cos_theta));
  const float cos_theta_minus_theta_u = cos_theta * cos_theta_u + sin_theta * sin_theta_u;

  float cos_theta_o, sin_theta_o;
  fast_sincosf(theta_o, &sin_theta_o, &cos_theta_o);

  float cos_theta_prime;
  if ((cos_theta > cos_theta_u) || (cos_theta_minus_theta_u > cos_theta_o)) {
    /* theta - theta_o - theta_u < 0 */
    kernel_assert((fast_acosf(cos_theta) - theta_o - fast_acosf(cos_theta_u)) < 1e-4f);
    cos_theta_prime = 1.0f;
  }
  else if ((cos_theta > cos_theta_u) || (theta_o + theta_e > M_PI_F) ||
           (cos_theta_minus_theta_u > cos(theta_o + theta_e))) {
    /* theta' = theta - theta_o - theta_u < theta_e */
    kernel_assert((fast_acosf(cos_theta) - theta_o - fast_acosf(cos_theta_u) - theta_e) < 1e-4f);
    const float sin_theta_minus_theta_u = safe_sqrtf(1.0f - sqr(cos_theta_minus_theta_u));
    cos_theta_prime = cos_theta_minus_theta_u * cos_theta_o +
                      sin_theta_minus_theta_u * sin_theta_o;
  }
  else {
    return 0.f;
  }

  /* cos_theta_i_prime = |cos(max{theta_i - theta_u, 0})| */
  float cos_theta_i_prime;
  if (cos_theta_i > cos_theta_u) {
    cos_theta_i_prime = 1.0f;
  }
  else {
    kernel_assert(fast_acosf(cos_theta_i) >= fast_acosf(cos_theta_u));
    const float sin_theta_i = safe_sqrtf(1.0f - sqr(cos_theta_i));
    cos_theta_i_prime = cos_theta_i * cos_theta_u + sin_theta_i * sin_theta_u;
  }

  /* TODO: find a good approximation for this value. */
  const float f_a = 1.0f;
  return fabsf(f_a * cos_theta_i_prime * energy / distance_squared * cos_theta_prime);
}

/* This is uniformly sampling the reservoir for now. */
ccl_device float light_tree_emitter_reservoir_weight(KernelGlobals kg,
                                                     const float3 P,
                                                     const float3 N,
                                                     int emitter_index)
{
  if (emitter_index < 0) {
    return 0.0f;
  }

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
                                               const bool has_transmission,
                                               int emitter_index)
{
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);

  const float3 bbox_min = make_float3(
      kemitter->bounding_box_min[0], kemitter->bounding_box_min[1], kemitter->bounding_box_min[2]);
  const float3 bbox_max = make_float3(
      kemitter->bounding_box_max[0], kemitter->bounding_box_max[1], kemitter->bounding_box_max[2]);
  const float3 bcone_axis = make_float3(kemitter->bounding_cone_axis[0],
                                        kemitter->bounding_cone_axis[1],
                                        kemitter->bounding_cone_axis[2]);

  return light_tree_node_importance(P,
                                    N,
                                    has_transmission,
                                    bbox_min,
                                    bbox_max,
                                    bcone_axis,
                                    kemitter->theta_o,
                                    kemitter->theta_e,
                                    kemitter->energy);
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
                                               const bool has_transmission,
                                               const ccl_global KernelLightTreeNode *knode)
{
  const float3 bbox_min = make_float3(
      knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
  const float3 bbox_max = make_float3(
      knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);
  const float3 bcone_axis = make_float3(
      knode->bounding_cone_axis[0], knode->bounding_cone_axis[1], knode->bounding_cone_axis[2]);

  return light_tree_node_importance(P,
                                    N,
                                    has_transmission,
                                    bbox_min,
                                    bbox_max,
                                    bcone_axis,
                                    knode->theta_o,
                                    knode->theta_e,
                                    knode->energy);
}

/* pick an emitter from a leaf node using resevoir sampling */
ccl_device int light_tree_cluster_select_emitter(KernelGlobals kg,
                                                 ccl_private float *randu,
                                                 const float3 P,
                                                 const float3 N,
                                                 const bool has_transmission,
                                                 const ccl_global KernelLightTreeNode *knode,
                                                 ccl_private float *pdf_factor)
{
  /* the first emitter in the cluster */
  int selected_prim_index = -knode->child_index;
  float total_weight = light_tree_emitter_importance(
      kg, P, N, has_transmission, selected_prim_index);
  float selected_weight = total_weight;

  for (int i = 1; i < knode->num_prims; i++) {
    int current_prim_index = -knode->child_index + i;
    float current_weight = light_tree_emitter_importance(
        kg, P, N, has_transmission, current_prim_index);
    total_weight += current_weight;

    if (total_weight == 0.0f) {
      continue;
    }

    float thresh = current_weight / total_weight;
    if (*randu < thresh) {
      selected_prim_index = current_prim_index;
      selected_weight = current_weight;
      *randu = *randu / thresh;
    }
    else {
      *randu = (*randu - thresh) / (1.0f - thresh);
    }
  }

  if (total_weight == 0.0f) {
    return -1;
  }

  *pdf_factor = selected_weight / total_weight;
  return selected_prim_index;
}

/* to-do: for now, we're not going to worry about being in a volume,
 * but this is how the other function determines whether we're in a volume or not. */
template<bool in_volume_segment>
ccl_device bool light_tree_sample(KernelGlobals kg,
                                  ccl_private const RNGState *rng_state,
                                  ccl_private float *randu,
                                  const float randv,
                                  const float time,
                                  const float3 P,
                                  const float3 N,
                                  const bool has_transmission,
                                  const int path_flag,
                                  const int bounce,
                                  ccl_private LightSample *ls,
                                  ccl_private float *pdf_factor)
{
  /* We keep track of the currently selected primitive and its weight,
   * as well as the total weight as part of the weighted reservoir sampling. */
  int current_light = -1;
  float current_reservoir_weight = -1.0f;
  float total_reservoir_weight = 0.0f;
  float pdf_node_emitter_selection = 1.0f;

  /* We need a stack to substitute for recursion. */
  const int stack_size = 32;
  int stack[stack_size];
  float pdfs_node_selection[stack_size];
  int stack_index = 0;
  stack[0] = 0;
  pdfs_node_selection[0] = 1.0f;

  /* First traverse the light tree until a leaf node is reached.
   * Also keep track of the probability of traversing to a given node,
   * so that we can scale our PDF accordingly later. */
  while (stack_index >= 0) {
    const float pdf_node_selection = pdfs_node_selection[stack_index];
    const int index = stack[stack_index];
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, index);

    /* If we're at a leaf node, we choose a primitive. Otherwise, we check if we should split
     * or traverse down the light tree. */
    if (knode->child_index <= 0) {
      float pdf_emitter_selection = 1.0f;
      const int selected_light = light_tree_cluster_select_emitter(
          kg, randu, P, N, has_transmission, knode, &pdf_emitter_selection);

      if (selected_light < 0) {
        stack_index--;
        continue;
      }

      const float light_reservoir_weight = light_tree_emitter_reservoir_weight(
          kg, P, N, selected_light);

      /* TODO: make pdf_node_emitter_selection part of the light_reservoir_weight, otherwise result
       * is suboptimal. */

      if (light_reservoir_weight == 0.0f) {
        stack_index--;
        continue;
      }
      total_reservoir_weight += light_reservoir_weight;

      /* We compute the probability of switching to the new candidate sample,
       * otherwise we stick with the old one. */
      const float selection_probability = light_reservoir_weight / total_reservoir_weight;
      if (*randu <= selection_probability) {
        *randu = *randu / selection_probability;
        current_light = selected_light;
        current_reservoir_weight = light_reservoir_weight;
        pdf_node_emitter_selection = pdf_node_selection * pdf_emitter_selection;
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
    if (light_tree_should_split(kg, P, knode) && stack_index < stack_size - 1) {
      stack[stack_index] = left_index;
      pdfs_node_selection[stack_index] = pdf_node_selection;
      stack[stack_index + 1] = right_index;
      pdfs_node_selection[stack_index + 1] = pdf_node_selection;
      stack_index++;
      continue;
    }

    /* If we don't split, then we need to choose sampling between the left or right child. */
    const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, left_index);
    const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes,
                                                                     right_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, has_transmission, left);
    const float right_importance = light_tree_cluster_importance(
        kg, P, N, has_transmission, right);
    const float total_importance = left_importance + right_importance;

    if (total_importance == 0.0f) {
      stack_index--;
      continue;
    }
    float left_probability = left_importance / (left_importance + right_importance);

    if (*randu <= left_probability) {
      stack[stack_index] = left_index;
      *randu = *randu / left_probability;
      pdfs_node_selection[stack_index] = pdf_node_selection * left_probability;
    }
    else {
      stack[stack_index] = right_index;
      *randu = (*randu - left_probability) / (1.0f - left_probability);
      pdfs_node_selection[stack_index] = pdf_node_selection * (1.0f - left_probability);
    }
  }

  if (total_reservoir_weight == 0.0f || current_light < 0) {
    return false;
  }

  const float pdf_reservoir = current_reservoir_weight / total_reservoir_weight;
  *pdf_factor *= pdf_node_emitter_selection * pdf_reservoir;

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

    const int mesh_shader_flag = kemitter->mesh_light.shader_flag;
    triangle_light_sample<in_volume_segment>(kg, prim, object, *randu, randv, time, ls, P);
    ls->shader |= mesh_shader_flag;

    return (ls->pdf > 0.0f);
  }

  const int lamp = -prim - 1;

  if (UNLIKELY(light_select_reached_max_bounces(kg, lamp, bounce))) {
    return false;
  }

  return light_sample<in_volume_segment>(kg, lamp, *randu, randv, P, path_flag, ls);
}

ccl_device float light_tree_distant_light_importance(KernelGlobals kg,
                                                     const float3 N,
                                                     const bool has_transmission,
                                                     const int index)
{
  ccl_global const KernelLightTreeDistantEmitter *kdistant = &kernel_data_fetch(
      light_tree_distant_group, index);

  if (kdistant->energy == 0.0f) {
    return 0.0f;
  }

  const float3 light_axis = make_float3(
      kdistant->direction[0], kdistant->direction[1], kdistant->direction[2]);
  float theta_i = fast_acosf(dot(N, light_axis));

  /* If the light is guaranteed to be behind the surface we're sampling, and the surface is opaque,
   * then we can give the light an importance of 0 as it contributes nothing to the surface. */
  if (!has_transmission && (theta_i - kdistant->bounding_radius > M_PI_2_F)) {
    return 0.0f;
  }

  if (theta_i > M_PI_2_F) {
    theta_i = M_PI_F - theta_i;
  }

  float cos_theta_i_prime = 1.0f;
  const float theta_i_prime = theta_i - kdistant->bounding_radius;
  if (theta_i_prime > 0.0f) {
    cos_theta_i_prime = fast_cosf(theta_i_prime);
  }

  /* to-do: find a good value for this. */
  const float f_a = 1.0f;
  float importance = f_a * cos_theta_i_prime * kdistant->energy;
  return importance;
}

template<bool in_volume_segment>
ccl_device bool light_tree_sample_distant_lights(KernelGlobals kg,
                                                 ccl_private const RNGState *rng_state,
                                                 ccl_private float *randu,
                                                 const float randv,
                                                 const float time,
                                                 const float3 P,
                                                 const float3 N,
                                                 const bool has_transmission,
                                                 const int path_flag,
                                                 const int bounce,
                                                 ccl_private LightSample *ls,
                                                 ccl_private float *pdf_factor)
{
  /* TODO: do single loop over lights to avoid computing importance twice? */

  const int num_distant_lights = kernel_data.integrator.num_distant_lights;
  float total_importance = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    total_importance += light_tree_distant_light_importance(kg, N, has_transmission, i);
  }

  if (total_importance == 0.0f) {
    return -1;
  }

  float light_cdf = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    const float light_pdf = light_tree_distant_light_importance(kg, N, has_transmission, i) /
                            total_importance;
    if (*randu <= light_cdf + light_pdf) {
      *randu = (*randu - light_cdf) / light_pdf;
      *pdf_factor *= light_pdf;
      ccl_global const KernelLightTreeDistantEmitter *kdistant = &kernel_data_fetch(
          light_tree_distant_group, i);

      const int lamp = kdistant->prim_id;

      if (UNLIKELY(light_select_reached_max_bounces(kg, lamp, bounce))) {
        return -1;
      }

      return light_sample<in_volume_segment>(kg, lamp, *randu, randv, P, path_flag, ls);
    }
    light_cdf += light_pdf;
  }

  /* TODO: change implementation so that even under precision issues, we always end
   * up selecting a light and this can't happen? */
  kernel_assert(false);
  return -1;
}

/* We need to be able to find the probability of selecting a given light for MIS.
 * Returns a conditioned pdf which is consistent with the pdf computed in `light_tree_sample()` */
ccl_device float light_tree_pdf(KernelGlobals kg,
                                ConstIntegratorState state,
                                const float3 P,
                                const float3 N,
                                const int path_flag,
                                const int prim)
{
  const bool has_transmission = (path_flag & PATH_RAY_MIS_HAD_TRANSMISSION);
  float distant_light_importance = light_tree_distant_light_importance(
      kg, N, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, has_transmission, kroot);
  }
  const float total_group_importance = light_tree_importance + distant_light_importance;
  kernel_assert(total_group_importance != 0.0f);
  float pdf = light_tree_importance / total_group_importance;

  const int emitter = (prim >= 0) ? kernel_data_fetch(triangle_to_tree, prim) :
                                    kernel_data_fetch(light_to_tree, ~prim);
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter);
  const int target_leaf = kemitter->parent_index;
  ccl_global const KernelLightTreeNode *kleaf = &kernel_data_fetch(light_tree_nodes, target_leaf);

  uint bit_trail = kleaf->bit_trail;

  /* We generate a random number to use for selecting a light. */
  RNGState rng_state;
  path_state_rng_load(state, &rng_state);
  /* to-do: is this the correct macro to use? */
  float randu = path_state_rng_1D(kg, &rng_state, PRNG_LIGHT);

  /* Keep track of the currently selected primitive and its weight,
   * as well as the total weight as part of the weighted reservoir sampling. */
  int selected_light = -1;
  float selected_reservoir_weight = -1.0f;
  float total_reservoir_weight = 0.0f;

  /* We need a stack to substitute for recursion. */
  const int stack_size = 32;
  int stack[stack_size];
  int stack_index = 0;
  stack[0] = 0;

  while (stack_index >= 0) {
    const bool traversing_target_branch = stack_index == 0;
    const int index = stack[stack_index];
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, index);

    /* TODO: interior nodes are processed before leaf nodes, but current implementation forces the
     * reader to think about the leaf node first. Maybe switch the two if-branches for better
     * understandability? */

    /* Leaf node */
    if (knode->child_index <= 0) {

      float pdf_emitter_selection = 1.0f;

      /* If the leaf node contains the target emitter, we are processing the last node.
       * We then iterate through the lights to find the target emitter.
       * Otherwise, we randomly select one. */
      if (index == target_leaf) {
        float target_emitter_importance = 0.0f;
        float total_emitter_importance = 0.0f;
        for (int i = 0; i < knode->num_prims; i++) {
          const int leaf_emitter = -knode->child_index + i;
          float light_importance = light_tree_emitter_importance(
              kg, P, N, has_transmission, leaf_emitter);
          if (leaf_emitter == emitter) {
            target_emitter_importance = light_importance;
            selected_light = leaf_emitter;
            selected_reservoir_weight = light_tree_emitter_reservoir_weight(
                kg, P, N, selected_light);
            total_reservoir_weight += selected_reservoir_weight;
          }
          total_emitter_importance += light_importance;
        }
        pdf_emitter_selection = target_emitter_importance / total_emitter_importance;
        const float pdf_reservoir = selected_reservoir_weight / total_reservoir_weight;
        pdf *= pdf_emitter_selection * pdf_reservoir;
      }
      else {
        selected_light = light_tree_cluster_select_emitter(
            kg, &randu, P, N, has_transmission, knode, &pdf_emitter_selection);

        if (selected_light < 0) {
          stack_index--;
          continue;
        }

        total_reservoir_weight += light_tree_emitter_reservoir_weight(kg, P, N, selected_light);
      }
      stack_index--;
      continue;
    }

    /* Interior node */
    const bool go_left = (bit_trail & 1) == 0;
    bit_trail >>= traversing_target_branch;
    const int left_index = index + 1;
    const int right_index = knode->child_index;
    if (light_tree_should_split(kg, P, knode) && stack_index < stack_size - 1) {
      /* Always store nodes on the target branch in stack[0] */
      stack[stack_index] = go_left ? left_index : right_index;
      stack_index++;
      stack[stack_index] = go_left ? right_index : left_index;
      continue;
    }

    /* No splitting, choose between the left or the right child */
    const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, left_index);
    const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes,
                                                                     right_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, has_transmission, left);
    const float right_importance = light_tree_cluster_importance(
        kg, P, N, has_transmission, right);
    const float total_importance = left_importance + right_importance;

    if (total_importance == 0.0f) {
      stack_index--;
      continue;
    }
    float left_probability = left_importance / (left_importance + right_importance);

    if (traversing_target_branch) {
      pdf *= go_left ? left_probability : (1.0f - left_probability);
      stack[stack_index] = go_left ? left_index : right_index;
    }
    else {
      if (randu <= left_probability) { /* traverse left */
        stack[stack_index] = left_index;
        randu = randu / left_probability;
      }
      else { /* traverse right */
        stack[stack_index] = right_index;
        randu = (randu - left_probability) / (1.0f - left_probability);
      }
    }
  }
  return pdf;
}

ccl_device float distant_lights_pdf(
    KernelGlobals kg, const float3 P, const float3 N, const int path_flag, const int prim)
{
  const bool has_transmission = (path_flag & PATH_RAY_MIS_HAD_TRANSMISSION);
  float distant_light_importance = light_tree_distant_light_importance(
      kg, N, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, has_transmission, kroot);
  }
  const float total_group_importance = light_tree_importance + distant_light_importance;
  kernel_assert(total_group_importance != 0.0f);
  float pdf = distant_light_importance / total_group_importance;

  /* The light_to_tree array doubles as a lookup table for
   * both the light tree as well as the distant lights group.*/
  const int distant_light = kernel_data_fetch(light_to_tree, prim);
  const int num_distant_lights = kernel_data.integrator.num_distant_lights;

  float emitter_importance = 0.0f;
  float total_importance = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    float importance = light_tree_distant_light_importance(kg, N, has_transmission, i);
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
                                                const int shader_flags,
                                                const int bounce,
                                                const uint32_t path_flag,
                                                ccl_private LightSample *ls)
{
  const bool has_transmission = (shader_flags & SD_BSDF_HAS_TRANSMISSION);
  /* to-do: with weighted reservoir sampling, we can also try picking a sample from the distant
   * light group and compare it to the sample from the light tree. */
  float distant_light_importance = light_tree_distant_light_importance(
      kg, N, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    light_tree_importance = light_tree_cluster_importance(kg, P, N, has_transmission, kroot);
  }
  const float total_importance = light_tree_importance + distant_light_importance;

  if (total_importance == 0.0f) {
    return false;
  }

  const float light_tree_probability = light_tree_importance / total_importance;

  float pdf_factor = 1.0f;
  bool ret;
  if (randu <= light_tree_probability) {
    randu = randu / light_tree_probability;
    pdf_factor *= light_tree_probability;
    ret = light_tree_sample<false>(kg,
                                   rng_state,
                                   &randu,
                                   randv,
                                   time,
                                   P,
                                   N,
                                   has_transmission,
                                   path_flag,
                                   bounce,
                                   ls,
                                   &pdf_factor);
  }
  else {
    randu = (randu - light_tree_probability) / (1.0f - light_tree_probability);
    pdf_factor *= (1.0f - light_tree_probability);
    ret = light_tree_sample_distant_lights<false>(kg,
                                                  rng_state,
                                                  &randu,
                                                  randv,
                                                  time,
                                                  P,
                                                  N,
                                                  has_transmission,
                                                  path_flag,
                                                  bounce,
                                                  ls,
                                                  &pdf_factor);
  }

  ls->pdf *= pdf_factor;
  return ret;
}

CCL_NAMESPACE_END
