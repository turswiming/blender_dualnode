#pragma once

#include "kernel/light/light.h"

CCL_NAMESPACE_BEGIN

ccl_device float light_tree_bounding_box_angle(const float3 bbox_min,
                                               const float3 bbox_max,
                                               const float3 P,
                                               const float3 point_to_centroid)
{
  /* Want to iterate through all 8 possible points of the bounding box. */
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
 * Both of the specialized functions obtain the necessary data before calling this function.
 * to-do: find a better way to handle this? or rename it to be more clear? */
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

  /* to-do: compare this with directly using fmaxf and cosf. */
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

ccl_device float light_tree_emitter_importance(KernelGlobals kg,
                                               const float3 P,
                                               const float3 N,
                                               int emitter_index)
{
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);

  /* Convert the data from the struct into float3 for calculations. */
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

ccl_device float light_tree_cluster_importance(KernelGlobals kg,
                                               const float3 P,
                                               const float3 N,
                                               const ccl_global KernelLightTreeNode *knode)
{
  /* Convert the data from the struct into float3 for calculations. */
  const float3 bbox_min = make_float3(
      knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
  const float3 bbox_max = make_float3(
      knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);
  const float3 bcone_axis = make_float3(
      knode->bounding_cone_axis[0], knode->bounding_cone_axis[1], knode->bounding_cone_axis[2]);

  return light_tree_node_importance(
      P, N, bbox_min, bbox_max, bcone_axis, knode->theta_o, knode->theta_e, knode->energy);
}

/* to-do: for now, we're not going to worry about being in a volume for now,
 * but this seems to be a good way to differentiate whether we're in a volume or not. */
template<bool in_volume_segment>
ccl_device int light_tree_sample(KernelGlobals kg,
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
  /* First traverse the light tree until a leaf node is reached. */
  /* Also keep track of the probability of traversing to a given node, */
  /* so that we can scale our PDF accordingly later. */
  int index = 0;

  /* to-do: is it better to generate a new random sample for each step of the traversal? */
  const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, index);
  while (knode->child_index > 0) {
    /* At an interior node, the left child is directly next to the parent,
     * while the right child is stored as the child index. */
    const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, index + 1);
    const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes, knode->child_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, left);
    const float right_importance = light_tree_cluster_importance(kg, P, N, right);
    float left_probability = left_importance / (left_importance + right_importance);

    if (*randu < left_probability) {
      index = index + 1;
      knode = left;
      *randu = *randu / left_probability;
      *pdf_factor *= left_probability;
    }
    else {
      index = knode->child_index;
      knode = right;
      *randu = (*randu - left_probability) / (1.0f - left_probability);
      *pdf_factor *= (1.0f - left_probability);
    }
  }

  /* Right now, sampling is done by incrementing the CDF by the PDF.
   * However, we first need to calculate the total importance so that we can normalize the CDF. */
  float total_emitter_importance = 0.0f;
  for (int i = 0; i < knode->num_prims; i++) {
    const int prim_index = -knode->child_index + i;
    total_emitter_importance += light_tree_emitter_importance(kg, P, N, prim_index);
  }

  /* to-do: need to handle a case when total importance is 0. */
  if (total_emitter_importance == 0.0f) {
    return false;
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
      ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(
          light_tree_emitters, prim_index);

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
    emitter_cdf += emitter_pdf;
  }

  /* We should never reach this point. */
  assert(false);
  return -1;
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
ccl_device int light_tree_sample_distant_lights(KernelGlobals kg,
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

/* We need to be able to find the probability of selecting a given light, for MIS. */
ccl_device float light_tree_pdf(KernelGlobals kg, const float3 P, const float3 N, const int prim)
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
  int parent = kemitter->parent_index;
  ccl_global const KernelLightTreeNode *kleaf = &kernel_data_fetch(light_tree_nodes, parent);

  /* First, we find the probability of selecting the primitive out of the leaf node. */
  float total_importance = 0.0f;
  float emitter_importance = 0.0f;
  for (int i = 0; i < kleaf->num_prims; i++) {
    /* At a leaf node, the negative value is the index into first prim. */
    int prim = i - kleaf->child_index;
    const float importance = light_tree_emitter_importance(kg, P, N, prim);
    if (prim == emitter) {
      emitter_importance = importance;
    }
    total_importance += importance;
  }
  pdf *= emitter_importance / total_importance;

  /* Next, we find the probability of traversing to that leaf node. */
  int child_index = parent;
  parent = kleaf->parent_index;
  while (parent != -1) {
    const ccl_global KernelLightTreeNode *kparent = &kernel_data_fetch(light_tree_nodes, parent);

    const int left_index = parent + 1;
    const int right_index = kparent->child_index;
    const ccl_global KernelLightTreeNode *kleft = &kernel_data_fetch(light_tree_nodes, left_index);
    const ccl_global KernelLightTreeNode *kright = &kernel_data_fetch(light_tree_nodes,
                                                                      right_index);

    const float left_importance = light_tree_cluster_importance(kg, P, N, kleft);
    const float right_importance = light_tree_cluster_importance(kg, P, N, kright);
    float left_probability = left_importance / (left_importance + right_importance);

    /* If the child index matches the left index, then we must've traversed left, otherwise right.
     */
    if (left_index == child_index) {
      pdf *= left_probability;
    }
    else {
      pdf *= (1.0f - left_probability);
    }

    child_index = parent;
    parent = kparent->parent_index;
  }

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
