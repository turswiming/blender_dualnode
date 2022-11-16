/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

#include "kernel/light/area.h"
#include "kernel/light/common.h"
#include "kernel/light/light.h"
#include "kernel/light/spot.h"
#include "kernel/light/triangle.h"

CCL_NAMESPACE_BEGIN

/* TODO: this seems like a relative expensive computation, and we can make it a lot cheaper
 * by using a bounding sphere instead of a bounding box. This will be more inaccurate, but it
 * might be fine when used along with the adaptive splitting. */
ccl_device float light_tree_cos_bounding_box_angle(const float3 bbox_min,
                                                   const float3 bbox_max,
                                                   const float3 P,
                                                   const float3 N,
                                                   const float3 point_to_centroid,
                                                   ccl_private bool &bbox_is_visible)
{
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
    bbox_is_visible |= dot(point_to_corner, N) > 0;
  }
  return cos_theta_u;
}

/* This is the general function for calculating the importance of either a cluster or an emitter.
 * Both of the specialized functions obtain the necessary data before calling this function. */
template<bool in_volume_segment>
ccl_device void light_tree_cluster_importance(const float3 N_or_D,
                                              const bool has_transmission,
                                              /* unnormalized if in_volume_segment */
                                              const float3 point_to_centroid,
                                              const float cos_theta_u,
                                              const float3 bcone_axis,
                                              const float max_distance,
                                              const float min_distance,
                                              const float t,
                                              const float theta_o,
                                              const float theta_e,
                                              const float energy,
                                              ccl_private float &max_importance,
                                              ccl_private float &min_importance)
{
  max_importance = 0.0f;
  min_importance = 0.0f;

  const float sin_theta_u = safe_sqrtf(1.0f - sqr(cos_theta_u));
  float cos_theta, cos_theta_i, sin_theta_i;
  /* cos(theta_i') in the paper, omitted for volume */
  float cos_min_incidence_angle = 1.0f;
  /* when sampling the light tree for the second time in `shade_volume.h` and when query the pdf in
   * `sample.h` */
  const bool in_volume = (dot(N_or_D, N_or_D) < 5e-4f);

  if (in_volume_segment) {
    const float3 D = N_or_D;
    const float3 v0 = -normalize(point_to_centroid);
    const float3 v1 = normalize(-point_to_centroid + D * fminf(t, 1e12f));

    const float3 o0 = v0;
    float3 o1, o2;
    make_orthonormals_tangent(o0, v1, &o1, &o2);

    const float dot_o0_a = dot(o0, bcone_axis);
    const float dot_o1_a = dot(o1, bcone_axis);
    const float cos_phi0 = dot_o0_a / sqrtf(sqr(dot_o0_a) + sqr(dot_o1_a));

    /* Eq. (6) */
    cos_theta = (dot_o1_a < 0 || dot(v0, v1) > cos_phi0) ?
                    fmaxf(dot_o0_a, dot(v1, bcone_axis)) : /* b_max */
                    dot(bcone_axis, cos_phi0 * o0 + safe_sqrtf(1.0f - sqr(cos_phi0)) * o1);
  }
  else {
    cos_theta = dot(bcone_axis, -point_to_centroid);
    if (!in_volume) {
      const float3 N = N_or_D;
      cos_theta_i = has_transmission ? fabsf(dot(point_to_centroid, N)) :
                                       dot(point_to_centroid, N);
      sin_theta_i = safe_sqrtf(1.0f - sqr(cos_theta_i));

      /* cos_min_incidence_angle = cos(max{theta_i - theta_u, 0}) = cos(theta_i') in the paper */
      cos_min_incidence_angle = cos_theta_i > cos_theta_u ?
                                    1.0f :
                                    cos_theta_i * cos_theta_u + sin_theta_i * sin_theta_u;
      /* If the node is guaranteed to be behind the surface we're sampling, and the surface is
       * opaque, then we can give the node an importance of 0 as it contributes nothing to the
       * surface. This is more accurate than the bbox test if we are calculating the importance of
       * an emitter with radius */
      if (!has_transmission && cos_min_incidence_angle < 0) {
        return;
      }
    }
  }

  /* minimum angle an emitter’s axis would form with the direction to the shading point,
   * cos(theta') in the paper */
  float cos_min_outgoing_angle;
  /* cos(theta - theta_u) */
  const float sin_theta = safe_sqrtf(1.0f - sqr(cos_theta));
  const float cos_theta_minus_theta_u = cos_theta * cos_theta_u + sin_theta * sin_theta_u;

  float cos_theta_o, sin_theta_o;
  fast_sincosf(theta_o, &sin_theta_o, &cos_theta_o);

  if ((cos_theta > cos_theta_u) || (cos_theta_minus_theta_u > cos_theta_o)) {
    /* theta - theta_o - theta_u < 0 */
    kernel_assert((fast_acosf(cos_theta) - theta_o - fast_acosf(cos_theta_u)) < 5e-4f);
    cos_min_outgoing_angle = 1.0f;
  }
  else if ((cos_theta > cos_theta_u) || (theta_o + theta_e > M_PI_F) ||
           (cos_theta_minus_theta_u > cos(theta_o + theta_e))) {
    /* theta' = theta - theta_o - theta_u < theta_e */
    kernel_assert((fast_acosf(cos_theta) - theta_o - fast_acosf(cos_theta_u) - theta_e) < 5e-4f);
    const float sin_theta_minus_theta_u = safe_sqrtf(1.0f - sqr(cos_theta_minus_theta_u));
    cos_min_outgoing_angle = cos_theta_minus_theta_u * cos_theta_o +
                             sin_theta_minus_theta_u * sin_theta_o;
  }
  else {
    /* cluster invisible */
    return;
  }

  /* TODO: find a good approximation for f_a. */
  const float f_a = 1.0f;
  /* TODO: also consider t (or theta_a, theta_b) for volume */
  max_importance = fabsf(f_a * cos_min_incidence_angle * energy * cos_min_outgoing_angle /
                         (in_volume_segment ? min_distance : sqr(min_distance)));

  /* TODO: also min importance for volume? */
  if (in_volume_segment || in_volume) {
    min_importance = max_importance;
    return;
  }

  /* compute mininum importance */
  /* cos_max_incidence_angle = cos(min{theta_i + theta_u, pi}) */
  const float cos_max_incidence_angle = fmaxf(
      cos_theta_i * cos_theta_u - sin_theta_i * sin_theta_u, 0.0f);

  /* cos(theta + theta_o + theta_u) if theta + theta_o + theta_u < theta_e, 0 otherwise */
  float cos_max_outgoing_angle;
  const float cos_theta_plus_theta_u = cos_theta * cos_theta_u - sin_theta * sin_theta_u;
  if (theta_e - theta_o < 0 || cos_theta < 0 || cos_theta_u < 0 ||
      cos_theta_plus_theta_u < cos(theta_e - theta_o)) {
    min_importance = 0.0f;
  }
  else {
    const float sin_theta_plus_theta_u = safe_sqrtf(1.0f - sqr(cos_theta_plus_theta_u));
    cos_max_outgoing_angle = cos_theta_plus_theta_u * cos_theta_o -
                             sin_theta_plus_theta_u * sin_theta_o;
    min_importance = fabsf(f_a * cos_max_incidence_angle * energy * cos_max_outgoing_angle /
                           sqr(max_distance));
  }
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

  /* TODO: reservoir is disabled for now */
  return 1.0f;

  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);
  const int prim = kemitter->prim_id;

  /* Triangles are handled normally for now. */
  if (prim < 0) {
    const int lamp = -prim - 1;
    const ccl_global KernelLight *klight = &kernel_data_fetch(lights, lamp);

    /* We use a special calculation to check if a light is
     * within the bounds of a spot or area light. */
    if (klight->type == LIGHT_SPOT) {
      return spot_light_tree_weight(klight, P, N);
    }
    else if (klight->type == LIGHT_AREA) {
      return area_light_tree_weight(klight, P, N);
    }
  }

  return 1.0f;
}

template<bool in_volume_segment>
ccl_device void light_tree_emitter_importance(KernelGlobals kg,
                                              const float3 P,
                                              const float3 N_or_D,
                                              const float t,
                                              const bool has_transmission,
                                              int emitter_index,
                                              ccl_private float &max_importance,
                                              ccl_private float &min_importance)
{
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         emitter_index);

  float3 bcone_axis = make_float3(kemitter->bounding_cone_axis[0],
                                  kemitter->bounding_cone_axis[1],
                                  kemitter->bounding_cone_axis[2]);
  float theta_o = kemitter->theta_o;
  float min_distance, distance;
  float max_distance = 0.0f;
  float cos_theta_u = 1.0f;
  float3 centroid = make_float3(
      kemitter->centroid[0], kemitter->centroid[1], kemitter->centroid[2]);
  float3 point_to_centroid = safe_normalize_len(centroid - P, &distance);
  bool bbox_is_visible = has_transmission;

  const int prim = kemitter->prim_id;
  if (prim < 0) {
    const int lamp = -prim - 1;
    const ccl_global KernelLight *klight = &kernel_data_fetch(lights, lamp);

    if (klight->type == LIGHT_SPOT || klight->type == LIGHT_POINT) {
      const float radius = klight->spot.radius;
      min_distance = distance;
      max_distance = sqrtf(sqr(radius) + sqr(distance));
      const float hypotenus = max_distance;
      cos_theta_u = distance / hypotenus;
      bbox_is_visible = true; /* will be tested later */
    }
    else { /* area light */
      const float3 extentu = make_float3(
          klight->area.extentu[0], klight->area.extentu[1], klight->area.extentu[2]);
      const float3 extentv = make_float3(
          klight->area.extentv[0], klight->area.extentv[1], klight->area.extentv[2]);
      for (int i = 0; i < 4; i++) {
        const float3 corner = ((i & 1) - 0.5f) * extentu + 0.5f * ((i & 2) - 1) * extentv +
                              centroid;
        float distance_point_to_corner;
        const float3 point_to_corner = safe_normalize_len(corner - P, &distance_point_to_corner);
        cos_theta_u = fminf(cos_theta_u, dot(point_to_centroid, point_to_corner));
        bbox_is_visible |= dot(point_to_corner, N_or_D) > 0;
        max_distance = fmaxf(max_distance, distance_point_to_corner);
      }
      /* TODO: a cheap substitute for minimal distance between point and primitive. Does it worth
       * the overhead to compute the accurate minimal distance? */
      min_distance = distance;
    }
    if (klight->type == LIGHT_POINT) {
      bcone_axis = -point_to_centroid; /* disk oriented normal */
      theta_o = 0.0f;
    }
  }
  else { /* mesh light */
    const int object = kemitter->mesh_light.object_id;
    float3 V[3];
    triangle_world_space_vertices(kg, object, prim, -1.0f, V);
    for (int i = 0; i < 3; i++) {
      const float3 corner = V[i];
      float distance_point_to_corner;
      const float3 point_to_corner = safe_normalize_len(corner - P, &distance_point_to_corner);
      cos_theta_u = fminf(cos_theta_u, dot(point_to_centroid, point_to_corner));
      bbox_is_visible |= dot(point_to_corner, N_or_D) > 0;
      max_distance = fmaxf(max_distance, distance_point_to_corner);
    }
    min_distance = distance;
  }

  if (!bbox_is_visible) {
    max_importance = 0.0f;
    min_importance = 0.0f;
    return;
  }

  /* TODO: better measure for single emitter */
  if (in_volume_segment) {
    const float3 D = N_or_D;
    const float3 closest_point = P + D * dot(point_to_centroid, D);
    /* minimal distance of the ray to the cluster */
    min_distance = len(centroid - closest_point);
    max_distance = min_distance;
    point_to_centroid = centroid - P;
  }

  light_tree_cluster_importance<in_volume_segment>(N_or_D,
                                                   has_transmission,
                                                   point_to_centroid,
                                                   cos_theta_u,
                                                   bcone_axis,
                                                   max_distance,
                                                   min_distance,
                                                   t,
                                                   theta_o,
                                                   kemitter->theta_e,
                                                   kemitter->energy,
                                                   max_importance,
                                                   min_importance);
}

ccl_device bool light_tree_should_split(KernelGlobals kg,
                                        const float3 P,
                                        const ccl_global KernelLightTreeNode *knode)
{
  /* TODO: don't split because it introduces variance. Maybe delete relevant field and functions
   * later. */
  return false;

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

  if (distance <= radius) {
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

template<bool in_volume_segment>
ccl_device void light_tree_node_importance(KernelGlobals kg,
                                           const float3 P,
                                           const float3 N_or_D,
                                           const float t,
                                           const bool has_transmission,
                                           const ccl_global KernelLightTreeNode *knode,
                                           ccl_private float &max_importance,
                                           ccl_private float &min_importance)
{
  if (knode->child_index <= 0 && knode->num_prims == 1) {
    /* at a leaf node and there is only one emitter */
    light_tree_emitter_importance<in_volume_segment>(
        kg, P, N_or_D, t, has_transmission, -knode->child_index, max_importance, min_importance);
  }
  else {
    const float3 bbox_min = make_float3(
        knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
    const float3 bbox_max = make_float3(
        knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);
    const float3 bcone_axis = make_float3(
        knode->bounding_cone_axis[0], knode->bounding_cone_axis[1], knode->bounding_cone_axis[2]);

    const float3 centroid = 0.5f * (bbox_min + bbox_max);
    float distance;
    float3 point_to_centroid = normalize_len(centroid - P, &distance);
    bool bbox_is_visible = has_transmission;
    float cos_theta_u = light_tree_cos_bounding_box_angle(
        bbox_min, bbox_max, P, N_or_D, point_to_centroid, bbox_is_visible);

    /* If the node is guaranteed to be behind the surface we're sampling, and the surface is
     * opaque, then we can give the node an importance of 0 as it contributes nothing to the
     * surface. */
    if (!bbox_is_visible) {
      max_importance = 0.0f;
      min_importance = 0.0f;
      return;
    }

    if (in_volume_segment) {
      const float3 D = N_or_D;
      const float3 closest_point = P + D * dot(point_to_centroid, D);
      /* minimal distance of the ray to the cluster */
      distance = len(centroid - closest_point);
      point_to_centroid = centroid - P;
    }
    /* clamp distance to half the radius of the cluster when splitting is disabled */
    distance = fmaxf(0.5f * len(centroid - bbox_max), distance);

    /* TODO: currently max_distance = min_distance, max_importance = min_importance for the nodes.
     * Do we need better weights for complex scenes? */
    light_tree_cluster_importance<in_volume_segment>(N_or_D,
                                                     has_transmission,
                                                     point_to_centroid,
                                                     cos_theta_u,
                                                     bcone_axis,
                                                     distance,
                                                     distance,
                                                     t,
                                                     knode->theta_o,
                                                     knode->theta_e,
                                                     knode->energy,
                                                     max_importance,
                                                     min_importance);
  }
}

ccl_device_inline void sample_resevoir(const int current_index,
                                       const float current_weight,
                                       ccl_private int &selected_index,
                                       ccl_private float &selected_weight,
                                       ccl_private float &total_weight,
                                       ccl_private float *rand)
{
  if (current_weight == 0.0f) {
    return;
  }
  total_weight += current_weight;
  float thresh = current_weight / total_weight;
  if (*rand < thresh) {
    selected_index = current_index;
    selected_weight = current_weight;
    *rand = *rand / thresh;
  }
  else {
    *rand = (*rand - thresh) / (1.0f - thresh);
  }
  return;
}

/* pick an emitter from a leaf node using resevoir sampling, keep two reservoirs for upper and
 * lower bounds */
template<bool in_volume_segment>
ccl_device int light_tree_cluster_select_emitter(KernelGlobals kg,
                                                 ccl_private float *rand,
                                                 const float3 P,
                                                 const float3 N_or_D,
                                                 const float t,
                                                 const bool has_transmission,
                                                 const ccl_global KernelLightTreeNode *knode,
                                                 ccl_private float *pdf_factor)
{
  float2 importance;
  float2 selected_importance = zero_float2();
  float2 total_importance = zero_float2();
  int selected_index = -1;

  /* Mark emitters with zero importance. Used for resevoir when total minimum importance = 0 */
  kernel_assert(knode->num_prims <= sizeof(uint) * 8);
  uint has_importance = 0;

  bool sample_max = (*rand > 0.5f); /* sampling using the maximum importance */
  *rand = *rand * 2.0f - float(sample_max);

  for (int i = 0; i < knode->num_prims; i++) {
    int current_index = -knode->child_index + i;
    /* maximum importance = importance[0], mininum importance = importance[1] */
    light_tree_emitter_importance<in_volume_segment>(
        kg, P, N_or_D, t, has_transmission, current_index, importance[0], importance[1]);

    sample_resevoir(current_index,
                    importance[!sample_max],
                    selected_index,
                    selected_importance[!sample_max],
                    total_importance[!sample_max],
                    rand);
    if (selected_index == current_index) {
      selected_importance[sample_max] = importance[sample_max];
    }
    total_importance[sample_max] += importance[sample_max];

    has_importance |= ((importance[0] > 0) << i);
  }

  if (total_importance[0] == 0.0f) {
    return -1;
  }

  if (total_importance[1] == 0.0f) {
    /* uniformly sample emitters with positive maximum importance */
    if (sample_max) {
      selected_importance[1] = 1.0f;
      total_importance[1] = float(popcount(has_importance));
    }
    else {
      selected_index = -1;
      for (int i = 0; i < knode->num_prims; i++) {
        int current_index = -knode->child_index + i;
        sample_resevoir(current_index,
                        float(has_importance & 1),
                        selected_index,
                        selected_importance[1],
                        total_importance[1],
                        rand);
        has_importance >>= 1;
      }
      light_tree_emitter_importance<in_volume_segment>(kg,
                                                       P,
                                                       N_or_D,
                                                       t,
                                                       has_transmission,
                                                       selected_index,
                                                       selected_importance[0],
                                                       importance[1]);
    }
  }

  *pdf_factor = average(selected_importance / total_importance);

  return selected_index;
}

template<bool in_volume_segment>
ccl_device_inline bool get_left_probability(KernelGlobals kg,
                                            const float3 P,
                                            const float3 N_or_D,
                                            const float t,
                                            const bool has_transmission,
                                            const int left_index,
                                            const int right_index,
                                            ccl_private float &left_probability)
{
  /* If we don't split, then we need to choose sampling between the left or right child. */
  const ccl_global KernelLightTreeNode *left = &kernel_data_fetch(light_tree_nodes, left_index);
  const ccl_global KernelLightTreeNode *right = &kernel_data_fetch(light_tree_nodes, right_index);

  float min_left_importance, max_left_importance, min_right_importance, max_right_importance;
  light_tree_node_importance<in_volume_segment>(
      kg, P, N_or_D, t, has_transmission, left, max_left_importance, min_left_importance);
  light_tree_node_importance<in_volume_segment>(
      kg, P, N_or_D, t, has_transmission, right, max_right_importance, min_right_importance);

  const float total_max_importance = max_left_importance + max_right_importance;
  if (total_max_importance == 0.0f) {
    return false;
  }
  const float total_min_importance = min_left_importance + min_right_importance;

  /* average two probabilities of picking the left child node using lower and upper bounds */
  const float probability_max = max_left_importance / total_max_importance;
  const float probability_min = total_min_importance > 0 ?
                                    min_left_importance / total_min_importance :
                                    0.5f * (float(max_left_importance > 0) +
                                            float(max_right_importance == 0.0f));
  left_probability = 0.5f * (probability_max + probability_min);
  return true;
}

template<bool in_volume_segment>
ccl_device bool light_tree_sample(KernelGlobals kg,
                                  ccl_private const RNGState *rng_state,
                                  ccl_private float *randu,
                                  const float randv,
                                  const float time,
                                  const float3 P,
                                  const float3 N_or_D,
                                  const float t,
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
      const int selected_light = light_tree_cluster_select_emitter<in_volume_segment>(
          kg, randu, P, N_or_D, t, has_transmission, knode, &pdf_emitter_selection);

      if (selected_light < 0) {
        stack_index--;
        continue;
      }

      const float light_reservoir_weight = light_tree_emitter_reservoir_weight(
          kg, P, N_or_D, selected_light);

      /* TODO: make pdf_node_emitter_selection part of the light_reservoir_weight, otherwise result
       * is suboptimal, or disable splitting and remove reservoir-related code . */
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

    float left_probability;
    if (!get_left_probability<in_volume_segment>(
            kg, P, N_or_D, t, has_transmission, left_index, right_index, left_probability)) {
      stack_index--;
      continue;
    }

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

template<bool in_volume_segment>
ccl_device float light_tree_distant_light_importance(KernelGlobals kg,
                                                     const float3 N,
                                                     const bool has_transmission,
                                                     const int index)
{
  if (in_volume_segment) {
    return 0.0f;
  }
  ccl_global const KernelLightTreeDistantEmitter *kdistant = &kernel_data_fetch(
      light_tree_distant_group, index);

  if (kdistant->energy == 0.0f) {
    return 0.0f;
  }

  const float3 light_axis = make_float3(
      kdistant->direction[0], kdistant->direction[1], kdistant->direction[2]);
  float theta_i = fast_acosf(dot(N, -light_axis));

  /* If the light is guaranteed to be behind the surface we're sampling, and the surface is
   * opaque, then we can give the light an importance of 0 as it contributes nothing to the
   * surface. */
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
                                                 const float t,
                                                 const bool has_transmission,
                                                 const int path_flag,
                                                 const int bounce,
                                                 ccl_private LightSample *ls,
                                                 ccl_private float *pdf_factor)
{
  kernel_assert(!in_volume_segment);
  /* TODO: do single loop over lights to avoid computing importance twice? */
  const int num_distant_lights = kernel_data.integrator.num_distant_lights;
  float total_importance = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    total_importance += light_tree_distant_light_importance<in_volume_segment>(
        kg, N, has_transmission, i);
  }

  if (total_importance == 0.0f) {
    return false;
  }

  float light_cdf = 0.0f;
  for (int i = 0; i < num_distant_lights; i++) {
    const float light_pdf = light_tree_distant_light_importance<in_volume_segment>(
                                kg, N, has_transmission, i) /
                            total_importance;
    if (*randu <= light_cdf + light_pdf) {
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

  /* TODO: change implementation so that even under precision issues, we always end
   * up selecting a light and this can't happen? */
  kernel_assert(false);
  return false;
}

/* We need to be able to find the probability of selecting a given light for MIS.
 * Returns a conditioned pdf which is consistent with the pdf computed in `light_tree_sample()`
 */
ccl_device float light_tree_pdf(KernelGlobals kg,
                                ConstIntegratorState state,
                                const float3 P,
                                const float3 N,
                                const int path_flag,
                                const int prim)
{
  const bool has_transmission = (path_flag & PATH_RAY_MIS_HAD_TRANSMISSION);
  float distant_light_importance = light_tree_distant_light_importance<false>(
      kg, N, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    float discard;
    light_tree_node_importance<false>(
        kg, P, N, 0, has_transmission, kroot, light_tree_importance, discard);
  }
  const float total_group_importance = light_tree_importance + distant_light_importance;
  kernel_assert(total_group_importance != 0.0f);
  float pdf = light_tree_importance / total_group_importance;

  const int target_emitter = (prim >= 0) ? kernel_data_fetch(triangle_to_tree, prim) :
                                           kernel_data_fetch(light_to_tree, ~prim);
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         target_emitter);
  const int target_leaf = kemitter->parent_index;
  ccl_global const KernelLightTreeNode *kleaf = &kernel_data_fetch(light_tree_nodes, target_leaf);

  uint bit_trail = kleaf->bit_trail;

  /* We generate a random number to use for selecting a light. */
  RNGState rng_state;
  path_state_rng_load(state, &rng_state);
  /* to-do: is this the correct macro to use? */
  float randu = path_state_rng_1D(kg, &rng_state, PRNG_LIGHT);

  /* Keep track of the currently selected primitive weight,
   * as well as the total weight as part of the weighted reservoir sampling. */
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

    /* TODO: interior nodes are processed before leaf nodes, but current implementation forces
     * the reader to think about the leaf node first. Maybe switch the two if-branches for better
     * understandability? */

    /* Leaf node */
    if (knode->child_index <= 0) {

      /* If the leaf node contains the target emitter, we are processing the last node.
       * We then iterate through the lights to find the target emitter.
       * Otherwise, we randomly select one. */
      if (index == target_leaf) {
        float target_max_importance = 0.0f;
        float target_min_importance = 0.0f;
        float total_max_importance = 0.0f;
        float total_min_importance = 0.0f;
        int num_has_importance = 0;
        for (int i = 0; i < knode->num_prims; i++) {
          const int emitter = -knode->child_index + i;
          float max_importance, min_importance;
          light_tree_emitter_importance<false>(
              kg, P, N, 0, has_transmission, emitter, max_importance, min_importance);
          num_has_importance += (max_importance > 0);
          if (emitter == target_emitter) {
            target_max_importance = max_importance;
            target_min_importance = min_importance;
            selected_reservoir_weight = light_tree_emitter_reservoir_weight(kg, P, N, emitter);
            total_reservoir_weight += selected_reservoir_weight;
          }
          total_max_importance += max_importance;
          total_min_importance += min_importance;
        }

        if (target_max_importance > 0.0f) {
          const float pdf_emitter_selection = 0.5f *
                                              (target_max_importance / total_max_importance +
                                               (total_min_importance > 0 ?
                                                    target_min_importance / total_min_importance :
                                                    1.0f / num_has_importance));
          const float pdf_reservoir = selected_reservoir_weight / total_reservoir_weight;
          pdf *= pdf_emitter_selection * pdf_reservoir;
        }
        else {
          pdf = 0.0f;
        }
      }
      else {
        int selected_light = -1;
        float pdf_emitter_selection = 1.0f;
        selected_light = light_tree_cluster_select_emitter<false>(
            kg, &randu, P, N, 0, has_transmission, knode, &pdf_emitter_selection);

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

    float left_probability;
    if (!get_left_probability<false>(
            kg, P, N, 0, has_transmission, left_index, right_index, left_probability)) {
      pdf = 0.0f;
      stack_index--;
      continue;
    }

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

ccl_device float light_tree_pdf_distant(
    KernelGlobals kg, const float3 P, const float3 N, const int path_flag, const int prim)
{
  const bool has_transmission = (path_flag & PATH_RAY_MIS_HAD_TRANSMISSION);
  float distant_light_importance = light_tree_distant_light_importance<false>(
      kg, N, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    float discard;
    light_tree_node_importance<false>(
        kg, P, N, 0, has_transmission, kroot, light_tree_importance, discard);
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
    float importance = light_tree_distant_light_importance<false>(kg, N, has_transmission, i);
    if (i == distant_light) {
      emitter_importance = importance;
    }
    total_importance += importance;
  }

  pdf *= emitter_importance / total_importance;
  return pdf;
}

template<bool in_volume_segment>
ccl_device bool light_tree_sample(KernelGlobals kg,
                                  ccl_private const RNGState *rng_state,
                                  float randu,
                                  const float randv,
                                  const float time,
                                  const float3 P,
                                  const float3 N_or_D,
                                  const float t,
                                  const int shader_flags,
                                  const int bounce,
                                  const uint32_t path_flag,
                                  ccl_private LightSample *ls)
{
  const bool has_transmission = (shader_flags & SD_BSDF_HAS_TRANSMISSION);
  kernel_assert(!in_volume_segment || (in_volume_segment && has_transmission));
  /* TODO: add distant lights to light tree (as one child node of the root node) */
  float distant_light_importance = light_tree_distant_light_importance<in_volume_segment>(
      kg, N_or_D, has_transmission, kernel_data.integrator.num_distant_lights);
  float light_tree_importance = 0.0f;
  if (kernel_data.integrator.num_distribution > kernel_data.integrator.num_distant_lights) {
    const ccl_global KernelLightTreeNode *kroot = &kernel_data_fetch(light_tree_nodes, 0);
    float discard;
    light_tree_node_importance<in_volume_segment>(
        kg, P, N_or_D, t, has_transmission, kroot, light_tree_importance, discard);
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
    ret = light_tree_sample<in_volume_segment>(kg,
                                               rng_state,
                                               &randu,
                                               randv,
                                               time,
                                               P,
                                               N_or_D,
                                               t,
                                               has_transmission,
                                               path_flag,
                                               bounce,
                                               ls,
                                               &pdf_factor);
  }
  else {
    randu = (randu - light_tree_probability) / (1.0f - light_tree_probability);
    pdf_factor *= (1.0f - light_tree_probability);
    ret = light_tree_sample_distant_lights<in_volume_segment>(kg,
                                                              rng_state,
                                                              &randu,
                                                              randv,
                                                              time,
                                                              P,
                                                              N_or_D,
                                                              t,
                                                              has_transmission,
                                                              path_flag,
                                                              bounce,
                                                              ls,
                                                              &pdf_factor);
  }
  ls->pdf_selection = pdf_factor;
  ls->pdf *= pdf_factor;
  return ret;
}

CCL_NAMESPACE_END
