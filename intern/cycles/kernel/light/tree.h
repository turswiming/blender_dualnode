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

  max_importance = 0.0f;
  min_importance = 0.0f;
  float theta_o = kemitter->theta_o;
  float min_distance, distance;
  float max_distance = 0.0f;
  float cos_theta_u = 1.0f;
  float3 bcone_axis, centroid, point_to_centroid;
  bool bbox_is_visible = has_transmission;

  const int prim = kemitter->prim_id;
  /* TODO: pack in functions and move to header files for respective light types */
  if (prim < 0) {
    const int lamp = ~prim;
    const ccl_global KernelLight *klight = &kernel_data_fetch(lights, lamp);
    centroid = make_float3(klight->co[0], klight->co[1], klight->co[2]);
    point_to_centroid = safe_normalize_len(centroid - P, &distance);

    if (klight->type == LIGHT_SPOT || klight->type == LIGHT_POINT) {
      bcone_axis = make_float3(klight->spot.dir[0], klight->spot.dir[1], klight->spot.dir[2]);
      const float radius = klight->spot.radius;
      min_distance = distance;
      max_distance = sqrtf(sqr(radius) + sqr(distance));
      const float hypotenus = max_distance;
      cos_theta_u = distance / hypotenus;
      bbox_is_visible = true; /* will be tested when computing the importance */
    }
    else if (klight->type == LIGHT_AREA) {
      bcone_axis = make_float3(klight->area.dir[0], klight->area.dir[1], klight->area.dir[2]);
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
    else { /* distant light */
      if (in_volume_segment) {
        return;
      }
      if (klight->type == LIGHT_DISTANT) {
        /* Treating it as a disk light 1 unit away */
        cos_theta_u = fast_cosf(theta_o);
        theta_o = 0.0f;
        max_distance = 1.0f / cos_theta_u;
      }
      else {
        /* Set an arbitrary direction for the background light. */
        centroid = -N_or_D;
        cos_theta_u = -1.0f;
        max_distance = 1.0f;
      }
      bcone_axis = centroid;
      point_to_centroid = -centroid;
      min_distance = 1.0f;
      bbox_is_visible = true; /* will be tested when computing the importance */
    }
    if (klight->type == LIGHT_POINT) {
      bcone_axis = -point_to_centroid; /* disk oriented normal */
      theta_o = 0.0f;
    }
  }
  else { /* mesh light */
    const int object = kemitter->mesh_light.object_id;
    float3 vertices[3];
    triangle_world_space_vertices(kg, object, prim, -1.0f, vertices);
    centroid = (vertices[0] + vertices[1] + vertices[2]) / 3.0f;
    point_to_centroid = safe_normalize_len(centroid - P, &distance);
    bcone_axis = safe_normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
    if (kemitter->mesh_light.emission_sampling == EMISSION_SAMPLING_BACK) {
      bcone_axis = -bcone_axis;
    }
    else if (kemitter->mesh_light.emission_sampling == EMISSION_SAMPLING_FRONT_BACK) {
      bcone_axis *= -signf(dot(bcone_axis, point_to_centroid));
    }
    theta_o = 0.0f;

    for (int i = 0; i < 3; i++) {
      const float3 corner = vertices[i];
      float distance_point_to_corner;
      const float3 point_to_corner = safe_normalize_len(corner - P, &distance_point_to_corner);
      cos_theta_u = fminf(cos_theta_u, dot(point_to_centroid, point_to_corner));
      bbox_is_visible |= dot(point_to_corner, N_or_D) > 0;
      max_distance = fmaxf(max_distance, distance_point_to_corner);
    }
    min_distance = distance;
  }

  /* TODO: rename `bbox_is_visible` to `is_visible`? */
  if (!bbox_is_visible) {
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
  max_importance = 0.0f;
  min_importance = 0.0f;
  if (knode->num_prims == 1) {
    /* At a leaf node with only one emitter */
    light_tree_emitter_importance<in_volume_segment>(
        kg, P, N_or_D, t, has_transmission, -knode->child_index, max_importance, min_importance);
  }
  else if (knode->num_prims != 0) {
    const float3 bcone_axis = make_float3(
        knode->bounding_cone_axis[0], knode->bounding_cone_axis[1], knode->bounding_cone_axis[2]);

    float3 point_to_centroid;
    float cos_theta_u;
    float distance;
    if (knode->bit_trail == 1) {
      /* distant light node */
      if (in_volume_segment) {
        return;
      }
      point_to_centroid = -bcone_axis;
      cos_theta_u = fmaxf(fast_cosf(knode->theta_o), 0.0f);
      distance = 1.0f;
    }
    else {
      const float3 bbox_min = make_float3(
          knode->bounding_box_min[0], knode->bounding_box_min[1], knode->bounding_box_min[2]);
      const float3 bbox_max = make_float3(
          knode->bounding_box_max[0], knode->bounding_box_max[1], knode->bounding_box_max[2]);

      const float3 centroid = 0.5f * (bbox_min + bbox_max);

      point_to_centroid = normalize_len(centroid - P, &distance);
      bool bbox_is_visible = has_transmission;
      cos_theta_u = light_tree_cos_bounding_box_angle(
          bbox_min, bbox_max, P, N_or_D, point_to_centroid, bbox_is_visible);

      /* If the node is guaranteed to be behind the surface we're sampling, and the surface is
       * opaque, then we can give the node an importance of 0 as it contributes nothing to the
       * surface. */
      if (!bbox_is_visible) {
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
    }
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
                                       ccl_private float &rand)
{
  if (current_weight == 0.0f) {
    return;
  }
  total_weight += current_weight;
  float thresh = current_weight / total_weight;
  if (rand <= thresh) {
    selected_index = current_index;
    selected_weight = current_weight;
    rand = rand / thresh;
  }
  else {
    rand = (rand - thresh) / (1.0f - thresh);
  }
  kernel_assert(rand >= 0.0f && rand <= 1.0f);
  return;
}

/* pick an emitter from a leaf node using resevoir sampling, keep two reservoirs for upper and
 * lower bounds */
template<bool in_volume_segment>
ccl_device int light_tree_cluster_select_emitter(KernelGlobals kg,
                                                 ccl_private float &rand,
                                                 const float3 P,
                                                 const float3 N_or_D,
                                                 const float t,
                                                 const bool has_transmission,
                                                 const ccl_global KernelLightTreeNode *knode,
                                                 ccl_private float *pdf_factor)
{
  float selected_importance[2] = {0.0f, 0.0f};
  float total_importance[2] = {0.0f, 0.0f};
  int selected_index = -1;

  /* Mark emitters with zero importance. Used for resevoir when total minimum importance = 0 */
  kernel_assert(knode->num_prims <= sizeof(uint) * 8);
  uint has_importance = 0;

  bool sample_max = (rand > 0.5f); /* sampling using the maximum importance */
  rand = rand * 2.0f - float(sample_max);

  for (int i = 0; i < knode->num_prims; i++) {
    int current_index = -knode->child_index + i;
    /* maximum importance = importance[0], mininum importance = importance[1] */
    float importance[2];
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

      float discard;
      light_tree_emitter_importance<in_volume_segment>(
          kg, P, N_or_D, t, has_transmission, selected_index, selected_importance[0], discard);
    }
  }

  *pdf_factor = 0.5f * (selected_importance[0] / total_importance[0] +
                        selected_importance[1] / total_importance[1]);

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
                                  float randu,
                                  float randv,
                                  const float time,
                                  const float3 P,
                                  const float3 N_or_D,
                                  const float t,
                                  const int shader_flags,
                                  const int bounce,
                                  const uint32_t path_flag,
                                  ccl_private LightSample *ls)
{
  if (!kernel_data.integrator.use_direct_light) {
    return false;
  }

  const bool has_transmission = (shader_flags & SD_BSDF_HAS_TRANSMISSION);
  float pdf_leaf = 1.0f;
  float pdf_emitter_from_leaf = 1.0f;
  int selected_light = -1;

  int node_index = 0; /* root node */

  /* Traverse the light tree until a leaf node is reached. */
  while (true) {
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, node_index);

    if (knode->child_index <= 0) {
      /* At a leaf node, we pick an emitter */
      selected_light = light_tree_cluster_select_emitter<in_volume_segment>(
          kg, randv, P, N_or_D, t, has_transmission, knode, &pdf_emitter_from_leaf);
      break;
    }

    /* At an interior node, the left child is directly after the parent,
     * while the right child is stored as the child index. */
    const int left_index = node_index + 1;
    const int right_index = knode->child_index;

    float left_probability;
    if (!get_left_probability<in_volume_segment>(
            kg, P, N_or_D, t, has_transmission, left_index, right_index, left_probability)) {
      return false; /* both child nodes have zero importance */
    }

    if (randu < left_probability) { /* go left */
      kernel_assert(left_probability > 0.0f);

      node_index = left_index;
      randu /= left_probability;
      pdf_leaf *= left_probability;
    }
    else {
      kernel_assert((1.0f - left_probability) > 0.0f);

      node_index = right_index;
      randu = (randu - left_probability) / (1.0f - left_probability);
      pdf_leaf *= (1.0f - left_probability);
    }
  }

  /* TODO: check `spot_light_tree_weight()` and `area_light_tree_weight()` */
  if (selected_light < 0) {
    return false;
  }

  /* Sample a point on the chosen emitter */
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         selected_light);

  /* TODO: this is the same code as light_distribution_sample, except the index is determined
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
    if (!triangle_light_sample<in_volume_segment>(kg, prim, object, randu, randv, time, ls, P)) {
      return false;
    }
    ls->shader |= mesh_shader_flag;
  }
  else {
    if (UNLIKELY(light_select_reached_max_bounces(kg, ~prim, bounce))) {
      return false;
    }

    if (!light_sample<in_volume_segment>(kg, ~prim, randu, randv, P, path_flag, ls)) {
      return false;
    }
  }

  ls->pdf_selection = pdf_leaf * pdf_emitter_from_leaf;
  ls->pdf *= ls->pdf_selection;
  return (ls->pdf > 0);
}

/* We need to be able to find the probability of selecting a given light for MIS. */
ccl_device float light_tree_pdf(
    KernelGlobals kg, const float3 P, const float3 N, const int path_flag, const int prim)
{
  const bool has_transmission = (path_flag & PATH_RAY_MIS_HAD_TRANSMISSION);
  /* Target emitter info */
  const int target_emitter = (prim >= 0) ? kernel_data_fetch(triangle_to_tree, prim) :
                                           kernel_data_fetch(light_to_tree, ~prim);
  ccl_global const KernelLightTreeEmitter *kemitter = &kernel_data_fetch(light_tree_emitters,
                                                                         target_emitter);
  const int target_leaf = kemitter->parent_index;
  ccl_global const KernelLightTreeNode *kleaf = &kernel_data_fetch(light_tree_nodes, target_leaf);
  uint bit_trail = kleaf->bit_trail;

  int node_index = 0; /* root node */

  float pdf = 1.0f;

  /* Traverse the light tree until we reach the target leaf node */
  while (true) {
    const ccl_global KernelLightTreeNode *knode = &kernel_data_fetch(light_tree_nodes, node_index);

    if (knode->child_index <= 0) {
      break;
    }

    /* Interior node */
    const int left_index = node_index + 1;
    const int right_index = knode->child_index;

    float left_probability;
    if (!get_left_probability<false>(
            kg, P, N, 0, has_transmission, left_index, right_index, left_probability)) {
      return 0.0f;
    }

    const bool go_left = (bit_trail & 1) == 0;
    bit_trail >>= 1;
    pdf *= go_left ? left_probability : (1.0f - left_probability);
    node_index = go_left ? left_index : right_index;

    if (pdf == 0) {
      return 0.0f;
    }
  }

  kernel_assert(node_index == target_leaf);

  /* Iterate through leaf node to find the probability of sampling the target emitter. */
  float target_max_importance = 0.0f;
  float target_min_importance = 0.0f;
  float total_max_importance = 0.0f;
  float total_min_importance = 0.0f;
  int num_has_importance = 0;
  for (int i = 0; i < kleaf->num_prims; i++) {
    const int emitter = -kleaf->child_index + i;
    float max_importance, min_importance;
    light_tree_emitter_importance<false>(
        kg, P, N, 0, has_transmission, emitter, max_importance, min_importance);
    num_has_importance += (max_importance > 0);
    if (emitter == target_emitter) {
      target_max_importance = max_importance;
      target_min_importance = min_importance;
    }
    total_max_importance += max_importance;
    total_min_importance += min_importance;
  }

  if (target_max_importance > 0.0f) {
    return pdf * 0.5f *
           (target_max_importance / total_max_importance +
            (total_min_importance > 0 ? target_min_importance / total_min_importance :
                                        1.0f / num_has_importance));
  }
  return 0.0f;
}

CCL_NAMESPACE_END
