/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#include "scene/light_tree.h"
#include "scene/mesh.h"
#include "scene/object.h"

CCL_NAMESPACE_BEGIN

float OrientationBounds::calculate_measure() const
{
  float theta_w = fminf(M_PI_F, theta_o + theta_e);
  float cos_theta_o = cosf(theta_o);
  float sin_theta_o = sinf(theta_o);

  return M_2PI_F * (1 - cos_theta_o) +
         M_PI_2_F * (2 * theta_w * sin_theta_o - cosf(theta_o - 2 * theta_w) -
                     2 * theta_o * sin_theta_o + cos_theta_o);
}

OrientationBounds merge(const OrientationBounds &cone_a, const OrientationBounds &cone_b)
{
  if (is_zero(cone_a.axis)) {
    return cone_b;
  }
  if (is_zero(cone_b.axis)) {
    return cone_a;
  }

  /* Set cone a to always have the greater theta_o. */
  const OrientationBounds *a = &cone_a;
  const OrientationBounds *b = &cone_b;
  if (cone_b.theta_o > cone_a.theta_o) {
    a = &cone_b;
    b = &cone_a;
  }

  float theta_d = safe_acosf(dot(a->axis, b->axis));
  float theta_e = fmaxf(a->theta_e, b->theta_e);

  /* Return axis and theta_o of a if it already contains b. */
  /* This should also be called when b is empty. */
  if (a->theta_o >= fminf(M_PI_F, theta_d + b->theta_o)) {
    return OrientationBounds({a->axis, a->theta_o, theta_e});
  }

  /* Compute new theta_o that contains both a and b. */
  float theta_o = (theta_d + a->theta_o + b->theta_o) * 0.5f;

  if (theta_o >= M_PI_F) {
    return OrientationBounds({a->axis, M_PI_F, theta_e});
  }

  /* Rotate new axis to be between a and b. */
  float theta_r = theta_o - a->theta_o;
  float3 new_axis = rotate_around_axis(a->axis, cross(a->axis, b->axis), theta_r);
  new_axis = normalize(new_axis);

  return OrientationBounds({new_axis, theta_o, theta_e});
}

LightTreePrimitive::LightTreePrimitive(Scene *scene,
                                       int prim_id,
                                       int object_id,
                                       bool is_double_sided)
    : prim_id(prim_id), object_id(object_id), is_double_sided(is_double_sided)
{
  if (is_triangle()) {
    calculate_triangle_vertices(scene);
  }
  calculate_centroid(scene);
  calculate_bbox(scene);
  calculate_bcone(scene);
  calculate_energy(scene);
}

void LightTreePrimitive::calculate_triangle_vertices(Scene *scene)
{
  assert(is_triangle());
  Object *object = scene->objects[object_id];
  Mesh *mesh = static_cast<Mesh *>(object->get_geometry());
  Mesh::Triangle triangle = mesh->get_triangle(prim_id);

  for (int i = 0; i < 3; i++) {
    vertices[i] = mesh->get_verts()[triangle.v[i]];
  }

  /* instanced mesh lights have not applied their transform at this point.
   * in this case, these points have to be transformed to get the proper
   * spatial bound. */
  if (!mesh->transform_applied) {
    const Transform &tfm = object->get_tfm();
    for (int i = 0; i < 3; i++) {
      vertices[i] = transform_point(&tfm, vertices[i]);
    }
  }
}

void LightTreePrimitive::calculate_centroid(Scene *scene)
{
  if (is_triangle()) {
    /* NOTE: the original implementation used the bounding box centroid, but primitive centroid
     * seems to work fine */
    centroid = (vertices[0] + vertices[1] + vertices[2]) / 3.0f;
  }
  else {
    centroid = scene->lights[object_id]->get_co();
  }
}

void LightTreePrimitive::calculate_bbox(Scene *scene)
{
  bbox = BoundBox::empty;

  if (is_triangle()) {
    for (int i = 0; i < 3; i++) {
      bbox.grow(vertices[i]);
    }
  }
  else {
    Light *lamp = scene->lights[object_id];
    LightType type = lamp->get_light_type();
    const float size = lamp->get_size();

    if (type == LIGHT_POINT || type == LIGHT_SPOT) {
      /* Point and spot lights can emit light from any point within its radius. */
      const float3 radius = make_float3(size);
      bbox.grow(centroid - radius);
      bbox.grow(centroid + radius);
    }
    else if (type == LIGHT_AREA) {
      /* For an area light, sizeu and sizev determine the 2 dimensions of the area light,
       * while axisu and axisv determine the orientation of the 2 dimensions.
       * We want to add all 4 corners to our bounding box. */
      const float3 half_extentu = 0.5 * lamp->get_sizeu() * lamp->get_axisu() * size;
      const float3 half_extentv = 0.5 * lamp->get_sizev() * lamp->get_axisv() * size;

      bbox.grow(centroid + half_extentu + half_extentv);
      bbox.grow(centroid + half_extentu - half_extentv);
      bbox.grow(centroid - half_extentu + half_extentv);
      bbox.grow(centroid - half_extentu - half_extentv);
    }
    else {
      /* No bounding box for distant lights */
    }
  }
}

void LightTreePrimitive::calculate_bcone(Scene *scene)
{
  bcone = OrientationBounds::empty;

  if (is_triangle()) {
    if (is_double_sided) {
      /* Any vector in the plane */
      bcone.axis = safe_normalize(vertices[0] - vertices[1]);
      bcone.theta_o = M_PI_2_F;
    }
    else {
      bcone.axis = safe_normalize(cross(vertices[1] - vertices[0], vertices[2] - vertices[0]));
      bcone.theta_o = 0;
    }

    bcone.theta_e = M_PI_2_F;
  }
  else {
    Light *lamp = scene->lights[object_id];
    LightType type = lamp->get_light_type();

    bcone.axis = normalize(lamp->get_dir());

    if (type == LIGHT_POINT) {
      bcone.theta_o = M_PI_F;
      bcone.theta_e = M_PI_2_F;
    }
    else if (type == LIGHT_SPOT) {
      bcone.theta_o = 0;
      bcone.theta_e = lamp->get_spot_angle() * 0.5f;
    }
    else if (type == LIGHT_AREA) {
      bcone.theta_o = 0;
      bcone.theta_e = lamp->get_spread() * 0.5f;
    }
    else if (type == LIGHT_DISTANT) {
      bcone.theta_o = tanf(0.5f * lamp->get_angle());
      bcone.theta_e = 0;
    }
    else if (type == LIGHT_BACKGROUND) {
      /* Set an arbitrary direction for the background light. */
      bcone.axis = make_float3(0.0f, 0.0f, 1.0f);
      /* TODO: this may depend on portal lights as well. */
      bcone.theta_o = M_PI_F;
      bcone.theta_e = 0;
    }
  }
}

void LightTreePrimitive::calculate_energy(Scene *scene)
{
  if (is_triangle()) {
    Object *object = scene->objects[object_id];
    Mesh *mesh = static_cast<Mesh *>(object->get_geometry());
    Shader *shader = static_cast<Shader *>(mesh->get_used_shaders()[mesh->get_shader()[prim_id]]);

    /* to-do: need a better way to handle this when textures are used. */
    float3 shader_estimate;
    shader->estimate_emission(shader_estimate);

    float area = triangle_area(vertices[0], vertices[1], vertices[2]);
    energy = area * scene->shader_manager->linear_rgb_to_gray(shader_estimate);
  }
  else {
    Light *lamp = scene->lights[object_id];
    LightType type = lamp->get_light_type();

    float3 strength = lamp->get_strength();
    if (type == LIGHT_AREA) {
      strength *= 0.25f; /* eval_fac scaling in `area.h` */
    }
    else if (type == LIGHT_SPOT || type == LIGHT_POINT) {
      strength *= 0.25f * M_1_PI_F; /* eval_fac scaling in `spot.h` and `point.h` */
    }
    else if (type == LIGHT_BACKGROUND) {
      /* integrate over cosine-weighted hemisphere */
      strength *= lamp->get_average_radiance() * M_PI_F;
    }

    if (lamp->get_shader()) {
      float3 shader_estimate;
      lamp->get_shader()->estimate_emission(shader_estimate);
      strength *= shader_estimate;
    }

    energy = scene->shader_manager->linear_rgb_to_gray(strength);
  }
}

void LightTreeBuildNode::init_leaf(const uint &offset,
                                   const uint &n,
                                   const BoundBox &b,
                                   const OrientationBounds &c,
                                   const float &e,
                                   const float &e_var,
                                   const uint &bits)
{
  bbox = b;
  bcone = c;
  energy = e;
  energy_variance = e_var;
  first_prim_index = offset;
  num_lights = n;

  children[0] = children[1] = nullptr;
  bit_trail = bits;
  is_leaf = true;
}

void LightTreeBuildNode::init_interior(LightTreeBuildNode *c0,
                                       LightTreeBuildNode *c1,
                                       const BoundBox &b,
                                       const OrientationBounds &c,
                                       const float &e,
                                       const float &e_var,
                                       const uint &bits)
{
  bbox = b;
  bcone = c;
  energy = e;
  energy_variance = e_var;
  first_prim_index = 0;
  num_lights = 0;

  children[0] = c0;
  children[1] = c1;
  bit_trail = bits;
  is_leaf = false;
}

LightTree::LightTree(const vector<LightTreePrimitive> &prims,
                     Scene *scene,
                     uint max_lights_in_leaf)
{
  if (prims.empty()) {
    return;
  }

  prims_ = prims;
  scene_ = scene;
  max_lights_in_leaf_ = max_lights_in_leaf;

  for (int i = 0; i < prims.size(); i++) {
    prims_[i].prim_num = i;
  }

  int total_nodes = 0;
  vector<LightTreePrimitive> ordered_prims;
  LightTreeBuildNode *root = recursive_build(0, prims.size(), total_nodes, ordered_prims, 0, 0);
  prims_ = ordered_prims;

  int offset = 0;
  nodes_.resize(total_nodes);
  flatten_tree(root, offset, -1);
}

const vector<LightTreePrimitive> &LightTree::get_prims() const
{
  return prims_;
}

const vector<PackedLightTreeNode> &LightTree::get_nodes() const
{
  return nodes_;
}

LightTreeBuildNode *LightTree::recursive_build(int start,
                                               int end,
                                               int &total_nodes,
                                               vector<LightTreePrimitive> &ordered_prims,
                                               uint bit_trail,
                                               int depth)
{
  LightTreeBuildNode *node = new LightTreeBuildNode();
  total_nodes++;
  BoundBox node_bbox = BoundBox::empty;
  OrientationBounds node_bcone = OrientationBounds::empty;
  BoundBox centroid_bounds = BoundBox::empty;
  float energy_total = 0.0;
  float energy_squared_total = 0.0;
  int num_prims = end - start;

  for (int i = start; i < end; i++) {
    const LightTreePrimitive &prim = prims_.at(i);
    node_bbox.grow(prim.bbox);
    node_bcone = merge(node_bcone, prim.bcone);
    centroid_bounds.grow(prim.centroid);

    energy_total += prim.energy;
    energy_squared_total += sqr(prim.energy);
  }

  /* Var(X) = E[X^2] - E[X]^2 */
  float energy_variance = (energy_squared_total / num_prims) - sqr(energy_total / num_prims);

  bool try_splitting = num_prims > 1 && len(centroid_bounds.size()) > 0.0f;
  int min_dim = -1, min_bucket = 0;
  bool should_split = false;
  if (try_splitting) {
    /* Find the best place to split the primitives into 2 nodes.
     * If the best split cost is no better than making a leaf node, make a leaf instead.*/
    float min_cost = min_split_saoh(
        centroid_bounds, start, end, node_bbox, node_bcone, min_dim, min_bucket);
    should_split = num_prims > max_lights_in_leaf_ || min_cost < energy_total;
  }
  if (should_split) {
    int middle;

    if (min_dim != -1) {
      /* Partition the primitives between start and end into the appropriate split,
       * based on the minimum dimension and minimum bucket returned from split_saoh.
       * This is an O(n) algorithm where we iterate from the left and right side,
       * and swaps the appropriate left and right elements until complete. */
      int left = start, right = end - 1;
      float bounding_dimension = (min_bucket + 1) * (centroid_bounds.size()[min_dim] /
                                                     LightTreeBucketInfo::num_buckets) +
                                 centroid_bounds.min[min_dim];
      while (left < right) {
        while (prims_[left].centroid[min_dim] <= bounding_dimension && left < end) {
          left++;
        }

        while (prims_[right].centroid[min_dim] > bounding_dimension && right >= start) {
          right--;
        }

        if (left < right) {
          std::swap(prims_[left].prim_num, prims_[right].prim_num);
          std::swap(prims_[left], prims_[right]);
        }
      }

      middle = left;
    }
    else {
      /* Degenerate case with many lights in the same place. */
      middle = (start + end) / 2;
    }

    /* At this point, we should expect right to be just past left,
     * where left points to the first element that belongs to the right node. */
    LightTreeBuildNode *left_node = recursive_build(
        start, middle, total_nodes, ordered_prims, bit_trail, depth + 1);
    LightTreeBuildNode *right_node = recursive_build(
        middle, end, total_nodes, ordered_prims, bit_trail | (1u << depth), depth + 1);
    node->init_interior(
        left_node, right_node, node_bbox, node_bcone, energy_total, energy_variance, bit_trail);
  }
  else {
    int first_prim_offset = ordered_prims.size();
    for (int i = start; i < end; i++) {
      int prim_num = prims_[i].prim_num;
      ordered_prims.push_back(prims_[prim_num]);
    }
    node->init_leaf(first_prim_offset,
                    num_prims,
                    node_bbox,
                    node_bcone,
                    energy_total,
                    energy_variance,
                    bit_trail);
  }
  return node;
}

float LightTree::min_split_saoh(const BoundBox &centroid_bbox,
                                int start,
                                int end,
                                const BoundBox &bbox,
                                const OrientationBounds &bcone,
                                int &min_dim,
                                int &min_bucket)
{
  /* Even though this factor is used for every bucket, we use it to compare
   * the min_cost and total_energy (when deciding between creating a leaf or interior node. */
  const float bbox_area = bbox.area();
  const bool has_area = bbox_area != 0.0f;
  const float total_area = has_area ? bbox_area : len(bbox.size());
  const float total_cost = total_area * bcone.calculate_measure();
  if (total_cost == 0.0f) {
    return FLT_MAX;
  }

  const float inv_total_cost = 1.0f / total_cost;
  const float3 extent = centroid_bbox.size();
  const float max_extent = max4(extent.x, extent.y, extent.z, 0.0f);

  /* Check each dimension to find the minimum splitting cost. */
  float min_cost = FLT_MAX;
  for (int dim = 0; dim < 3; dim++) {
    /* If the centroid bounding box is 0 along a given dimension, skip it. */
    if (centroid_bbox.size()[dim] == 0.0f) {
      continue;
    }

    const float inv_extent = 1 / (centroid_bbox.size()[dim]);

    /* Fill in buckets with primitives. */
    vector<LightTreeBucketInfo> buckets(LightTreeBucketInfo::num_buckets);
    for (int i = start; i < end; i++) {
      const LightTreePrimitive &prim = prims_[i];

      /* Place primitive into the appropriate bucket,
       * where the centroid box is split into equal partitions. */
      int bucket_idx = LightTreeBucketInfo::num_buckets *
                       (prim.centroid[dim] - centroid_bbox.min[dim]) * inv_extent;
      if (bucket_idx == LightTreeBucketInfo::num_buckets) {
        bucket_idx = LightTreeBucketInfo::num_buckets - 1;
      }

      buckets[bucket_idx].count++;
      buckets[bucket_idx].energy += prim.energy;
      buckets[bucket_idx].bbox.grow(prim.bbox);
      buckets[bucket_idx].bcone = merge(buckets[bucket_idx].bcone, prim.bcone);
    }

    /* Calculate the cost of splitting at each point between partitions. */
    vector<float> bucket_costs(LightTreeBucketInfo::num_buckets - 1);
    float energy_L, energy_R;
    BoundBox bbox_L, bbox_R;
    OrientationBounds bcone_L, bcone_R;
    for (int split = 0; split < LightTreeBucketInfo::num_buckets - 1; split++) {
      energy_L = 0;
      energy_R = 0;
      bbox_L = BoundBox::empty;
      bbox_R = BoundBox::empty;
      bcone_L = OrientationBounds::empty;
      bcone_R = OrientationBounds::empty;

      for (int left = 0; left <= split; left++) {
        if (buckets[left].bbox.valid()) {
          energy_L += buckets[left].energy;
          bbox_L.grow(buckets[left].bbox);
          bcone_L = merge(bcone_L, buckets[left].bcone);
        }
      }

      for (int right = split + 1; right < LightTreeBucketInfo::num_buckets; right++) {
        if (buckets[right].bbox.valid()) {
          energy_R += buckets[right].energy;
          bbox_R.grow(buckets[right].bbox);
          bcone_R = merge(bcone_R, buckets[right].bcone);
        }
      }

      /* Calculate the cost of splitting using the heuristic as described in the paper. */
      const float area_L = has_area ? bbox_L.area() : len(bbox_L.size());
      const float area_R = has_area ? bbox_R.area() : len(bbox_R.size());
      float left = (bbox_L.valid()) ? energy_L * area_L * bcone_L.calculate_measure() : 0.0f;
      float right = (bbox_R.valid()) ? energy_R * area_R * bcone_R.calculate_measure() : 0.0f;
      float regularization = max_extent * inv_extent;
      bucket_costs[split] = regularization * (left + right) * inv_total_cost;

      if (bucket_costs[split] < min_cost) {
        min_cost = bucket_costs[split];
        min_dim = dim;
        min_bucket = split;
      }
    }
  }
  return min_cost;
}

int LightTree::flatten_tree(const LightTreeBuildNode *node, int &offset, int parent)
{
  PackedLightTreeNode *current_node = &nodes_[offset];
  current_node->bbox = node->bbox;
  current_node->bcone = node->bcone;
  current_node->energy = node->energy;
  current_node->energy_variance = node->energy_variance;
  current_node->bit_trail = node->bit_trail;
  int current_index = offset;
  offset++;

  /* If current node contains lights, then it is a leaf node.
   * Otherwise, create interior node and children recursively. */
  if (node->num_lights > 0) {
    current_node->first_prim_index = node->first_prim_index;
    current_node->num_lights = node->num_lights;
    current_node->is_leaf_node = true;
  }
  else {
    current_node->num_lights = 0;
    current_node->is_leaf_node = false;

    /* The first child is located directly to the right of the parent. */
    flatten_tree(node->children[0], offset, current_index);
    current_node->second_child_index = flatten_tree(node->children[1], offset, current_index);
  }

  return current_index;
}

CCL_NAMESPACE_END
