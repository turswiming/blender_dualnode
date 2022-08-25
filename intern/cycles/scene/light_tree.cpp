/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#include "scene/mesh.h"
#include "scene/object.h"
#include "scene/light_tree.h"

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

OrientationBounds merge(const OrientationBounds& cone_a, 
						const OrientationBounds& cone_b)
{
  if (is_zero(cone_a.axis)) {
    return cone_b;
  }
  else if (is_zero(cone_b.axis)) {
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
  else {
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
}

BoundBox LightTreePrimitive::calculate_bbox(Scene *scene) const
{
  BoundBox bbox = BoundBox::empty;

  if (prim_id >= 0) {
    Object *object = scene->objects[object_id];
    Mesh *mesh = static_cast<Mesh*>(object->get_geometry());
    Mesh::Triangle triangle = mesh->get_triangle(prim_id);

    float3 p[3] = {mesh->get_verts()[triangle.v[0]],
                   mesh->get_verts()[triangle.v[1]],
                   mesh->get_verts()[triangle.v[2]]};

    /* instanced mesh lights have not applied their transform at this point.
     * in this case, these points have to be transformed to get the proper
     * spatial bound. */
    if (!mesh->transform_applied) {
      const Transform &tfm = object->get_tfm();
      for (int i = 0; i < 3; i++) {
        p[i] = transform_point(&tfm, p[i]);
      }
    }

    for (int i = 0; i < 3; i++) {
      bbox.grow(p[i]);
    }
  }
  else {
    Light *lamp = scene->lights[lamp_id];
    LightType type = lamp->get_light_type();
    const float3 center = lamp->get_co();
    const float size = lamp->get_size();

    if (type == LIGHT_POINT || type == LIGHT_SPOT) {
      /* Point and spot lights can emit light from any point within its radius. */
      const float3 radius = make_float3(size);
      bbox.grow(center - radius);
      bbox.grow(center + radius);
    }
    else if (type == LIGHT_AREA) {
      /* For an area light, sizeu and sizev determine the 2 dimensions of the area light,
       * while axisu and axisv determine the orientation of the 2 dimensions. 
       * We want to add all 4 corners to our bounding box. */
      const float3 half_extentu = 0.5 * lamp->get_sizeu() * lamp->get_axisu() * size;
      const float3 half_extentv = 0.5 * lamp->get_sizev() * lamp->get_axisv() * size;

      bbox.grow(center + half_extentu + half_extentv);
      bbox.grow(center + half_extentu - half_extentv);
      bbox.grow(center - half_extentu + half_extentv);
      bbox.grow(center - half_extentu - half_extentv);
    }
    else {
      /* This should never be reached during construction. */
      assert(false);
    }
  }

  return bbox;
}

OrientationBounds LightTreePrimitive::calculate_bcone(Scene *scene) const
{
  OrientationBounds bcone = OrientationBounds::empty;

  if (prim_id >= 0) {
    Object *object = scene->objects[object_id];
    Mesh *mesh = static_cast<Mesh *>(object->get_geometry());
    Mesh::Triangle triangle = mesh->get_triangle(prim_id);

    float3 p[3] = {mesh->get_verts()[triangle.v[0]],
                   mesh->get_verts()[triangle.v[1]],
                   mesh->get_verts()[triangle.v[2]]};

    /* instanced mesh lights have not applied their transform at this point.
     * in this case, these points have to be transformed to get the proper
     * spatial bound. */
    if (!mesh->transform_applied) {
      const Transform &tfm = object->get_tfm();
      for (int i = 0; i < 3; i++) {
        p[i] = transform_point(&tfm, p[i]);
      }
    }

    float3 normal = cross(p[1] - p[0], p[2] - p[0]);
    const float normlen = len(normal);
    if (normlen == 0.0f) {
      normal = make_float3(1.0f, 0.0f, 0.0f);
    }
    normal = normal / normlen;

    bcone.axis = normal;

    /* to-do: is there a better way to handle this case where both sides of the triangle are visible? 
     * Right now, we assume that the normal axis is within pi radians of the triangle normal. */
    bcone.theta_o = M_PI_F;
    bcone.theta_e = M_PI_2_F;
  }
  else {
    Light *lamp = scene->lights[lamp_id];
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
    else {
      /* This should never be reached during construction. */
      assert(false);
    }
  }
  
  return bcone;
}

float LightTreePrimitive::calculate_energy(Scene *scene) const
{
  float3 strength = make_float3(0.0f);

  if (prim_id >= 0) {
    Object *object = scene->objects[object_id];
    Mesh *mesh = static_cast<Mesh *>(object->get_geometry());
    Mesh::Triangle triangle = mesh->get_triangle(prim_id);
    Shader *shader = static_cast<Shader*>(mesh->get_used_shaders()[mesh->get_shader()[prim_id]]);

    /* to-do: need a better way to handle this when textures are used. */
    if (!shader->is_constant_emission(&strength)) {
      strength = make_float3(1.0f);
    }

     float3 p[3] = {mesh->get_verts()[triangle.v[0]],
                   mesh->get_verts()[triangle.v[1]],
                   mesh->get_verts()[triangle.v[2]]};

    /* instanced mesh lights have not applied their transform at this point.
     * in this case, these points have to be transformed to get the proper
     * spatial bound. */
    if (!mesh->transform_applied) {
      const Transform &tfm = object->get_tfm();
      for (int i = 0; i < 3; i++) {
        p[i] = transform_point(&tfm, p[i]);
      }
    }

    float area = triangle_area(p[0], p[1], p[2]);
    strength *= area;
  }
  else {
    Light *lamp = scene->lights[lamp_id];
    strength = lamp->get_strength();
  }

  return fabsf(scene->shader_manager->linear_rgb_to_gray(strength));
}

void LightTreeBuildNode::init_leaf(
    uint offset, uint n, const BoundBox &b, const OrientationBounds &c, float e, float e_var, uint bits)
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

void LightTreeBuildNode::init_interior(LightTreeBuildNode *c0, LightTreeBuildNode *c1)
{
  bbox = merge(c0->bbox, c1->bbox);
  bcone = merge(c0->bcone, c1->bcone);
  energy = c0->energy + c1->energy;
  energy_variance = c0->energy_variance + c1->energy_variance;
  first_prim_index = 0;
  num_lights = 0;

  children[0] = c0;
  children[1] = c1;
  is_leaf = false;
}

LightTree::LightTree(const vector<LightTreePrimitive> &prims, Scene *scene, uint max_lights_in_leaf)
{
  if (prims.size() == 0) {
    return;
  }

  prims_ = prims;
  scene_ = scene;
  max_lights_in_leaf_ = max_lights_in_leaf;

  vector<LightTreePrimitiveInfo> build_data;
  for (int i = 0; i < prims.size(); i++) {
    LightTreePrimitiveInfo prim_info;
    prim_info.bbox = prims[i].calculate_bbox(scene);
    prim_info.bcone = prims[i].calculate_bcone(scene);
    prim_info.energy = prims[i].calculate_energy(scene);
    prim_info.centroid = prim_info.bbox.center();
    prim_info.prim_num = i;
    build_data.push_back(prim_info);
  }

  int total_nodes = 0;
  vector<LightTreePrimitive> ordered_prims;
  LightTreeBuildNode *root = recursive_build(build_data, 0, prims.size(), total_nodes, ordered_prims, 0, 0);
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

LightTreeBuildNode *LightTree::recursive_build(vector<LightTreePrimitiveInfo> &primitive_info,
                                               int start,
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
    const LightTreePrimitiveInfo &prim = primitive_info.at(i);
    node_bbox.grow(prim.bbox);
    node_bcone = merge(node_bcone, prim.bcone);
    centroid_bounds.grow(prim.centroid);

    energy_total += prim.energy;
    energy_squared_total += prim.energy * prim.energy;
  }

  /* Var(X) = E[X^2] - E[X]^2 */
  float energy_variance = (energy_squared_total / num_prims) - (energy_total / num_prims) * (energy_total / num_prims);
  if (num_prims == 1 || len(centroid_bounds.size()) == 0.0f) {
    int first_prim_offset = ordered_prims.size();
    for (int i = start; i < end; i++) {
      int prim_num = primitive_info[i].prim_num;
      ordered_prims.push_back(prims_[prim_num]);
    }
    node->init_leaf(first_prim_offset, num_prims, node_bbox, node_bcone, energy_total, energy_variance, bit_trail);
  }
  else {
    /* Find the best place to split the primitives into 2 nodes. 
     * If the best split cost is no better than making a leaf node, make a leaf instead.*/
    float min_cost;
    int min_dim, min_bucket;
    split_saoh(centroid_bounds,
               primitive_info,
               start,
               end,
               node_bbox,
               node_bcone,
               min_cost,
               min_dim,
               min_bucket);
    if (num_prims > max_lights_in_leaf_ || min_cost < energy_total) {
      /* Partition the primitives between start and end into the appropriate split,
       * based on the minimum dimension and minimum bucket returned from split_saoh. 
       * This is an O(n) algorithm where we iterate from the left and right side,
       * and swaps the appropriate left and right elements until complete. */
      int left = start, right = end - 1;
      float bounding_dimension = (min_bucket + 1) * (centroid_bounds.size()[min_dim] /
                                               LightTreeBucketInfo::num_buckets) +
                                 centroid_bounds.min[min_dim];
      while (left < right) {
        while (primitive_info[left].centroid[min_dim] <= bounding_dimension && left < end) {
          left++;
        }

        while (primitive_info[right].centroid[min_dim] > bounding_dimension && right >= start) {
          right--;
        }

        if (left < right) {
          LightTreePrimitiveInfo temp = primitive_info[left];
          primitive_info[left] = primitive_info[right];
          primitive_info[right] = temp;
        }
      }

      /* At this point, we should expect right to be just past left,
       * where left points to the first element that belongs to the right node. */
      LightTreeBuildNode *left_node = recursive_build(
          primitive_info, start, left, total_nodes, ordered_prims, bit_trail, depth + 1);
      LightTreeBuildNode *right_node = recursive_build(
          primitive_info, left, end, total_nodes, ordered_prims, bit_trail | (1u << bit_trail), depth + 1);
      node->init_interior(left_node, right_node);
    }
    else {
      int first_prim_offset = ordered_prims.size();
      for (int i = start; i < end; i++) {
        int prim_num = primitive_info[i].prim_num;
        ordered_prims.push_back(prims_[prim_num]);
      }
      node->init_leaf(first_prim_offset, num_prims, node_bbox, node_bcone, energy_total, energy_variance, bit_trail);
    }
  }

  return node;
}

void LightTree::split_saoh(const BoundBox &centroid_bbox,
                           const vector<LightTreePrimitiveInfo> &primitive_info,
                           int start,
                           int end,
                           const BoundBox &bbox,
                           const OrientationBounds &bcone,
                           float& min_cost,
                           int& min_dim,
                           int& min_bucket)
{
  /* Even though this factor is used for every bucket, we use it to compare
   * the min_cost and total_energy (when deciding between creating a leaf or interior node. */
  const float inv_total_cost = 1 / (bbox.area() * bcone.calculate_measure());
  const float3 extent = centroid_bbox.size();
  const float max_extent = max4(extent.x, extent.y, extent.z, 0.0f);

  /* Check each dimension to find the minimum splitting cost. */
  min_cost = FLT_MAX;
  for (int dim = 0; dim < 3; dim++) {
    /* If the centroid bounding box is 0 along a given dimension, skip it. */
    if (centroid_bbox.size()[dim] == 0.0f) {
      continue;
    }

    const float inv_extent = 1 / (centroid_bbox.size()[dim]);
    
    /* Fill in buckets with primitives. */
    vector<LightTreeBucketInfo> buckets(LightTreeBucketInfo::num_buckets);
    for (int i = start; i < end; i++) {
      const LightTreePrimitiveInfo &primitive = primitive_info[i];

      /* Place primitive into the appropriate bucket,
       * where the centroid box is split into equal partitions. */
      int bucket_idx = LightTreeBucketInfo::num_buckets *
                       (primitive.centroid[dim] - centroid_bbox.min[dim]) * inv_extent;
      if (bucket_idx == LightTreeBucketInfo::num_buckets)
      {
        bucket_idx = LightTreeBucketInfo::num_buckets - 1;
      }

      buckets[bucket_idx].count++;
      buckets[bucket_idx].energy += primitive.energy;
      buckets[bucket_idx].bbox.grow(primitive.bbox);
      buckets[bucket_idx].bcone = merge(buckets[bucket_idx].bcone, primitive.bcone);
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
      float left = (bbox_L.valid()) ? energy_L * bbox_L.area() * bcone_L.calculate_measure() : 0.0f;
      float right = (bbox_R.valid()) ? energy_R * bbox_R.area() * bcone_R.calculate_measure() : 0.0f;
      float regularization = max_extent * inv_extent;
      bucket_costs[split] = regularization * (left + right) * inv_total_cost;

      if (bucket_costs[split] < min_cost) {
        min_cost = bucket_costs[split];
        min_dim = dim;
        min_bucket = split;
      }
    }
  }
}

int LightTree::flatten_tree(const LightTreeBuildNode *node, int &offset, int parent)
{
  PackedLightTreeNode *current_node = &nodes_[offset];
  current_node->bbox = node->bbox;
  current_node->bcone = node->bcone;
  current_node->energy = node->energy;
  current_node->energy_variance = node->energy_variance;
  int current_index = offset;
  offset++;

  /* If current node contains lights, then it is a leaf node.
   * Otherwise, create interior node and children recursively. */
  if (node->num_lights > 0) {
    current_node->first_prim_index = node->first_prim_index;
    current_node->num_lights = node->num_lights;
    current_node->bit_trail = node->bit_trail;
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