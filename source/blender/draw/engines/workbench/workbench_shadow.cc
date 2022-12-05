/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 *
 * Shadow:
 *
 * Use stencil shadow buffer to cast a sharp shadow over opaque surfaces.
 *
 * After the main pre-pass we render shadow volumes using custom depth & stencil states to
 * set the stencil of shadowed area to anything but 0.
 *
 * Then the shading pass will shade the areas with stencil not equal 0 differently.
 */

#include "BKE_object.h"
#include "BLI_math.h"
#include "DRW_render.h"

#include "workbench_private.hh"

#define DEBUG_SHADOW_VOLUME 0

namespace blender::workbench {

static void compute_parallel_lines_nor_and_dist(const float2 v1,
                                                const float2 v2,
                                                const float2 v3,
                                                float2 r_line[2])
{
  r_line[0] = v2 - v1;
  /* Find orthogonal vector. */
  SWAP(float, r_line[0].x, r_line[0].y);
  r_line[0].x *= -1;
  /* Edge distances. */
  r_line[1].x = math::dot(r_line[0], v1);
  r_line[1].y = math::dot(r_line[0], v3);
  /* Make sure r_line[1].x is the minimum. */
  if (r_line[1].x > r_line[1].y) {
    SWAP(float, r_line[1].x, r_line[1].y);
  }
}

PassMain::Sub *&ShadowPass::get_pass_ptr(bool depth_pass, bool manifold, bool cap /*= false*/)
{
  return passes_[depth_pass][manifold][cap];
}

GPUShader *ShadowPass::get_shader(bool depth_pass, bool manifold, bool cap /*= false*/)
{
  GPUShader *&shader = shaders[depth_pass][manifold][cap];

  if (shader == nullptr) {
    std::string create_info_name = "workbench_next_shadow";
    create_info_name += (depth_pass) ? "_pass" : "_fail";
    create_info_name += (manifold) ? "_manifold" : "_no_manifold";
    create_info_name += (cap) ? "_caps" : "_no_caps";
#if DEBUG_SHADOW_VOLUME
    create_info_name += "_debug";
#endif
    shader = GPU_shader_create_from_info_name(create_info_name.c_str());
  }
  return shader;
}

void ShadowPass::init(const SceneState &scene_state, SceneResources &resources)
{
  enabled_ = scene_state.draw_shadows;
  if (!enabled_) {
    resources.world_buf.shadow_mul = 0.0f;
    resources.world_buf.shadow_add = 1.0f;
    return;
  }
  const Scene &scene = *scene_state.scene;

  direction_ws = scene.display.light_direction;
  /* Turn the light in a way where it's more user friendly to control. */
  SWAP(float, direction_ws.y, direction_ws.z);
  direction_ws *= float3(-1, 1, -1);

  /* Shadow direction. */
  float4x4 view_matrix;
  DRW_view_viewmat_get(NULL, view_matrix.ptr(), false);
  resources.world_buf.shadow_direction_vs = float4(view_matrix.ref_3x3() * direction_ws);

  /* Clamp to avoid overshadowing and shading errors. */
  float focus = clamp_f(scene.display.shadow_focus, 0.0001f, 0.99999f);
  resources.world_buf.shadow_shift = scene.display.shadow_shift;
  resources.world_buf.shadow_focus = 1.0f - focus * (1.0f - resources.world_buf.shadow_shift);
  resources.world_buf.shadow_mul = scene_state.shading.shadow_intensity;
  resources.world_buf.shadow_add = 1.0f - resources.world_buf.shadow_mul;
}

void ShadowPass::update()
{
  changed = !compare_v3v3(direction_ws, cached_direction, 1e-5f);

  if (changed) {
    const float3 up = {0.0f, 0.0f, 1.0f};
    matrix = float4x4::identity();

    /* TODO: fix singularity. */
    copy_v3_v3(matrix[2], direction_ws);
    cross_v3_v3v3(matrix[0], matrix[2], up);
    normalize_v3(matrix[0]);
    cross_v3_v3v3(matrix[1], matrix[2], matrix[0]);

    matrix_inv = matrix.inverted();

    cached_direction = direction_ws;
  }

  float planes[6][4];
  DRW_culling_frustum_planes_get(nullptr, planes);
  /* we only need the far plane. */
  far_plane = planes[2];

  BoundBox frustum_corners;
  DRW_culling_frustum_corners_get(nullptr, &frustum_corners);

  float3 shadow_near_corners[4];
  shadow_near_corners[0] = matrix_inv * float3(frustum_corners.vec[0]);
  shadow_near_corners[1] = matrix_inv * float3(frustum_corners.vec[3]);
  shadow_near_corners[2] = matrix_inv * float3(frustum_corners.vec[7]);
  shadow_near_corners[3] = matrix_inv * float3(frustum_corners.vec[4]);

  INIT_MINMAX(near_min, near_max);
  for (int i = 0; i < 4; i++) {
    minmax_v3v3_v3(near_min, near_max, shadow_near_corners[i]);
  }

  compute_parallel_lines_nor_and_dist(float2(shadow_near_corners[0]),
                                      float2(shadow_near_corners[1]),
                                      float2(shadow_near_corners[2]),
                                      near_sides[0]);
  compute_parallel_lines_nor_and_dist(float2(shadow_near_corners[1]),
                                      float2(shadow_near_corners[2]),
                                      float2(shadow_near_corners[0]),
                                      near_sides[1]);
}

void ShadowPass::sync()
{
  if (!enabled_) {
    return;
  }

  update();

#if DEBUG_SHADOW_VOLUME
  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_BLEND_ADD_FULL;
  DRWState depth_pass_state = state | DRW_STATE_DEPTH_LESS;
  DRWState depth_fail_state = state | DRW_STATE_DEPTH_GREATER_EQUAL;
#else
  DRWState state = DRW_STATE_DEPTH_LESS | DRW_STATE_STENCIL_ALWAYS;
  DRWState depth_pass_state = state | DRW_STATE_WRITE_STENCIL_SHADOW_PASS;
  DRWState depth_fail_state = state | DRW_STATE_WRITE_STENCIL_SHADOW_FAIL;
#endif

  /* TODO(fclem): Merge into one pass with sub-passes. */
  pass_ps.init();
  pass_ps.state_set(depth_pass_state);
  pass_ps.state_stencil(0xFF, 0xFF, 0xFF);

  /* TODO(Miguel Pozo) */
  pass_ps.clear_stencil(0);

  fail_ps.init();
  fail_ps.state_set(depth_fail_state);
  fail_ps.state_stencil(0xFF, 0xFF, 0xFF);

  /* Stencil Shadow passes. */
  for (bool manifold : {false, true}) {
    {
      PassMain::Sub *&ps = get_pass_ptr(true, manifold);
      ps = &pass_ps.sub(manifold ? "manifold" : "non_manifold");
      ps->shader_set(get_shader(true, manifold));
    }
    {
      PassMain::Sub *&ps = get_pass_ptr(false, manifold, false);
      ps = &fail_ps.sub(manifold ? "NoCaps.manifold" : "NoCaps.non_manifold");
      ps->shader_set(get_shader(false, manifold, false));
    }
    {
      PassMain::Sub *&ps = get_pass_ptr(false, manifold, true);
      ps = &fail_ps.sub(manifold ? "Caps.manifold" : "Caps.non_manifold");
      ps->shader_set(get_shader(false, manifold, true));
    }
  }
}

void ShadowPass::object_sync(Manager &manager,
                             ObjectRef &ob_ref,
                             SceneState &scene_state,
                             const bool has_transp_mat)
{
  if (!enabled_) {
    return;
  }

  Object *ob = ob_ref.object;
  bool is_manifold;
  GPUBatch *geom_shadow = DRW_cache_object_edge_detection_get(ob, &is_manifold);
  if (geom_shadow == nullptr) {
    return;
  }

  ObjectData *engine_object_data = (ObjectData *)DRW_drawdata_ensure(
      &ob->id, &draw_engine_workbench_next, sizeof(ObjectData), &ObjectData::init, nullptr);

  ObjectShadowData &object_data = engine_object_data->shadow_data;

  if (object_data.cast_visible_shadow(ob, *this)) {
    mul_v3_mat3_m4v3(object_data.direction, ob->world_to_object, direction_ws);

    bool use_shadow_pass_technique = !object_data.camera_in_object_shadow(ob, *this);

    /* Shadow pass technique needs object to be have all its surface opaque. */
    if (has_transp_mat) {
      use_shadow_pass_technique = false;
    }

    /* TODO (Miguel Pozo):
     * Disable use_shadow_pass_tecnique when there are "in front" objects in the scene.
     */

    /* We cannot use Shadow Pass technique on non-manifold object (see T76168). */
    if (use_shadow_pass_technique && !is_manifold && (scene_state.cull_state != 0)) {
      use_shadow_pass_technique = false;
    }

    if (use_shadow_pass_technique) {
      PassMain::Sub &ps = *get_pass_ptr(true, is_manifold);
      ps.push_constant("lightDirection", object_data.direction);
      ps.push_constant("lightDistance", 1e5f);
      ResourceHandle handle = manager.resource_handle(ob_ref);
      ps.draw(geom_shadow, handle);
#if DEBUG_SHADOW_VOLUME
      DRW_debug_bbox(&object_data.bbox, float4(1.0f, 0.0f, 0.0f, 1.0f));
#endif
    }
    else {
      float extrude_distance = object_data.shadow_distance(ob, *this);

      /* TODO(fclem): only use caps if they are in the view frustum. */
      const bool need_caps = true;
      if (need_caps) {
        PassMain::Sub &ps = *get_pass_ptr(false, is_manifold, true);
        ps.push_constant("lightDirection", object_data.direction);
        ps.push_constant("lightDistance", extrude_distance);
        ResourceHandle handle = manager.resource_handle(ob_ref);
        ps.draw(DRW_cache_object_surface_get(ob), handle);
      }

      PassMain::Sub &ps = *get_pass_ptr(false, is_manifold, false);
      ps.push_constant("lightDirection", object_data.direction);
      ps.push_constant("lightDistance", extrude_distance);
      ResourceHandle handle = manager.resource_handle(ob_ref);
      ps.draw(geom_shadow, handle);
#if DEBUG_SHADOW_VOLUME
      DRW_debug_bbox(&object_data.bbox, float4(1.0f, 0.0f, 0.0f, 1.0f));
#endif
    }
  }
}

void ShadowPass::draw(Manager &manager, View &view, SceneResources &resources, int2 resolution)
{
  if (!enabled_) {
    return;
  }

  Framebuffer fb = {"FBShadows"};
  fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx),
            GPU_ATTACHMENT_TEXTURE(resources.color_tx));
  fb.bind();

  manager.submit(pass_ps, view);
  manager.submit(fail_ps, view);
}

void ObjectShadowData::init()
{
  bbox_dirty = true;
}

const BoundBox *ObjectShadowData::get_bbox(Object *ob, ShadowPass &shadow_pass)
{
  if (bbox_dirty || shadow_pass.changed) {
    float4x4 matrix = shadow_pass.matrix_inv * ob->object_to_world;

    /* Get AABB in shadow space. */
    INIT_MINMAX(min, max);

    /* From object space to shadow space */
    const BoundBox *_bbox = BKE_object_boundbox_get(ob);
    for (int i : IndexRange(8)) {
      float3 corner = matrix * float3(_bbox->vec[i]);
      math::min_max(corner, min, max);
    }
    depth = max.z - min.z;
    /* Extend towards infinity. */
    max.z += 1e4f;

    /* Get extended AABB in world space. */
    BKE_boundbox_init_from_minmax(&bbox, min, max);
    for (int i = 0; i < 8; i++) {
      mul_m4_v3(shadow_pass.matrix.ptr(), bbox.vec[i]);
    }
    bbox_dirty = false;
  }

  return &bbox;
}

bool ObjectShadowData::cast_visible_shadow(Object *ob, ShadowPass &shadow_pass)
{
  const BoundBox *_bbox = get_bbox(ob, shadow_pass);
  const DRWView *default_view = DRW_view_default_get();
  return DRW_culling_box_test(default_view, _bbox);
}

float ObjectShadowData::shadow_distance(Object *ob, ShadowPass &shadow_pass)
{
  const BoundBox *_bbox = get_bbox(ob, shadow_pass);
  float dist = 1e4f;

  for (int corner : {0, 3, 4, 7}) {
    float dist_isect;
    bool isect = isect_ray_plane_v3(_bbox->vec[corner],
                                    shadow_pass.cached_direction,
                                    shadow_pass.far_plane,
                                    &dist_isect,
                                    true);
    if (!isect) {
      /* All rays are parallels. If one fails, the other will too. */
      break;
    }
    if (dist_isect < dist) {
      dist = dist_isect;
    }
  }
  return max_ii(dist - depth, 0);
}

bool ObjectShadowData::camera_in_object_shadow(Object *ob, ShadowPass &shadow_pass)
{
  /* Just to be sure the min, max are updated. */
  get_bbox(ob, shadow_pass);
  /* Test if near plane is in front of the shadow. */
  if (min.z > shadow_pass.near_max.z) {
    return false;
  }

  /* Separation Axis Theorem test */

  /* Test bbox sides first (faster) */
  if ((min.x > shadow_pass.near_max.x) || (max.x < shadow_pass.near_min.x) ||
      (min.y > shadow_pass.near_max.y) || (max.y < shadow_pass.near_min.y)) {
    return false;
  }
  /* Test projected near rectangle sides */
  const float2 pts[4] = {
      {min.x, min.y},
      {min.x, max.y},
      {max.x, min.y},
      {max.x, max.y},
  };

  for (int i : IndexRange(2)) {
    float min_dst = FLT_MAX, max_dst = -FLT_MAX;
    for (int j = 0; j < 4; j++) {
      float dst = math::dot(shadow_pass.near_sides[i][0], pts[j]);
      /* Do min max */
      if (min_dst > dst) {
        min_dst = dst;
      }
      if (max_dst < dst) {
        max_dst = dst;
      }
    }

    if ((shadow_pass.near_sides[i][1].x > max_dst) || (shadow_pass.near_sides[i][1].y < min_dst)) {
      return false;
    }
  }
  /* No separation axis found. Both shape intersect. */
  return true;
}

}  // namespace blender::workbench
