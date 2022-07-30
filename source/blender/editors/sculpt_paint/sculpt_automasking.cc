/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup edsculpt
 */

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_hash.h"
#include "BLI_index_range.hh"
#include "BLI_math.h"
#include "BLI_task.h"

#include "DNA_brush_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_brush.h"
#include "BKE_context.h"
#include "BKE_mesh.h"
#include "BKE_mesh_mapping.h"
#include "BKE_object.h"
#include "BKE_paint.h"
#include "BKE_pbvh.h"
#include "BKE_scene.h"

#include "DEG_depsgraph.h"

#include "WM_api.h"
#include "WM_message.h"
#include "WM_toolsystem.h"
#include "WM_types.h"

#include "ED_object.h"
#include "ED_screen.h"
#include "ED_sculpt.h"
#include "paint_intern.h"
#include "sculpt_intern.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "bmesh.h"

#include <cmath>
#include <cstdlib>

using blender::IndexRange;

AutomaskingCache *SCULPT_automasking_active_cache_get(SculptSession *ss)
{
  if (ss->cache) {
    return ss->cache->automasking;
  }
  if (ss->filter_cache) {
    return ss->filter_cache->automasking;
  }
  return nullptr;
}

bool SCULPT_is_automasking_mode_enabled(const Sculpt *sd,
                                        const Brush *br,
                                        const eAutomasking_flag mode)
{
  if (br) {
    return br->automasking_flags & mode || sd->automasking_flags & mode;
  }
  return sd->automasking_flags & mode;
}

bool SCULPT_is_automasking_enabled(const Sculpt *sd, const SculptSession *ss, const Brush *br)
{
  if (br && SCULPT_stroke_is_dynamic_topology(ss, br)) {
    return false;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_TOPOLOGY)) {
    return true;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_FACE_SETS)) {
    return true;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_BOUNDARY_EDGES)) {
    return true;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_BOUNDARY_FACE_SETS)) {
    return true;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_BRUSH_NORMAL)) {
    return true;
  }
  if (SCULPT_is_automasking_mode_enabled(sd, br, BRUSH_AUTOMASKING_VIEW_NORMAL)) {
    return true;
  }
  return false;
}

static int sculpt_automasking_mode_effective_bits(const Sculpt *sculpt, const Brush *brush)
{
  if (brush) {
    return sculpt->automasking_flags | brush->automasking_flags;
  }
  return sculpt->automasking_flags;
}

bool SCULPT_automasking_needs_normal(const SculptSession *ss,
                                     const Sculpt *sculpt,
                                     const Brush *brush)
{
  int flags = sculpt_automasking_mode_effective_bits(sculpt, brush);

  return flags & (BRUSH_AUTOMASKING_BRUSH_NORMAL | BRUSH_AUTOMASKING_VIEW_NORMAL);
}

static float sculpt_automasking_normal_calc(AutomaskingCache *automasking,
                                            SculptSession *ss,
                                            PBVHVertRef vertex,
                                            const float normal[3],
                                            float limit_lower,
                                            float limit_upper,
                                            AutomaskingNodeData *automask_data)
{
  float normal_v[3];

  if (automask_data->have_orig_data) {
    copy_v3_v3(normal_v, automask_data->orig_data.no);
  }
  else {
    SCULPT_vertex_normal_get(ss, vertex, normal_v);
  }

  float angle = saacos(dot_v3v3(normal, normal_v));

  /* note that limit is pre-divided by M_PI */

  if (angle > limit_lower && angle < limit_upper) {
    float t = 1.0f - (angle - limit_lower) / (limit_upper - limit_lower);

    /* smoothstep */
    t = t * t * (3.0 - 2.0 * t);

    return t;
  }
  else if (angle > limit_upper) {
    return 0.0f;
  }

  return 1.0f;
}

static bool SCULPT_automasking_needs_factors_cache(const Sculpt *sd, const Brush *brush)
{

  const int automasking_flags = sculpt_automasking_mode_effective_bits(sd, brush);
  if (automasking_flags & BRUSH_AUTOMASKING_TOPOLOGY) {
    return true;
  }
  if (automasking_flags & (BRUSH_AUTOMASKING_BOUNDARY_EDGES,
                           BRUSH_AUTOMASKING_BOUNDARY_FACE_SETS,
                           BRUSH_AUTOMASKING_BRUSH_NORMAL,
                           BRUSH_AUTOMASKING_VIEW_NORMAL)) {
    return brush && brush->automasking_boundary_edges_propagation_steps != 1;
  }
  return false;
}

static float automasking_brush_normal_factor(AutomaskingCache *automasking,
                                             SculptSession *ss,
                                             PBVHVertRef vertex,
                                             AutomaskingNodeData *automask_data)
{
  float falloff = automasking->settings.start_normal_falloff * M_PI;
  float initial_normal[3];

  if (ss->cache) {
    copy_v3_v3(initial_normal, ss->cache->initial_normal);
  }
  else {
    copy_v3_v3(initial_normal, ss->filter_cache->initial_normal);
  }

  return sculpt_automasking_normal_calc(automasking,
                                        ss,
                                        vertex,
                                        initial_normal,
                                        automasking->settings.start_normal_limit - falloff * 0.5f,
                                        automasking->settings.start_normal_limit + falloff * 0.5f,
                                        automask_data);
}

static float automasking_view_normal_factor(AutomaskingCache *automasking,
                                            SculptSession *ss,
                                            PBVHVertRef vertex,
                                            AutomaskingNodeData *automask_data)
{
  float falloff = automasking->settings.view_normal_falloff * M_PI;

  float view_normal[3];

  if (ss->cache) {
    copy_v3_v3(view_normal, ss->cache->view_normal);
  }
  else {
    copy_v3_v3(view_normal, ss->filter_cache->view_normal);
  }

  return sculpt_automasking_normal_calc(automasking,
                                        ss,
                                        vertex,
                                        view_normal,
                                        automasking->settings.view_normal_limit,
                                        automasking->settings.view_normal_limit + falloff,
                                        automask_data);
}

static float automasking_view_occlusion_factor(AutomaskingCache *automasking,
                                               SculptSession *ss,
                                               PBVHVertRef vertex,
                                               AutomaskingNodeData *automask_data)
{
  int index = BKE_pbvh_vertex_to_index(ss->pbvh, vertex);

  if (!automasking->occluded[index]) {
    automasking->occluded[index] = SCULPT_vertex_is_occluded(ss, vertex, true) ? 2 : 1;
  }

  return automasking->occluded[index] == 2;
}

float SCULPT_automasking_factor_get(AutomaskingCache *automasking,
                                    SculptSession *ss,
                                    PBVHVertRef vert,
                                    AutomaskingNodeData *automask_data)
{
  if (!automasking) {
    return 1.0f;
  }

  int index = BKE_pbvh_vertex_to_index(ss->pbvh, vert);

  /* If the cache is initialized with valid info, use the cache. This is used when the
   * automasking information can't be computed in real time per vertex and needs to be
   * initialized for the whole mesh when the stroke starts. */
  if (automasking->factor) {
    return automasking->factor[index];
  }

  bool do_occlusion = (automasking->settings.flags &
                       (BRUSH_AUTOMASKING_VIEW_OCCLUSION | BRUSH_AUTOMASKING_VIEW_NORMAL)) ==
                      (BRUSH_AUTOMASKING_VIEW_OCCLUSION | BRUSH_AUTOMASKING_VIEW_NORMAL);
  if (do_occlusion && automasking_view_occlusion_factor(automasking, ss, vert, automask_data)) {
    return 0.0f;
  }

  if (automasking->settings.flags & BRUSH_AUTOMASKING_FACE_SETS) {
    if (!SCULPT_vertex_has_face_set(ss, vert, automasking->settings.initial_face_set)) {
      return 0.0f;
    }
  }

  if (automasking->settings.flags & BRUSH_AUTOMASKING_BOUNDARY_EDGES) {
    if (SCULPT_vertex_is_boundary(ss, vert)) {
      return 0.0f;
    }
  }

  if (automasking->settings.flags & BRUSH_AUTOMASKING_BOUNDARY_FACE_SETS) {
    if (!SCULPT_vertex_has_unique_face_set(ss, vert)) {
      return 0.0f;
    }
  }

  float mask = 1.0f;

  if ((ss->cache || ss->filter_cache) &&
      (automasking->settings.flags & BRUSH_AUTOMASKING_BRUSH_NORMAL)) {
    mask *= automasking_brush_normal_factor(automasking, ss, vert, automask_data);
  }

  if ((ss->cache || ss->filter_cache) &&
      (automasking->settings.flags & BRUSH_AUTOMASKING_VIEW_NORMAL)) {
    mask *= automasking_view_normal_factor(automasking, ss, vert, automask_data);
  }

  return mask;
}

void SCULPT_automasking_cache_free(AutomaskingCache *automasking)
{
  if (!automasking) {
    return;
  }

  MEM_SAFE_FREE(automasking->factor);
  MEM_SAFE_FREE(automasking->occluded);
  MEM_SAFE_FREE(automasking);
}

static bool sculpt_automasking_is_constrained_by_radius(Brush *br)
{
  /* 2D falloff is not constrained by radius. */
  if (br->falloff_shape == PAINT_FALLOFF_SHAPE_TUBE) {
    return false;
  }

  if (ELEM(br->sculpt_tool, SCULPT_TOOL_GRAB, SCULPT_TOOL_THUMB, SCULPT_TOOL_ROTATE)) {
    return true;
  }
  return false;
}

struct AutomaskFloodFillData {
  float *automask_factor;
  float radius;
  bool use_radius;
  float location[3];
  char symm;
};

static bool automask_floodfill_cb(SculptSession *ss,
                                  PBVHVertRef from_v,
                                  PBVHVertRef to_v,
                                  bool UNUSED(is_duplicate),
                                  void *userdata)
{
  AutomaskFloodFillData *data = (AutomaskFloodFillData *)userdata;
  int from_v_i = BKE_pbvh_vertex_to_index(ss->pbvh, from_v);
  int to_v_i = BKE_pbvh_vertex_to_index(ss->pbvh, to_v);

  data->automask_factor[to_v_i] = 1.0f;
  data->automask_factor[from_v_i] = 1.0f;
  return (!data->use_radius ||
          SCULPT_is_vertex_inside_brush_radius_symm(
              SCULPT_vertex_co_get(ss, to_v), data->location, data->radius, data->symm));
}

static float *SCULPT_topology_automasking_init(Sculpt *sd, Object *ob, float *automask_factor)
{
  SculptSession *ss = ob->sculpt;
  Brush *brush = BKE_paint_brush(&sd->paint);

  if (BKE_pbvh_type(ss->pbvh) == PBVH_FACES && !ss->pmap) {
    BLI_assert_msg(0, "Topology masking: pmap missing");
    return nullptr;
  }

  const int totvert = SCULPT_vertex_count_get(ss);
  for (int i : IndexRange(totvert)) {
    automask_factor[i] = 0.0f;
  }

  /* Flood fill automask to connected vertices. Limited to vertices inside
   * the brush radius if the tool requires it. */
  SculptFloodFill flood;
  SCULPT_floodfill_init(ss, &flood);
  const float radius = ss->cache ? ss->cache->radius : FLT_MAX;
  SCULPT_floodfill_add_active(sd, ob, ss, &flood, radius);

  AutomaskFloodFillData fdata = {nullptr};

  fdata.automask_factor = automask_factor;
  fdata.radius = radius;
  fdata.use_radius = ss->cache && sculpt_automasking_is_constrained_by_radius(brush);
  fdata.symm = SCULPT_mesh_symmetry_xyz_get(ob);

  copy_v3_v3(fdata.location, SCULPT_active_vertex_co_get(ss));
  SCULPT_floodfill_execute(ss, &flood, automask_floodfill_cb, &fdata);
  SCULPT_floodfill_free(&flood);

  return automask_factor;
}

static float *sculpt_face_sets_automasking_init(Sculpt *sd, Object *ob, float *automask_factor)
{
  SculptSession *ss = ob->sculpt;
  Brush *brush = BKE_paint_brush(&sd->paint);

  if (!SCULPT_is_automasking_enabled(sd, ss, brush)) {
    return nullptr;
  }

  if (BKE_pbvh_type(ss->pbvh) == PBVH_FACES && !ss->pmap) {
    BLI_assert_msg(0, "Face Sets automasking: pmap missing");
    return nullptr;
  }

  int tot_vert = SCULPT_vertex_count_get(ss);
  int active_face_set = SCULPT_active_face_set_get(ss);
  for (int i : IndexRange(tot_vert)) {
    PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);

    if (!SCULPT_vertex_has_face_set(ss, vertex, active_face_set)) {
      automask_factor[i] *= 0.0f;
    }
  }

  return automask_factor;
}

#define EDGE_DISTANCE_INF -1

float *SCULPT_boundary_automasking_init(Object *ob,
                                        eBoundaryAutomaskMode mode,
                                        int propagation_steps,
                                        float *automask_factor)
{
  SculptSession *ss = ob->sculpt;

  if (!ss->pmap) {
    BLI_assert_msg(0, "Boundary Edges masking: pmap missing");
    return nullptr;
  }

  const int totvert = SCULPT_vertex_count_get(ss);
  int *edge_distance = (int *)MEM_callocN(sizeof(int) * totvert, "automask_factor");

  for (int i : IndexRange(totvert)) {
    PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);

    edge_distance[i] = EDGE_DISTANCE_INF;
    switch (mode) {
      case AUTOMASK_INIT_BOUNDARY_EDGES:
        if (SCULPT_vertex_is_boundary(ss, vertex)) {
          edge_distance[i] = 0;
        }
        break;
      case AUTOMASK_INIT_BOUNDARY_FACE_SETS:
        if (!SCULPT_vertex_has_unique_face_set(ss, vertex)) {
          edge_distance[i] = 0;
        }
        break;
    }
  }

  for (int propagation_it : IndexRange(propagation_steps)) {
    for (int i : IndexRange(totvert)) {
      PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);

      if (edge_distance[i] != EDGE_DISTANCE_INF) {
        continue;
      }
      SculptVertexNeighborIter ni;
      SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
        if (edge_distance[ni.index] == propagation_it) {
          edge_distance[i] = propagation_it + 1;
        }
      }
      SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);
    }
  }

  for (int i : IndexRange(totvert)) {
    if (edge_distance[i] == EDGE_DISTANCE_INF) {
      continue;
    }
    const float p = 1.0f - ((float)edge_distance[i] / (float)propagation_steps);
    const float edge_boundary_automask = pow2f(p);
    automask_factor[i] *= (1.0f - edge_boundary_automask);
  }

  MEM_SAFE_FREE(edge_distance);
  return automask_factor;
}

static void SCULPT_automasking_cache_settings_update(AutomaskingCache *automasking,
                                                     SculptSession *ss,
                                                     Sculpt *sd,
                                                     Brush *brush)
{
  automasking->settings.flags = sculpt_automasking_mode_effective_bits(sd, brush);
  automasking->settings.initial_face_set = SCULPT_active_face_set_get(ss);

  automasking->settings.view_normal_limit = sd->automasking_view_normal_limit;
  automasking->settings.view_normal_falloff = sd->automasking_view_normal_falloff;
  automasking->settings.start_normal_limit = sd->automasking_start_normal_limit;
  automasking->settings.start_normal_falloff = sd->automasking_start_normal_falloff;
}

void sculpt_normal_occlusion_automasking_fill(AutomaskingCache *automasking,
                                              Object *ob,
                                              float *factor,
                                              eAutomasking_flag mode)
{
  SculptSession *ss = ob->sculpt;
  const int totvert = SCULPT_vertex_count_get(ss);

  /* No need to build original data since this is only called at the beginning of strokes.*/
  AutomaskingNodeData nodedata;
  nodedata.have_orig_data = false;

  for (int i = 0; i < totvert; i++) {
    PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);

    if ((int)mode & BRUSH_AUTOMASKING_BRUSH_NORMAL) {
      factor[i] *= automasking_brush_normal_factor(automasking, ss, vertex, &nodedata);
    }
    if ((int)mode & BRUSH_AUTOMASKING_VIEW_NORMAL) {
      if ((int)mode & BRUSH_AUTOMASKING_VIEW_OCCLUSION) {
        factor[i] *= automasking_view_occlusion_factor(automasking, ss, vertex, &nodedata);
      }

      factor[i] *= automasking_view_normal_factor(automasking, ss, vertex, &nodedata);
    }
  }
}

AutomaskingCache *SCULPT_automasking_cache_init(Sculpt *sd, Brush *brush, Object *ob)
{
  SculptSession *ss = ob->sculpt;
  const int totvert = SCULPT_vertex_count_get(ss);

  if (!SCULPT_is_automasking_enabled(sd, ss, brush)) {
    return nullptr;
  }

  AutomaskingCache *automasking = (AutomaskingCache *)MEM_callocN(sizeof(AutomaskingCache),
                                                                  "automasking cache");
  SCULPT_automasking_cache_settings_update(automasking, ss, sd, brush);
  SCULPT_boundary_info_ensure(ob);

  if (sculpt_automasking_mode_effective_bits(sd, brush) & BRUSH_AUTOMASKING_VIEW_OCCLUSION) {
    automasking->occluded = (char *)MEM_callocN(totvert, "automasking->occluded");
  }

  if (!SCULPT_automasking_needs_factors_cache(sd, brush)) {
    return automasking;
  }

  automasking->factor = (float *)MEM_malloc_arrayN(totvert, sizeof(float), "automask_factor");
  for (int i : IndexRange(totvert)) {
    automasking->factor[i] = 1.0f;
  }

  const int boundary_propagation_steps = brush ?
                                             brush->automasking_boundary_edges_propagation_steps :
                                             1;

  if (SCULPT_is_automasking_mode_enabled(sd, brush, BRUSH_AUTOMASKING_TOPOLOGY)) {
    SCULPT_vertex_random_access_ensure(ss);
    SCULPT_topology_automasking_init(sd, ob, automasking->factor);
  }
  if (SCULPT_is_automasking_mode_enabled(sd, brush, BRUSH_AUTOMASKING_FACE_SETS)) {
    SCULPT_vertex_random_access_ensure(ss);
    sculpt_face_sets_automasking_init(sd, ob, automasking->factor);
  }

  int normal_bits = sculpt_automasking_mode_effective_bits(sd, brush) &
                    (BRUSH_AUTOMASKING_BRUSH_NORMAL | BRUSH_AUTOMASKING_VIEW_NORMAL |
                     BRUSH_AUTOMASKING_VIEW_OCCLUSION);

  if (normal_bits) {
    sculpt_normal_occlusion_automasking_fill(
        automasking, ob, automasking->factor, (eAutomasking_flag)normal_bits);
  }

  if (SCULPT_is_automasking_mode_enabled(sd, brush, BRUSH_AUTOMASKING_BOUNDARY_EDGES)) {
    SCULPT_vertex_random_access_ensure(ss);
    SCULPT_boundary_automasking_init(
        ob, AUTOMASK_INIT_BOUNDARY_EDGES, boundary_propagation_steps, automasking->factor);
  }
  if (SCULPT_is_automasking_mode_enabled(sd, brush, BRUSH_AUTOMASKING_BOUNDARY_FACE_SETS)) {
    SCULPT_vertex_random_access_ensure(ss);
    SCULPT_boundary_automasking_init(
        ob, AUTOMASK_INIT_BOUNDARY_FACE_SETS, boundary_propagation_steps, automasking->factor);
  }

  return automasking;
}

bool SCULPT_automasking_needs_origco(const SculptSession *ss, const Sculpt *sd, const Brush *br)
{
  if (br &&
      br->automasking_flags & (BRUSH_AUTOMASKING_BRUSH_NORMAL | BRUSH_AUTOMASKING_VIEW_NORMAL)) {
    return true;
  }

  if (sd &&
      sd->automasking_flags & (BRUSH_AUTOMASKING_BRUSH_NORMAL | BRUSH_AUTOMASKING_VIEW_NORMAL)) {
    return true;
  }

  return false;
}
