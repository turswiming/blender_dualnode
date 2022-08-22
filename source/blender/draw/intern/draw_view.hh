/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 */

#include "DRW_gpu_wrapper.hh"

#include "draw_shader_shared.h"

namespace blender::draw {

class Manager;

/* TODO deduplicate. */
using ObjectBoundsBuf = StorageArrayBuffer<ObjectBounds, 128>;

class View {
  friend Manager;

 private:
  UniformBuffer<ViewInfos> data_;
  /** Result of the visibility computation. 1 bit per resource ID. */
  StorageArrayBuffer<uint4, 1, true> object_visibility_buf;

  const char *debug_name;

  bool do_visibility_ = true;
  bool dirty_ = true;

 public:
  View(const char *name) : object_visibility_buf(name), debug_name(name){};

  void set_clip_planes(Span<float4> planes);

  void sync(const float4x4 &view_mat, const float4x4 &win_mat);

 private:
  /** Called from draw manager. */
  void bind();
  void compute_visibility(ObjectBoundsBuf &bounds, uint resource_len);

  void update_view_vectors();
  void update_viewport_size();

  void frustum_boundbox_calc(BoundBox &bbox);
  void frustum_culling_planes_calc();
  void frustum_culling_sphere_calc(const BoundBox &bbox, BoundSphere &bsphere);
};

}  // namespace blender::draw
