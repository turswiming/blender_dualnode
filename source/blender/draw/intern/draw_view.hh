/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 */

#include "DRW_gpu_wrapper.hh"

#include "draw_shader_shared.h"

namespace blender::draw {

class View {
 public:
  UniformBuffer<ViewInfos> data_;

 private:
  /** Result of the visibility computation. 1 bit per resource ID. */
  StorageArrayBuffer<uint4, 1, true> object_visibility_buf;

  const char *debug_name;

  View *parent = nullptr;
  bool do_visibility = true;

 public:
  View(const char *name) : debug_name(name), object_visibility_buf(name){};

  void sync();

 private:
  void compute_visibility(Manager::ObjectBoundsBuf &bounds, uint resource_len)
  {
    object_visibility_buf.resize(divide_ceil_u(resource_len, 128));

    if (do_visibility == false) {
      object_visibility_buf.clear(0xFFFFFFFFu);
      return;
    }

    uint thread_groups = divide_ceil_u(resource_len, DRW_VISIBILITY_GROUP_SIZE);
    GPUShader *shader = draw_shader_visibility_get();
    GPU_shader_bind(shader);
    GPU_storagebuf_bind(bounds, 0);
    GPU_compute_dispatch(shader, thread_groups, 1, 1);
  }
};

}  // namespace blender::draw
