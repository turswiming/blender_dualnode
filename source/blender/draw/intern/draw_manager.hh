/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 */

#include "BLI_sys_types.h"

#include "draw_pass.hh"

#include <string>

namespace blender::draw {

struct ObjectRef {
  Object *object;
  /** Dupli object that corresponds to the current object. */
  DupliObject *dupli_object;
  /** Object that created the dupli-list the current object is part of. */
  Object *dupli_parent;
};

class Manager {
  using ObjectMatricesBuf = StorageArrayBuffer<ObjectMatrices, 128>;
  using ObjectBoundsBuf = StorageArrayBuffer<ObjectBounds, 128>;
  using ObjectInfosBuf = StorageArrayBuffer<ObjectInfos, 128>;

 public:
  /**
   * Buffers containing all object data. Referenced by resource index.
   * Exposed as public members for shader access after sync.
   */
  ObjectMatricesBuf matrix_buf;
  ObjectBoundsBuf bounds_buf;
  ObjectInfosBuf infos_buf;

 private:
  uint resource_len = 0;
  Object *object = nullptr;

  Object *object_active = nullptr;

 public:
  /**
   * Create a new resource handle for the given object. Can be called multiple time with the same
   * object **successively** without duplicating the data.
   */
  ResourceHandle resource_handle(const ObjectRef ref);
  /**
   * Get resource id for a loose matrix. The draw-calls for this resource handle won't be culled
   * and there won't be any associated object info / bounds. Assumes correct handedness / winding.
   */
  ResourceHandle resource_handle(const float4x4 &model_matrix);

  /**
   * Populate additional per resource data on demand.
   */
  void object_attributes(ResourceHandle handle, Object &object, Span<GPUMaterial *> materials);

  /**
   * Submit a pass for drawing. All resource reference will be dereferenced and commands will be
   * sent to GPU.
   */
  void submit(const PassSimple &pass);
  void submit(const PassMain &pass);

 private:
  /**
   * Reset all buffers to be refilled.
   */
  void begin_sync();

  /**
   * Finalize the object data on GPU.
   */
  void end_sync();
};

inline ResourceHandle Manager::resource_handle(const ObjectRef ref)
{
  bool is_active_object = (ref.dupli_object ? ref.dupli_parent : ref.object) == object_active;
  matrix_buf.get_or_resize(resource_len).init(*ref.object);
  bounds_buf.get_or_resize(resource_len).init(*ref.object);
  infos_buf.get_or_resize(resource_len).init(ref, is_active_object);
  return ResourceHandle(resource_len++, (ref.object->transflag & OB_NEG_SCALE) != 0);
}

inline ResourceHandle Manager::resource_handle(const float4x4 &model_matrix)
{
  matrix_buf.get_or_resize(resource_len).init(model_matrix);
  return ResourceHandle(resource_len++, false);
}

inline void Manager::object_attributes(ResourceHandle handle,
                                       Object &object,
                                       Span<GPUMaterial *> materials)
{
  /* TODO */
  (void)handle;
  (void)object;
  (void)materials;
}

}  // namespace blender::draw

/* TODO(@fclem): This is for testing. The manager should be passed to the engine through the
 * callbacks. */
blender::draw::Manager *DRW_manager_get();
blender::draw::ObjectRef DRW_object_ref_get(Object *object);
