/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#ifndef GPU_SHADER
#  include "BLI_span.hh"
#  include "GPU_shader_shared_utils.h"

namespace blender::draw::command {

struct RecordingState;

#endif

/* -------------------------------------------------------------------- */
/** \name Multi Draw
 * \{ */

/**
 * A DrawGroup allow to split the command stream into batch-able chunks of commands with
 * the same render state.
 */
struct DrawGroup {
  /** Index of next DrawGroup from the same header. */
  uint next;

  /** Index of the first command after sorting. */
  uint command_start;
  /** Total number of commands (including inverted facing). Needed to issue the draw call. */
  uint command_len;
  /** Number of non inverted scaling commands in this Group. */
  uint front_facing_len;

  /** GPUBatch values to be copied to DrawCommand after sorting (if not overriden). */
  int vertex_len; /** NOTE: Negative if using indexed draw. */
  uint _pad0;

#ifdef GPU_SHADER
  uint _pad1 _pad2;
#else
  /* NOTE: Union just to make sure the struct has always the same size on all platform. */
  union {
    /** Needed to create the correct draw call. */
    GPUBatch *gpu_batch;
    uint _pad1[2];
  };
#endif
};
BLI_STATIC_ASSERT(sizeof(DrawGroup) == 32, "DrawGroup might not have the same size on GPU and CPU")

/**
 * Representation of a future draw call inside a DrawGroup. This #DrawPrototype is then
 * converted into #DrawCommand on GPU after visibility and compaction. Multiple
 * #DrawPrototype might get merged into the same final #DrawCommand.
 */
struct DrawPrototype {
  /* Reference to parent DrawGroup to get the GPUBatch vertex / instance count. */
  uint group_id;
  /* Resource handle associated with this call. Also reference visibility. */
  uint resource_handle;
  /* Override of GPUBatch values. (uint)-1 otherwise. */
  uint vertex_len;
  uint instance_len;
};

/** \} */

#ifndef GPU_SHADER
};  // namespace blender::draw::command
#endif
