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
 * A Command::MultiDraw allow to split the command stream into batch-able chunks of commands with
 * the same render state.
 */
struct MultiDraw {
  /** Index of next MultiDraw from the same Command::Header. */
  uint next;

  /** Index of the first command after sorting. */
  uint command_start;
#if defined(GPU_SHADER) && !defined(GPU_METAL)
  /* No support for ushort. */
  uint cmd_len_packed;
#  define _cmd_len (cmd_len_packed & 0xFFFFu)
#  define _inverted_len (cmd_len_packed >> 16u)
#else
  /**
   * NOTE(fclem): We have to make room to be able to stick the GPUBatch pointer at the end.
   */
  /** Number of commands. Needed to issue the draw call. */
  ushort command_len;
  /** Number of non inverted scaling commands in this Group. */
  ushort front_facing_len;
#endif
  /** GPUBatch values to be copied to DrawCommand after sorting (if not overriden). */
  uint vertex_count;
  uint instance_count;
  uint base_vertex; /** NOTE: (uint)-1 if non indexed draw. */
#ifdef GPU_SHADER
  uint _pad0, _pad1;
#else
  /* NOTE: Union just to make sure the struct has always the same size on all platform. */
  union {
    /** Needed to create the correct draw call. */
    GPUBatch *gpu_batch;
    uint _pad0[2];
  };

  void execute(RecordingState &state, Span<MultiDraw> multi_draw_buf, uint command_id) const;
#endif
};
BLI_STATIC_ASSERT(sizeof(MultiDraw) == 32, "MultiDraw might not have the same size on GPU and CPU")

/**
 * Representation of a future draw call inside a MultiDraw. This #DrawDescription is then converted
 * into #DrawCommand on GPU after visibility and compaction. Multiple #DrawDescription might get
 * merged into the same final #DrawCommand.
 */
struct DrawDescription {
  /* Override of GPUBatch values. (uint)-1 otherwise. */
  uint vertex_first;
  uint vertex_count;
  uint instance_first;
  uint instance_count;
  /* Resource ID associated with this call. */
  uint resource_id;
  /* Reference to parent MultiDraw to get the GPUBatch vertex / instance count. */
  uint multi_draw_id;
};

/** \} */

#ifndef GPU_SHADER
};  // namespace blender::draw::command
#endif
