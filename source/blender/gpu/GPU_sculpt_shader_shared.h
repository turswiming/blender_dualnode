/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#ifndef USE_GPU_SHADER_CREATE_INFO
#  include "GPU_shader_shared_utils.h"
#endif

struct TrianglePaintInput {
  int4 vert_indices;
  /**
   * Delta barycentric coordinates between 2 neighboring UV's in the U direction.
   *
   * Only the first two coordinates are stored. The third should be recalculated on the fly.
   */
  float2 delta_barycentric_coord_u;
  float2 _pad;
};
BLI_STATIC_ASSERT_ALIGN(TrianglePaintInput, 16)

struct PackedPixelRow {
  /** Barycentric coordinate of the first pixel. */
  float2 start_barycentric_coord;
  /** Image coordinate starting of the first pixel. First 16 bits is u, other 16 bits is v. */
  uint start_image_coordinate;

  /**
   * 16 bits: Number of sequential pixels encoded in this package.
   * 16 bits: Reference to the pbvh triangle index.
   */
  uint encoded;
};

#define PIXEL_ROW_START_IMAGE_COORD(row) \
  ivec2((row.start_image_coordinate & 0xffff0000) >> 16, row.start_image_coordinate & 0xffff)
#define PIXEL_ROW_LEN(row) uint((row.encoded & 0xffff0000) >> 16)
#define PIXEL_ROW_PRIM_INDEX(row) uint(row.encoded & 0xffff);

BLI_STATIC_ASSERT_ALIGN(TrianglePaintInput, 16)
