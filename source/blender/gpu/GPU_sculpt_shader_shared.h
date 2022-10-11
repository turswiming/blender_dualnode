/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#ifndef USE_GPU_SHADER_CREATE_INFO
#  include "GPU_shader_shared_utils.h"
#endif

struct TrianglePaintInput {
  int3 vert_indices;
  float _pad0;
  /**
   * Delta barycentric coordinates between 2 neighboring UV's in the U direction.
   *
   * Only the first two coordinates are stored. The third should be recalculated on the fly.
   */
  float2 delta_barycentric_coord;
  float2 _pad1;
};
BLI_STATIC_ASSERT_ALIGN(TrianglePaintInput, 16)

struct PackedPixelRow {
  /** Barycentric coordinate of the first pixel. */
  float2 start_barycentric_coord;
  /** Image coordinate starting of the first pixel. First 16 bits is u, other 16 bits is v. */
  uint start_image_coordinate;

  /**
   * 16 bits: Reference to the pbvh triangle index.
   * 16 bits: Number of sequential pixels encoded in this package.
   */
  uint encoded;
};
BLI_STATIC_ASSERT_ALIGN(TrianglePaintInput, 16)

#define PIXEL_ROW_START_IMAGE_COORD(row) \
  ivec2(row.start_image_coordinate & 0xffff, (row.start_image_coordinate & 0xffff0000) >> 16)
#define PIXEL_ROW_LEN(row) uint(row.encoded & 0xffff);
#define PIXEL_ROW_PRIM_INDEX(row) uint((row.encoded & 0xffff0000) >> 16)

struct PaintBrushTestData {
  /* world to local matrix for clipping plane tests. */
  float4x4 symm_rot_mat_inv;
};
BLI_STATIC_ASSERT_ALIGN(PaintBrushTestData, 16)

struct PaintBrushData {
  float4 color;
  PaintBrushTestData test;
  float alpha;
  int falloff_shape;
  float _pad0[2];
};
BLI_STATIC_ASSERT_ALIGN(PaintBrushData, 16)

struct PaintStepData {
  float3 location;
  float radius;
  /* Circle falloff. */
  float4 plane_view;
  float hardness;
  float strength;
  int mirror_symmetry_pass;
  int _pad0[1];
};
BLI_STATIC_ASSERT_ALIGN(PaintStepData, 16);

struct PaintTileData {
  int tile_number;
  int layer_id;
  int2 sub_tile_id;
};
BLI_STATIC_ASSERT_ALIGN(PaintTileData, 16);
