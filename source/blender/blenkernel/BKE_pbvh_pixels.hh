/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#pragma once

#include "BLI_math.h"
#include "BLI_math_vec_types.hh"
#include "BLI_rect.h"
#include "BLI_vector.hh"

#include "DNA_image_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_image.h"
#include "BKE_image_wrappers.hh"

#include "IMB_imbuf_types.h"

#include "GPU_sculpt_shader_shared.h"
#include "GPU_storage_buffer.h"

namespace blender::bke::pbvh::pixels {

/**
 * Data shared between pixels that belong to the same triangle.
 *
 * Data is stored as a list of structs, grouped by usage to improve performance (improves CPU
 * cache prefetching).
 */
struct Triangles {
  /** Data accessed by the inner loop of the painting brush. */
  Vector<TrianglePaintInput> paint_input;
  GPUStorageBuf *gpu_buffer = nullptr;

 public:
  void append(const int3 vert_indices)
  {
    TrianglePaintInput triangle;
    triangle.vert_indices = int4(vert_indices.x, vert_indices.y, vert_indices.z, 0.0f);
    triangle.delta_barycentric_coord_u = float2(0.0f);
    this->paint_input.append(triangle);
  }

  TrianglePaintInput &get_paint_input(const int index)
  {
    return paint_input[index];
  }

  const TrianglePaintInput &get_paint_input(const int index) const
  {
    return paint_input[index];
  }

  ~Triangles()
  {
    clear();
  }

  /** Clear data associated with self. */
  void clear();
  void ensure_gpu_buffer();

  uint64_t size() const
  {
    return paint_input.size();
  }

  uint64_t mem_size() const
  {
    return paint_input.size() * sizeof(TrianglePaintInput);
  }
};

/**
 * Encode sequential pixels to reduce memory footprint.
 */
struct PackedPixelRow {
  /** Barycentric coordinate of the first pixel. */
  float2 start_barycentric_coord;
  /** Image coordinate starting of the first pixel. */
  ushort2 start_image_coordinate;
  /** Number of sequential pixels encoded in this package. */
  ushort num_pixels;
  /** Reference to the pbvh triangle index. */
  ushort triangle_index;
};

/**
 * Node pixel data containing the pixels for a single UDIM tile.
 */
struct UDIMTilePixels {
  /** UDIM Tile number. */
  short tile_number;

  struct {
    bool dirty : 1;
  } flags;

  /* Dirty region of the tile in image space. */
  rcti dirty_region;

  Vector<PackedPixelRow> pixel_rows;
  int64_t gpu_buffer_offset;

  UDIMTilePixels()
  {
    flags.dirty = false;
    BLI_rcti_init_minmax(&dirty_region);
  }

  void mark_dirty(const PackedPixelRow &pixel_row)
  {
    int2 start_image_coord(pixel_row.start_image_coordinate.x, pixel_row.start_image_coordinate.y);
    BLI_rcti_do_minmax_v(&dirty_region, start_image_coord);
    BLI_rcti_do_minmax_v(&dirty_region, start_image_coord + int2(pixel_row.num_pixels + 1, 0));
    flags.dirty = true;
  }

  void clear_dirty()
  {
    BLI_rcti_init_minmax(&dirty_region);
    flags.dirty = false;
  }
};

struct UDIMTileUndo {
  short tile_number;
  rcti region;

  UDIMTileUndo(short tile_number, rcti &region) : tile_number(tile_number), region(region)
  {
  }
};

struct NodeData {
  struct {
    bool dirty : 1;
  } flags;

  Vector<UDIMTilePixels> tiles;
  Vector<UDIMTileUndo> undo_regions;
  Triangles triangles;

  struct {
    /** Contains GPU buffer for all associated pixels. Tiles have a range inside this buffer
     * (#UDIMTilePixels.start_index, #UDIMTilePixels.end_index). */
    GPUStorageBuf *pixels = nullptr;
  } gpu_buffers;

  NodeData()
  {
    flags.dirty = false;
  }

  ~NodeData()
  {
    clear_data();
  }

  UDIMTilePixels *find_tile_data(const image::ImageTileWrapper &image_tile)
  {
    for (UDIMTilePixels &tile : tiles) {
      if (tile.tile_number == image_tile.get_tile_number()) {
        return &tile;
      }
    }
    return nullptr;
  }

  void rebuild_undo_regions()
  {
    undo_regions.clear();
    for (UDIMTilePixels &tile : tiles) {
      rcti region;
      BLI_rcti_init_minmax(&region);
      for (PackedPixelRow &pixel_row : tile.pixel_rows) {
        BLI_rcti_do_minmax_v(
            &region, int2(pixel_row.start_image_coordinate.x, pixel_row.start_image_coordinate.y));
        BLI_rcti_do_minmax_v(&region,
                             int2(pixel_row.start_image_coordinate.x + pixel_row.num_pixels + 1,
                                  pixel_row.start_image_coordinate.y + 1));
      }
      undo_regions.append(UDIMTileUndo(tile.tile_number, region));
    }
  }

  void mark_region(Image &image, const image::ImageTileWrapper &image_tile, ImBuf &image_buffer)
  {
    UDIMTilePixels *tile = find_tile_data(image_tile);
    if (tile && tile->flags.dirty) {
      if (image_buffer.planes == 8) {
        image_buffer.planes = 32;
        BKE_image_partial_update_mark_full_update(&image);
      }
      else {
        BKE_image_partial_update_mark_region(
            &image, image_tile.image_tile, &image_buffer, &tile->dirty_region);
      }
      tile->clear_dirty();
    }
  }

  void clear_data()
  {
    tiles.clear();
    triangles.clear();
    if (gpu_buffers.pixels) {
      GPU_storagebuf_free(gpu_buffers.pixels);
      gpu_buffers.pixels = nullptr;
    }
  }

  void ensure_gpu_buffers()
  {
    triangles.ensure_gpu_buffer();
    if (gpu_buffers.pixels == nullptr) {
      build_pixels_gpu_buffer();
    }
  }

  static void free_func(void *instance)
  {
    NodeData *node_data = static_cast<NodeData *>(instance);
    MEM_delete(node_data);
  }

 private:
  void build_pixels_gpu_buffer();
};

NodeData &BKE_pbvh_pixels_node_data_get(PBVHNode &node);
void BKE_pbvh_pixels_mark_image_dirty(PBVHNode &node, Image &image, ImageUser &image_user);

}  // namespace blender::bke::pbvh::pixels
