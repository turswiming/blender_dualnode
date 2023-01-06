/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#include "BLI_array.hh"
#include "BLI_bit_vector.hh"
#include "BLI_math.h"
#include "BLI_math_vec_types.hh"
#include "BLI_vector.hh"

#include "IMB_imbuf.h"
#include "IMB_imbuf_types.h"

#include "BKE_image_wrappers.hh"
#include "BKE_pbvh.h"
#include "BKE_pbvh_pixels.hh"

#include "pbvh_intern.h"
#include "pbvh_pixels_copy.hh"
#include "pbvh_uv_islands.hh"

namespace blender::bke::pbvh::pixels {

enum class CoordSpace {
  UV,
  Tile,
};

template<CoordSpace Space> struct Vertex {
  float2 co;
};

template<CoordSpace Space> struct Edge {
  Vertex<Space> v1;
  Vertex<Space> v2;
};

class NonManifoldTileEdges : public Vector<Edge<CoordSpace::Tile>> {};

class NonManifoldUVEdges : public Vector<Edge<CoordSpace::UV>> {
 public:
  NonManifoldUVEdges(const uv_islands::MeshData &mesh_data)
  {
    reserve(count_non_manifold_edges(mesh_data));

    for (const uv_islands::MeshPrimitive &mesh_primitive : mesh_data.primitives) {
      for (int i = 0; i < 3; i++) {
        const uv_islands::MeshEdge &mesh_edge = *mesh_primitive.edges[i];
        if (is_manifold(mesh_edge)) {
          continue;
        }
        Edge<CoordSpace::UV> edge;
        edge.v1.co = find_uv_vert(mesh_primitive, mesh_edge.vert1).uv;
        edge.v2.co = find_uv_vert(mesh_primitive, mesh_edge.vert1).uv;
        append(edge);
      }
    }
  }

  NonManifoldTileEdges extract_tile_edges(const image::ImageTileWrapper image_tile,
                                          const int2 tile_resolution)
  {
    NonManifoldTileEdges result;
    // TODO add edges that intersects with the given tile.
    // Convert the space from uv to tile.
    return result;
  }

 private:
  static int64_t count_non_manifold_edges(const uv_islands::MeshData &mesh_data)
  {
    int64_t result = 0;
    for (const uv_islands::MeshPrimitive &mesh_primitive : mesh_data.primitives) {
      for (int i = 0; i < 3; i++) {
        const uv_islands::MeshEdge &mesh_edge = *mesh_primitive.edges[i];
        if (is_manifold(mesh_edge)) {
          continue;
        }
        result += 1;
      }
    }
    return result;
  }

  static const uv_islands::MeshUVVert &find_uv_vert(
      const uv_islands::MeshPrimitive &mesh_primitive, const uv_islands::MeshVertex *mesh_vertex)
  {
    for (const uv_islands::MeshUVVert &uv_vertex : mesh_primitive.vertices) {
      if (uv_vertex.vertex == mesh_vertex) {
        return uv_vertex;
      }
    }
    // TODO: Use cleaner interface.
    BLI_assert_unreachable();
    static uv_islands::MeshUVVert dummy;
    return dummy;
  }
  static bool is_manifold(const uv_islands::MeshEdge mesh_edge)
  {
    return mesh_edge.primitives.size() == 2;
  }
};

class PixelNodesTileData : public Vector<std::reference_wrapper<UDIMTilePixels>> {
 public:
  PixelNodesTileData(PBVH &pbvh, const image::ImageTileWrapper &image_tile)
  {
    reserve(count_nodes(pbvh, image_tile));

    for (PBVHNode &node : MutableSpan(pbvh.nodes, pbvh.totnode)) {
      if (should_add_node(node, image_tile)) {
        NodeData &node_data = *static_cast<NodeData *>(node.pixels.node_data);
        UDIMTilePixels &tile_pixels = *node_data.find_tile_data(image_tile);
        append(tile_pixels);
      }
    }
  }

 private:
  static bool should_add_node(PBVHNode &node, const image::ImageTileWrapper &image_tile)
  {
    if ((node.flag & PBVH_Leaf) == 0) {
      return false;
    }
    if (node.pixels.node_data == nullptr) {
      return false;
    }
    NodeData &node_data = *static_cast<NodeData *>(node.pixels.node_data);
    if (node_data.find_tile_data(image_tile) == nullptr) {
      return false;
    }
    return true;
  }

  static int64_t count_nodes(PBVH &pbvh, const image::ImageTileWrapper &image_tile)
  {
    int64_t result = 0;
    for (PBVHNode &node : MutableSpan(pbvh.nodes, pbvh.totnode)) {
      if (should_add_node(node, image_tile)) {
        result++;
      }
    }
    return result;
  }
};

/**
 * Row contains intermediate data per pixel for a single image row. It is used during updating to
 * encode pixels.
 */

struct Row {
  enum class PixelType {
    Undecided,
    /** This pixel is directly affected by a brush and doesn't need to be solved. */
    Brush,
    /** This pixel will be copid from another pixel to solve non-manifold edge bleeding. */
    CopyFromClosestEdge,
  };

  struct Elem {
    PixelType type = PixelType::Undecided;
    /**
     * Distance to the closest edge that can be sourced to fix an edge bleed.
     * A distance of 0.0 means that the pixel is being drawn on directly and
     * doesn't need to be checked.
     */
    float distance = 0.0f;
    PixelCopyCommand copy_command;

    Elem() = default;

    Elem(int2 co)
    {
      copy_command.destination = co;
      copy_command.source_1 = co;
      copy_command.source_2 = co;
      copy_command.mix_factor = 0.0f;
    }
  };

  int row_number = 0;
  Array<Elem> pixels;

  Row(int64_t width) : pixels(width)
  {
  }

  void reinit(int y)
  {
    row_number = y;
    for (int x = 0; x < pixels.size(); x++) {
      pixels[x] = Elem(int2(x, y));
    }
  }

  void mask_brush_pixels(const PixelNodesTileData &nodes_tile_pixels)
  {
    for (const UDIMTilePixels &tile_pixels : nodes_tile_pixels) {
      for (const PackedPixelRow &encoded_pixels : tile_pixels.pixel_rows) {
        if (encoded_pixels.start_image_coordinate.y != row_number) {
          continue;
        }
        for (int x = encoded_pixels.start_image_coordinate.x;
             x < encoded_pixels.start_image_coordinate.x + encoded_pixels.num_pixels;
             x++) {
          pixels[x].type = PixelType::Brush;
        }
      }
    }
  }

  void print_debug() const
  {
    for (const Elem &pixel : pixels) {
      printf("%d", pixel.type);
    }
    printf("\n");
  }
};

static void copy_pixels_reinit(PixelCopyTiles &tiles)
{
  tiles.clear();
}

void BKE_pbvh_pixels_copy_update(PBVH &pbvh,
                                 Image &image,
                                 ImageUser &image_user,
                                 const uv_islands::MeshData &mesh_data)
{
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  copy_pixels_reinit(pbvh_data.tiles_copy_pixels);
  NonManifoldUVEdges non_manifold_edges(mesh_data);
  if (non_manifold_edges.is_empty()) {
    printf("Early exit: No non manifold edges detected\n");
    return;
  }

  ImageUser tile_user = image_user;
  LISTBASE_FOREACH (ImageTile *, tile, &image.tiles) {
    const image::ImageTileWrapper image_tile = image::ImageTileWrapper(tile);
    tile_user.tile = image_tile.get_tile_number();

    ImBuf *tile_buffer = BKE_image_acquire_ibuf(&image, &tile_user, nullptr);
    if (tile_buffer == nullptr) {
      continue;
    }
    const PixelNodesTileData nodes_tile_pixels(pbvh, image_tile);

    ushort2 tile_resolution(tile_buffer->x, tile_buffer->y);
    BKE_image_release_ibuf(&image, tile_buffer, nullptr);

    PixelCopyTile copy_tile(image_tile.get_tile_number());
    Row per_pixel_solution(tile_resolution.x);

    for (int y = 0; y < tile_resolution.y; y++) {
      per_pixel_solution.reinit(y);
      per_pixel_solution.mask_brush_pixels(nodes_tile_pixels);
      per_pixel_solution.print_debug();
    }

    pbvh_data.tiles_copy_pixels.tiles.append(copy_tile);
  }
}

void BKE_pbvh_pixels_copy_pixels(PBVH &pbvh,
                                 Image &image,
                                 ImageUser &image_user,
                                 image::TileNumber tile_number)
{
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  std::optional<std::reference_wrapper<PixelCopyTile>> pixel_tile =
      pbvh_data.tiles_copy_pixels.find_tile(tile_number);
  if (!pixel_tile.has_value()) {
    return;
  }

  ImageUser tile_user = image_user;
  tile_user.tile = tile_number;
  ImBuf *tile_buffer = BKE_image_acquire_ibuf(&image, &tile_user, nullptr);
  if (tile_buffer == nullptr) {
    return;
  }
  pixel_tile->get().copy_pixels(*tile_buffer);

  BKE_image_release_ibuf(&image, tile_buffer, nullptr);
}

}  // namespace blender::bke::pbvh::pixels
