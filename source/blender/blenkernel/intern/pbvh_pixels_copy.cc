/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#include "BLI_array.hh"
#include "BLI_bit_vector.hh"
#include "BLI_math.h"
#include "BLI_math_vector.hh"
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

rcti get_bounds(const Edge<CoordSpace::Tile> &tile_edge)
{
  rcti bounds;
  BLI_rcti_init_minmax(&bounds);
  BLI_rcti_do_minmax_v(&bounds, int2(tile_edge.v1.co));
  BLI_rcti_do_minmax_v(&bounds, int2(tile_edge.v2.co));
  return bounds;
}

void add_margin(rcti &bounds, int margin)
{
  bounds.xmin -= margin;
  bounds.xmax += margin;
  bounds.ymin -= margin;
  bounds.ymax += margin;
}

void clamp(rcti &bounds, int2 resolution)
{
  rcti clamping_bounds;
  int2 xy;
  BLI_rcti_init(&clamping_bounds, 0, resolution.x - 1, 0, resolution.y - 1);
  BLI_rcti_clamp(&bounds, &clamping_bounds, xy);
}

const Vertex<CoordSpace::Tile> convert_coord_space(const Vertex<CoordSpace::UV> &uv_vertex,
                                                   const image::ImageTileWrapper image_tile,
                                                   const int2 tile_resolution)
{
  return Vertex<CoordSpace::Tile>{(uv_vertex.co - float2(image_tile.get_tile_offset())) *
                                  float2(tile_resolution)};
}

const Edge<CoordSpace::Tile> convert_coord_space(const Edge<CoordSpace::UV> &uv_edge,
                                                 const image::ImageTileWrapper image_tile,
                                                 const int2 tile_resolution)
{
  return Edge<CoordSpace::Tile>{
      convert_coord_space(uv_edge.v1, image_tile, tile_resolution),
      convert_coord_space(uv_edge.v2, image_tile, tile_resolution),
  };
}

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
        edge.v2.co = find_uv_vert(mesh_primitive, mesh_edge.vert2).uv;
        append(edge);
      }
    }
  }

  NonManifoldTileEdges extract_tile_edges(const image::ImageTileWrapper image_tile,
                                          const int2 tile_resolution) const
  {
    NonManifoldTileEdges result;
    // TODO: Only add edges that intersects with the given tile.
    // TODO: Clamp edges to tile bounds.

    for (const Edge<CoordSpace::UV> &uv_edge : *this) {
      const Edge<CoordSpace::Tile> tile_edge = convert_coord_space(
          uv_edge, image_tile, tile_resolution);
      result.append(tile_edge);
    }
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

struct Rows {
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
      float distance = std::numeric_limits<float>::max();
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
    Row() = delete;
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
            pixels[x].distance = 0.0f;
          }
        }
      }
    }

    void determine_copy_pixels(const NonManifoldTileEdges &tile_edges,
                               int search_margin,
                               const int2 tile_resolution)
    {
      for (const Edge<CoordSpace::Tile> &tile_edge : tile_edges) {
        rcti edge_bounds = get_bounds(tile_edge);
        add_margin(edge_bounds, search_margin);
        clamp(edge_bounds, tile_resolution);

        if (edge_bounds.ymax < row_number) {
          continue;
        }
        if (edge_bounds.ymin > row_number) {
          continue;
        }

        for (const int x : IndexRange(edge_bounds.xmin, edge_bounds.xmax - edge_bounds.xmin)) {
          Elem &pixel = pixels[x];
          switch (pixel.type) {
            case PixelType::Brush: {
              break;
            }
            case PixelType::Undecided:
            case PixelType::CopyFromClosestEdge: {
              const float2 point(pixel.copy_command.destination);
              float2 closest_edge_point;
              closest_to_line_v2(closest_edge_point, point, tile_edge.v1.co, tile_edge.v2.co);
              float distance_to_edge = blender::math::distance_squared(closest_edge_point, point);
              if (distance_to_edge > pixel.distance) {
                break;
              }

              // TODO:
              //  Find upto 2 valid pixels to copy from.
              //  Determine the mix factor between the two pixels
              //  Store result in the pixel.command.

              pixel.distance = distance_to_edge;
              pixel.type = PixelType::CopyFromClosestEdge;

              break;
            }
          }
        }
      }
    }

    void solution2(Rows &rows)
    {
      for (int x : pixels.index_range()) {
        Elem &elem = pixels[x];
        /* Skip pixels that already used directly via a brush. */
        if (elem.type == PixelType::Brush) {
          continue;
        }

        rcti bounds;
        BLI_rcti_init(&bounds, x, x, row_number, row_number);
        add_margin(bounds, 8);
        clamp(bounds, rows.resolution);

        float found_distance = std::numeric_limits<float>().max();
        int2 found_source(0);

        for (int sy : IndexRange(bounds.ymin, BLI_rcti_size_y(&bounds))) {
          Row &row = rows.rows[sy];
          for (int sx : IndexRange(bounds.xmin, BLI_rcti_size_x(&bounds))) {
            Elem &source = row.pixels[sx];
            if (source.type != PixelType::Brush) {
              continue;
            }
            float new_distance = blender::math::distance(float2(sx, sy), float2(x, row_number));
            if (found_distance > new_distance) {
              found_source = int2(sx, sy);
              found_distance = new_distance;
            }
          }
        }

        if (found_distance == std::numeric_limits<float>().max()) {
          continue;
        }
        elem.type = PixelType::CopyFromClosestEdge;
        elem.distance = found_distance;
        elem.copy_command.source_1 = found_source;
        // TODO: find second source by looking at neighbouring pixels of source_1.
        elem.copy_command.source_2 = found_source;
        elem.copy_command.mix_factor = 0.0f;
      }
    }

    static bool can_be_extended_with(const PixelCopyGroup &group, const PixelCopyCommand &command)
    {
      PixelCopyCommand last_command = last_copy_command(group);
      /* Can only extend when pushing the next pixel. */
      if (last_command.destination.x != command.destination.x - 1 ||
          last_command.destination.y != command.destination.y) {
        return false;
      }
      /* Can only extend when */
      int2 delta_source_1 = last_command.source_1 - command.source_1;
      if (max_ii(UNPACK2(blender::math::abs(delta_source_1))) > 127) {
        return false;
      }
      return true;
    }

    static void extend_with(PixelCopyGroup &group, const PixelCopyCommand &command)
    {
      PixelCopyCommand last_command = last_copy_command(group);
      PixelCopyItem new_item = {char2(command.source_1 - last_command.source_1),
                                char2(command.source_2 - command.source_1),
                                uint8_t(command.mix_factor * 255)};
      group.items.append(new_item);
    }

    static PixelCopyCommand last_copy_command(const PixelCopyGroup &group)
    {
      PixelCopyCommand last_command(group);
      for (const PixelCopyItem &item : group.items) {
        last_command.apply(item);
      }
      return last_command;
    }

    void pack_into(Vector<PixelCopyGroup> groups) const
    {
      for (const Elem &elem : pixels) {
        if (elem.type == PixelType::CopyFromClosestEdge) {
          if (groups.is_empty() || !can_be_extended_with(groups.last(), elem.copy_command)) {
            PixelCopyGroup new_group = {
                elem.copy_command.destination - int2(1, 0), elem.copy_command.source_1, {}};
            groups.append(new_group);
          }
          extend_with(groups.last(), elem.copy_command);
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

  int2 resolution;
  int margin;
  Vector<Row> rows;

  Rows(int2 resolution, int margin, const PixelNodesTileData &node_tile_pixels)
      : resolution(resolution), margin(margin)
  {
    Row row_template(resolution.x);
    rows.resize(resolution.y, row_template);
    for (int row_number : rows.index_range()) {
      rows[row_number].reinit(row_number);
      rows[row_number].mask_brush_pixels(node_tile_pixels);
    }
  }

  void find_copy_source()
  {
    for (Row &row : rows) {
      row.solution2(*this);
    }
  }
};  // namespace blender::bke::pbvh::pixels

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
  const NonManifoldUVEdges non_manifold_edges(mesh_data);
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

    int2 tile_resolution(tile_buffer->x, tile_buffer->y);
    BKE_image_release_ibuf(&image, tile_buffer, nullptr);

    NonManifoldTileEdges tile_edges = non_manifold_edges.extract_tile_edges(image_tile,
                                                                            tile_resolution);
    PixelCopyTile copy_tile(image_tile.get_tile_number());

    Rows rows(tile_resolution, image.seam_margin, nodes_tile_pixels);
    rows.find_copy_source();

    for (int y = 0; y < tile_resolution.y; y++) {
      Rows::Row &row = rows.rows[y];
      // row.determine_copy_pixels(tile_edges, image.seam_margin, tile_resolution);
      row.print_debug();
      row.pack_into(copy_tile.groups);
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
