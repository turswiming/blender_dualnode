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

#include "PIL_time_utildefines.h"

namespace blender::bke::pbvh::pixels {

/** Coordinate space of a coordinate. */
enum class CoordSpace {
  /**
   * Coordinate is in UV coordinate space. As in unmodified from mesh data.
   */
  UV,

  /**
   * Coordinate is in Tile coordinate space.
   *
   * With tile coordinate space each unit is a single pixel of the tile.
   * Range is [0..buffer width].
   */
  Tile,
};

template<CoordSpace Space> struct Vertex {
  float2 coordinate;
};

template<CoordSpace Space> struct Edge {
  Vertex<Space> vertex_1;
  Vertex<Space> vertex_2;
};

/** Calculate the bounds of the given edge. */
rcti get_bounds(const Edge<CoordSpace::Tile> &tile_edge)
{
  rcti bounds;
  BLI_rcti_init_minmax(&bounds);
  BLI_rcti_do_minmax_v(&bounds, int2(tile_edge.vertex_1.coordinate));
  BLI_rcti_do_minmax_v(&bounds, int2(tile_edge.vertex_2.coordinate));
  return bounds;
}

/** Add a margin to the given bounds. */
void add_margin(rcti &bounds, int margin)
{
  bounds.xmin -= margin;
  bounds.xmax += margin;
  bounds.ymin -= margin;
  bounds.ymax += margin;
}

/** Clamp bounds to be between 0,0 and the given resolution. */
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
  return Vertex<CoordSpace::Tile>{(uv_vertex.coordinate - float2(image_tile.get_tile_offset())) *
                                  float2(tile_resolution)};
}

const Edge<CoordSpace::Tile> convert_coord_space(const Edge<CoordSpace::UV> &uv_edge,
                                                 const image::ImageTileWrapper image_tile,
                                                 const int2 tile_resolution)
{
  return Edge<CoordSpace::Tile>{
      convert_coord_space(uv_edge.vertex_1, image_tile, tile_resolution),
      convert_coord_space(uv_edge.vertex_2, image_tile, tile_resolution),
  };
}

class NonManifoldTileEdges : public Vector<Edge<CoordSpace::Tile>> {};

class NonManifoldUVEdges : public Vector<Edge<CoordSpace::UV>> {
 public:
  NonManifoldUVEdges(const uv_islands::MeshData &mesh_data)
  {
    int num_non_manifold_edges = count_non_manifold_edges(mesh_data);
    reserve(num_non_manifold_edges);
    for (const int primitive_id : mesh_data.looptris.index_range()) {
      for (const int edge_id : mesh_data.primitive_to_edge_map[primitive_id]) {
        if (is_manifold(mesh_data, edge_id)) {
          continue;
        }
        const MLoopTri &loop_tri = mesh_data.looptris[primitive_id];
        const uv_islands::MeshEdge &mesh_edge = mesh_data.edges[edge_id];
        Edge<CoordSpace::UV> edge;

        edge.vertex_1.coordinate = find_uv(mesh_data, loop_tri, mesh_edge.vert1);
        edge.vertex_2.coordinate = find_uv(mesh_data, loop_tri, mesh_edge.vert2);
        append(edge);
      }
    }
    BLI_assert_msg(size() == num_non_manifold_edges,
                   "Incorrect number of non manifold edges added. ");
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
    for (const int primitive_id : mesh_data.looptris.index_range()) {
      for (const int edge_id : mesh_data.primitive_to_edge_map[primitive_id]) {
        if (is_manifold(mesh_data, edge_id)) {
          continue;
        }
        result += 1;
      }
    }
    return result;
  }

  static bool is_manifold(const uv_islands::MeshData &mesh_data, const int edge_id)
  {
    return mesh_data.edge_to_primitive_map[edge_id].size() == 2;
  }

  static float2 find_uv(const uv_islands::MeshData &mesh_data,
                        const MLoopTri &loop_tri,
                        int vertex_i)
  {
    for (int i = 0; i < 3; i++) {
      int loop_i = loop_tri.tri[i];
      const MLoop &loop = mesh_data.loops[loop_i];
      if (loop.v == vertex_i) {
        return mesh_data.uv_map[loop_i];
      }
    }
    BLI_assert_unreachable();
    return float2(0.0f);
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
      Selected,
      /** This pixel will be copid from another pixel to solve non-manifold edge bleeding. */
      CopyFromClosestEdge,
    };

    struct Pixel {
      PixelType type = PixelType::Undecided;
      float distance = std::numeric_limits<float>::max();
      CopyPixelCommand copy_command;

      Pixel() = default;

      Pixel(int2 coordinate)
      {
        copy_command.destination = coordinate;
        copy_command.source_1 = coordinate;
        copy_command.source_2 = coordinate;
        copy_command.mix_factor = 0.0f;
      }
    };

    int row_number = 0;
    Array<Pixel> pixels;
    Row() = delete;
    Row(int64_t width) : pixels(width)
    {
    }

    void reinit(int y)
    {
      row_number = y;
      for (int x = 0; x < pixels.size(); x++) {
        pixels[x] = Pixel(int2(x, y));
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

    void mark_for_evaluation(const Rows &rows, const NonManifoldTileEdges &tile_edges)
    {
      for (const Edge<CoordSpace::Tile> &tile_edge : tile_edges) {
        rcti edge_bounds = get_bounds(tile_edge);
        add_margin(edge_bounds, rows.margin);
        clamp(edge_bounds, rows.resolution);

        if (edge_bounds.ymax < row_number) {
          continue;
        }
        if (edge_bounds.ymin > row_number) {
          continue;
        }

        for (const int x : IndexRange(edge_bounds.xmin, BLI_rcti_size_x(&edge_bounds))) {
          Pixel &pixel = pixels[x];
          if (pixel.type != PixelType::Undecided) {
            continue;
          }

          const float2 point(pixel.copy_command.destination);
          float2 closest_edge_point;
          closest_to_line_segment_v2(closest_edge_point,
                                     point,
                                     tile_edge.vertex_1.coordinate,
                                     tile_edge.vertex_2.coordinate);
          float distance_to_edge = blender::math::distance(closest_edge_point, point);
          if (distance_to_edge < rows.margin) {
            pixel.type = PixelType::Selected;
          }
        }
      }
    }

    int2 find_second_source(Rows &rows, int2 destination, int2 first_source)
    {
      rcti search_bounds;
      BLI_rcti_init(&search_bounds,
                    max_ii(destination.x - 1, 0),
                    min_ii(destination.x + 1, rows.resolution.x - 1),
                    max_ii(destination.y - 1, 0),
                    min_ii(destination.y + 1, rows.resolution.y - 1));
      /* Initialize to the first source, so when no other source could be found it will use the
       * first_source. */
      int2 found_source = first_source;
      float found_distance = std::numeric_limits<float>().max();
      for (int sy : IndexRange(search_bounds.ymin, BLI_rcti_size_y(&search_bounds) + 1)) {
        for (int sx : IndexRange(search_bounds.xmin, BLI_rcti_size_x(&search_bounds) + 1)) {
          int2 source(sx, sy);
          /* Skip first source as it should be the closest and already selected. */
          if (source == first_source) {
            continue;
          }
          if (rows.rows[sy].pixels[sx].type != PixelType::Brush) {
            continue;
          }

          float new_distance = blender::math::distance(float2(destination), float2(source));
          if (new_distance < found_distance) {
            found_distance = new_distance;
            found_source = source;
          }
        }
      }
      return found_source;
    }

    float determine_mix_factor(int2 destination, int2 source_1, int2 source_2)
    {
      if (source_1 == source_2) {
        return 0.0f;
      }
      return dist_to_line_segment_v2(float2(destination), float2(source_1), float2(source_2));
    }

    void find_copy_source(Rows &rows)
    {
      for (int x : pixels.index_range()) {
        Pixel &elem = pixels[x];
        /* Skip pixels that are not selected for evaluation. */
        if (elem.type != PixelType::Selected) {
          continue;
        }

        rcti bounds;
        BLI_rcti_init(&bounds, x, x, row_number, row_number);
        add_margin(bounds, rows.margin);
        clamp(bounds, rows.resolution);

        float found_distance = std::numeric_limits<float>().max();
        int2 found_source(0);

        for (int sy : IndexRange(bounds.ymin, BLI_rcti_size_y(&bounds))) {
          Row &row = rows.rows[sy];
          for (int sx : IndexRange(bounds.xmin, BLI_rcti_size_x(&bounds))) {
            Pixel &source = row.pixels[sx];
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
        elem.copy_command.source_2 = find_second_source(
            rows, elem.copy_command.destination, found_source);
        elem.copy_command.mix_factor = determine_mix_factor(
            elem.copy_command.destination, elem.copy_command.source_1, elem.copy_command.source_2);
      }
    }

    static bool can_be_extended_with(const CopyPixelGroup &group, const CopyPixelCommand &command)
    {
      CopyPixelCommand last_command = last_copy_command(group);
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

    static void extend_with(CopyPixelGroup &group, const CopyPixelCommand &command)
    {
      CopyPixelCommand last_command = last_copy_command(group);
      DeltaCopyPixelCommand delta_command = last_command.encode_delta(command);
      group.deltas.append(delta_command);
    }

    static CopyPixelCommand last_copy_command(const CopyPixelGroup &group)
    {
      CopyPixelCommand last_command(group);
      for (const DeltaCopyPixelCommand &item : group.deltas) {
        last_command.apply(item);
      }
      return last_command;
    }

    void pack_into(Vector<CopyPixelGroup> &groups) const
    {
      for (const Pixel &elem : pixels) {
        if (elem.type == PixelType::CopyFromClosestEdge) {
          if (groups.is_empty() || !can_be_extended_with(groups.last(), elem.copy_command)) {
            CopyPixelGroup new_group = {
                elem.copy_command.destination - int2(1, 0), elem.copy_command.source_1, {}};
            groups.append(new_group);
          }
          extend_with(groups.last(), elem.copy_command);
        }
      }
    }

    void print_debug() const
    {
      for (const Pixel &elem : pixels) {
        printf("%d", elem.type);
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
      row.find_copy_source(*this);
    }
  }

  void mark_for_evaluation(const NonManifoldTileEdges &tile_edges)
  {
    for (Row &row : rows) {
      row.mark_for_evaluation(*this, tile_edges);
    }
  }

  void pack_into(Vector<CopyPixelGroup> &groups) const
  {
    for (const Row &row : rows) {
      row.pack_into(groups);
    }
  }

  void print_debug() const
  {
    for (const Row &row : rows) {
      row.print_debug();
    }
    printf("\n");
  }

};  // namespace blender::bke::pbvh::pixels

static void copy_pixels_reinit(CopyPixelTiles &tiles)
{
  tiles.clear();
}

void BKE_pbvh_pixels_copy_update(PBVH &pbvh,
                                 Image &image,
                                 ImageUser &image_user,
                                 const uv_islands::MeshData &mesh_data)
{
  TIMEIT_START(pbvh_pixels_copy_update);
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  copy_pixels_reinit(pbvh_data.tiles_copy_pixels);
  const NonManifoldUVEdges non_manifold_edges(mesh_data);
  if (non_manifold_edges.is_empty()) {
    /* Early exit: No non manifold edges detected. */
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
    CopyPixelTile copy_tile(image_tile.get_tile_number());

    Rows rows(tile_resolution, image.seam_margin, nodes_tile_pixels);
    rows.mark_for_evaluation(tile_edges);
    rows.find_copy_source();
    rows.pack_into(copy_tile.groups);
    pbvh_data.tiles_copy_pixels.tiles.append(copy_tile);
  }
  TIMEIT_END(pbvh_pixels_copy_update);
}

void BKE_pbvh_pixels_copy_pixels(PBVH &pbvh,
                                 Image &image,
                                 ImageUser &image_user,
                                 image::TileNumber tile_number)
{
  // TIMEIT_START(pbvh_pixels_copy_pixels);
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  std::optional<std::reference_wrapper<CopyPixelTile>> pixel_tile =
      pbvh_data.tiles_copy_pixels.find_tile(tile_number);
  if (!pixel_tile.has_value()) {
    /* No pixels need to be copied. */
    return;
  }

  ImageUser tile_user = image_user;
  tile_user.tile = tile_number;
  ImBuf *tile_buffer = BKE_image_acquire_ibuf(&image, &tile_user, nullptr);
  if (tile_buffer == nullptr) {
    /* No tile buffer found to copy. */
    return;
  }
  pixel_tile->get().copy_pixels(*tile_buffer);

  BKE_image_release_ibuf(&image, tile_buffer, nullptr);
  // TIMEIT_END(pbvh_pixels_copy_pixels);
}

}  // namespace blender::bke::pbvh::pixels
