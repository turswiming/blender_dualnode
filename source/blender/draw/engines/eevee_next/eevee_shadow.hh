/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation.
 */

/** \file
 * \ingroup eevee
 *
 * The shadow module manages shadow update tagging & shadow rendering.
 */

#pragma once

#include "BLI_pool.hh"
#include "BLI_vector.hh"

#include "GPU_batch.h"

#include "eevee_material.hh"
#include "eevee_shader.hh"
#include "eevee_shader_shared.hh"

namespace blender::eevee {

class Instance;
class ShadowModule;
class ShadowPipeline;
struct Light;

enum eCubeFace {
  /* Ordering by culling order. If cone aperture is shallow, we cull the later view. */
  Z_NEG = 0,
  X_POS,
  X_NEG,
  Y_POS,
  Y_NEG,
  Z_POS,
};

/* To be applied after view matrix. Follow same order as eCubeFace. */
constexpr static const float shadow_face_mat[6][4][4] = {
    {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}},   /* Z_NEG */
    {{0, 0, -1, 0}, {-1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}}, /* X_POS */
    {{0, 0, 1, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},   /* X_NEG */
    {{1, 0, 0, 0}, {0, 0, -1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},  /* Y_POS */
    {{-1, 0, 0, 0}, {0, 0, 1, 0}, {0, 1, 0, 0}, {0, 0, 0, 1}},  /* Y_NEG */
    {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, 1}}, /* Z_POS */
};

/* Converts to [-SHADOW_TILEMAP_RES / 2..SHADOW_TILEMAP_RES / 2] for XY and [0..1] for Z. */
constexpr static const float shadow_clipmap_scale_mat[4][4] = {{SHADOW_TILEMAP_RES / 2, 0, 0, 0},
                                                               {0, SHADOW_TILEMAP_RES / 2, 0, 0},
                                                               {0, 0, 0.5, 0},
                                                               {0, 0, 0.5, 1}};

/* -------------------------------------------------------------------- */
/** \name Tile-Map
 *
 * Stores indirection table and states of each tile of a virtual shadow-map.
 * One tile-map has the effective resolution of `pagesize * tile_map_resolution`.
 * Each tile-map overhead is quite small if they do not have any pages allocated.
 *
 * \{ */

struct ShadowTileMap : public ShadowTileMapData {
  static constexpr int64_t tile_map_resolution = SHADOW_TILEMAP_RES;
  static constexpr int64_t tiles_count = tile_map_resolution * tile_map_resolution;
  /**
   * Maximum "bounding" angle of a tile inside a cube-map.
   * Half the diagonal of tile since we test using the tile center.
   */
  static float tile_cone_half_angle;

  /** Level of detail for clipmap. */
  int level = INT_MAX;
  /** Cube face index. */
  eCubeFace cubeface = Z_NEG;
  /** Cached, used for detecting updates. */
  float4x4 object_mat;
  /** Near and far clip distances. For clip-map they are updated after sync. */
  float near, far;

 public:
  ShadowTileMap(int tiles_index_)
  {
    tiles_index = tiles_index_;
    this->set_dirty();
  }

  void sync_clipmap(const float3 &camera_position,
                    const float4x4 &object_mat_,
                    float near_,
                    float far_,
                    int2 origin_offset,
                    int clipmap_level);

  void sync_cubeface(
      const float4x4 &object_mat, float near, float far, float cone_aperture, eCubeFace face);

  static float tilemap_coverage_get(int lvl)
  {
    /* This function should be kept in sync with shadow_directional_clipmap_level(). */
    /* \note: If we would to introduce a global scaling option it would be here. */
    return powf(2.0f, lvl);
  }

  static float tile_size_get(int lvl)
  {
    return tilemap_coverage_get(lvl) / tile_map_resolution;
  }

  float4x4 winmat_get() const;

  void debug_draw() const;

  /* For external callers. Use this in order to not miss an update. */
  void set_level(int clipmap_level)
  {
    if (assign_if_different(level, clipmap_level)) {
      set_dirty();
    }
  }

  void set_is_cubemap(bool1 is_cubemap_)
  {
    if (assign_if_different(is_cubeface, is_cubemap_)) {
      set_dirty();
    }
  }

  void set_dirty()
  {
    grid_shift = int2(SHADOW_TILEMAP_RES);
  }

  void set_updated()
  {
    grid_shift = int2(0);
  }
};

/**
 * The tile-maps are managed on CPU and associated with each light shadow object.
 *
 * The number of tile-maps & tiles is unbounded (to the limit of SSBOs), but the actual number
 * used for rendering is caped to 4096. This is to simplify tile-maps management on CPU.
 *
 * At sync end, all tile-maps are grouped by light inside the ShadowTileMapDataBuf so that each
 * light has a contiguous range of tile-maps to refer to.
 *
 * The tile-map atlas has a fixed 64x64 size. So it can contain 4096 tile-map of 16x16 pixels each.
 */
struct ShadowTileMapPool {
 public:
  /** Limit the width of the texture. */
  static constexpr int64_t maps_per_row = SHADOW_TILEMAP_PER_ROW;

  /** Vector containing available tile range indices in the ShadowTileDataBuf. */
  Vector<uint> free_indices;
  /** Pool containing shadow tile structure on CPU. */
  Pool<ShadowTileMap> tilemap_pool;
  /** Sorted descriptions for each tilemap in the pool. Updated each frame. */
  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  /** Previously used tile-maps that needs to release their tiles/pages. Updated each frame. */
  ShadowTileMapDataBuf tilemaps_unused = {"tilemaps_unused"};
  /** All possible tiles. A range of tiles tile is referenced by a tile-map. */
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  /** Texture equivalent of ShadowTileDataBuf but grouped by light. */
  Texture tilemap_tx = {"tilemap_tx"};
  /** If false a clear pass will be issue to init the page heaps and clear the tiles ownership. */
  bool tilemap_initialized = false;
  /** Number of free tile-maps at the end of the previous sync. */
  int64_t last_free_len = 0;

 public:
  ShadowTileMapPool();

  ShadowTileMap *acquire();

  /**
   * Push the given list of ShadowTileMap onto the free stack. Their pages will be free.
   */
  void release(Span<ShadowTileMap *> free_list);

  void end_sync(ShadowModule &module);
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Shadow Casters & Receivers
 *
 * \{ */

/* Can be either a shadow caster or a shadow receiver. */
struct ShadowObject {
  ResourceHandle resource_handle = {0};
  bool used = true;
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name ShadowModule
 *
 * Manages shadow atlas and shadow region data.
 * \{ */

class ShadowModule {
  friend ShadowPunctual;
  friend ShadowDirectional;
  friend ShadowPipeline;
  friend ShadowTileMapPool;

 public:
  /** Need to be first because of destructor order. */
  ShadowTileMapPool tilemap_pool;

  Pool<ShadowPunctual> punctual_pool;
  Pool<ShadowDirectional> directional_pool;

 private:
  Instance &inst_;

  /** Map of shadow casters to track deletion & update of intersected shadows. */
  Map<ObjectKey, ShadowObject> objects_;

  /* -------------------------------------------------------------------- */
  /** \name Tilemap Management
   * \{ */

  PassSimple tilemap_setup_ps_ = {"TilemapSetup"};
  PassMain tilemap_usage_ps_ = {"TagUsage"};
  PassSimple tilemap_update_ps_ = {"TilemapUpdate"};

  PassMain::Sub *tilemap_usage_transparent_ps_ = nullptr;
  GPUBatch *box_batch_ = nullptr;

  Framebuffer usage_tag_fb;

  /** List of Resource IDs (to get bounds) for tagging passes. */
  StorageVectorBuffer<uint, 128> past_casters_updated_ = {"PastCastersUpdated"};
  StorageVectorBuffer<uint, 128> curr_casters_updated_ = {"CurrCastersUpdated"};

  int3 dispatch_depth_scan_size_;
  float tilemap_pixel_radius_;
  float screen_pixel_radius_inv_;

  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name Page Management
   * \{ */

  Texture atlas_tx_ = {"shadow_atlas_tx_"};

  /** Pool of unallocated pages waiting to be assigned to specific tiles in the tilemap atlas. */
  ShadowPageHeapBuf pages_free_data_ = {"PagesFreeBuf"};
  /** Pool of cached tiles waiting to be reused. */
  ShadowPageCacheBuf pages_cached_data_ = {"PagesCachedBuf"};
  /** Infos for book keeping and debug. */
  ShadowPagesInfoDataBuf pages_infos_data_ = {"PagesInfosBuf"};

  int3 copy_dispatch_size_;
  int3 scan_dispatch_size_;
  int rendering_tilemap_;
  int rendering_lod_;
  bool do_full_update = true;

  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name Rendering
   * \{ */

  /** Multi-View containing a maximum of 64 view to be rendered with the shadow pipeline. */
  View shadow_multi_view_ = {"ShadowMultiView", 64};

  /** An array mapping view index to tilemap index. */
  UniformArrayBuffer<int, 64> view_to_tilemap_buf_;

  /** \} */

  /* -------------------------------------------------------------------- */
  /** \name Debugging
   * \{ */

  /** Display informations about the virtual shadows. */
  PassSimple debug_draw_ps_ = {"Shadow.Debug"};

  /** \} */

  /** Scene immutable parameters. */
  int shadow_page_size_ = 256;
  /** Maximum number of allocated pages. Maximum value is SHADOW_MAX_TILEMAP. */
  int shadow_page_len_ = SHADOW_MAX_TILEMAP;

 public:
  ShadowModule(Instance &inst) : inst_(inst){};
  ~ShadowModule(){};

  void init();

  void begin_sync();
  /** Register a shadow caster or receiver. */
  void sync_object(const ObjectHandle &handle,
                   const ResourceHandle &resource_handle,
                   bool is_shadow_caster,
                   bool is_alpha_blend);
  void end_sync();

  void set_lights_data();

  void set_view(View &view);

  void debug_end_sync();
  void debug_draw(View &view, GPUFrameBuffer *view_fb);

  template<typename T> void bind_resources(draw::detail::PassBase<T> *pass)
  {
    pass->bind_texture(SHADOW_ATLAS_TEX_SLOT, &atlas_tx_);
    pass->bind_texture(SHADOW_TILEMAPS_TEX_SLOT, &tilemap_pool.tilemap_tx);
  }

 private:
  void remove_unused();
  void debug_page_map_call(DRWPass *pass);
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Shadow
 *
 * A shadow component is associated to a `eevee::Light` and manages its associated Tile-maps.
 * \{ */

class ShadowPunctual : public NonCopyable, NonMovable {
 private:
  ShadowModule &shadows_;
  /** Tile-map for each cube-face needed (in eCubeFace order). */
  Vector<ShadowTileMap *> tilemaps_;
  /** Area light size. */
  float size_x_, size_y_;
  /** Shape type. */
  eLightType light_type_;
  /** Random position on the light. In world space. */
  float3 random_offset_;
  /** Light position. */
  float3 position_;
  /** Near and far clip distances. */
  float far_, near_;
  /** View space offset to apply to the shadow. */
  float bias_;
  /** Number of tile-maps needed to cover the light angular extents. */
  int tilemaps_needed_;
  /** Visibility cone angle from the light source. */
  int cone_aperture_;

 public:
  ShadowPunctual(ShadowModule &module) : shadows_(module){};
  ShadowPunctual(ShadowPunctual &&other)
      : shadows_(other.shadows_), tilemaps_(std::move(other.tilemaps_)){};

  ~ShadowPunctual()
  {
    shadows_.tilemap_pool.release(tilemaps_);
  }

  /**
   * Sync shadow parameters but do not allocate any shadow tile-maps.
   */
  void sync(eLightType light_type,
            const float4x4 &object_mat,
            float cone_aperture,
            float near_clip,
            float far_clip,
            float bias);

  /**
   * Release the tile-maps that will not be used in the current frame.
   */
  void release_excess_tilemaps();

  /**
   * Allocate shadow tile-maps and setup views for rendering.
   */
  void end_sync(Light &light);
};

class ShadowDirectional : public NonCopyable, NonMovable {
 private:
  ShadowModule &shadows_;
  /** Tile-map for each clip-map level. */
  Vector<ShadowTileMap *> tilemaps_;
  /** User minimum resolution. */
  float min_resolution_;
  /** View space offset to apply to the shadow. */
  float bias_;
  /** Near and far clip distances. For clip-map, when they are updated after sync. */
  float near_, far_;
  /** Offset of the lowest clip-map relative to the highest one. */
  int2 base_offset_;
  /** Copy of object matrix. Normalized. */
  float4x4 object_mat_;
  /** Current range of clip-map levels covered by this shadow. */
  IndexRange lods_range;

 public:
  ShadowDirectional(ShadowModule &module) : shadows_(module){};
  ShadowDirectional(ShadowDirectional &&other)
      : shadows_(other.shadows_), tilemaps_(std::move(other.tilemaps_)){};

  ~ShadowDirectional()
  {
    shadows_.tilemap_pool.release(tilemaps_);
  }

  /**
   * Sync shadow parameters but do not allocate any shadow tile-maps.
   */
  void sync(const float4x4 &object_mat, float bias, float min_resolution);

  /**
   * Release the tile-maps that will not be used in the current frame.
   */
  void release_excess_tilemaps(const Camera &camera);

  /**
   * Allocate shadow tile-maps and setup views for rendering.
   */
  void end_sync(Light &light, const Camera &camera);

 private:
  IndexRange clipmap_level_range(const Camera &camera);
};

/** \} */

}  // namespace blender::eevee
