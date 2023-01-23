
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(common_shape_lib.glsl)

/* ---------------------------------------------------------------------- */
/** \name Tilemap data
 * \{ */

int shadow_tile_index(ivec2 tile)
{
  return tile.x + tile.y * SHADOW_TILEMAP_RES;
}

ivec2 shadow_tile_coord(int tile_index)
{
  return ivec2(tile_index % SHADOW_TILEMAP_RES, tile_index / SHADOW_TILEMAP_RES);
}

/* Return bottom left pixel position of the tilemap inside the tilemap atlas. */
ivec2 shadow_tilemap_start(int tilemap_index)
{
  return SHADOW_TILEMAP_RES *
         ivec2(tilemap_index % SHADOW_TILEMAP_PER_ROW, tilemap_index / SHADOW_TILEMAP_PER_ROW);
}

ivec2 shadow_tile_coord_in_atlas(ivec2 tile, int tilemap_index)
{
  return shadow_tilemap_start(tilemap_index) + tile;
}

/**
 * Return tile index inside `tiles_buf` for a given tile coordinate inside a specific LOD.
 * `tiles_index` should be `ShadowTileMapData.tiles_index`.
 */
int shadow_tile_offset(ivec2 tile, int tiles_index, int lod)
{
  const int lod0_width = SHADOW_TILEMAP_RES / 1;
  const int lod1_width = SHADOW_TILEMAP_RES / 2;
  const int lod2_width = SHADOW_TILEMAP_RES / 4;
  const int lod3_width = SHADOW_TILEMAP_RES / 8;
  const int lod4_width = SHADOW_TILEMAP_RES / 16;
  const int lod0_size = lod0_width * lod0_width;
  const int lod1_size = lod1_width * lod1_width;
  const int lod2_size = lod2_width * lod2_width;
  const int lod3_size = lod3_width * lod3_width;
  const int lod4_size = lod4_width * lod4_width;

  int offset = tiles_index;
  switch (lod) {
    case 4:
      offset += lod0_size + lod1_size + lod2_size + lod3_size;
      offset += tile.y * lod4_width;
      break;
    case 3:
      offset += lod0_size + lod1_size + lod2_size;
      offset += tile.y * lod3_width;
      break;
    case 2:
      offset += lod0_size + lod1_size;
      offset += tile.y * lod2_width;
      break;
    case 1:
      offset += lod0_size;
      offset += tile.y * lod1_width;
      break;
    case 0:
    default:
      offset += tile.y * lod0_width;
      break;
  }
  offset += tile.x;
  return offset;
}

/** \} */

/* ---------------------------------------------------------------------- */
/** \name Load / Store functions.
 * \{ */

/** \note: Will clamp if out of bounds. */
ShadowTileData shadow_tile_load(usampler2D tilemaps_tx, ivec2 tile_co, int tilemap_index)
{
  /* NOTE(@fclem): This clamp can hide some small imprecision at clipmap transition.
   * Can be disabled to check if the clipmap is well centered. */
  tile_co = clamp(tile_co, ivec2(0), ivec2(SHADOW_TILEMAP_RES - 1));
  uint tile_data =
      texelFetch(tilemaps_tx, shadow_tile_coord_in_atlas(tile_co, tilemap_index), 0).x;
  return shadow_tile_unpack(tile_data);
}

/* This function should be the inverse of ShadowTileMap::tilemap_coverage_get. */
int shadow_directional_clipmap_level(LightData light, vec3 lP)
{
  /* WORKAROUND(fclem): We need to hide one tile worth of data to hide the moving transition.
   * This could be avoided by preselecting 2 levels and checking which one is available based on
   * tile boundaries. But this is expensive. */
  const float narrowing = float(SHADOW_TILEMAP_RES) / (float(SHADOW_TILEMAP_RES) - 1.0001);
  int clipmap_lod = int(ceil(log2(length(lP) * narrowing)));
  /* Since the distance is centered around the camera (and thus by extension the tilemap),
   * we need to multiply by 2 to get the lod level which covers the following range:
   * [-tilemap_coverage_get(lod)/2..tilemap_coverage_get(lod)/2] */
  clipmap_lod += 1;
  return clamp(clipmap_lod, light.clipmap_lod_min, light.clipmap_lod_max);
}

struct ShadowClipmapCoordinates {
  /* Index of the tilemap to containing the tile. */
  int tilemap_index;
  /* LOD of the tile to load relative to the min level. Always positive */
  int clipmap_lod_relative;
  /* Tile coordinate inside the tilemap. */
  ivec2 tile_coord;
  /* UV coordinates in [0..SHADOW_TILEMAP_RES) range. */
  vec2 uv;
};

/* Retain sign bit and avoid costly int division. */
ivec2 shadow_decompress_grid_offset(ivec2 ofs, int relative_lod)
{
  return ((ofs & 0xFFFF) >> relative_lod) - ((ofs >> 16) >> relative_lod);
}

/**
 * \a lP shading point position in light space (world unit) and translated to camera position
 * snapped to smallest clipmap level.
 */
ShadowClipmapCoordinates shadow_directional_coordinates(LightData light, vec3 lP)
{
  ShadowClipmapCoordinates ret;

  int clipmap_lod = shadow_directional_clipmap_level(light, lP - light._position);
  /* This difference needs to be less than 32 for the later shift to be valid.
   * This is ensured by ShadowDirectional::clipmap_level_range(). */
  ret.clipmap_lod_relative = clipmap_lod - light.clipmap_lod_min;

  ret.tilemap_index = light.tilemap_index + ret.clipmap_lod_relative;

  /* Compute offset in tile. */
  ivec2 clipmap_offset = shadow_decompress_grid_offset(light.clipmap_base_offset,
                                                       ret.clipmap_lod_relative);

  ret.uv = lP.xy - vec2(light._clipmap_origin_x, light._clipmap_origin_y);
  ret.uv /= exp2(clipmap_lod);
  ret.uv = ret.uv * float(SHADOW_TILEMAP_RES) + float(SHADOW_TILEMAP_RES / 2);
  ret.uv -= vec2(clipmap_offset);

  /* Clamp to avoid out of tilemap access. */
  ret.tile_coord = clamp(ivec2(ret.uv), ivec2(0.0), ivec2(SHADOW_TILEMAP_RES - 1));
  return ret;
}

/** \} */

/* ---------------------------------------------------------------------- */
/** \name Frustum shapes.
 * \{ */

vec3 shadow_tile_corner_persp(ShadowTileMapData tilemap, ivec2 tile)
{
  return tilemap.corners[1].xyz + tilemap.corners[2].xyz * float(tile.x) +
         tilemap.corners[3].xyz * float(tile.y);
}

Pyramid shadow_tilemap_cubeface_bounds(ShadowTileMapData tilemap,
                                       ivec2 tile_start,
                                       const ivec2 extent)
{
  Pyramid shape;
  shape.corners[0] = tilemap.corners[0].xyz;
  shape.corners[1] = shadow_tile_corner_persp(tilemap, tile_start + ivec2(0, 0));
  shape.corners[2] = shadow_tile_corner_persp(tilemap, tile_start + ivec2(extent.x, 0));
  shape.corners[3] = shadow_tile_corner_persp(tilemap, tile_start + extent);
  shape.corners[4] = shadow_tile_corner_persp(tilemap, tile_start + ivec2(0, extent.y));
  return shape;
}

/** \} */
