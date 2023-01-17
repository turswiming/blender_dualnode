
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

  int offset = tiles_index * SHADOW_TILEDATA_PER_TILEMAP;
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
int shadow_directional_clipmap_level(LightData light, float distance_to_camera)
{
  /* Since the distance is centered around the camera (and thus by extension the tilemap),
   * we need to multiply by 2 to get the lod level which covers the following range:
   * [-tilemap_coverage_get(lod)/2..tilemap_coverage_get(lod)/2] */
  int clipmap_lod = int(ceil(log2(distance_to_camera))) + 1;
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
ivec2 divide_by_two_n(ivec2 val, int exponent)
{
  return (abs(val) >> exponent) * sign(val);
}

/**
 * \a lP shading point position in light space (world unit).
 */
ShadowClipmapCoordinates shadow_directional_coordinates(LightData light,
                                                        vec3 lP,
                                                        float distance_to_camera)
{
  ShadowClipmapCoordinates ret;

  int clipmap_lod = shadow_directional_clipmap_level(light, distance_to_camera);
  /* This difference needs to be less than 32 for the later shift to be valid.
   * This is ensured by ShadowDirectional::clipmap_level_range(). */
  ret.clipmap_lod_relative = clipmap_lod - light.clipmap_lod_min;

  ret.tilemap_index = light.tilemap_index + ret.clipmap_lod_relative;

  /* Compute offset of the clipmap from the largest LOD. */
  ivec2 clipmap_offset = divide_by_two_n(light.clipmap_base_offset, ret.clipmap_lod_relative);
  // clipmap_offset += float(SHADOW_TILEMAP_RES / 2);

  /* Compute how many time we need to subdivide. */
  // float clipmap_res_mul = float(1 << (light.clipmap_lod_max - clipmap_lod));
  /* Scale to [-SHADOW_TILEMAP_RES/2..SHADOW_TILEMAP_RES/2] range for largest LOD. */
  // clipmap_res_mul *= light._clipmap_scale; /* TODO !!!! DOESNT MATCH OFFSET SCALING */

  /* TODO(fclem): This could be optimized. */
  float level_size = pow(2.0, float(clipmap_lod));
  /* [0..SHADOW_TILEMAP_RES] range for target LOD. */
  ret.uv = ((lP.xy / level_size) + 0.5) * float(SHADOW_TILEMAP_RES) - vec2(clipmap_offset);

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

vec3 shadow_tile_corner_ortho(ShadowTileMapData tilemap, ivec2 tile, const bool far)
{
  return tilemap.corners[0].xyz + tilemap.corners[1].xyz * float(tile.x) +
         tilemap.corners[2].xyz * float(tile.y) + tilemap.corners[3].xyz * float(far);
}

Box shadow_tilemap_clipmap_bounds(ShadowTileMapData tilemap, ivec2 tile_start, const ivec2 extent)
{
  Box shape;
  shape.corners[0] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(0, 0), false);
  shape.corners[1] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(extent.x, 0), false);
  shape.corners[2] = shadow_tile_corner_ortho(tilemap, tile_start + extent, false);
  shape.corners[3] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(0, extent.y), false);
  shape.corners[4] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(0, 0), true);
  shape.corners[5] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(extent.x, 0), true);
  shape.corners[6] = shadow_tile_corner_ortho(tilemap, tile_start + extent, true);
  shape.corners[7] = shadow_tile_corner_ortho(tilemap, tile_start + ivec2(0, extent.y), true);
  return shape;
}

/** \} */
