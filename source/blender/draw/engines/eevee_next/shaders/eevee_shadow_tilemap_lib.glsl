
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
 */
int shadow_tile_offset(ivec2 tile, int tilemap_index, int lod)
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

  int offset = tilemap_index * SHADOW_TILEDATA_PER_TILEMAP;
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
  /* Bias to avoid sampling outside of the clipmap level. This leaves some padding between each
   * level because of the rounding of the camera tile position, there might be cases where the
   * camera will see further than the clipmap allows. Tweak if needed. */
  float bias = 0.18;
  /* Why do we need to bias by 2 here? I don't know... */
  int clipmap_lod = int(ceil(log2(distance_to_camera) + bias)) + 2;
  return clamp(clipmap_lod, light.clipmap_lod_min, light.clipmap_lod_max);
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
