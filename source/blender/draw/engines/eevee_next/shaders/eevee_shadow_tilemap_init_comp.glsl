
/**
 * Virtual shadowmapping: Setup phase for tilemaps.
 *
 * Clear the usage flag.
 * Also tag for update shifted tiles for directional shadow clipmaps.
 * Dispatched with one local thread per LOD0 tile and one workgroup per tilemap.
 */

#pragma BLENDER_REQUIRE(gpu_shader_utildefines_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

shared int directional_range_changed;

void main()
{
  uint tilemap_index = gl_GlobalInvocationID.z;
  ShadowTileMapData tilemap = tilemaps_buf[tilemap_index];

  barrier();

  if (all(equal(gl_LocalInvocationID, uvec3(0)))) {
    /* Reset shift to not tag for update more than once per sync cycle. */
    tilemaps_buf[tilemap_index].grid_shift = ivec2(0);

    if (!tilemap.is_cubeface) {
      directional_range_changed = int(
          (orderedIntBitsToFloat(tilemap.clip_near) != tilemap._clip_near_new) ||
          (orderedIntBitsToFloat(tilemap.clip_far) != tilemap._clip_far_new));
      if (directional_range_changed != 0) {
        /* NOTE(fclem): This assumes clip near/far are computed each time the init phase runs. */
        tilemaps_buf[tilemap_index]._clip_far_stored = tilemaps_buf[tilemap_index]._clip_far_new;
        tilemaps_buf[tilemap_index]._clip_near_stored = tilemaps_buf[tilemap_index]._clip_near_new;
        tilemaps_buf[tilemap_index]._clip_far_new = orderedIntBitsToFloat(tilemap.clip_far);
        tilemaps_buf[tilemap_index]._clip_near_new = orderedIntBitsToFloat(tilemap.clip_near);
        tilemaps_buf[tilemap_index].clip_far = floatBitsToOrderedInt(FLT_MAX);
        tilemaps_buf[tilemap_index].clip_near = floatBitsToOrderedInt(-FLT_MAX);
      }
    }
  }

  barrier();

  ivec2 tile_co = ivec2(gl_GlobalInvocationID.xy);
  ivec2 tile_shifted = ivec2(uvec2(tile_co + tilemap.grid_offset) % SHADOW_TILEMAP_RES);
  tile_shifted += tilemap.grid_shift;

  /* This tile was shifted in and contains old information.
   * Note that cubemap always shift all tiles on update. */
  bool do_update = !in_range_inclusive(tile_shifted, ivec2(0), ivec2(SHADOW_TILEMAP_RES - 1));

  /* TODO(fclem): Might be better to resize the depth stored instead of a full render update. */
  if (directional_range_changed != 0) {
    do_update = true;
  }

  int lod_max = (tilemap.is_cubeface) ? SHADOW_TILEMAP_LOD : 0;
  uint lod_size = uint(SHADOW_TILEMAP_RES);
  for (int lod = 0; lod <= lod_max; lod++, lod_size >>= 1u) {
    if (all(lessThan(tile_co, ivec2(lod_size)))) {
      int tile_index = shadow_tile_offset(tile_co, tilemap.tiles_index, lod);
      ShadowTileDataPacked tile = tiles_buf[tile_index];

      if (do_update) {
        tile |= SHADOW_DO_UPDATE;
      }
      else if (flag_test(tile, SHADOW_IS_RENDERED)) {
        tile &= ~SHADOW_DO_UPDATE;
        tile &= ~SHADOW_IS_RENDERED;
      }
      tile &= ~SHADOW_IS_USED;

      tiles_buf[tile_index] = tile;
    }
  }
}