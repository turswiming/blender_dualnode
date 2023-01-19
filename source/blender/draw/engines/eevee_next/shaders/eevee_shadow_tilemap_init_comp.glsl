
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

ShadowTileDataPacked init_tile_data(ShadowTileDataPacked tile, bool do_update)
{
  if (flag_test(tile, SHADOW_IS_RENDERED)) {
    tile &= ~(SHADOW_DO_UPDATE | SHADOW_IS_RENDERED);
  }
  if (do_update) {
    tile |= SHADOW_DO_UPDATE;
  }
  tile &= ~SHADOW_IS_USED;
  return tile;
}

void main()
{
  uint tilemap_index = gl_GlobalInvocationID.z;
  ShadowTileMapData tilemap = tilemaps_buf[tilemap_index];

  barrier();

  if (all(equal(gl_LocalInvocationID, uvec3(0)))) {
    /* Reset shift to not tag for update more than once per sync cycle. */
    tilemaps_buf[tilemap_index].grid_shift = ivec2(0);

    if (!tilemap.is_cubeface) {
      ShadowTileMapClip clip_data = tilemaps_clip_buf[tilemap_index];
      float clip_near_new = orderedIntBitsToFloat(clip_data.clip_near);
      float clip_far_new = orderedIntBitsToFloat(clip_data.clip_far);
      directional_range_changed = int((clip_near_new != clip_data.clip_near_stored) ||
                                      (clip_far_new != clip_data.clip_far_stored));
      if (directional_range_changed != 0) {
        /* NOTE(fclem): This assumes clip near/far are computed each time the init phase runs. */
        tilemaps_clip_buf[tilemap_index].clip_near_stored = clip_near_new;
        tilemaps_clip_buf[tilemap_index].clip_far_stored = clip_far_new;
        /* Reset for next update. */
        tilemaps_clip_buf[tilemap_index].clip_near = floatBitsToOrderedInt(-FLT_MAX);
        tilemaps_clip_buf[tilemap_index].clip_far = floatBitsToOrderedInt(FLT_MAX);
      }
    }
  }

  barrier();

  ivec2 tile_co = ivec2(gl_GlobalInvocationID.xy);
  ivec2 tile_shifted = tile_co + tilemap.grid_shift;
  ivec2 tile_wrapped = ivec2(tile_shifted % SHADOW_TILEMAP_RES);

  /* If this tile was shifted in and contains old information, update it.
   * Note that cubemap always shift all tiles on update. */
  bool do_update = !in_range_inclusive(tile_shifted, ivec2(0), ivec2(SHADOW_TILEMAP_RES - 1));

  /* TODO(fclem): Might be better to resize the depth stored instead of a full render update. */
  if (!tilemap.is_cubeface && directional_range_changed != 0) {
    do_update = true;
  }

  int lod_max = (tilemap.is_cubeface) ? SHADOW_TILEMAP_LOD : 0;
  uint lod_size = uint(SHADOW_TILEMAP_RES);
  for (int lod = 0; lod <= lod_max; lod++, lod_size >>= 1u) {
    bool thread_active = all(lessThan(tile_co, ivec2(lod_size)));
    ShadowTileDataPacked tile;
    if (thread_active) {
      int tile_load = shadow_tile_offset(tile_wrapped, tilemap.tiles_index, lod);
      tile = init_tile_data(tiles_buf[tile_load], do_update);
    }

    /* Uniform control flow for barrier. Needed to avoid race condition on shifted loads. */
    barrier();

    if (thread_active) {
      int tile_store = shadow_tile_offset(tile_co, tilemap.tiles_index, lod);
      tiles_buf[tile_store] = tile;
    }
  }
}