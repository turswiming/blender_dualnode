
/**
 * Virtual shadowmapping: Tilemap to texture conversion.
 *
 * For all visible light tilemaps, copy page coordinate to a texture.
 * This avoids one level of indirection when evaluating shadows and allows
 * to use a sampler instead of a SSBO bind.
 */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

shared uvec2 update_min;
shared uvec2 update_max;

void main()
{
  if (all(equal(gl_LocalInvocationID, uvec3(0)))) {
    update_min = uvec2(9999);
    update_max = uvec2(0);
  }
  barrier();

  int tilemap_index = int(gl_GlobalInvocationID.z);
  ivec2 tile_co = ivec2(gl_GlobalInvocationID.xy);

  ivec2 atlas_texel = shadow_tile_coord_in_atlas(tile_co, tilemap_index);

  uint lod_size = uint(SHADOW_TILEMAP_RES);

  ShadowTileMapData tilemap_data = tilemaps_buf[tilemap_index];
  int lod_max = tilemap_data.is_cubeface ? SHADOW_TILEMAP_LOD : 0;

  int lod_valid = 0;
  bool do_update = false;
  uvec2 page_valid;
  for (int lod = lod_max; lod >= 0; lod--) {
    int tile_index = shadow_tile_offset(tile_co >> lod, tilemap_index, lod);

    ShadowTileData tile = shadow_tile_unpack(tiles_buf[tile_index]);

    if (tile.is_used && tile.do_update) {
      do_update = true;
    }

    /* Save highest lod for this thread. */
    if (tile.is_used && lod > 0) {
      /* Reload the page in case there was an allocation in the valid thread. */
      page_valid = tile.page;
      lod_valid = lod;
    }
    else if (lod == 0 && lod_valid != 0 && !tile.is_allocated) {
      /* If the tile is not used, store the valid LOD level in LOD0. */
      tile.page = page_valid;
      tile.lod = lod_valid;
      /* This is not a real ownership. It is just a tag so that the shadowing is deemed correct. */
      tile.is_allocated = true;
    }

    if (lod == 0) {
      imageStore(tilemaps_img, atlas_texel, uvec4(shadow_tile_pack(tile)));
    }
  }

  if (do_update) {
    atomicMin(update_min.x, gl_LocalInvocationID.x);
    atomicMin(update_min.y, gl_LocalInvocationID.y);
    atomicMax(update_max.x, gl_LocalInvocationID.x);
    atomicMax(update_max.y, gl_LocalInvocationID.y);
  }

  barrier();

  if (all(equal(gl_LocalInvocationID, uvec3(0)))) {
    if (update_min.x != 9999) {
      int view_index = atomicAdd(pages_infos_buf.view_count, 4);
      view_infos_buf[view_index].persmat = tilemap_data.winmat * tilemap_data.viewmat;
      view_infos_buf[view_index].persinv = inverse(view_infos_buf[view_index].persmat);
      view_infos_buf[view_index].viewmat = tilemap_data.viewmat;
      view_infos_buf[view_index].viewinv = inverse(tilemap_data.viewmat);
      view_infos_buf[view_index].winmat = tilemap_data.winmat;
      view_infos_buf[view_index].wininv = inverse(tilemap_data.winmat);
      // view_infos_buf[view_index].clip_planes;
      // view_infos_buf[view_index].viewvecs;
      // view_infos_buf[view_index].viewcamtexcofac;
      // view_infos_buf[view_index].viewport_size;
      // view_infos_buf[view_index].viewport_size_inverse;
      // view_infos_buf[view_index].frustum_corners;
      // view_infos_buf[view_index].frustum_planes;
      // view_infos_buf[view_index].frustum_bound_sphere;
    }
  }

  if (all(equal(gl_GlobalInvocationID, uvec3(0)))) {
    /* Clamp it as it can underflow if there is too much tile present on screen. */
    pages_infos_buf.page_free_count = max(pages_infos_buf.page_free_count, 0);
  }
}