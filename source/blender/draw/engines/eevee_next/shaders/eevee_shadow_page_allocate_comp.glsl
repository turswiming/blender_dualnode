
/**
 * Virtual shadowmapping: Allocation.
 *
 * Allocates pages to tiles needing them.
 * Note that allocation can fail, in this case the tile is left with no page.
 */

#pragma BLENDER_REQUIRE(eevee_shadow_page_ops_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

void main()
{
  ShadowTileMapData tilemap_data = tilemaps_buf[gl_GlobalInvocationID.z];

  int tile_index = tilemap_data.tiles_index * SHADOW_TILEDATA_PER_TILEMAP +
                   int(gl_LocalInvocationID.x);

  ShadowTileData tile = shadow_tile_unpack(tiles_buf[tile_index]);
  if (tile.is_used && !tile.is_allocated) {
    shadow_page_alloc(tile);
    tiles_buf[tile_index] = shadow_tile_pack(tile);
  }

  /* Gather pages to clear here because we have a one thread per page dispatch for this compute. */
  if (tile.is_used && tile.do_update) {
    uint clear_page_index = atomicAdd(clear_dispatch_buf.num_groups_z, 1u);
    clear_page_buf[clear_page_index] = packUvec2x16(tile.page);
  }
}