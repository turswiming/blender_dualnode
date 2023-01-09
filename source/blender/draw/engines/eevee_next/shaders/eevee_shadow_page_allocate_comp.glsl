
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

  if (tile.is_used) {
    atomicAdd(pages_infos_buf.page_used_count, 1);
  }
  if (tile.do_update) {
    atomicAdd(pages_infos_buf.page_update_count, 1);
  }
  if (tile.is_allocated) {
    atomicAdd(pages_infos_buf.page_allocated_count, 1);
  }
  if (tile.is_rendered) {
    atomicAdd(pages_infos_buf.page_rendered_count, 1);
  }
  if (tile.is_cached) {
    atomicAdd(pages_infos_buf.page_cached_count, 1);
  }
}