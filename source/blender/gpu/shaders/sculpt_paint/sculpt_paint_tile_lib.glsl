ivec2 paint_tile_coord_to_sub_tile_id(ivec2 coord)
{
  return coord / ivec2(SUB_TILE_SIZE);
}

bool paint_tile_search(int tile_number, ivec2 coord, out PaintTileData r_paint_tile)
{
  int2 sub_tile_id = paint_tile_coord_to_sub_tile_id(coord);
  for (int i = 0; i < paint_tile_buf_len; i++) {
    if (paint_tile_buf[i].tile_number == tile_number &&
        paint_tile_buf[i].sub_tile_id == sub_tile_id) {
      r_paint_tile = paint_tile_buf[i];
      return true;
    }
  }
  return false;
}

void paint_tile_get_layer(int layer_id, out PaintTileData r_paint_tile)
{
  r_paint_tile = paint_tile_buf[layer_id];
}

ivec3 paint_tile_coord_from_udim(int tile_number, ivec2 coord, out PaintTileData r_paint_tile)
{
  if (paint_tile_search(tile_number, coord, r_paint_tile)) {
    return ivec3(coord - r_paint_tile.sub_tile_id * ivec2(SUB_TILE_SIZE), r_paint_tile.layer_id);
  }

  return ivec3(0);
}

ivec2 paint_tile_coord_to_udim(PaintTileData paint_tile, ivec2 coord)
{
  return paint_tile.sub_tile_id * ivec2(SUB_TILE_SIZE) + coord;
}
