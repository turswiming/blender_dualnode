#pragma BLENDER_REQUIRE(sculpt_paint_tile_lib.glsl)

void main()
{
  PaintTileData paint_tile;
  paint_tile_get_layer(layer_id, paint_tile);
  if (!paint_tile.in_use_frame) {
    return;
  }

  ivec3 coord_in = ivec3(gl_GlobalInvocationID.xy, layer_id);
  vec4 paint_color = imageLoad(paint_tiles_img, coord_in);
  paint_color.a = 1.0;

  ivec2 coord_out = paint_tile_coord_to_udim(paint_tile, coord_in.xy);
  imageStore(texture_img, coord_out, paint_color);
}