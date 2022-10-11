#pragma BLENDER_REQUIRE(sculpt_paint_tile_lib.glsl)

void main()
{
  ivec2 coord_out = ivec2(gl_GlobalInvocationID.xy);
  PaintTileData paint_tile;
  ivec3 coord_in = paint_tile_coord_from_udim(1001, coord_out, paint_tile);
  vec4 paint_color = imageLoad(paint_tiles_img, coord_in);
  paint_color.a = 1.0;
  imageStore(texture_img, coord_out, paint_color);
}