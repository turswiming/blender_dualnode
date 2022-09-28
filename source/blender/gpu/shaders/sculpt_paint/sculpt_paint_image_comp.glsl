#pragma BLENDER_REQUIRE(sculpt_paint_image_lib.glsl)

void main()
{
  PackedPixelRow row = pixel_row_buf[gl_GlobalInvocationID.x + pixel_row_offset];
  TrianglePaintInput triangle = paint_input[PIXEL_ROW_PRIM_INDEX(row)];
  ivec2 image_coord = PIXEL_ROW_START_IMAGE_COORD(row);

  uint row_len = PIXEL_ROW_LEN(row);

  vec3 pos;
  vec3 delta;
  vec3 co1 = vec3(vert_coord_buf[triangle.vert_indices.x]);
  vec3 co2 = vec3(vert_coord_buf[triangle.vert_indices.y]);
  vec3 co3 = vec3(vert_coord_buf[triangle.vert_indices.z]);

  SCULPT_get_row_pos_and_delta(co1, co2, co3, triangle, row, pos, delta);

  for (int i = 0; i < row_len; i++) {
    float distance;
    bool test_result = SCULPT_brush_test_sphere(paint_brush_buf.test, pos, distance);
    // pos += delta;
    if (!test_result) {
      continue;
    }

    imageStore(out_img, image_coord + ivec2(i, 0), paint_brush_buf.color);
  }
}