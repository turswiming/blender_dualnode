#pragma BLENDER_REQUIRE(sculpt_paint_image_lib.glsl)

bool SCULPT_brush_test(PaintBrushTestData test_data,
                       PaintStepData step_data,
                       vec3 co,
                       out float dist)
{
#if defined(BRUSH_TEST_SPHERE)
  return SCULPT_brush_test_sphere(test_data, step_data, co, dist);
#elif defined(BRUSH_TEST_CIRCLE)
  return SCULPT_brush_test_circle(test_data, step_data, co, dist);
#else
  dist = 0.0;
  return true;
#endif
}

void main()
{
  PackedPixelRow row = pixel_row_buf[gl_GlobalInvocationID.x + pixel_row_offset];
  TrianglePaintInput triangle = paint_input[PIXEL_ROW_PRIM_INDEX(row)];
  ivec2 image_coord = PIXEL_ROW_START_IMAGE_COORD(row);

  uint row_len = PIXEL_ROW_LEN(row);

  vec3 co1 = vec3(vert_coord_buf[triangle.vert_indices.x]);
  vec3 co2 = vec3(vert_coord_buf[triangle.vert_indices.y]);
  vec3 co3 = vec3(vert_coord_buf[triangle.vert_indices.z]);

  vec3 pos;
  vec3 delta;
  SCULPT_get_row_pos_and_delta(co1, co2, co3, triangle, row, pos, delta);

  for (int x = 0; x < row_len; x++) {
    /* TODO: Do clipping test. */
    vec4 color;
    bool color_read = false;

    for (int step_index = paint_step_range[0]; step_index < paint_step_range[1]; step_index++) {
      PaintStepData step_data = paint_step_buf[step_index];

      float distance;
      bool test_result = SCULPT_brush_test(paint_brush_buf.test, step_data, pos, distance);
      if (test_result) {
        if (!color_read) {
          color = imageLoad(out_img, image_coord);
          color_read = true;
        }
        // TODO: blend with color...
        float factor = SCULPT_hardness_factor(distance, step_data.hardness, step_data.radius);
        float curve_factor = SCULPT_curve_strength(factor, paint_brush_buf.falloff_shape);
        vec4 final_paint_color = SCULPT_blend_color(
            color, paint_brush_buf.color * curve_factor * step_data.strength);
        final_paint_color *= paint_brush_buf.alpha;
        color = SCULPT_blend_color(color, final_paint_color);
      }
    }
    if (color_read) {
      imageStore(out_img, image_coord, color);
    }
    image_coord.x += 1;
    pos += delta;
  }
}