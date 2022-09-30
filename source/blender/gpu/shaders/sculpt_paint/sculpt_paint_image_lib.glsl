bool SCULPT_brush_test_sphere(PaintBrushTestData test_data,
                              PaintStepData step_data,
                              vec3 co,
                              out float dist)
{
  dist = distance(vec3(step_data.location), co);
  return dist <= step_data.radius;
}

void SCULPT_get_row_pos_and_delta(vec3 co1,
                                  vec3 co2,
                                  vec3 co3,
                                  TrianglePaintInput triangle,
                                  PackedPixelRow row,
                                  out vec3 pos,
                                  out vec3 delta)
{

  vec3 barycentric_weights = vec3(row.start_barycentric_coord,
                                  1.0 - row.start_barycentric_coord.x -
                                      row.start_barycentric_coord.y);
  mat3 coords = mat3(co1, co2, co3);
  pos = coords * barycentric_weights;

  vec2 next_barycentric = row.start_barycentric_coord + triangle.delta_barycentric_coord;
  barycentric_weights = vec3(next_barycentric, 1.0 - next_barycentric.x - next_barycentric.y);
  vec3 next_pos = coords * barycentric_weights;
  delta = next_pos - pos;
}
