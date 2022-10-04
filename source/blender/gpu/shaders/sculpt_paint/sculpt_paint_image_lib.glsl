bool SCULPT_brush_test_sphere(PaintBrushTestData test_data,
                              PaintStepData step_data,
                              vec3 co,
                              out float dist)
{
  dist = distance(vec3(step_data.location), co);
  return dist <= step_data.radius;
}

float plane_point_side_v3(vec4 plane, vec3 co)
{
  return dot(co, vec3(plane)) + plane.w;
}

vec3 closest_to_plane_normalized_v3(vec4 plane, vec3 co)
{
  float side = plane_point_side_v3(plane, co);
  return (vec3(plane) * -side) + co;
}

bool SCULPT_brush_test_circle(PaintBrushTestData test_data,
                              PaintStepData step_data,
                              vec3 co,
                              out float dist)
{
  vec3 proj = closest_to_plane_normalized_v3(step_data.plane_view, co);
  return SCULPT_brush_test_sphere(test_data, step_data, proj, dist);
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
