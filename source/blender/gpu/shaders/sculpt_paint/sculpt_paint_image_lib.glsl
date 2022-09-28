bool SCULPT_brush_test_sphere(PaintBrushTestData test_data, vec3 co, out float dist)
{
  dist = distance(vec3(test_data.location), co);
  if (dist > test_data.radius) {
    return false;
  }

  // TODO(jbakker): do clipping test.
  return true;
}

vec3 interp_v3_v3v3v3(vec3 co1, vec3 co2, vec3 co3, vec3 weights)
{
  /* TODO(jbakker): there should be a glsl function for this. mat3*vec3?*/
  return vec3(co1.x * weights.x + co2.x * weights.y + co3.x * weights.z,
              co1.y * weights.x + co2.y * weights.y + co3.y * weights.z,
              co1.z * weights.x + co2.z * weights.y + co3.z * weights.z);
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
  pos = interp_v3_v3v3v3(co1, co2, co3, barycentric_weights);

  vec2 next_barycentric = row.start_barycentric_coord + triangle.delta_barycentric_coord;
  barycentric_weights = vec3(next_barycentric, 1.0 - next_barycentric - next_barycentric);
  vec3 next_pos = interp_v3_v3v3v3(co1, co2, co3, barycentric_weights);
  delta = next_pos - pos;
}
