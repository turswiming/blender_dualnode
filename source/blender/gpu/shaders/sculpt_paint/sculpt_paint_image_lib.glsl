
/* -------------------------------------------------------------------- */
/** \name Brush testing
 * \{ */

vec4 SCULPT_blend_color(vec4 src1, vec4 src2)
{
  return src1 * (1.0 - src2.a) + src2;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Brush testing
 * \{ */

float SCULPT_curve_strength(float factor, int curve_type)
{
  if (factor > 1.0) {
    return 0.0;
  }
  float p = 1.0 - factor;

  switch (curve_type) {
    case 0 /*BRUSH_CURVE_CUSTOM*/:
      return factor;
    case 1 /*BRUSH_CURVE_SMOOTH*/:
      return 3.0 * p * p - 2.0 * p * p * p;
    case 2 /*BRUSH_CURVE_SPHERE*/:
      return sqrt(2.0 * p - p * p);
    case 3 /*BRUSH_CURVE_ROOT*/:
      return sqrt(p);
    case 4 /*BRUSH_CURVE_SHARP*/:
      return p * p;
    case 5 /*BRUSH_CURVE_LIN*/:
      return p;
    case 6 /*BRUSH_CURVE_POW4*/:
      return p * p * p * p;
    case 7 /*BRUSH_CURVE_INVSQUARE*/:
      return p * (2.0 - p);
    case 8 /*BRUSH_CURVE_CONSTANT*/:
      return 1.0;
    case 9 /*BRUSH_CURVE_SMOOTHER*/:
      return p * p * p * (p * (p * 6.0 - 15.0) + 10.0);
    default:
      return factor;
  }

  return factor;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Brush testing
 * \{ */
float SCULPT_hardness_factor(float dist, float hardness, float radius)
{
  float p = dist / radius;
  if (p < hardness) {
    return 0.0;
  }
  else if (hardness >= 1.0) {
    return 1.0;
  }
  return (p - hardness) / (1.0 - hardness);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Brush testing
 * \{ */

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

/** \} */

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
