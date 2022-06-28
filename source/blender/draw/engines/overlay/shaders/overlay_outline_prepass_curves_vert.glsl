
#pragma BLENDER_REQUIRE(common_view_clipping_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_hair_lib.glsl)

uint outline_colorid_get(void)
{
  int flag = int(abs(ObjectInfo.w));
  bool is_active = (flag & DRW_BASE_ACTIVE) != 0;

  if (isTransform) {
    return 0u; /* colorTransform */
  }
  else if (is_active) {
    return 3u; /* colorActive */
  }
  else {
    return 1u; /* colorSelect */
  }

  return 0u;
}

/* Replace top 2 bits (of the 16bit output) by outlineId.
 * This leaves 16K different IDs to create outlines between objects.
  vec3 world_pos = point_object_to_world(pos);
 * SHIFT = (32 - (16 - 2)) */
#define SHIFT 18u

void main()
{
  bool is_persp = (ProjectionMatrix[3][3] == 0.0);
  float time, thick_time, thickness;
  vec3 world_pos, tan, binor;
  hair_get_pos_tan_binor_time(is_persp,
                              ModelMatrixInverse,
                              ViewMatrixInverse[3].xyz,
                              ViewMatrixInverse[2].xyz,
                              world_pos,
                              tan,
                              binor,
                              time,
                              thickness,
                              thick_time);

  vec4 pos_ndc = point_world_to_ndc(world_pos);

  if (hairThicknessRes > 1) {
    vec3 orig_pos;
    orig_pos = world_pos + binor * -thick_time;
    vec4 orig_pos_ndc = point_world_to_ndc(orig_pos);
    vec3 orig_pos_view = point_world_to_view(orig_pos);
    gl_Position.xyz = orig_pos_view;
    return;
    vec4 d = pos_ndc - orig_pos_ndc;
    float distance = length(d.xy);
    /*if (distance < 0.01) {
      distance = 0.01;
    }*/
    pos_ndc = orig_pos_ndc + d * distance;
  }

  gl_Position = pos_ndc;

#ifdef USE_GEOM
  vert.pos = point_world_to_view(world_pos);
#endif

  /* Small bias to always be on top of the geom. */
  gl_Position.z -= 1e-3;

  /* ID 0 is nothing (background) */
  interp.ob_id = uint(resource_handle + 1);

  /* Should be 2 bits only [0..3]. */
  uint outline_id = outline_colorid_get();

  /* Combine for 16bit uint target. */
  interp.ob_id = (outline_id << 14u) | ((interp.ob_id << SHIFT) >> SHIFT);

  view_clipping_distances(world_pos);
}
