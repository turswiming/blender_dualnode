
#pragma BLENDER_REQUIRE(common_view_clipping_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  GPU_INTEL_VERTEX_SHADER_WORKAROUND

  /* Use local variable to workaround macro unrolling issue in shaderc. */
  vec3 pos = pPosition;
  gl_Position = point_world_to_ndc(pos);
  finalColor = pColor;
  gl_PointSize = pSize;

  view_clipping_distances(pos);
}
