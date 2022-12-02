
#pragma BLENDER_REQUIRE(common_view_clipping_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  GPU_INTEL_VERTEX_SHADER_WORKAROUND

  /* Spir-V GLSL pre-processors chokes when passed directly. */
  vec3 pos = pPosition;
  gl_Position = point_world_to_ndc(pos);
  finalColor = pColor;
  gl_PointSize = pSize;

  view_clipping_distances(pos);
}
