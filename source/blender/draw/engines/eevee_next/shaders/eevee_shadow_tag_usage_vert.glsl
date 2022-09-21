
/**
 * Virtual shadowmapping: Usage tagging
 *
 * Shadow pages are only allocated if they are visible.
 * This pass scan the depth buffer and tag all tiles that are needed for light shadowing as
 * needed.
 */

#pragma BLENDER_REQUIRE(common_view_lib.glsl)

void main()
{
  interp.P = point_object_to_world(pos);
  interp.vP = point_world_to_view(interp.P);

  gl_Position = point_world_to_ndc(interp.P);
}