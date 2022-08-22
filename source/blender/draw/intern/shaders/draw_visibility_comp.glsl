
/**
 * Compute visibility of each resource bounds for a given view.
 */
/* TODO(fclem): This could be augmented by a 2 pass occlusion culling system. */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)

shared uint shared_result[4];

void main()
{
  if (gl_LocalInvocationID.x == 0) {
    shared_result[0] = 0u;
    shared_result[1] = 0u;
    shared_result[2] = 0u;
    shared_result[3] = 0u;
  }

  barrier();

  uint resource_id = gl_GlobalInvocationID.x;

  barrier();

  if (gl_LocalInvocationID.x == 0) {
    visibility_buf[gl_WorkGroupID.x] = uvec4(
        shared_result[0], shared_result[1], shared_result[2], shared_result[3]);
  }
}