
/**
 * Compute visibility of each resource bounds for a given view.
 */
/* TODO(fclem): This could be augmented by a 2 pass occlusion culling system. */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(common_intersect_lib.glsl)
#pragma BLENDER_REQUIRE(common_debug_print_lib.glsl)

shared uint shared_result;

void main()
{
  if (gl_GlobalInvocationID.x >= resource_len) {
    return;
  }

  ObjectBounds bounds = bounds_buf[gl_GlobalInvocationID.x];

  if (bounds.bounding_sphere.w != -1.0) {
    IsectBox box = isect_data_setup(bounds.bounding_corners[0].xyz,
                                    bounds.bounding_corners[1].xyz,
                                    bounds.bounding_corners[2].xyz,
                                    bounds.bounding_corners[3].xyz);
    if (intersect_view(box) == false) {
      uint bit = 1u << gl_LocalInvocationID.x;
      atomicAnd(visibility_buf[gl_WorkGroupID.x], ~bit);
    }
  }
}