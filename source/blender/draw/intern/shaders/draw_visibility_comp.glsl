
/**
 * Compute visibility of each resource bounds for a given view.
 */
/* TODO(fclem): This could be augmented by a 2 pass occlusion culling system. */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(common_intersect_lib.glsl)

const uint uint_len_per_group = gl_WorkGroupSize.x / 32u;
shared uint shared_result[uint_len_per_group];

void main()
{
  if (gl_LocalInvocationID.x == 0) {
    for (uint i = 0; i < uint_len_per_group; i++) {
      shared_result[i] = 0u;
    }
  }

  barrier();

  uint resource_id = min(gl_GlobalInvocationID.x, resource_len - 1u);

  ObjectBounds bounds = bounds_buf[resource_id];
  if (bounds.bounding_sphere.w != -1.0) {
    IsectBox box = isect_data_setup(bounds.bounding_corners[0].xyz,
                                    bounds.bounding_corners[1].xyz,
                                    bounds.bounding_corners[2].xyz,
                                    bounds.bounding_corners[3].xyz);
    if (intersect_view(box)) {
      uint result = 1u << (gl_LocalInvocationID.x % 32u);
      atomicOr(shared_result[gl_LocalInvocationID.x / 32u], result);
    }
  }

  barrier();

  if (gl_LocalInvocationID.x == 0) {
    for (uint i = 0; i < uint_len_per_group; i++) {
      visibility_buf[gl_WorkGroupID.x * uint_len_per_group + i] = shared_result[i];
    }
  }
}