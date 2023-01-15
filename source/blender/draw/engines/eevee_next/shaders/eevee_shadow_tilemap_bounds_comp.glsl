
/**
 * Virtual shadowmapping: Bounds computation for directional shadows.
 *
 * Iterate through all shadow casters and extract min/max per directional shadow.
 * This needs to happen first in the pipeline to allow tagging all relevant tilemap as dirty if
 * their range changes.
 */

#pragma BLENDER_REQUIRE(gpu_shader_utildefines_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(common_intersect_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_light_iter_lib.glsl)

shared int global_min;
shared int global_max;

void main()
{
  uint index = gl_GlobalInvocationID.x;
  /* Keep uniform control flow. Do not return. */
  index = min(index, uint(resource_len) - 1);

  uint resource_id = casters_id_buf[index];
  ObjectBounds bounds = bounds_buf[resource_id];
  IsectBox box = isect_data_setup(bounds.bounding_corners[0].xyz,
                                  bounds.bounding_corners[1].xyz,
                                  bounds.bounding_corners[2].xyz,
                                  bounds.bounding_corners[3].xyz);

  LIGHT_FOREACH_BEGIN_DIRECTIONAL(light_cull_buf, l_idx)
  {
    LightData light = light_buf[l_idx];

    if (gl_LocalInvocationID.x == 0) {
      global_min = floatBitsToOrderedInt(FLT_MAX);
      global_max = floatBitsToOrderedInt(-FLT_MAX);
    }

    barrier();

    float local_min = FLT_MAX;
    float local_max = -FLT_MAX;
    for (int i = 0; i < 8; i++) {
      const float epsilon = 1e-16;
      float z = dot(box.corners[i].xyz, light._back);
      local_min = min(local_min, z - epsilon);
      local_max = max(local_max, z + epsilon);
    }

    /* Intermediate result. Min/Max of a compute group. */
    atomicMin(global_min, floatBitsToOrderedInt(local_min));
    atomicMax(global_max, floatBitsToOrderedInt(local_max));

    barrier();

    if (gl_LocalInvocationID.x == 0) {
      /* Final result. Min/Max of the whole dispatch. */
      atomicMin(light_buf[l_idx].clip_far, global_min);
      atomicMax(light_buf[l_idx].clip_near, global_max);
      /* TODO(fclem): This feel unecessary but we currently have no indexing from
       * tilemap to lights. This is because the lights are selected by culling phase. */
      for (int i = light.tilemap_index; i <= light.tilemap_last; i++) {
        atomicMin(tilemaps_buf[i].clip_far, global_min);
        atomicMax(tilemaps_buf[i].clip_near, global_max);
      }
    }
  }
  LIGHT_FOREACH_END
}