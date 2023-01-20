
/**
 * Compute visibility of each resource bounds for a given view.
 */
/* TODO(fclem): This could be augmented by a 2 pass occlusion culling system. */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(common_intersect_lib.glsl)

shared uint shared_result;

void mask_visibility_bit(uint view_id)
{
  if (view_len > 1) {
    uint index = gl_GlobalInvocationID.x * uint(visibility_word_per_draw) + (view_id / 32u);
    visibility_buf[index] &= ~(1u << view_id);
  }
  else {
    atomicAnd(visibility_buf[gl_WorkGroupID.x], ~(1u << gl_LocalInvocationID.x));
  }
}

void main()
{
  if (gl_GlobalInvocationID.x >= resource_len) {
    return;
  }

  mat4 model_mat = matrix_buf[gl_GlobalInvocationID.x].model;
  ObjectBounds bounds = bounds_buf[gl_GlobalInvocationID.x];

  if (bounds.test_enabled) {
    vec3 origin = transform_point(model_mat, bounds.center - bounds.size);
    vec3 side_x = transform_point(model_mat, bounds.center + bounds.size * vec3(1, -1, -1)) -
                  origin;
    vec3 side_y = transform_point(model_mat, bounds.center + bounds.size * vec3(-1, 1, -1)) -
                  origin;
    vec3 side_z = transform_point(model_mat, bounds.center + bounds.size * vec3(-1, -1, 1)) -
                  origin;

    vec3 center = origin + ((side_x + side_y + side_z) / 2.0f);
    float radius = distance(origin, center);
    float inscribed_radius = min(min(bounds.size.x, bounds.size.y), bounds.size.z);

    IsectBox box = isect_data_setup(origin, side_x, side_y, side_z);
    Sphere bounding_sphere = Sphere(center, radius);
    Sphere inscribed_sphere = Sphere(center, inscribed_radius);

    for (drw_view_id = 0; drw_view_id < view_len; drw_view_id++) {
      if (intersect_view(inscribed_sphere)) {
        /* Visible */
      }
      else if (intersect_view(bounding_sphere) == false || intersect_view(box) == false) {
        /* Not visible. */
        mask_visibility_bit(drw_view_id);
      }
    }
  }
}
