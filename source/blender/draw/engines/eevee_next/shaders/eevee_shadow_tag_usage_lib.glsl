
/**
 * Virtual shadowmapping: Usage tagging
 *
 * Shadow pages are only allocated if they are visible.
 * This pass scan the depth buffer and tag all tiles that are needed for light shadowing as
 * needed.
 */

#pragma BLENDER_REQUIRE(common_intersect_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_geom_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_light_iter_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_light_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_lib.glsl)

void shadow_tag_usage_tilemap(uint l_idx, vec3 P, float dist_to_cam, const bool is_directional)
{
  LightData light = light_buf[l_idx];

  if (light.tilemap_index == LIGHT_NO_SHADOW) {
    return;
  }

  int lod = 0;
  ivec2 tile_co;
  int tilemap_index = light.tilemap_index;
  if (is_directional) {
    int clipmap_lod = shadow_directional_clipmap_level(light, dist_to_cam);
    int clipmap_lod_relative = clipmap_lod - light.clipmap_lod_min;
    /* Compute how many time we need to subdivide. */
    float clipmap_res_mul = float(1 << (light.clipmap_lod_max - clipmap_lod));
    /* Compute offset of the clipmap from the largest LOD. */
    vec2 clipmap_offset = vec2(abs(light.clipmap_base_offset) >> clipmap_lod_relative) *
                          sign(light.clipmap_base_offset);

    /* [-SHADOW_TILEMAP_RES/2..SHADOW_TILEMAP_RES/2] range for highest LOD. */
    /* TODO mat only rotate to light space. Need to apply offset + scale to remap to
     * [-SHADOW_TILEMAP_RES/2..SHADOW_TILEMAP_RES/2] */
    vec3 lP = transform_point(light.object_mat, P);
    tile_co = ivec2(floor(lP.xy * clipmap_res_mul - clipmap_offset)) + SHADOW_TILEMAP_RES / 2;
    tile_co = clamp(tile_co, ivec2(0), ivec2(SHADOW_TILEMAP_RES - 1));
    tilemap_index += clipmap_lod_relative;
    tilemap_index = clamp(tilemap_index, light.tilemap_index, light.tilemap_last);
  }
  else {
    vec3 lL = light_world_to_local(light, P - light._position);
    float dist_to_light = length(lL);
    if (dist_to_light > light.influence_radius_max) {
      return;
    }
    /* How much a shadow map pixel covers a final image pixel. */
    float footprint_ratio = dist_to_light * (tilemap_pixel_radius * screen_pixel_radius_inv);
    /* Project the radius to the screen. 1 unit away from the camera the same way
     * pixel_world_radius_inv was computed. Not needed in orthographic mode. */
    bool is_persp = (ProjectionMatrix[3][3] == 0.0);
    if (is_persp) {
      footprint_ratio /= dist_to_cam;
    }
    lod = int(ceil(-log2(footprint_ratio)));
    lod = clamp(lod, 0, SHADOW_TILEMAP_LOD);

    int face_id = shadow_punctual_face_index_get(lL);
    lL = shadow_punctual_local_position_to_face_local(face_id, lL);

    uint lod_res = uint(SHADOW_TILEMAP_RES) >> uint(lod);
    tile_co = ivec2(((lL.xy / abs(lL.z)) * 0.5 + 0.5) * float(lod_res));
    tile_co = clamp(tile_co, ivec2(0), ivec2(lod_res - 1));
    tilemap_index += face_id;
  }

  int tile_index = shadow_tile_offset(tile_co, tilemap_index, lod);
  atomicOr(tiles_buf[tile_index], SHADOW_IS_USED);
}

void shadow_tag_usage(vec3 vP, vec3 P, vec2 pixel)
{
  float dist_to_cam = length(vP);

  LIGHT_FOREACH_BEGIN_DIRECTIONAL(light_cull_buf, l_idx)
  {
    shadow_tag_usage_tilemap(l_idx, P, dist_to_cam, true);
  }
  LIGHT_FOREACH_END

  LIGHT_FOREACH_BEGIN_LOCAL(light_cull_buf, light_zbin_buf, light_tile_buf, pixel, vP.z, l_idx)
  {
    shadow_tag_usage_tilemap(l_idx, P, dist_to_cam, false);
  }
  LIGHT_FOREACH_END
}