
/**
 * Debug drawing for virtual shadowmaps.
 * See eShadowDebug for more information.
 */

#pragma BLENDER_REQUIRE(common_debug_print_lib.glsl)
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_light_iter_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_light_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_sampling_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

/** Control the scaling of the tilemap splat. */
const float pixel_scale = 4.0;

vec3 debug_random_color(ivec2 v)
{
  float r = interlieved_gradient_noise(vec2(v), 0.0, 0.0);
  return hue_gradient(r);
}

vec3 debug_random_color(int v)
{
  return debug_random_color(ivec2(v, 0));
}

void debug_tile_print(ShadowTileData tile, ivec4 tile_coord)
{
  drw_print("Tile (", tile_coord.x, ",", tile_coord.y, ") in Tilemap ", tile_coord.z, " : ");
  drw_print(tile.lod);
  drw_print(tile.page);
  drw_print(tile.cache_index);
}

vec3 debug_tile_state_color(ShadowTileData tile)
{
  if (tile.lod > 0) {
    return vec3(1, 0.5, 0) * float(tile.lod) / float(SHADOW_TILEMAP_LOD);
  }
  if (tile.do_update && tile.is_used) {
    /* Error color because at this point (at the end of the frame), there should be no visible
     * update remaining. */
    return vec3(1, 0, 0);
  }
  if (tile.is_used) {
    return vec3(0, 1, 0);
  }
  vec3 col = vec3(0);
  if (tile.is_cached) {
    col += vec3(0.2, 0, 0.5);
    if (tile.do_update) {
      col += vec3(0.8, 0, 0);
    }
  }
  return col;
}

ShadowTileData debug_tile_get(vec3 P, LightData light)
{
  vec2 uv;
  if (light.type == LIGHT_SUN) {
    /* [-SHADOW_TILEMAP_RES/2..SHADOW_TILEMAP_RES/2] range for highest LOD. */
    vec3 lP = transform_point(light.object_mat, P);
    return shadow_directional_tile_get(shadow_tilemaps_tx, light, cameraPos, lP, P, uv);
  }
  else {
    vec3 lL = light_world_to_local(light, P - light._position);
    return shadow_punctual_tile_get(shadow_tilemaps_tx, light, lL, uv);
  }
}

LightData debug_light_get()
{
  LIGHT_FOREACH_BEGIN_LOCAL_NO_CULL(light_cull_buf, l_idx)
  {
    LightData light = light_buf[l_idx];
    if (light.tilemap_index == debug_tilemap_index) {
      return light;
    }
  }
  LIGHT_FOREACH_END

  LIGHT_FOREACH_BEGIN_DIRECTIONAL(light_cull_buf, l_idx)
  {
    LightData light = light_buf[l_idx];
    if (light.tilemap_index == debug_tilemap_index) {
      return light;
    }
  }
  LIGHT_FOREACH_END
}

void debug_tile_state(vec3 P, LightData light)
{
  ShadowTileData tile = debug_tile_get(P, light);
  out_color_add = vec4(debug_tile_state_color(tile), 0) * 0.5;
  out_color_mul = vec4(0.5);
}

void main()
{
  /* Default to no output. */
  gl_FragDepth = 1.0;
  out_color_add = vec4(0.0);
  out_color_mul = vec4(1.0);

  float depth = texelFetch(hiz_tx, ivec2(gl_FragCoord.xy), 0).r;
  vec3 P = get_world_space_from_depth(uvcoordsvar.xy, depth);
  /* Make it pass the depth test. */
  gl_FragDepth = depth - 1e-6;

  LightData light = debug_light_get();

  if (depth != 1.0) {
    switch (eDebugMode(debug_mode)) {
      case DEBUG_SHADOW_TILEMAPS:
        debug_tile_state(P, light);
        break;
    }
  }
}