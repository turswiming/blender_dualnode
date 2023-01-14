
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
    /* Uses data from another LOD. */
    return neon_gradient(float(tile.lod) / float(SHADOW_TILEMAP_LOD));
  }
  if (tile.do_update && tile.is_used) {
    /* Updated. */
    return vec3(0.5, 1, 0);
  }
  if (tile.is_used) {
    /* Used but was cached. */
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

ShadowTileData debug_tile_get(vec3 P, LightData light, out vec2 uv)
{
  vec3 lNg = vec3(1.0, 0.0, 0.0);
  if (light.type == LIGHT_SUN) {
    vec3 lP = shadow_world_to_local(light, P);
    float bias;
    return shadow_directional_tile_get(shadow_tilemaps_tx, light, cameraPos, lP, P, lNg, uv, bias);
  }
  else {
    vec3 lL = light_world_to_local(light, P - light._position);
    float bias;
    return shadow_punctual_tile_get(shadow_tilemaps_tx, light, lL, lNg, uv, bias);
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

/** Return true if a pixel was written. */
bool debug_tilemaps(vec3 P, LightData light)
{
  const int debug_tile_size_px = 4;
  ivec2 px = ivec2(gl_FragCoord.xy) / debug_tile_size_px;
  int tilemap = px.x / SHADOW_TILEMAP_RES;
  int tilemap_index = light.tilemap_index + tilemap;
  if ((px.y < SHADOW_TILEMAP_RES) && (tilemap_index <= light.tilemap_last)) {
    /* Debug actual values in the tilemap buffer. */
    ShadowTileMapData tilemap = tilemaps_buf[tilemap_index];
    int tile_index = shadow_tile_offset(px % SHADOW_TILEMAP_RES, tilemap.tiles_index, 0);
    ShadowTileData tile = shadow_tile_unpack(tiles_buf[tile_index]);
    /* Leave 1 px border between tilemaps. */
    if (!any(
            equal(ivec2(gl_FragCoord.xy) % (SHADOW_TILEMAP_RES * debug_tile_size_px), ivec2(0)))) {
      gl_FragDepth = 0.0;
      out_color_add = vec4(debug_tile_state_color(tile), 0.0);
      out_color_mul = vec4(0.0);

      if (ivec2(gl_FragCoord.xy) == ivec2(0)) {
        drw_print(light.object_mat);
      }
      return true;
    }
  }
  return false;
}

void debug_tile_state(vec3 P, LightData light)
{
  vec2 unused_uv;
  ShadowTileData tile = debug_tile_get(P, light, unused_uv);
  out_color_add = vec4(debug_tile_state_color(tile), 0) * 0.5;
  out_color_mul = vec4(0.5);
}

void debug_atlas_values(vec3 P, LightData light)
{
  vec2 uv;
  ShadowTileData tile = debug_tile_get(P, light, uv);
  float depth = shadow_tile_depth_get(shadow_atlas_tx, tile, uv);

  out_color_add = vec4(vec3(1.0 / (1.0 + depth)), 0);
  out_color_mul = vec4(0.0);
}

void debug_random_tile_color(vec3 P, LightData light)
{
  vec2 unused_uv;
  ShadowTileData tile = debug_tile_get(P, light, unused_uv);
  out_color_add = vec4(debug_random_color(ivec2(tile.page)), 0) * 0.9;
  out_color_mul = vec4(0.1);
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

  if (debug_tilemaps(P, light)) {
    return;
  }

  if (depth != 1.0) {
    switch (eDebugMode(debug_mode)) {
      case DEBUG_SHADOW_TILEMAPS:
        debug_tile_state(P, light);
        break;
      case DEBUG_SHADOW_VALUES:
        debug_atlas_values(P, light);
        break;
      case DEBUG_SHADOW_TILE_RANDOM_COLOR:
        debug_random_tile_color(P, light);
        break;
    }
  }
}