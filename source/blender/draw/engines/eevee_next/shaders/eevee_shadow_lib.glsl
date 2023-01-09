
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

/** \a unormalized_uv is the uv coordinates for the whole tilemap [0..SHADOW_TILEMAP_RES]. */
vec2 shadow_page_uv_transform(uvec2 page, uint lod, vec2 unormalized_uv)
{
  vec2 page_texel = fract(unormalized_uv / float(1u << lod));
  /* Fix float imprecision that can make some pixel sample the wrong page. */
  page_texel *= 0.999999;
  /* Assumes atlas is squared. */
  return (vec2(page) + page_texel) / vec2(SHADOW_PAGE_PER_ROW);
}

/* ---------------------------------------------------------------------- */
/** \name Shadow Sampling Functions
 * \{ */

/* Turns local light coordinate into shadow region index. Matches eCubeFace order. */
int shadow_punctual_face_index_get(vec3 lL)
{
  vec3 aP = abs(lL);
  if (all(greaterThan(aP.xx, aP.yz))) {
    return (lL.x > 0.0) ? 1 : 2;
  }
  else if (all(greaterThan(aP.yy, aP.xz))) {
    return (lL.y > 0.0) ? 3 : 4;
  }
  else {
    return (lL.z > 0.0) ? 5 : 0;
  }
}

/* Transform vector to face local coordinate. */
vec3 shadow_punctual_local_position_to_face_local(int face_id, vec3 lL)
{
  switch (face_id) {
    case 1:
      return vec3(-lL.y, lL.z, -lL.x);
    case 2:
      return vec3(lL.y, lL.z, lL.x);
    case 3:
      return vec3(lL.x, lL.z, -lL.y);
    case 4:
      return vec3(-lL.x, lL.z, lL.y);
    case 5:
      return vec3(lL.x, -lL.y, -lL.z);
    default:
      return lL;
  }
}

/* Returns minimum bias needed for a given geometry birlak abd a shadowmap page. */
float shadow_slope_bias_get(LightData light, vec3 lNg, int lod)
{
  if (is_punctual) {
    /* Assume 90° aperture projection frustum for punctual shadows. Note near/far do not mater. */
    /* TODO(fclem): Optimize. This can be easily precomputed. */
    mat4 normal_mat = invert(transpose(projection__perspective(-light.clip_near,
                                                               light.clip_near,
                                                               -light.clip_near,
                                                               light.clip_near,
                                                               light.clip_near,
                                                               light.clip_far)));
    vec3 ndc_Ng = transform_direction(normal_mat, lNg);
  }

  /* Get slope from normal vector. */
  vec2 ndc_slope = ndc_Ng.xy / -ndc_Ng.z;
  /* Slope bias definition from fixed pipeline. */
  float bias = dot(vec2(1.0), abs(ndc_slope));
  /* Bias for 1 pixel of LOD 0. */
  bias *= 1.0 / (SHADOW_TILEMAP_RES * shadow_page_size_);
  /* Compensate for each increasing lod level as the space between pixels increases. */
  bias *= float(1u << lod);
  return bias;
}

ShadowTileData shadow_punctual_tile_get(
    usampler2D tilemaps_tx, LightData light, vec3 lL, vec3 lNg, out vec2 uv, out float bias)
{
  int face_id = shadow_punctual_face_index_get(lL);
  lL = shadow_punctual_local_position_to_face_local(face_id, lL);
  lNg = shadow_punctual_local_position_to_face_local(face_id, lNg);
  /* UVs in [-1..+1] range. */
  uv = lL.xy / abs(lL.z);
  /* UVs in [0..SHADOW_TILEMAP_RES] range. */
  const float lod0_res = float(SHADOW_TILEMAP_RES / 2);
  uv = uv * lod0_res + lod0_res;
  ivec2 tile_co = ivec2(floor(uv));
  int tilemap_index = light.tilemap_index + face_id;
  ShadowTileData tile = shadow_tile_load(tilemaps_tx, tile_co, tilemap_index);
  bias = shadow_slope_bias_get(light, lNg, tile.lod);
  return tile;
}

ShadowTileData shadow_directional_tile_get(
    usampler2D tilemaps_tx, LightData light, vec3 camera_P, vec3 lP, vec3 P, vec3 lNg, out vec2 uv)
{
  int clipmap_lod = shadow_directional_clipmap_level(light, distance(P, camera_P));
  int clipmap_lod_relative = clipmap_lod - light.clipmap_lod_min;
  int tilemap_index = clamp(
      light.tilemap_index + clipmap_lod_relative, light.tilemap_index, light.tilemap_last);
  /* Compute how many time we need to subdivide. */
  float clipmap_res_mul = float(1 << (light.clipmap_lod_max - clipmap_lod));
  /* Compute offset of the clipmap from the largest LOD. */
  vec2 clipmap_offset = vec2(abs(light.clipmap_base_offset) >> clipmap_lod_relative) *
                        sign(light.clipmap_base_offset);

  uv = (lP.xy * clipmap_res_mul - clipmap_offset) + float(SHADOW_TILEMAP_RES / 2);
  ivec2 tile_co = ivec2(floor(uv));
  ShadowTileData tile = shadow_tile_load(tilemaps_tx, tile_co, tilemap_index);
  bias = shadow_slope_bias_get(light, lNg, tile.lod);
  return tile;
}

float shadow_tile_depth_get(sampler2D atlas_tx, ShadowTileData tile, vec2 uv)
{
  float depth = 1.0;
  if (tile.is_allocated) {
    vec2 shadow_uv = shadow_page_uv_transform(tile.page, tile.lod, uv);
    depth = texture(atlas_tx, shadow_uv).r;
  }
  return depth;
}

/* Return world distance delta from light between shading point and first occluder. */
float shadow_sample(sampler2D atlas_tx,
                    usampler2D tilemaps_tx,
                    LightData light,
                    vec3 lL,
                    vec3 lNg,
                    float receiver_dist,
                    vec3 P)
{

  if (light.type == LIGHT_SUN) {
    /* [-SHADOW_TILEMAP_RES/2..SHADOW_TILEMAP_RES/2] range for highest LOD. */
    vec3 lP = transform_point(light.object_mat, P);
    vec2 uv;
    ShadowTileData tile = shadow_directional_tile_get(
        tilemaps_tx, light, cameraPos, lP, P, lNg, uv, bias);
    float occluder_z = shadow_tile_depth_get(atlas_tx, tile, uv);
    /* Transform to world space distance. */
    return (lP.z - occluder_z) * abs(light.clip_far - light.clip_near);
  }
  else {
    vec2 uv;
    ShadowTileData tile = shadow_punctual_tile_get(tilemaps_tx, light, lL, uv, bias);
    float occluder_z = shadow_tile_depth_get(atlas_tx, tile, uv);
    /* TODO(fclem): Store linear depth. */
    occluder_z = linear_depth(true, occluder_z, light.clip_far, light.clip_near);
    /* Take into account the cubemap projection. We want the radial distance. */
    float occluder_dist = receiver_dist * occluder_z / max_v3(abs(lL));
    return receiver_dist - occluder_dist;
  }
}

/** \} */
