
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

/** \a unormalized_uv is the uv coordinates for the whole tilemap [0..SHADOW_TILEMAP_RES]. */
vec2 shadow_page_uv_transform(uvec2 page, uint lod, vec2 unormalized_uv)
{
  /* TODO(fclem): It should be possible to just saturate(unormalized_uv - tile_co << lod). */
  vec2 page_texel = fract(unormalized_uv / float(1u << lod));
  /* Fix float imprecision that can make some pixel sample the wrong page. */
  page_texel *= 0.999999;
  /* Assumes atlas is squared. */
  return (vec2(page) + page_texel) / vec2(SHADOW_PAGE_PER_ROW);
}

/* Rotate vector to light's local space. Used for directional shadows. */
vec3 shadow_world_to_local(LightData ld, vec3 L)
{
  /* Avoid relying on compiler to optimize this.
   * vec3 lL = transpose(mat3(ld.object_mat)) * L; */
  vec3 lL;
  lL.x = dot(ld.object_mat[0].xyz, L);
  lL.y = dot(ld.object_mat[1].xyz, L);
  lL.z = dot(ld.object_mat[2].xyz, L);
  return lL;
}

/* ---------------------------------------------------------------------- */
/** \name Shadow Sampling Functions
 * \{ */

/* Turns local light coordinate into shadow region index. Matches eCubeFace order.
 * \note lL does not need to be normalized. */
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
float shadow_slope_bias_get(LightData light, vec3 lNg, uint lod)
{
#if 0 /* TODO(fclem): finish */
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
  float bias = length_manhattan(ndc_slope);
  /* Bias for 1 pixel of LOD 0. */
  bias *= 1.0 / (SHADOW_TILEMAP_RES * shadow_page_size_);
  /* Compensate for each increasing lod level as the space between pixels increases. */
  bias *= float(1u << lod);
  return bias;
#endif
  const float quantization_bias = 1e-20;
  return light.shadow_bias * float(1u << lod) + quantization_bias;
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

ShadowTileData shadow_directional_tile_get(usampler2D tilemaps_tx,
                                           LightData light,
                                           vec3 camera_P,
                                           vec3 lP,
                                           vec3 P,
                                           vec3 lNg,
                                           out vec2 uv,
                                           out float bias)
{
  ShadowClipmapCoordinates coord = shadow_directional_coordinates(
      light, lP, distance(camera_P, P));
  uv = coord.uv;

  ShadowTileData tile = shadow_tile_load(tilemaps_tx, coord.tile_coord, coord.tilemap_index);
  bias = shadow_slope_bias_get(light, lNg, tile.lod);
  return tile;
}

float shadow_tile_depth_get(sampler2D atlas_tx, ShadowTileData tile, vec2 uv)
{
  float depth = FLT_MAX;
  if (tile.is_allocated) {
    vec2 shadow_uv = shadow_page_uv_transform(tile.page, tile.lod, uv);
    depth = texture(atlas_tx, shadow_uv).r;
  }
  return depth;
}

struct ShadowSample {
  /* Signed delta in world units from the shading point to the occluder. Negative if occluded. */
  float occluder_delta;
  float bias;
};

/* TODO(fclem) use utildef version. */
float shadow_orderedIntBitsToFloat(int int_value)
{
  return intBitsToFloat((int_value < 0) ? (int_value ^ 0x7FFFFFFF) : int_value);
}

ShadowSample shadow_sample(sampler2D atlas_tx,
                           usampler2D tilemaps_tx,
                           LightData light,
                           vec3 lL,
                           vec3 lNg,
                           float receiver_dist,
                           vec3 P,
                           vec3 camera_P)
{
  ShadowSample samp;
  float occluder_dist;
  if (light.type == LIGHT_SUN) {
    vec3 lP = shadow_world_to_local(light, P);
    vec2 uv;
    ShadowTileData tile = shadow_directional_tile_get(
        tilemaps_tx, light, camera_P, lP, P, lNg, uv, samp.bias);

    occluder_dist = shadow_tile_depth_get(atlas_tx, tile, uv);
    /* Shadow is stored positive only for atomic operation.
     * So the encoded distance is positive and increasing from the near plane.
     * Bias back to get world distance. */
    occluder_dist = occluder_dist - shadow_orderedIntBitsToFloat(light.clip_near);
    /* Receiver distance needs to also be increasing.
     * Negate since Z distance follows opengl convention of neg Z as forward. */
    receiver_dist = -lP.z;
  }
  else {
    vec2 uv;
    ShadowTileData tile = shadow_punctual_tile_get(tilemaps_tx, light, lL, lNg, uv, samp.bias);
    occluder_dist = shadow_tile_depth_get(atlas_tx, tile, uv);
  }
  samp.occluder_delta = occluder_dist - receiver_dist;
  return samp;
}

/** \} */
