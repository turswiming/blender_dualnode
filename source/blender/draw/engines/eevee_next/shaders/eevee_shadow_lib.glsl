
#pragma BLENDER_REQUIRE(eevee_shadow_tilemap_lib.glsl)

/** \a unormalized_uv is the uv coordinates for the whole tilemap [0..SHADOW_TILEMAP_RES]. */
vec2 shadow_page_uv_transform(uvec2 page, uint lod, vec2 unormalized_uv, ivec2 tile_coord)
{
  /* TODO(fclem): It should be possible to just saturate(unormalized_uv - tile_co << lod). */
  vec2 page_texel = unormalized_uv / float(1u << lod);
  /* Fix float imprecision that can make some pixel sample the wrong page. */
  page_texel = saturate(page_texel - vec2(tile_coord >> lod));
  /* Assumes atlas is squared. */
  return (vec2(page) + page_texel) / vec2(SHADOW_PAGE_PER_ROW);
}

/* Rotate vector to light's local space and apply location. Used for directional shadows. */
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

/* TODO(fclem) use utildef version. */
float shadow_orderedIntBitsToFloat(int int_value)
{
  return intBitsToFloat((int_value < 0) ? (int_value ^ 0x7FFFFFFF) : int_value);
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

mat4x4 shadow_load_normal_matrix(LightData light)
{
  if (light.type != LIGHT_SUN) {
    return mat4x4(vec4(1.0, 0.0, 0.0, 0.0),
                  vec4(0.0, 1.0, 0.0, 0.0),
                  vec4(0.0, 0.0, 0.0, -1.0),
                  vec4(0.0, 0.0, light.normal_mat_packed.x, light.normal_mat_packed.y));
  }
  else {
    return mat4x4(vec4(light.normal_mat_packed.x, 0.0, 0.0, 0.0),
                  vec4(0.0, light.normal_mat_packed.x, 0.0, 0.0),
                  vec4(0.0, 0.0, 1.0, 0.0),
                  vec4(0.0, 0.0, 0.0, 1.0));
  }
}

/* Returns minimum bias (in world space unit) needed for a given geometry normal and a shadowmap
 * page to avoid self shadowing artifacts. */
float shadow_slope_bias_get(LightData light, vec3 lNg, vec3 lP, uint lod)
{
  /* Create a normal plane equation and go through the normal projection matrix. */
  vec4 lNg_plane = vec4(lNg, -dot(lNg, lP));
  vec4 ndc_Ng = shadow_load_normal_matrix(light) * lNg_plane;
  /* Get slope from normal vector. */
  vec2 ndc_slope = ndc_Ng.xy / ndc_Ng.z;
  /* Slope bias definition from fixed pipeline. */
  float bias = abs(ndc_slope.x) + abs(ndc_slope.y);
  /* Clamp out to avoid the bias going to infinity. Remember this is in NDC space. */
  bias = clamp(bias, 0.0, 100.0);
  /* Bias for 1 pixel of LOD 0. */
  bias *= 1.0 / (SHADOW_TILEMAP_RES * SHADOW_PAGE_RES);
  /* Compensate for each increasing lod level as the space between pixels increases. */
  bias *= float(1u << lod);

  return bias;
}

struct ShadowTileSample {
  /* Signed delta in world units from the shading point to the occluder. Negative if occluded. */
  ShadowTileData tile;
  /* Tile coordinate inside the tilemap [0..SHADOW_TILEMAP_RES). */
  ivec2 tile_coord;
  /* UV coordinate inside the tilemap [0..SHADOW_TILEMAP_RES). */
  vec2 uv;
  /* Minimum slope bias to apply during comparison. */
  float bias;
};

ShadowTileSample shadow_punctual_tile_get(usampler2D tilemaps_tx,
                                          LightData light,
                                          vec3 lP,
                                          vec3 lNg)
{
  ShadowTileSample samp;

  int face_id = shadow_punctual_face_index_get(lP);
  lP = shadow_punctual_local_position_to_face_local(face_id, lP);
  lNg = shadow_punctual_local_position_to_face_local(face_id, lNg);
  /* UVs in [-1..+1] range. */
  samp.uv = lP.xy / abs(lP.z);
  /* UVs in [0..SHADOW_TILEMAP_RES] range. */
  const float lod0_res = float(SHADOW_TILEMAP_RES / 2);
  samp.uv = samp.uv * lod0_res + lod0_res;
  samp.tile_coord = ivec2(floor(samp.uv));

  int tilemap_index = light.tilemap_index + face_id;

  samp.tile = shadow_tile_load(tilemaps_tx, samp.tile_coord, tilemap_index);
  samp.bias = shadow_slope_bias_get(light, lNg, lP, samp.tile.lod);
  /* Bias uv sample for LODs since custom raster aligns LOD0 and LOD1 pixels
   * instead of centering them. */
  if (samp.tile.lod != 0.0) {
    samp.uv += 0.5 / float(SHADOW_TILEMAP_RES * SHADOW_PAGE_RES);
  }
  return samp;
}

ShadowTileSample shadow_directional_tile_get(usampler2D tilemaps_tx,
                                             LightData light,
                                             vec3 lP,
                                             vec3 lNg)
{
  ShadowTileSample samp;

  ShadowClipmapCoordinates coord = shadow_directional_coordinates(light, lP);
  samp.uv = coord.uv;
  samp.tile_coord = coord.tile_coord;

  samp.tile = shadow_tile_load(tilemaps_tx, coord.tile_coord, coord.tilemap_index);
  samp.bias = shadow_slope_bias_get(light, lNg, lP, coord.clipmap_lod_relative + samp.tile.lod);
  return samp;
}

float shadow_tile_depth_get(sampler2D atlas_tx, ShadowTileSample samp)
{
  /* Far plane distance but with a bias to make sure there will be no shadowing.
   * But also not FLT_MAX since it can cause issue with projection. */
  float depth = 1.1;

  if (samp.tile.is_allocated) {
    vec2 shadow_uv = shadow_page_uv_transform(
        samp.tile.page, samp.tile.lod, samp.uv, samp.tile_coord);
    depth = texture(atlas_tx, shadow_uv).r;
  }
  return depth;
}

struct ShadowSample {
  /* Signed delta in world units from the shading point to the occluder. Negative if occluded. */
  float occluder_delta;
  float bias;
};

float shadow_punctual_linear_depth(float z, float near, float far)
{
  return (near * far) / (z * (near - far) + far);
}

float shadow_directional_linear_depth(float z, float near, float far)
{
  return z * (near - far) - near;
}

ShadowSample shadow_sample(sampler2D atlas_tx,
                           usampler2D tilemaps_tx,
                           LightData light,
                           vec3 lL,
                           vec3 lNg,
                           float receiver_dist,
                           vec3 P)
{
  ShadowSample samp;
  float occluder_dist;
  if (light.type == LIGHT_SUN) {
    vec3 lP = shadow_world_to_local(light, P);
    ShadowTileSample tile = shadow_directional_tile_get(tilemaps_tx, light, lP, lNg);

    float occluder_ndc = shadow_tile_depth_get(atlas_tx, tile);

    float near = shadow_orderedIntBitsToFloat(light.clip_near);
    float far = shadow_orderedIntBitsToFloat(light.clip_far);
    occluder_dist = shadow_directional_linear_depth(occluder_ndc, near, far);
    /* Receiver distance needs to also be increasing.
     * Negate since Z distance follows opengl convention of neg Z as forward. */
    receiver_dist = -lP.z;
    samp.bias = 0.0;
    // tile.bias *(near - far);
  }
  else {
    vec3 lP = lL;
    ShadowTileSample tile = shadow_punctual_tile_get(tilemaps_tx, light, lP, lNg);
    float occluder_ndc = shadow_tile_depth_get(atlas_tx, tile);

    float near = shadow_orderedIntBitsToFloat(light.clip_near);
    float far = shadow_orderedIntBitsToFloat(light.clip_far);
    /* Shadow is stored as gl_FragCoord.z. Convert to radial distance along with the bias. */
    float occluder_z = shadow_punctual_linear_depth(occluder_ndc, near, far);
    float occluder_z_bias = shadow_punctual_linear_depth(occluder_ndc + tile.bias, near, far);
    float radius_divisor = receiver_dist / max_v3(abs(lL));
    occluder_dist = occluder_z * radius_divisor;
    samp.bias = (occluder_z_bias - occluder_z) * radius_divisor;
  }
  samp.occluder_delta = occluder_dist - receiver_dist;
  return samp;
}

/** \} */