/* Directive for resetting the line numbering so the failing tests lines can be printed.
 * This conflict with the shader compiler error logging scheme.
 * Comment out for correct compilation error line. */
#line 5

#pragma BLENDER_REQUIRE(gpu_shader_math_vector_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_shadow_lib.glsl)
#pragma BLENDER_REQUIRE(gpu_shader_test_lib.glsl)

#define TEST(a, b) if (true)

void main()
{
  TEST(eevee_shadow, DirectionalClipmapLevel)
  {
    LightData light;
    light.clipmap_lod_min = -5;
    light.clipmap_lod_max = 8;
    EXPECT_EQ(shadow_directional_clipmap_level(light, 0.0), light.clipmap_lod_min);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 0.49), 0);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 0.5), 0);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 0.51), 1);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 0.99), 1);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 1.0), 1);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 1.01), 2);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 12.5), 5);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 12.51), 5);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 15.9999), 5);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 16.0), 5);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 16.00001), 6);
    EXPECT_EQ(shadow_directional_clipmap_level(light, 5000.0), light.clipmap_lod_max);
    /* Produces NaN / Inf, Undefined behavior. */
    // EXPECT_EQ(shadow_directional_clipmap_level(light, FLT_MAX), light.clipmap_lod_max);
  }

  TEST(eevee_shadow, DirectionalClipmapCoordinates)
  {
    ShadowClipmapCoordinates coords;
    vec3 lP, camera_lP;

    LightData light;
    light.clipmap_lod_min = 0; /* Range [-0.5..0.5]. */
    light.clipmap_lod_max = 2; /* Range [-2..2]. */
    light.tilemap_index = light.clipmap_lod_min;
    light.tilemap_last = light.clipmap_lod_max;
    float lod_min_tile_size = pow(2.0, float(light.clipmap_lod_min)) / float(SHADOW_TILEMAP_RES);
    float lod_max_half_size = pow(2.0, float(light.clipmap_lod_max)) / 2.0;
    light._clipmap_scale = float(SHADOW_TILEMAP_RES / 2) / lod_max_half_size;

    camera_lP = vec3(0.0, 0.0, 0.0);
    /* Follows ShadowDirectional::end_sync(). */
    light.clipmap_base_offset = ivec2(round(camera_lP.xy / lod_min_tile_size));
    EXPECT_EQ(light.clipmap_base_offset, ivec2(0));

    /* Test UVs and tile mapping. */

    lP = vec3(1e-5, 1e-5, 0.0);
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 0);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES / 2), 1e-3);

    lP = vec3(-1e-5, -1e-5, 0.0);
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 0);
    EXPECT_EQ(coords.tile_coord, ivec2((SHADOW_TILEMAP_RES / 2) - 1));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES / 2), 1e-3);

    lP = vec3(-0.5, -0.5, 0.0); /* Min of first LOD. */
    coords = shadow_directional_coordinates(light, lP, 0.0 /* Force LOD 0. */);
    EXPECT_EQ(coords.tilemap_index, 0);
    EXPECT_EQ(coords.tile_coord, ivec2(0));
    EXPECT_NEAR(coords.uv, vec2(0), 1e-3);

    lP = vec3(0.5, 0.5, 0.0); /* Max of first LOD. */
    coords = shadow_directional_coordinates(light, lP, 0.0 /* Force LOD 0. */);
    EXPECT_EQ(coords.tilemap_index, 0);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES - 1));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES), 1e-3);

    /* Test clipmap level selection. */

    camera_lP = vec3(2.0, 2.0, 0.0);
    /* Follows ShadowDirectional::end_sync(). */
    light.clipmap_base_offset = ivec2(round(camera_lP.xy / lod_min_tile_size));
    EXPECT_EQ(light.clipmap_base_offset, ivec2(32));

    lP = vec3(2.00001, 2.00001, 0.0);
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 0);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES / 2), 1e-3);

    lP = vec3(1.50001, 1.50001, 0.0);
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 1);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 4));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES / 4), 1e-3);

    lP = vec3(1.00001, 1.00001, 0.0);
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 2);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 4));
    EXPECT_NEAR(coords.uv, vec2(SHADOW_TILEMAP_RES / 4), 1e-3);

    lP = vec3(-0.0001, -0.0001, 0.0); /* Out of bounds. */
    coords = shadow_directional_coordinates(light, lP, distance(lP, camera_lP));
    EXPECT_EQ(coords.tilemap_index, 2);
    EXPECT_EQ(coords.tile_coord, ivec2(0));
    EXPECT_NEAR(coords.uv, vec2(0), 1e-3);

    /* Test clipmap offset. */

    light.clipmap_base_offset = ivec2(31, 1);
    lP = vec3(2.0001, 0.0001, 0.0);

    coords = shadow_directional_coordinates(light, lP, 0.0 /* Force LOD 0. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2) + ivec2(1, -1));

    coords = shadow_directional_coordinates(light, lP, 2.0 /* Force LOD 1. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2) + ivec2(1, 0));

    coords = shadow_directional_coordinates(light, lP, 4.0 /* Force LOD 2. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2) + ivec2(1, 0));

    coords = shadow_directional_coordinates(light, lP, 8.0 /* Force LOD 3. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2) + ivec2(1, 0));

    /* Test clipmap negative offsets. */

    light.clipmap_base_offset = ivec2(-31, -1);
    lP = vec3(-2.0001, -0.0001, 0.0);

    coords = shadow_directional_coordinates(light, lP, 0.0 /* Force LOD 0. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2 - 1) + ivec2(-1, 1));

    coords = shadow_directional_coordinates(light, lP, 2.0 /* Force LOD 1. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2 - 1) + ivec2(-1, 0));

    coords = shadow_directional_coordinates(light, lP, 4.0 /* Force LOD 2. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2 - 1) + ivec2(-1, 0));

    coords = shadow_directional_coordinates(light, lP, 8.0 /* Force LOD 3. */);
    EXPECT_EQ(coords.tile_coord, ivec2(SHADOW_TILEMAP_RES / 2 - 1) + ivec2(-1, 0));
  }
}