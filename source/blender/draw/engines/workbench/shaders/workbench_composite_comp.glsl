
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_common_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_matcap_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_world_light_lib.glsl)

void main()
{
  ivec2 texel = ivec2(gl_GlobalInvocationID.xy);
  vec2 uv = (vec2(texel) + 0.5) / vec2(textureSize(normal_tx, 0));
  /* Normal and Incident vector are in viewspace. Lighting is evaluated in viewspace. */
  vec3 V = get_view_vector_from_screen_uv(uv);
  vec3 N = workbench_normal_decode(texture(normal_tx, uv));
  vec4 mat_data = texture(material_tx, uv);
  float depth = texture(depth_tx, uv).r;
  depth = min(depth, texture(depth_in_front_tx, uv).r);

  vec3 base_color = mat_data.rgb;
  vec4 color = world_data.background_color;

  /* Background pixels. */
  if (depth != 1.0) {
#ifdef WORKBENCH_LIGHTING_MATCAP
    /* When using matcaps, mat_data.a is the back-face sign. */
    N = (mat_data.a > 0.0) ? N : -N;
    color.rgb = get_matcap_lighting(matcap_tx, base_color, N, V);
#endif

#ifdef WORKBENCH_LIGHTING_STUDIO
    float roughness, metallic;
    workbench_float_pair_decode(mat_data.a, roughness, metallic);
    color.rgb = get_world_lighting(base_color, roughness, metallic, N, V);
#endif

#ifdef WORKBENCH_LIGHTING_FLAT
    color.rgb = base_color;
#endif
  }

  /* TODO(fclem): Port the TAA shader that does this tranformation. */
  /* Use log2 space to avoid highlights creating too much aliasing. */
  color = log2(color + 0.5);

  imageStore(out_color_img, texel, color);
}
