
#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_common_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_matcap_lib.glsl)
#pragma BLENDER_REQUIRE(workbench_world_light_lib.glsl)

void main()
{
  vec2 uv = vec2(gl_GlobalInvocationID.xy) / vec2(textureSize(normal_tx, 0));

  /* Normal and Incident vector are in viewspace. Lighting is evaluated in viewspace. */
  vec3 I = get_view_vector_from_screen_uv(uv);
  vec3 N = workbench_normal_decode(texture(normal_tx, uv));
  vec4 mat_data = texture(material_tx, uv);

  vec3 base_color = mat_data.rgb;

  float roughness, metallic;
  workbench_float_pair_decode(mat_data.a, roughness, metallic);

  vec4 color = vec4(0.0, 0.0, 0.0, 1.0);

#ifdef WORKBENCH_LIGHTING_MATCAP
  /* When using matcaps, mat_data.a is the back-face sign. */
  N = (mat_data.a > 0.0) ? N : -N;

  color.rgb = get_matcap_lighting(matcap_tx, base_color, N, I);
#endif

#ifdef WORKBENCH_LIGHTING_STUDIO
  color.rgb = get_world_lighting(base_color, roughness, metallic, N, I);
#endif

#ifdef WORKBENCH_LIGHTING_FLAT
  color.rgb = base_color;
#endif

  imageStore(out_color_img, ivec2(gl_GlobalInvocationID.xy), color);
}
