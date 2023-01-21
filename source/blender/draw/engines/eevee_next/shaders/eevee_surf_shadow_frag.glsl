
/**
 * Virtual Shadow map output.
 *
 * Meshes are rasterize onto an empty framebuffer. Each generated fragment then checks which
 * virtual page it is supposed to go and load the physical page adress.
 * If a physical page exists, we then use atomicMin to mimic a less-than depth test and write to
 * the destination texel.
 **/

#pragma BLENDER_REQUIRE(common_view_lib.glsl)
#pragma BLENDER_REQUIRE(common_math_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_attributes_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_surf_lib.glsl)
#pragma BLENDER_REQUIRE(eevee_nodetree_lib.glsl)

void write_depth(ivec2 texel_co, const int lod, ivec2 tile_co, float depth)
{
  ivec2 texel_co_lod = texel_co >> lod;
  if (!all(equal(texel_co_lod << lod, texel_co))) {
    return;
  }

  ivec3 render_map_coord = ivec3(tile_co >> lod, shadow_interp.view_id);
  uint page_packed = texelFetch(shadow_render_map_tx, render_map_coord, lod).r;
  /* Return if no valid page. */
  if (page_packed == 0xFFFFFFFFu) {
    return;
  }
  ivec2 page = ivec2(unpackUvec2x16(page_packed));
  ivec2 texel_in_page = texel_co_lod % pages_infos_buf.page_size;
  ivec2 out_texel = page * pages_infos_buf.page_size + texel_in_page;

  uint u_depth = floatBitsToUint(depth);
  /* Quantization bias. Equivalent to nextafter in C without all the safety. 1 is not enough. */
  u_depth += 2;

  imageAtomicMin(shadow_atlas_img, out_texel, u_depth);
}

void main()
{
  drw_view_id = shadow_interp.view_id;

  ivec2 texel_co = ivec2(gl_FragCoord.xy);
  ivec2 tile_co = texel_co / pages_infos_buf.page_size;

  float depth = gl_FragCoord.z;
  float slope_bias = 0.0;
  fwidth(depth);
  write_depth(texel_co, 0, tile_co, depth + slope_bias);

  /* Only needed for local lights. */
  bool is_persp = (drw_view.winmat[3][3] == 0.0);
  if (is_persp) {
    /* Note that even if texel center is offset, we store unmodified depth.
     * We increase bias instead at sampling time. */
    write_depth(texel_co, 1, tile_co, depth + slope_bias * 2.0);
    write_depth(texel_co, 2, tile_co, depth + slope_bias * 4.0);
    write_depth(texel_co, 3, tile_co, depth + slope_bias * 8.0);
    write_depth(texel_co, 4, tile_co, depth + slope_bias * 16.0);
  }
}
