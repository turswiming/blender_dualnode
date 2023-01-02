
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
  /* We need to select the lod0 pixel closest to the lod pixel which is
   * located at the center of a lod0 pixel quad.
   * From the 4 possible pixels we choose the top left corner. */
  const int lod_center_offset = (lod > 0) ? (1 << (lod - 1)) : 0;
  if (!all(equal(((texel_co >> lod) << lod) + lod_center_offset, texel_co))) {
    return;
  }

  ivec3 render_map_coord = ivec3(tile_co >> lod, shadow_interp.view_id);
  uint page_packed = texelFetch(shadow_render_map_tx, render_map_coord, lod).r;
  if (page_packed != 0xFFFFFFFFu) {
    ivec2 page = ivec2(unpackUvec2x16(page_packed));
    ivec2 texel_in_page = (texel_co >> lod) % pages_infos_buf.page_size;
    ivec2 out_texel = page * pages_infos_buf.page_size + texel_in_page;

    imageAtomicMin(shadow_atlas_img, out_texel, floatBitsToUint(depth));
  }
}

void main()
{
  ivec2 texel_co = ivec2(gl_FragCoord.xy);
  ivec2 tile_co = texel_co / pages_infos_buf.page_size;
  float depth = gl_FragCoord.z;

  write_depth(texel_co, 0, tile_co, depth);

  /* We have to compensate the output pixel position being different than the input pixel's.
   * This is only half a pixel since we chose one pixel inside the quad. */
  depth -= DFDX_SIGN * 0.5 * dFdx(gl_FragCoord.z) + DFDY_SIGN * 0.5 * dFdy(gl_FragCoord.z);

  write_depth(texel_co, 1, tile_co, depth);
  write_depth(texel_co, 2, tile_co, depth);
  write_depth(texel_co, 3, tile_co, depth);
  write_depth(texel_co, 4, tile_co, depth);
}
