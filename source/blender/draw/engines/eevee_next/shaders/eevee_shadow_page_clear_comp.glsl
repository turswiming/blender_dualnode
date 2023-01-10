
/**
 * Virtual shadowmapping: Page Clear.
 *
 * Equivalent to a framebuffer depth clear but only for pages pushed to the clear_page_buf.
 */

#pragma BLENDER_REQUIRE(common_math_lib.glsl)

void main()
{
  if (gl_GlobalInvocationID == uvec3(0)) {
    drw_print(pages_infos_buf.page_used_count);
    drw_print(pages_infos_buf.page_update_count);
    drw_print(pages_infos_buf.page_allocated_count);
    drw_print(pages_infos_buf.page_rendered_count);
    drw_print(pages_infos_buf.page_cached_count);
  }

  uvec2 page_co = unpackUvec2x16(clear_page_buf[gl_GlobalInvocationID.z]);
  uvec2 page_texel = page_co * pages_infos_buf.page_size + gl_GlobalInvocationID.xy;

  imageStore(atlas_img, ivec2(page_texel), uvec4(floatBitsToUint(FLT_MAX)));
}