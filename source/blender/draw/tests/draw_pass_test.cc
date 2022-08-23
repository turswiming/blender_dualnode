/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "draw_manager.hh"
#include "draw_shader.h"
#include "draw_testing.hh"

namespace blender::draw {

static void test_draw_pass_all_commands()
{
  Texture tex;
  tex.ensure_2d(GPU_RGBA16, int2(1));

  UniformBuffer<uint4> ubo;
  ubo.push_update();

  StorageBuffer<uint4> ssbo;
  ssbo.push_update();

  float alpha = 0.0f;
  int3 dispatch_size(1);

  PassSimple pass = {"test.all_commands"};
  pass.init();
  pass.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_STENCIL);
  pass.clear_color_depth_stencil(float4(0.25f, 0.5f, 100.0f, -2000.0f), 0.5f, 0xF0);
  pass.state_stencil(0x80, 0x0F, 0x8F);
  pass.shader_set(GPU_shader_get_builtin_shader(GPU_SHADER_3D_IMAGE_MODULATE_ALPHA));
  pass.bind("image", tex);
  pass.bind("image", &tex);
  pass.bind("missing_image", as_image(tex));  /* Should not crash. */
  pass.bind("missing_image", as_image(&tex)); /* Should not crash. */
  pass.bind("missing_ubo", ubo);              /* Should not crash. */
  pass.bind("missing_ubo", &ubo);             /* Should not crash. */
  pass.bind("missing_ssbo", ssbo);            /* Should not crash. */
  pass.bind("missing_ssbo", &ssbo);           /* Should not crash. */
  pass.push_constant("alpha", alpha);
  pass.push_constant("alpha", &alpha);
  pass.push_constant("ModelViewProjectionMatrix", float4x4::identity());
  pass.draw_procedural(GPU_PRIM_TRIS, 1, 3);

  /* Should not crash even if shader is not a compute. This is because we only serialize. */
  /* TODO(fclem): Use real compute shader. */
  pass.shader_set(GPU_shader_get_builtin_shader(GPU_SHADER_3D_IMAGE_MODULATE_ALPHA));
  pass.dispatch(dispatch_size);
  pass.dispatch(&dispatch_size);
  pass.barrier(GPU_BARRIER_SHADER_IMAGE_ACCESS);

  /* Change references. */
  alpha = 1.0f;
  dispatch_size = int3(2);

  std::string result = pass.serialize();
  std::stringstream expected;
  expected << ".test.all_commands" << std::endl;
  expected << "  .state_set(6)" << std::endl;
  expected << "  .clear(color=(0.25, 0.5, 100, -2000), depth=0.5, stencil=0b11110000))"
           << std::endl;
  expected << "  .stencil_set(write_mask=0b10000000, compare_mask=0b00001111, reference=0b10001111"
           << std::endl;
  expected << "  .shader_bind(gpu_shader_3D_image_modulate_alpha)" << std::endl;
  expected << "  .bind_texture(0)" << std::endl;
  expected << "  .bind_texture_ref(0)" << std::endl;
  expected << "  .bind_image(-1)" << std::endl;
  expected << "  .bind_image_ref(-1)" << std::endl;
  expected << "  .bind_uniform_buf(-1)" << std::endl;
  expected << "  .bind_uniform_buf_ref(-1)" << std::endl;
  expected << "  .bind_storage_buf(-1)" << std::endl;
  expected << "  .bind_storage_buf_ref(-1)" << std::endl;
  expected << "  .push_constant(2, data=0)" << std::endl;
  expected << "  .push_constant(2, data=1)" << std::endl;
  expected << "  .push_constant(0, data=(" << std::endl;
  expected << "(   1.000000,    0.000000,    0.000000,    0.000000)" << std::endl;
  expected << "(   0.000000,    1.000000,    0.000000,    0.000000)" << std::endl;
  expected << "(   0.000000,    0.000000,    1.000000,    0.000000)" << std::endl;
  expected << "(   0.000000,    0.000000,    0.000000,    1.000000)" << std::endl;
  expected << ")" << std::endl;
  expected << ")" << std::endl;
  expected << "  .draw(inst_len=1, vert_len=3, vert_first=0, res_id=0)" << std::endl;
  expected << "  .shader_bind(gpu_shader_3D_image_modulate_alpha)" << std::endl;
  expected << "  .dispatch(1, 1, 1)" << std::endl;
  expected << "  .dispatch_ref(2, 2, 2)" << std::endl;
  expected << "  .barrier(4)" << std::endl;

  EXPECT_EQ(result, expected.str());

  DRW_shape_cache_free();
}
DRAW_TEST(draw_pass_all_commands)

static void test_draw_pass_sub_ordering()
{
  PassSimple pass = {"test.sub_ordering"};
  pass.init();
  pass.shader_set(GPU_shader_get_builtin_shader(GPU_SHADER_3D_IMAGE_MODULATE_ALPHA));
  pass.push_constant("test_pass", 1);

  PassSimple::Sub &sub1 = pass.sub("Sub1");
  sub1.push_constant("test_sub1", 11);

  PassSimple::Sub &sub2 = pass.sub("Sub2");
  sub2.push_constant("test_sub2", 21);

  /* Will execute after both sub. */
  pass.push_constant("test_pass", 2);

  /* Will execute after sub1. */
  sub2.push_constant("test_sub2", 22);

  /* Will execute before sub2. */
  sub1.push_constant("test_sub1", 12);

  /* Will execute before end of pass. */
  sub2.push_constant("test_sub2", 23);

  std::string result = pass.serialize();
  std::stringstream expected;
  expected << ".test.sub_ordering" << std::endl;
  expected << "  .shader_bind(gpu_shader_3D_image_modulate_alpha)" << std::endl;
  expected << "  .push_constant(-1, data=1)" << std::endl;
  expected << "  .Sub1" << std::endl;
  expected << "    .push_constant(-1, data=11)" << std::endl;
  expected << "  .Sub2" << std::endl;
  expected << "    .push_constant(-1, data=21)" << std::endl;
  expected << "    .push_constant(-1, data=22)" << std::endl;
  expected << "    .push_constant(-1, data=23)" << std::endl;
  expected << "  .push_constant(-1, data=2)" << std::endl;

  EXPECT_EQ(result, expected.str());
}
DRAW_TEST(draw_pass_sub_ordering)

}  // namespace blender::draw