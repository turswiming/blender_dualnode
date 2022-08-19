/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "draw_manager.hh"
#include "draw_shader.h"
#include "draw_testing.hh"

namespace blender::draw {

static void test_draw_pass_all_commands()
{
  Manager drw;

  Texture tex;
  tex.ensure_2d(GPU_RGBA16, int2(1));

  UniformBuffer<uint4> ubo;
  ubo.push_update();

  StorageBuffer<uint4> ssbo;
  ssbo.push_update();

  float alpha = 0.0f;
  int3 dispatch_size(1);

  Pass pass = {"test.all_commands"};
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
  pass.draw(GPU_PRIM_TRIS, 1, 3);

  /* Should not crash even if shader is not a compute. This is because we only serialize. */
  /* TODO(fclem): Use real compute shader. */
  pass.shader_set(GPU_shader_get_builtin_shader(GPU_SHADER_3D_IMAGE_MODULATE_ALPHA));
  pass.dispatch(dispatch_size);
  pass.dispatch(&dispatch_size);
  pass.barrier(GPU_BARRIER_SHADER_IMAGE_ACCESS);

  /* Change references. */
  alpha = 1.0f;
  dispatch_size = int3(2);

  std::string result = drw.serialize(pass);
  StringRefNull expected = "";

  EXPECT_EQ(result, expected);
}
DRAW_TEST(draw_pass_all_commands)

}  // namespace blender::draw