/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

#include "smaa_textures.h"

namespace blender::workbench {

AntiAliasingPass::AntiAliasingPass()
{
  smaa_edge_detect_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_0");
  smaa_aa_weight_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_1");
  smaa_resolve_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_2");

  smaa_search_tx.ensure_2d(GPU_R8, {SEARCHTEX_WIDTH, SEARCHTEX_HEIGHT});
  GPU_texture_update(smaa_search_tx, GPU_DATA_UBYTE, searchTexBytes);
  GPU_texture_filter_mode(smaa_search_tx, true);

  smaa_area_tx.ensure_2d(GPU_RG8, {AREATEX_WIDTH, AREATEX_HEIGHT});
  GPU_texture_update(smaa_area_tx, GPU_DATA_UBYTE, areaTexBytes);
  GPU_texture_filter_mode(smaa_area_tx, true);
}

AntiAliasingPass::~AntiAliasingPass()
{
  if (smaa_edge_detect_sh) {
    GPU_shader_free(smaa_edge_detect_sh);
  }
  if (smaa_aa_weight_sh) {
    GPU_shader_free(smaa_aa_weight_sh);
  }
  if (smaa_resolve_sh) {
    GPU_shader_free(smaa_resolve_sh);
  }
}

void AntiAliasingPass::init(bool reset_taa)
{
  is_playback = DRW_state_is_playback();
  is_navigating = DRW_state_is_navigating();

  if (reset_taa || is_playback || is_navigating) {
    taa_sample = 0;
  }
}

void AntiAliasingPass::sync(SceneResources &resources)
{
  {
    smaa_edge_detect_ps_.init();
    smaa_edge_detect_ps_.state_set(DRW_STATE_WRITE_COLOR);
    smaa_edge_detect_ps_.shader_set(smaa_edge_detect_sh);
    smaa_edge_detect_ps_.bind_texture("colorTex", &resources.color_tx);
    smaa_edge_detect_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
    smaa_edge_detect_ps_.clear_color(float4(0.0f));
    smaa_edge_detect_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
  }
  {
    smaa_aa_weight_ps_.init();
    smaa_aa_weight_ps_.state_set(DRW_STATE_WRITE_COLOR);
    smaa_aa_weight_ps_.shader_set(smaa_aa_weight_sh);
    smaa_aa_weight_ps_.bind_texture("edgesTex", &smaa_edge_tx);
    smaa_aa_weight_ps_.bind_texture("areaTex", smaa_area_tx);
    smaa_aa_weight_ps_.bind_texture("searchTex", smaa_search_tx);
    smaa_aa_weight_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
    smaa_aa_weight_ps_.clear_color(float4(0.0f));
    smaa_aa_weight_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
  }
  {
    smaa_resolve_ps_.init();
    smaa_resolve_ps_.state_set(DRW_STATE_WRITE_COLOR);
    smaa_resolve_ps_.shader_set(smaa_resolve_sh);
    smaa_resolve_ps_.bind_texture("blendTex", &smaa_weight_tx);
    smaa_resolve_ps_.bind_texture("colorTex", &resources.color_tx);
    smaa_resolve_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
    smaa_resolve_ps_.push_constant("mixFactor", &smaa_mix_factor, 1);
    smaa_resolve_ps_.push_constant("taaAccumulatedWeight", &taa_weight_accum, 1);
    smaa_resolve_ps_.clear_color(float4(0.0f));
    smaa_resolve_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
  }
}

void AntiAliasingPass::draw(Manager &manager,
                            View &view,
                            GPUTexture *depth_tx,
                            GPUTexture *color_tx)
{
  int2 size = {GPU_texture_width(depth_tx), GPU_texture_height(depth_tx)};

  taa_weight_accum = 1.0f; /* TODO */

  smaa_viewport_metrics = float4(1.0f / size.x, 1.0f / size.y, size.x, size.y);
  smaa_mix_factor = 1.0f; /* TODO */

  smaa_edge_tx.acquire(size, GPU_RG8);
  smaa_edge_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(smaa_edge_tx));
  smaa_edge_fb.bind();
  manager.submit(smaa_edge_detect_ps_, view);

  smaa_weight_tx.acquire(size, GPU_RGBA8);
  smaa_weight_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(smaa_weight_tx));
  smaa_weight_fb.bind();
  manager.submit(smaa_aa_weight_ps_, view);
  smaa_edge_tx.release();

  smaa_resolve_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(color_tx));
  smaa_resolve_fb.bind();
  manager.submit(smaa_resolve_ps_, view);
  smaa_weight_tx.release();
}

}  // namespace blender::workbench
