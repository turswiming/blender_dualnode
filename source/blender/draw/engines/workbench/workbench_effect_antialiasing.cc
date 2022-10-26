/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

#include "BLI_jitter_2d.h"
#include "smaa_textures.h"

namespace blender::workbench {

class TaaSamples {
  void init_samples(blender::Array<float2> &samples, const int size)
  {
    samples = blender::Array<float2>(size);
    BLI_jitter_init((float(*)[2])samples.begin(), size);

    /* find closest element to center */
    int closest_index = 0;
    float closest_squared_distance = 1.0f;

    for (int i : samples.index_range()) {
      float2 sample = samples[i];
      printf("%f : %f", sample.x, sample.y);
      const float squared_dist = len_squared_v2(sample);
      if (squared_dist < closest_squared_distance) {
        closest_squared_distance = squared_dist;
        closest_index = i;
      }
    }

    float2 closest_sample = samples[closest_index];

    for (float2 &sample : samples) {
      /* move jitter samples so that closest sample is in center */
      sample -= closest_sample;
      /* Avoid samples outside range (wrap around). */
      sample = {fmodf(sample.x + 0.5f, 1.0f), fmodf(sample.y + 0.5f, 1.0f)};
      /* Recenter the distribution[-1..1]. */
      sample = (sample * 2.0f) - 1.0f;
    }

    /* swap center sample to the start of the array */
    if (closest_index != 0) {
      swap_v2_v2(samples[0], samples[closest_index]);
    }

    /* Sort list based on farthest distance with previous. */
    for (int i = 0; i < size - 2; i++) {
      float squared_dist = 0.0;
      int index = i;
      for (int j = i + 1; j < size; j++) {
        const float _squared_dist = len_squared_v2(samples[i] - samples[j]);
        if (_squared_dist > squared_dist) {
          squared_dist = _squared_dist;
          index = j;
        }
      }
      swap_v2_v2(samples[i + 1], samples[index]);
    }
  }

 public:
  blender::Array<float2> x5;
  blender::Array<float2> x8;
  blender::Array<float2> x11;
  blender::Array<float2> x16;
  blender::Array<float2> x32;

  TaaSamples()
  {
    init_samples(x5, 5);
    init_samples(x8, 8);
    init_samples(x11, 11);
    init_samples(x16, 16);
    init_samples(x32, 32);
  }
};

static TaaSamples TAA_SAMPLES = TaaSamples();

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
