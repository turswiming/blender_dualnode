/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 *
 * Cavity Effect:
 *
 * We use Screen Space Ambient Occlusion (SSAO) to enhance geometric details of the surfaces.
 * We also use a Curvature effect computed only using the surface normals.
 *
 * This is done after the opaque pass. It only affects the opaque surfaces.
 */

#include "../eevee/eevee_lut.h" /* TODO: find somewhere to share blue noise Table. */
#include "BLI_rand.h"
#include "workbench_private.hh"

namespace blender::workbench {

void CavityEffect::init(const SceneState &scene_state, UniformBuffer<WorldData> &world_buf)
{
  cavity_enabled_ = scene_state.draw_cavity;
  curvature_enabled_ = scene_state.draw_curvature;

  const int ssao_samples = scene_state.scene->display.matcap_ssao_samples;
  int sample_count = min_ii(max_ii(1, scene_state.samples_len) * ssao_samples, max_samples_);
  const int max_iter_count = sample_count / ssao_samples;

  sample_ = scene_state.sample % max_iter_count;

  world_buf.cavity_sample_start = ssao_samples * sample_;
  world_buf.cavity_sample_end = ssao_samples * (sample_ + 1);

  world_buf.cavity_sample_count_inv = 1.0f / (world_buf.cavity_sample_end -
                                              world_buf.cavity_sample_start);
  world_buf.cavity_jitter_scale = 1.0f / 64.0f;

  world_buf.cavity_valley_factor = scene_state.shading.cavity_valley_factor;
  world_buf.cavity_ridge_factor = scene_state.shading.cavity_ridge_factor;
  world_buf.cavity_attenuation = scene_state.scene->display.matcap_ssao_attenuation;
  world_buf.cavity_distance = scene_state.scene->display.matcap_ssao_distance;

  world_buf.curvature_ridge = 0.5f /
                              max_ff(square_f(scene_state.shading.curvature_ridge_factor), 1e-4f);
  world_buf.curvature_valley = 0.7f / max_ff(square_f(scene_state.shading.curvature_valley_factor),
                                             1e-4f);

  if (cavity_enabled_ || scene_state.draw_dof) {
    setup_resources(ssao_samples, sample_count);
  }
}

void CavityEffect::setup_resources(int iteration_samples, int total_samples)
{
  if (sample_count_ != total_samples) {
    sample_count_ = total_samples;
    const float iteration_samples_inv = 1.0f / iteration_samples;

    /* Create disk samples using Hammersley distribution */
    for (int i : IndexRange(sample_count_)) {
      float it_add = (i / iteration_samples) * 0.499f;
      float r = fmodf((i + 0.5f + it_add) * iteration_samples_inv, 1.0f);
      double dphi;
      BLI_hammersley_1d(i, &dphi);

      float phi = (float)dphi * 2.0f * M_PI + it_add;
      samples_buf[i].x = cosf(phi);
      samples_buf[i].y = sinf(phi);
      /* This deliberately distribute more samples
       * at the center of the disk (and thus the shadow). */
      samples_buf[i].z = r;
    }

    samples_buf.push_update();

    const float total_samples_inv = 1.0f / iteration_samples;

    /* Create blue noise jitter texture */
    const int jitter_texel_count = jitter_tx_size_ * jitter_tx_size_;
    static float4 jitter[jitter_texel_count];
    for (int i = 0; i < jitter_texel_count; i++) {
      float phi = blue_noise[i][0] * 2.0f * M_PI;
      /* This rotate the sample per pixels */
      jitter[i].x = cosf(phi);
      jitter[i].y = sinf(phi);
      /* This offset the sample along its direction axis (reduce banding) */
      float bn = blue_noise[i][1] - 0.5f;
      bn = clamp_f(bn, -0.499f, 0.499f); /* fix fireflies */
      jitter[i].z = bn * total_samples_inv;
      jitter[i].w = blue_noise[i][1];
    }

    jitter_tx.ensure_2d(GPU_RGBA16F, int2(jitter_tx_size_), jitter[0]);
  }
}

void CavityEffect::setup_resolve_pass(PassSimple &pass, Texture &object_id_tx)
{
  if (cavity_enabled_) {
    pass.bind_ubo("cavity_samples", samples_buf);
    pass.bind_texture("jitter_tx", &jitter_tx, eGPUSamplerState::GPU_SAMPLER_REPEAT);
  }
  if (curvature_enabled_) {
    pass.bind_texture("object_id_tx", &object_id_tx);
  }
}

}  // namespace blender::workbench
