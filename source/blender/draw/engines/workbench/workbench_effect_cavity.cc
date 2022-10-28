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

void CavityEffect::init(const DrawConfig &config, UniformBuffer<WorldData> &world_buf)
{
  cavity_enabled = config.draw_cavity;
  curvature_enabled = config.draw_curvature;

  const int ssao_samples = config.scene->display.matcap_ssao_samples;
  const int sample_count = min_ii(max_ii(1, config.aa_samples) * ssao_samples, MAX_SAMPLES);
  const int max_iter_count = max_ii(1, sample_count / ssao_samples);

  if (config.reset_taa) {
    sample = 0;
  }
  sample = sample++ % max_iter_count;

  world_buf.cavity_sample_start = ssao_samples * sample;
  world_buf.cavity_sample_end = ssao_samples * (sample + 1);

  world_buf.cavity_sample_count_inv = 1.0f / (world_buf.cavity_sample_end -
                                              world_buf.cavity_sample_start);
  world_buf.cavity_jitter_scale = 1.0f / 64.0f;

  world_buf.cavity_valley_factor = config.shading.cavity_valley_factor;
  world_buf.cavity_ridge_factor = config.shading.cavity_ridge_factor;
  world_buf.cavity_attenuation = config.scene->display.matcap_ssao_attenuation;
  world_buf.cavity_distance = config.scene->display.matcap_ssao_distance;

  world_buf.curvature_ridge = 0.5f /
                              max_ff(square_f(config.shading.curvature_ridge_factor), 1e-4f);
  world_buf.curvature_valley = 0.7f /
                               max_ff(square_f(config.shading.curvature_valley_factor), 1e-4f);

  if (cavity_enabled || config.draw_dof) {
    setup_resources(sample_count);
  }
}

void CavityEffect::setup_resources(int sample_count)
{
  if (this->sample_count != sample_count) {
    this->sample_count = sample_count;
    const float sample_count_inv = 1.0f / sample_count;

    /* Create disk samples using Hammersley distribution */
    for (int i : IndexRange(sample_count)) {
      float it_add = (i / sample_count) * 0.499f;
      float r = fmodf((i + 0.5f + it_add) * sample_count_inv, 1.0f);
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

    /* Create blue noise jitter texture */
    const int jitter_texel_count = JITTER_TEX_SIZE * JITTER_TEX_SIZE;
    static float4 jitter[jitter_texel_count];
    for (int i = 0; i < jitter_texel_count; i++) {
      float phi = blue_noise[i][0] * 2.0f * M_PI;
      /* This rotate the sample per pixels */
      jitter[i].x = cosf(phi);
      jitter[i].y = sinf(phi);
      /* This offset the sample along its direction axis (reduce banding) */
      float bn = blue_noise[i][1] - 0.5f;
      bn = clamp_f(bn, -0.499f, 0.499f); /* fix fireflies */
      jitter[i].z = bn * sample_count_inv;
      jitter[i].w = blue_noise[i][1];
    }

    jitter_tx.ensure_2d(GPU_RGBA16F, int2(JITTER_TEX_SIZE), jitter[0]);
  }
}

void CavityEffect::setup_resolve_pass(PassSimple &pass, Texture &object_id_tx)
{
  if (cavity_enabled) {
    pass.bind_ubo("cavity_samples", samples_buf);
    pass.bind_texture("jitter_tx", &jitter_tx, eGPUSamplerState::GPU_SAMPLER_REPEAT);
  }
  if (curvature_enabled) {
    pass.bind_texture("object_id_tx", &object_id_tx);
  }
}

}  // namespace blender::workbench
