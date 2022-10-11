/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#pragma once

#include "draw_manager.hh"

#include "overlay_grid.hh"

namespace blender::draw::overlay {

class ShaderCache {
  Map<StringRefNull, std::array<GPUShader *, 2>> cache;

  int clipping_enabled = 0;
};

class SceneResources {
  ShaderCache shaders;

  // UniformBuffer<ThemeColorData> theme_colors;
  // Texture color_ramp = {"color_ramp"};

  void weight_ramp_init()
  {
    /* Weight Painting color ramp texture */
    // bool user_weight_ramp = (U.flag & USER_CUSTOM_RANGE) != 0;

    // if (weight_ramp_custom != user_weight_ramp ||
    //     (user_weight_ramp && memcmp(&weight_ramp_copy, &U.coba_weight, sizeof(ColorBand)) != 0))
    //     {
    //   DRW_TEXTURE_FREE_SAFE(G_draw.weight_ramp);
    // }

    // if (G_draw.weight_ramp == NULL) {
    //   weight_ramp_custom = user_weight_ramp;
    //   memcpy(&weight_ramp_copy, &U.coba_weight, sizeof(ColorBand));

    //   G_draw.weight_ramp = DRW_create_weight_colorramp_texture();
    // }
  }
};

class Instance {
 public:
  ShaderCache shaders;

  /* WORKAROUND: Legacy. Move to grid pass. */
  GPUUniformBuf *grid_ubo = nullptr;

  Framebuffer overlay_fb = {"overlay_fb"};

  Grid grid;

  ~Instance()
  {
    DRW_UBO_FREE_SAFE(grid_ubo);
  }

  void init();
  void begin_sync();
  void object_sync(ObjectRef &ob_ref);
  void end_sync();
  void draw(Manager &manager);
};

}  // namespace blender::draw::overlay
