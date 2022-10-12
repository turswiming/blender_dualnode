/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#pragma once

#include "draw_manager.hh"

#include "overlay_background.hh"
#include "overlay_grid.hh"
#include "overlay_metaball.hh"

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

  /** Global types. */
  Resources resources;
  State state;

  /** Overlay types. */
  Background background;
  Metaballs metaballs;
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

 private:
  bool object_is_edit_mode(const Object *ob)
  {
    if (DRW_object_is_in_edit_mode(ob)) {
      /* Also check for context mode as the object mode is not 100% reliable. (see T72490) */
      switch (ob->type) {
        case OB_MESH:
          return state.ctx_mode == CTX_MODE_EDIT_MESH;
        case OB_ARMATURE:
          return state.ctx_mode == CTX_MODE_EDIT_ARMATURE;
        case OB_CURVES_LEGACY:
          return state.ctx_mode == CTX_MODE_EDIT_CURVE;
        case OB_SURF:
          return state.ctx_mode == CTX_MODE_EDIT_SURFACE;
        case OB_LATTICE:
          return state.ctx_mode == CTX_MODE_EDIT_LATTICE;
        case OB_MBALL:
          return state.ctx_mode == CTX_MODE_EDIT_METABALL;
        case OB_FONT:
          return state.ctx_mode == CTX_MODE_EDIT_TEXT;
        case OB_CURVES:
          return state.ctx_mode == CTX_MODE_EDIT_CURVES;
        case OB_POINTCLOUD:
        case OB_VOLUME:
          /* No edit mode yet. */
          return false;
      }
    }
    return false;
  }
};

}  // namespace blender::draw::overlay
