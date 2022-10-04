/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 */

#include "GPU_shader.h"

#include "sculpt_intern.h"

#include <sstream>

static struct SCULPT_Shaders {
  GPUShader *paint_image_comp_sh[BRUSH_MAX_VARIATIONS];
  GPUShader *paint_image_merge_comp_sh;
} sh_data;

extern "C" {
GPUShader *SCULPT_shader_paint_image_get(BrushVariationFlags variation_flags)
{
  int index = static_cast<int>(variation_flags);

  if (sh_data.paint_image_comp_sh[index] == nullptr) {

    std::stringstream info_name;
    info_name << "sculpt_paint_image";
    info_name << (variation_flags & BRUSH_TEST_CIRCLE ? "_circle" : "_sphere");
    BrushVariationFlags curve_falloff = static_cast<BrushVariationFlags>(
        variation_flags & BRUSH_VARIATION_FALLOFF_MASK);
    switch (curve_falloff) {
      case BRUSH_VARIATION_FALLOFF_CUSTOM:
        info_name << "_custom";
        break;
      case BRUSH_VARIATION_FALLOFF_SMOOTH:
        info_name << "_smooth";
        break;
      case BRUSH_VARIATION_FALLOFF_SPHERE:
        info_name << "_sphere";
        break;
      case BRUSH_VARIATION_FALLOFF_ROOT:
        info_name << "_root";
        break;
      case BRUSH_VARIATION_FALLOFF_SHARP:
        info_name << "_sharp";
        break;
      case BRUSH_VARIATION_FALLOFF_LIN:
        info_name << "_lin";
        break;
      case BRUSH_VARIATION_FALLOFF_POW4:
        info_name << "_pow4";
        break;
      case BRUSH_VARIATION_FALLOFF_INVSQUARE:
        info_name << "_invsquare";
        break;
      case BRUSH_VARIATION_FALLOFF_SMOOTHER:
        info_name << "_smoother";
        break;
      default:
        BLI_assert_unreachable();
    }

    printf("%s create shader %s\n", __func__, info_name.str().c_str());

    sh_data.paint_image_comp_sh[index] = GPU_shader_create_from_info_name(info_name.str().c_str());
  }
  return sh_data.paint_image_comp_sh[index];
}

GPUShader *SCULPT_shader_paint_image_merge_get(void)
{
  if (sh_data.paint_image_merge_comp_sh == nullptr) {
    sh_data.paint_image_merge_comp_sh = GPU_shader_create_from_info_name(
        "sculpt_paint_image_merge_compute");
  }
  return sh_data.paint_image_merge_comp_sh;
}

void SCULPT_shader_free(void)
{
  GPUShader **sh_data_as_array = (GPUShader **)&sh_data;
  for (int i = 0; i < (sizeof(SCULPT_Shaders) / sizeof(GPUShader *)); i++) {
    if (sh_data_as_array[i]) {
      GPU_shader_free(sh_data_as_array[i]);
      sh_data_as_array[i] = nullptr;
    }
  }
}
}
