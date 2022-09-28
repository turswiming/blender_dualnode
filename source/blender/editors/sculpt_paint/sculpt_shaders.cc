/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. */

/** \file
 * \ingroup draw_engine
 */

#include "GPU_shader.h"

#include "sculpt_intern.h"

static struct SCULPT_Shaders {
  GPUShader *paint_image_comp_sh;
} sh_data;

extern "C" {
GPUShader *SCULPT_shader_paint_image_get(void)
{
  if (sh_data.paint_image_comp_sh == nullptr) {
    sh_data.paint_image_comp_sh = GPU_shader_create_from_info_name("sculpt_paint_image_compute");
  }
  return sh_data.paint_image_comp_sh;
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
