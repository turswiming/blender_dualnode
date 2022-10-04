/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "gpu_shader_create_info.hh"

GPU_SHADER_CREATE_INFO(sculpt_paint_image_compute)
    .local_group_size(1, 1, 1)
    .image(0, GPU_RGBA16F, Qualifier::READ_WRITE, ImageType::FLOAT_2D, "out_img")
    .storage_buf(0, Qualifier::READ, "PackedPixelRow", "pixel_row_buf[]")
    .storage_buf(1, Qualifier::READ, "TrianglePaintInput", "paint_input[]")
    .storage_buf(2, Qualifier::READ, "vec3", "vert_coord_buf[]")
    .storage_buf(3, Qualifier::READ, "PaintStepData", "paint_step_buf[]")
    .uniform_buf(0, "PaintBrushData", "paint_brush_buf")
    .push_constant(Type::INT, "pixel_row_offset")
    .push_constant(Type::IVEC2, "paint_step_range")
    .compute_source("sculpt_paint_image_comp.glsl")
    .typedef_source("GPU_sculpt_shader_shared.h");

GPU_SHADER_CREATE_INFO(sculpt_paint_image_merge_compute)
    .local_group_size(1, 1, 1)
    .image(0, GPU_RGBA16F, Qualifier::READ, ImageType::FLOAT_2D, "in_paint_img")
    .image(1, GPU_RGBA16F, Qualifier::READ_WRITE, ImageType::FLOAT_2D, "out_img")
    .compute_source("sculpt_paint_image_merge_comp.glsl")
    .typedef_source("GPU_sculpt_shader_shared.h")
    .do_static_compilation(true);

/* -------------------------------------------------------------------- */
/** \name Brush variations
 * \{ */

GPU_SHADER_CREATE_INFO(sculpt_paint_test_sphere).define("BRUSH_TEST_SPHERE");
GPU_SHADER_CREATE_INFO(sculpt_paint_test_circle).define("BRUSH_TEST_CIRCLE");

#define SCULPT_PAINT_FINAL_VARIATION(name, ...) \
  GPU_SHADER_CREATE_INFO(name).additional_info(__VA_ARGS__).do_static_compilation(true);

#define SCULPT_PAINT_TEST_VARIATIONS(name, ...) \
  SCULPT_PAINT_FINAL_VARIATION(name##_sphere, "sculpt_paint_test_sphere", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_circle, "sculpt_paint_test_circle", __VA_ARGS__)

SCULPT_PAINT_TEST_VARIATIONS(sculpt_paint_image, "sculpt_paint_image_compute")

/** \} */
