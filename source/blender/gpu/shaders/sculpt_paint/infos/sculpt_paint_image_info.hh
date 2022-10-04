/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "gpu_shader_create_info.hh"

GPU_SHADER_CREATE_INFO(sculpt_paint_image_compute)
    .local_group_size(1, 1, 1)
    .image(0, GPU_RGBA32F, Qualifier::READ_WRITE, ImageType::FLOAT_2D, "out_img")
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
    .image(0, GPU_RGBA32F, Qualifier::READ, ImageType::FLOAT_2D, "in_paint_img")
    .image(1, GPU_RGBA16F, Qualifier::READ_WRITE, ImageType::FLOAT_2D, "out_img")
    .compute_source("sculpt_paint_image_merge_comp.glsl")
    .typedef_source("GPU_sculpt_shader_shared.h")
    .do_static_compilation(true);

/* -------------------------------------------------------------------- */
/** \name Brush variations
 * \{ */

GPU_SHADER_CREATE_INFO(sculpt_paint_test_sphere).define("BRUSH_TEST_SPHERE");
GPU_SHADER_CREATE_INFO(sculpt_paint_test_circle).define("BRUSH_TEST_CIRCLE");

GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_custom).define("BRUSH_CURVE_PRESET 0");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_smooth).define("BRUSH_CURVE_PRESET 1");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_sphere).define("BRUSH_CURVE_PRESET 2");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_root).define("BRUSH_CURVE_PRESET 3");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_sharp).define("BRUSH_CURVE_PRESET 4");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_lin).define("BRUSH_CURVE_PRESET 5");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_pow4).define("BRUSH_CURVE_PRESET 6");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_invsquare).define("BRUSH_CURVE_PRESET 7");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_constant).define("BRUSH_CURVE_PRESET 8");
GPU_SHADER_CREATE_INFO(sculpt_paint_falloff_smoother).define("BRUSH_CURVE_PRESET 9");

#define SCULPT_PAINT_FINAL_VARIATION(name, ...) \
  GPU_SHADER_CREATE_INFO(name).additional_info(__VA_ARGS__).do_static_compilation(true);

#define SCULPT_PAINT_CURVE_VARIATION(name, ...) \
  SCULPT_PAINT_FINAL_VARIATION(name##_custom, "sculpt_paint_falloff_custom", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_smooth, "sculpt_paint_falloff_smooth", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_sphere, "sculpt_paint_falloff_sphere", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_root, "sculpt_paint_falloff_root", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_sharp, "sculpt_paint_falloff_sharp", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_lin, "sculpt_paint_falloff_lin", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_pow4, "sculpt_paint_falloff_pow4", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_invsquare, "sculpt_paint_falloff_invsquare", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_constant, "sculpt_paint_falloff_constant", __VA_ARGS__) \
  SCULPT_PAINT_FINAL_VARIATION(name##_smoother, "sculpt_paint_falloff_smoother", __VA_ARGS__)

#define SCULPT_PAINT_TEST_VARIATIONS(name, ...) \
  SCULPT_PAINT_CURVE_VARIATION(name##_sphere, "sculpt_paint_test_sphere", __VA_ARGS__) \
  SCULPT_PAINT_CURVE_VARIATION(name##_circle, "sculpt_paint_test_circle", __VA_ARGS__)

SCULPT_PAINT_TEST_VARIATIONS(sculpt_paint_image, "sculpt_paint_image_compute")

/** \} */
