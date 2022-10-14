/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "gpu_shader_create_info.hh"
#include "workbench_defines.hh"

/* -------------------------------------------------------------------- */
/** \name Base Composite
 * \{ */

GPU_SHADER_CREATE_INFO(workbench_composite)
    .sampler(0, ImageType::FLOAT_2D, "normalBuffer", Frequency::PASS)
    .sampler(1, ImageType::FLOAT_2D, "materialBuffer", Frequency::PASS)
    .uniform_buf(WB_WORLD_SLOT, "WorldData", "world_data", Frequency::PASS)
    .push_constant(Type::BOOL, "forceShadowing")
    .fragment_out(0, Type::VEC4, "fragColor")
    .typedef_source("workbench_shader_shared.h")
    .fragment_source("workbench_composite_frag.glsl")
    .additional_info("draw_fullscreen", "draw_view");

GPU_SHADER_CREATE_INFO(workbench_next_composite)
    .local_group_size(8, 8)
    .sampler(3, ImageType::FLOAT_2D, "normal_tx")
    .sampler(4, ImageType::FLOAT_2D, "material_tx")
    .sampler(5, ImageType::DEPTH_2D, "depth_tx")
    .uniform_buf(WB_WORLD_SLOT, "WorldData", "world_data")
    .push_constant(Type::BOOL, "forceShadowing")
    .image(0, GPU_RGBA16F, Qualifier::WRITE, ImageType::FLOAT_2D, "out_color_img")
    .typedef_source("workbench_shader_shared.h")
    .compute_source("workbench_composite_comp.glsl")
    .additional_info("draw_view");

/** \} */

/* -------------------------------------------------------------------- */
/** \name Lighting Type
 * \{ */

GPU_SHADER_CREATE_INFO(workbench_composite_studio)
    .define("WORKBENCH_LIGHTING_STUDIO")
    .additional_info("workbench_composite")
    .do_static_compilation(true);

GPU_SHADER_CREATE_INFO(workbench_composite_matcap)
    .define("WORKBENCH_LIGHTING_MATCAP")
    .sampler(2, ImageType::FLOAT_2D, "matcap_diffuse_tx", Frequency::PASS)
    .sampler(3, ImageType::FLOAT_2D, "matcap_specular_tx", Frequency::PASS)
    .additional_info("workbench_composite")
    .do_static_compilation(true);

GPU_SHADER_CREATE_INFO(workbench_composite_flat)
    .define("WORKBENCH_LIGHTING_FLAT")
    .additional_info("workbench_composite")
    .do_static_compilation(true);

/** \} */

/* -------------------------------------------------------------------- */
/** \name Lighting Type
 * \{ */

GPU_SHADER_CREATE_INFO(workbench_next_resolve_opaque_studio)
    .define("WORKBENCH_LIGHTING_STUDIO")
    .additional_info("workbench_next_composite")
    .do_static_compilation(true);

GPU_SHADER_CREATE_INFO(workbench_next_resolve_opaque_matcap)
    .define("WORKBENCH_LIGHTING_MATCAP")
    .sampler(WB_MATCAP_SLOT, ImageType::FLOAT_2D_ARRAY, "matcap_tx")
    .additional_info("workbench_next_composite")
    .do_static_compilation(true);

GPU_SHADER_CREATE_INFO(workbench_next_resolve_opaque_flat)
    .define("WORKBENCH_LIGHTING_FLAT")
    .additional_info("workbench_next_composite")
    .do_static_compilation(true);

/** \} */
