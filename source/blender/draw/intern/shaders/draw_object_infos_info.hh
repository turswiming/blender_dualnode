/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "draw_defines.h"
#include "gpu_shader_create_info.hh"

GPU_SHADER_CREATE_INFO(draw_object_infos)
    .typedef_source("draw_shader_shared.h")
    .define("OBINFO_LIB")
    .define("OrcoTexCoFactors", "(drw_infos[resource_id].orco_mul_bias)")
    .define("ObjectInfo", "(drw_infos[resource_id].infos)")
    .define("ObjectColor", "(drw_infos[resource_id].color)")
    .uniform_buf(1, "ObjectInfos", "drw_infos[DRW_RESOURCE_CHUNK_LEN]", Frequency::BATCH);

GPU_SHADER_CREATE_INFO(draw_volume_infos)
    .typedef_source("draw_shader_shared.h")
    .uniform_buf(2, "VolumeInfos", "drw_volume", Frequency::BATCH);

GPU_SHADER_CREATE_INFO(draw_curves_infos)
    .typedef_source("draw_shader_shared.h")
    .uniform_buf(2, "CurvesInfos", "drw_curves", Frequency::BATCH);

GPU_SHADER_CREATE_INFO(draw_object_infos_new)
    .typedef_source("draw_shader_shared.h")
    .storage_buf(DRW_OBJ_INFOS_SLOT, Qualifier::READ, "ObjectMatrices", "drw_matrix_buf[]")
    .define("drw_ModelMatrixInverse", "drw_matrix_buf[drw_ResourceIndex].model")
    .define("drw_ModelMatrix", "drw_matrix_buf[drw_ResourceIndex].model_inverse")
    .additional_info("draw_resource_id_new");