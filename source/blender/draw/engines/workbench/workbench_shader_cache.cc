/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

namespace blender::workbench {

ShaderCache::~ShaderCache()
{
  for (auto i : IndexRange(shading_type_len)) {
    for (auto j : IndexRange(color_type_len)) {
      for (auto k : IndexRange(geometry_type_len)) {
        for (auto l : IndexRange(pipeline_type_len)) {
          DRW_SHADER_FREE_SAFE(prepass_shader_cache_[i][j][k][l]);
        }
      }
    }
  }
  for (auto i : IndexRange(shading_type_len)) {
    for (auto j : IndexRange(pipeline_type_len)) {
      DRW_SHADER_FREE_SAFE(resolve_shader_cache_[i][j]);
    }
  }
}

GPUShader *ShaderCache::prepass_shader_get(ePipelineType pipeline_type,
                                           eGeometryType geometry_type,
                                           eColorType color_type,
                                           eShadingType shading_type)
{
  GPUShader *&shader_ptr = prepass_shader_cache_[static_cast<int>(pipeline_type)][static_cast<int>(
      geometry_type)][static_cast<int>(color_type)][static_cast<int>(shading_type)];

  if (shader_ptr != nullptr) {
    return shader_ptr;
  }
  std::string info_name = "workbench_next_prepass_";
  switch (geometry_type) {
    case eGeometryType::MESH:
      info_name += "mesh_";
      break;
    case eGeometryType::CURVES:
      info_name += "curves_";
      break;
    case eGeometryType::POINTCLOUD:
      info_name += "ptcloud_";
      break;
  }
  switch (pipeline_type) {
    case ePipelineType::OPAQUE:
      info_name += "opaque_";
      break;
    case ePipelineType::TRANSPARENT:
      info_name += "transparent_";
      break;
    case ePipelineType::SHADOW:
      info_name += "shadow_";
      break;
  }
  switch (shading_type) {
    case eShadingType::FLAT:
      info_name += "flat_";
      break;
    case eShadingType::STUDIO:
      info_name += "studio_";
      break;
    case eShadingType::MATCAP:
      info_name += "matcap_";
      break;
  }
  switch (color_type) {
    case eColorType::MATERIAL:
      info_name += "material";
      break;
    case eColorType::TEXTURE:
      info_name += "texture";
      break;
  }
  /* TODO Clipping */
  info_name += "_no_clip";
  shader_ptr = GPU_shader_create_from_info_name(info_name.c_str());
  prepass_shader_cache_[static_cast<int>(pipeline_type)][static_cast<int>(geometry_type)]
                       [static_cast<int>(color_type)][static_cast<int>(shading_type)] = shader_ptr;
  return shader_ptr;
}

GPUShader *ShaderCache::resolve_shader_get(ePipelineType pipeline_type, eShadingType shading_type)
{
  GPUShader *&shader_ptr =
      resolve_shader_cache_[static_cast<int>(shading_type)][static_cast<int>(pipeline_type)];

  if (shader_ptr != nullptr) {
    return shader_ptr;
  }
  std::string info_name = "workbench_next_resolve_";
  switch (pipeline_type) {
    case ePipelineType::OPAQUE:
      info_name += "opaque_";
      break;
    case ePipelineType::TRANSPARENT:
      info_name += "transparent_";
      break;
    case ePipelineType::SHADOW:
      BLI_assert_unreachable();
      break;
  }
  switch (shading_type) {
    case eShadingType::FLAT:
      info_name += "flat";
      break;
    case eShadingType::STUDIO:
      info_name += "studio";
      break;
    case eShadingType::MATCAP:
      info_name += "matcap";
      break;
  }
  shader_ptr = GPU_shader_create_from_info_name(info_name.c_str());
  return shader_ptr;
}

}  // namespace blender::workbench
