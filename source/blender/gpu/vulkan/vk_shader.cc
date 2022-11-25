/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "vk_shader.hh"

#include "vk_backend.hh"
#include "vk_shader_log.hh"

#include "BLI_string_utils.h"
#include "BLI_vector.hh"

namespace blender::gpu {

static const std::string to_stage_name(shaderc_shader_kind stage)
{
  switch (stage) {
    case shaderc_compute_shader:
      return std::string("compute");

    default:
      BLI_assert_msg(false, "Do not know how to convert shaderc_shader_kind to stage name.");
      break;
  }
  return std::string("unknown stage");
}

static std::string combine_sources(Span<const char *> sources)
{
  char *sources_combined = BLI_string_join_arrayN((const char **)sources.data(), sources.size());
  return std::string(sources_combined);
}

static char *glsl_patch_get()
{
  static char patch[512] = "\0";
  if (patch[0] != '\0') {
    return patch;
  }

  size_t slen = 0;
  /* Version need to go first. */
  STR_CONCAT(patch, slen, "#version 430\n");

  BLI_assert(slen < sizeof(patch));
  return patch;
}

Vector<uint32_t> VKShader::compile_glsl_to_spirv(Span<const char *> sources,
                                                 shaderc_shader_kind stage)
{
  std::string combined_sources = combine_sources(sources);
  VKBackend &backend = static_cast<VKBackend &>(*VKBackend::get());
  shaderc::Compiler &compiler = backend.get_shaderc_compiler();
  shaderc::CompileOptions options;

  shaderc::SpvCompilationResult module = compiler.CompileGlslToSpv(
      combined_sources, stage, name, options);
  if (module.GetNumErrors() != 0 || module.GetNumWarnings() != 0) {
    std::string log = module.GetErrorMessage();
    Vector<char> logcstr(log.c_str(), log.c_str() + log.size() + 1);

    VKLogParser parser;
    print_log(sources,
              logcstr.data(),
              to_stage_name(stage).c_str(),
              module.GetCompilationStatus() != shaderc_compilation_status_success,
              &parser);
  }

  if (module.GetCompilationStatus() != shaderc_compilation_status_success) {
    return Vector<uint32_t>();
  }

  return Vector<uint32_t>(module.cbegin(), module.cend());
}

void VKShader::build_shader_module(Span<uint32_t> spirv_module, VkShaderModule *r_shader_module)
{
  VkShaderModuleCreateInfo create_info = {};
  create_info.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
  create_info.codeSize = spirv_module.size() * sizeof(uint32_t);
  create_info.pCode = spirv_module.data();

  VKContext &context = *static_cast<VKContext *>(VKContext::get());

  VkResult result = vkCreateShaderModule(
      context.device_get(), &create_info, nullptr, r_shader_module);
  if (result != VK_SUCCESS) {
    *r_shader_module = VK_NULL_HANDLE;
  }
}

VKShader::VKShader(const char *name) : Shader(name)
{
  context_ = VKContext::get();
}
VKShader::~VKShader()
{
  VkDevice device = context_->device_get();
  if (compute_module_ != VK_NULL_HANDLE) {
    vkDestroyShaderModule(device, compute_module_, nullptr);
    compute_module_ = VK_NULL_HANDLE;
  }
}

void VKShader::build_shader_module(MutableSpan<const char *> sources,
                                   shaderc_shader_kind stage,
                                   VkShaderModule *r_shader_module)
{
  BLI_assert_msg(ELEM(stage, shaderc_compute_shader),
                 "Only forced ShaderC shader kinds are supported.");
  sources[0] = glsl_patch_get();
  Vector<uint32_t> spirv_module = compile_glsl_to_spirv(sources, shaderc_compute_shader);
  build_shader_module(spirv_module, &compute_module_);
}

void VKShader::vertex_shader_from_glsl(MutableSpan<const char *> /*sources*/)
{
}

void VKShader::geometry_shader_from_glsl(MutableSpan<const char *> /*sources*/)
{
}

void VKShader::fragment_shader_from_glsl(MutableSpan<const char *> /*sources*/)
{
}

void VKShader::compute_shader_from_glsl(MutableSpan<const char *> sources)
{
  build_shader_module(sources, shaderc_compute_shader, &compute_module_);
}

bool VKShader::finalize(const shader::ShaderCreateInfo * /*info*/)
{
  return false;
}

void VKShader::transform_feedback_names_set(Span<const char *> /*name_list*/,
                                            eGPUShaderTFBType /*geom_type*/)
{
}

bool VKShader::transform_feedback_enable(GPUVertBuf *)
{
  return false;
}

void VKShader::transform_feedback_disable()
{
}

void VKShader::bind()
{
}

void VKShader::unbind()
{
}

void VKShader::uniform_float(int /*location*/,
                             int /*comp_len*/,
                             int /*array_size*/,
                             const float * /*data*/)
{
}
void VKShader::uniform_int(int /*location*/,
                           int /*comp_len*/,
                           int /*array_size*/,
                           const int * /*data*/)
{
}

std::string VKShader::resources_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

std::string VKShader::vertex_interface_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

std::string VKShader::fragment_interface_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

std::string VKShader::geometry_interface_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

std::string VKShader::geometry_layout_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

std::string VKShader::compute_layout_declare(const shader::ShaderCreateInfo & /*info*/) const
{
  return std::string();
}

int VKShader::program_handle_get() const
{
  return -1;
}

}  // namespace blender::gpu