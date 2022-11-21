/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "vk_shader.hh"

#include "vk_backend.hh"

#include "BLI_string_utils.h"
#include "BLI_vector.hh"

namespace blender::gpu {

static std::string combine_sources(MutableSpan<const char *> sources)
{
  char *sources_combined = BLI_string_join_arrayN((const char **)sources.data(), sources.size());
  return std::string(sources_combined);
}

Vector<uint32_t> VKShader::compile_glsl_to_spirv(StringRef source, shaderc_shader_kind kind)
{
  VKBackend &backend = static_cast<VKBackend &>(*VKBackend::get());
  shaderc::Compiler &compiler = backend.get_shaderc_compiler();
  shaderc::CompileOptions options;

  shaderc::SpvCompilationResult module = compiler.CompileGlslToSpv(source, kind, name, options);

  if (module.GetCompilationStatus() != shaderc_compilation_status_success) {
    // TODO(jbakker): error handling.
  }

  return Vector<uint32_t>(module.cbegin(), module.cend());
}

void VKShader::build_shader_module(Span<uint32_t> spirv_module, VkShaderModule *r_shader_module)
{
  VkShaderModuleCreateInfo createInfo = {};
  createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
  createInfo.codeSize = spirv_module.size() * sizeof(uint32_t);
  createInfo.pCode = spirv_module.data();

  // TODO(jbakker): retrieve allocator and device
  // VkResult result = vkCreateShaderModule(*device, &create_info, nullptr, r_shader_module);
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
  std::string source = combine_sources(sources);

  Vector<uint32_t> spirv_module = compile_glsl_to_spirv(StringRef(source), shaderc_compute_shader);
  build_shader_module(spirv_module, &compute_module_);
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