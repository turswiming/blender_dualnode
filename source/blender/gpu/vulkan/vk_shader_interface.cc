/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "vk_shader_interface.hh"

namespace blender::gpu {

void VKShaderInterface::init(const shader::ShaderCreateInfo &info)
{
  using namespace blender::gpu::shader;

  ssbo_len_ = 0;

  Vector<ShaderCreateInfo::Resource> all_resources;
  all_resources.extend(info.pass_resources_);
  all_resources.extend(info.batch_resources_);

  for (ShaderCreateInfo::Resource &res : all_resources) {
    switch (res.bind_type) {
      case ShaderCreateInfo::Resource::BindType::IMAGE:
      case ShaderCreateInfo::Resource::BindType::SAMPLER:
      case ShaderCreateInfo::Resource::BindType::UNIFORM_BUFFER:
        // BLI_assert_msg(false, "not implemented yet");
        break;
      case ShaderCreateInfo::Resource::BindType::STORAGE_BUFFER:
        ssbo_len_++;
        break;
    }
  }

  int32_t input_tot_len = ssbo_len_;
  inputs_ = static_cast<ShaderInput *>(
      MEM_calloc_arrayN(input_tot_len, sizeof(ShaderInput), __func__));
  ShaderInput *input = inputs_;

  name_buffer_ = (char *)MEM_mallocN(info.interface_names_size_, "name_buffer");
  uint32_t name_buffer_offset = 0;

  for (const ShaderCreateInfo::Resource &res : all_resources) {
    if (res.bind_type == ShaderCreateInfo::Resource::BindType::STORAGE_BUFFER) {
      copy_input_name(input, res.storagebuf.name, name_buffer_, name_buffer_offset);
      input->location = input->binding = res.slot;
      enabled_ssbo_mask_ |= (1 << input->binding);
      input++;
    }
  }

  sort_inputs();
}

}  // namespace blender::gpu
