/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#pragma once

#include "gpu_shader_create_info.hh"
#include "gpu_shader_interface.hh"

namespace blender::gpu {
class VKShaderInterface : public ShaderInterface {
 public:
  VKShaderInterface() = default;

  void init(const shader::ShaderCreateInfo &info);
};
}  // namespace blender::gpu
