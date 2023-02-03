/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2023 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#pragma once

#include "gpu_context_private.hh"
#include "vk_context.hh"

#ifdef __APPLE__
#  include <MoltenVK/vk_mvk_moltenvk.h>
#else
#  include <vulkan/vulkan.h>
#endif

namespace blender::gpu {

/**
 * Class for handing vulkan buffers (allocation/updating/binding).
 */
class VKBuffer {
  int64_t size_in_bytes_;
  VkBuffer vk_buffer_ = VK_NULL_HANDLE;
  VmaAllocation allocation_ = VK_NULL_HANDLE;

 public:
  VKBuffer() = default;
  virtual ~VKBuffer();

  /** Has this buffer been allocated? */
  bool is_allocated();

  bool create(VKContext &context,
              int64_t size,
              GPUUsageType usage,
              VkBufferUsageFlagBits buffer_usage);
  bool update(VKContext &context, const void *data);
  bool free(VKContext &context);
  bool map(VKContext &context, void **r_mapped_memory);
  void unmap(VKContext &context);
};
}  // namespace blender::gpu
