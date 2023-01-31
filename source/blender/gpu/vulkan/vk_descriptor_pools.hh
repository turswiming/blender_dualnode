/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#pragma once

#include "BLI_vector.hh"

#ifdef __APPLE__
#  include <MoltenVK/vk_mvk_moltenvk.h>
#else
#  include <vulkan/vulkan.h>
#endif

namespace blender::gpu {

/**
 * List of VkDescriptorPools.
 *
 * In Vulkan a pool is constructed with a certain size. When more is needed it is adviced to
 * construct a second pool. VKDescriptorPools will keep track of those pools and construct
 * new pools when the previous one is exhausted.
 */
class VKDescriptorPools {
  /** Number of pool sizes */
  static constexpr uint32_t POOL_SIZE_UNIFORM_BUFFER = 1000;
  static constexpr uint32_t POOL_SIZE_DESCRIPTOR_SETS = 1000;

  Vector<VkDescriptorPool> pools_;

 public:
  VKDescriptorPools();
  ~VKDescriptorPools();

 private:
  void new_pool(VkDevice vk_device);
};
}  // namespace blender::gpu