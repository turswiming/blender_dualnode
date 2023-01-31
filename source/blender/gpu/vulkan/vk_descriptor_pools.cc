/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup gpu
 */

#include "vk_descriptor_pools.hh"

namespace blender::gpu {
VKDescriptorPools::VKDescriptorPools()
{
}

VKDescriptorPools::~VKDescriptorPools()
{
}

void VKDescriptorPools::new_pool(VkDevice vk_device)
{
  Vector<VkDescriptorPoolSize> pool_sizes = {
      {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, POOL_SIZE_UNIFORM_BUFFER},
  };
  VkDescriptorPoolCreateInfo pool_info = {};
  pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
  pool_info.flags = 0;
  pool_info.maxSets = POOL_SIZE_DESCRIPTOR_SETS;
  pool_info.poolSizeCount = pool_sizes.size();
  pool_info.pPoolSizes = pool_sizes.data();
  VkDescriptorPool descriptor_pool = VK_NULL_HANDLE;
  VkResult result = vkCreateDescriptorPool(vk_device, &pool_info, nullptr, &descriptor_pool);
  UNUSED_VARS(result);
  pools_.append(descriptor_pool);
}

}  // namespace blender::gpu