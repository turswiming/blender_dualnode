/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "GPU_storage_buffer.h"

#include "BLI_vector.hh"

#include "gpu_testing.hh"

namespace blender::gpu::tests {

constexpr size_t SIZE = 128;
constexpr size_t SIZE_IN_BYTES = SIZE * sizeof(int);

static void test_gpu_storage_buffer_create_update_read()
{
  GPUStorageBuf *ssbo = GPU_storagebuf_create_ex(
      SIZE_IN_BYTES, nullptr, GPU_USAGE_STATIC, __func__);
  EXPECT_NE(ssbo, nullptr);

  /* Upload some dummy data. */
  int data[SIZE];
  for (int i : IndexRange(SIZE)) {
    data[i] = i;
  }
  GPU_storagebuf_update(ssbo, data);

  GPU_storagebuf_free(ssbo);
}

GPU_TEST(gpu_storage_buffer_create_update_read);

}  // namespace blender::gpu::tests