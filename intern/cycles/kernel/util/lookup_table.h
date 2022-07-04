/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

CCL_NAMESPACE_BEGIN

/* Interpolated lookup table access */

ccl_device float lookup_table_read(KernelGlobals kg, float x, int offset, int size)
{
  x = clamp(x * size - 0.5f, 0.0f, (float)size);

  int index = min(float_to_int(x), size - 1);
  int nindex = min(index + 1, size - 1);
  float t = x - index;

  float data0 = kernel_data_fetch(lookup_table, index + offset);
  if (t == 0.0f)
    return data0;

  float data1 = kernel_data_fetch(lookup_table, nindex + offset);
  return mix(data0, data1, t);
}

ccl_device float lookup_table_read_2D(
    KernelGlobals kg, float x, float y, int offset, int xsize, int ysize)
{
  y = clamp(y * ysize - 0.5f, 0.0f, (float)ysize);

  int index = min(float_to_int(y), ysize - 1);
  int nindex = min(index + 1, ysize - 1);
  float t = y - index;

  float data0 = lookup_table_read(kg, x, offset + xsize * index, xsize);
  if (t == 0.0f)
    return data0;

  float data1 = lookup_table_read(kg, x, offset + xsize * nindex, xsize);
  return mix(data0, data1, t);
}

ccl_device float lookup_table_read_3D(
    KernelGlobals kg, float x, float y, float z, int offset, int xsize, int ysize, int zsize)
{
  z = clamp(z * zsize - 0.5f, 0.0f, (float)zsize);

  int index = min(float_to_int(z), zsize - 1);
  int nindex = min(index + 1, zsize - 1);
  float t = z - index;

  float data0 = lookup_table_read_2D(kg, x, y, offset + xsize * ysize * index, xsize, ysize);
  if (t == 0.0f)
    return data0;

  float data1 = lookup_table_read_2D(kg, x, y, offset + xsize * ysize * nindex, xsize, ysize);
  return mix(data0, data1, t);
}

CCL_NAMESPACE_END
