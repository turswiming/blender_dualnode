/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#include "BLI_math_vec_types.hh"
#include "BLI_vector.hh"
#include "BLI_math.h"

#include "IMB_imbuf.h"
#include "IMB_imbuf_types.h"

namespace blender::kernel::pbvh::pixels {
template<typename T, int Channels = 4> struct ImageBufferAccessor {
  ImBuf &image_buffer;

  ImageBufferAccessor(ImBuf &image_buffer) : image_buffer(image_buffer)
  {
  }

  float4 read_pixel(const int2 coordinate)
  {
    int offset = (coordinate.y * image_buffer.x + coordinate.x) * Channels;
    return float4(&image_buffer.rect_float[offset]);
  }

  void write_pixel(const int2 coordinate, float4 new_value)
  {
    int offset = (coordinate.y * image_buffer.x + coordinate.x) * Channels;
    copy_v4_v4(&image_buffer.rect_float[offset], new_value);
  }
};

struct PixelCopyItem {
  uint8_t delta_destination_x;
  char2 delta_source_1;
  char2 delta_source_2;
  uint8_t mix_factor;
};

struct PixelCopyGroup {
  int2 destination;
  Vector<PixelCopyItem> items;
};

/** Pixel copy command to mix 2 source pixels and write to a destination pixel. */
struct PixelCopyCommand {
  /** Pixel coordinate to write to. */
  int2 destination;
  /** Pixel coordinate to read first source from. */
  int2 source_1;
  /** Pixel coordinate to read second source from. */
  int2 source_2;
  /** Factor to mix between first and second source. */
  float mix_factor;

  PixelCopyCommand(const PixelCopyGroup &group)
      : destination(group.destination),
        source_1(group.destination),
        source_2(group.destination),
        mix_factor(0.0f)
  {
  }

  template<typename T> void mix_source_and_write_destination(ImageBufferAccessor<T> &tile_buffer) const
  {
    float4 source_color_1 = tile_buffer.read_pixel(source_1);
    float4 source_color_2 = tile_buffer.read_pixel(source_2);
    float4 destination_color = source_color_1 * (1.0f - mix_factor) + source_color_2 * mix_factor;
    tile_buffer.write_pixel(destination, destination_color);
  }
};

PixelCopyCommand &operator+=(PixelCopyCommand &command, const PixelCopyItem &item)
{
  command.destination.x += int(item.delta_destination_x);
  command.source_1 += int2(item.delta_source_1);
  command.source_2 += int2(item.delta_source_2);
  command.mix_factor = float(item.mix_factor) / 255.0f;
  return command;
}

struct PixelCopyGroups {
  Vector<PixelCopyGroup> groups;

  void copy_pixels(ImBuf &tile_buffer) const
  {
    if (tile_buffer.rect_float) {
      ImageBufferAccessor<float4> accessor(tile_buffer);
      copy_pixels<float4>(accessor);
    }
    else {
      ImageBufferAccessor<int> accessor(tile_buffer);
      copy_pixels<int>(accessor);
    }
  }

 private:
  template<typename T> void copy_pixels(ImageBufferAccessor<T> &image_buffer) const
  {
    for (const PixelCopyGroup &group : groups) {
      PixelCopyCommand copy_command(group);
      for (const PixelCopyItem &item : group.items) {
        copy_command += item;
        copy_command.mix_source_and_write_destination<T>(image_buffer);
      }
    }
  }
};

}  // namespace blender::kernel::pbvh::pixels
