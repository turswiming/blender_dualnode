/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#include "BLI_math.h"
#include "BLI_math_vec_types.hh"
#include "BLI_vector.hh"

#include "IMB_imbuf.h"
#include "IMB_imbuf_types.h"

#include "BKE_image_wrappers.hh"
#include "BKE_pbvh.h"
#include "BKE_pbvh_pixels.hh"

namespace blender::bke::pbvh::pixels {

static void copy_pixels_reinit(PixelCopyTiles &tiles)
{
  tiles.clear();
}

void BKE_pbvh_pixels_copy_update(PBVH &pbvh, Image &image, ImageUser &image_user)
{
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  copy_pixels_reinit(pbvh_data.tiles_copy_pixels);

  ImageUser tile_user = image_user;
  LISTBASE_FOREACH (ImageTile *, tile, &image.tiles) {
    image::ImageTileWrapper image_tile = image::ImageTileWrapper(tile);
    tile_user.tile = image_tile.get_tile_number();

    ImBuf *tile_buffer = BKE_image_acquire_ibuf(&image, &tile_user, nullptr);
    if (tile_buffer == nullptr) {
      continue;
    }

    ushort2 tile_resolution(tile_buffer->x, tile_buffer->y);

    BKE_image_release_ibuf(&image, tile_buffer, nullptr);
  }
}

void BKE_pbvh_pixels_copy_pixels(PBVH &pbvh,
                                 Image &image,
                                 ImageUser &image_user,
                                 image::TileNumber tile_number)
{
  PBVHData &pbvh_data = BKE_pbvh_pixels_data_get(pbvh);
  std::optional<std::reference_wrapper<PixelCopyTile>> pixel_tile =
      pbvh_data.tiles_copy_pixels.find_tile(tile_number);
  if (!pixel_tile.has_value()) {
    return;
  }

  ImageUser tile_user = image_user;
  tile_user.tile = tile_number;
  ImBuf *tile_buffer = BKE_image_acquire_ibuf(&image, &tile_user, nullptr);
  if (tile_buffer == nullptr) {
    return;
  }
  pixel_tile->get().copy_pixels(*tile_buffer);

  BKE_image_release_ibuf(&image, tile_buffer, nullptr);
}

}  // namespace blender::bke::pbvh::pixels
