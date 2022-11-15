/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021  NVIDIA Corporation. All rights reserved. */
#pragma once

#include <pxr/usd/usd/stage.h>

#include <string>

namespace blender::io::usd {

  /* Invoke the the USD asset resolver to copy assets. */
  bool copy_usd_asset(const char *src, const char *dst, bool overwrite);

  /* Return the path to a 'textures' directory which is a sibling to the given
   * stage's root layer path. */
  std::string get_textures_dir(const pxr::UsdStageRefPtr stage);

  bool usd_path_exists(const char *src);

  bool usd_paths_equal(const char *p1, const char *p2);

}  // namespace blender::io::usd
