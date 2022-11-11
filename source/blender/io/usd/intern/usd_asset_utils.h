/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021  NVIDIA Corporation. All rights reserved. */
#pragma once

namespace blender::io::usd {

  /* Invoke the the USD asset resolver to copy assets. */
  bool copy_usd_asset(const char *src, const char *dst, bool overwrite);

}  // namespace blender::io::usd
