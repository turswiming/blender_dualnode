/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 */

#include "BLI_sys_types.h"

#include "draw_pass.hh"

#include <string>

namespace blender::draw {

class Manager {

 public:
  /**
   * Submit a pass for drawing. All resource reference will be dereferenced and commands will be
   * sent to GPU.
   */
  void submit(const PassSimple &pass);
  void submit(const PassMain &pass);
};

}  // namespace blender::draw

/* TODO(@fclem): This is for testing. The manager should be passed to the engine through the
 * callbacks. */
blender::draw::Manager *DRW_manager_get();
