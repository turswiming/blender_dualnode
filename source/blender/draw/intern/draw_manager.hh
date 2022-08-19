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
  /* TODO(fclem): Why doesn't this work? */
  // void submit(const PassCommon &pass);
  void submit(const Pass &pass);
  void submit(const PassSimple &pass);
  template<int64_t len> void submit(const PassImmutable<len> &pass);

  /**
   * Will turn the command stream into a debug string.
   */
  /* TODO(fclem): Why doesn't this work? */
  // std::string serialize(const PassCommon &pass);
  std::string serialize(const Pass &pass);
  std::string serialize(const PassSub &pass);
  std::string serialize(const PassSimple &pass);
  template<int64_t len> std::string serialize(const PassImmutable<len> &pass);
};

}  // namespace blender::draw

/* TODO(@fclem): This is for testing. The manager should be passed to the engine through the
 * callbacks. */
blender::draw::Manager *DRW_manager_get();
