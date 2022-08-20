/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#include "draw_manager.hh"

namespace blender::draw {

void Manager::submit(const Pass &pass)
{
  command::RecordingState state;
  pass.submit(state);
}

std::string Manager::serialize(const Pass &pass)
{
  return pass.serialize();
}

}  // namespace blender::draw
