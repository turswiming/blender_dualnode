/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#include "draw_manager.hh"

namespace blender::draw {

void Manager::submit(const PassSimple &pass)
{
  command::RecordingState state;
  pass.submit(state);
}

void Manager::submit(const PassMain &pass)
{
  command::RecordingState state;
  pass.submit(state);
}

}  // namespace blender::draw
