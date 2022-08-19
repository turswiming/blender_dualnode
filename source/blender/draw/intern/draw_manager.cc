/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 *
 */

#include "draw_manager.hh"

namespace blender::draw {

void Manager::submit(const Pass &pass)
{
  command::RecordingState state;
  //   pass.submit(&state);
}

void Manager::submit(const PassSimple &pass)
{
  command::RecordingState state;
  //   pass.submit(&state);
}

template<int64_t len> void Manager::submit(const PassImmutable<len> &pass)
{
  command::RecordingState state;
  pass.submit(&state);
}

std::string Manager::serialize(const Pass &pass)
{
  //   return pass.serialize();
  return "";
}

std::string Manager::serialize(const PassSub &pass)
{
  //   return pass.serialize();
  return "";
}

std::string Manager::serialize(const PassSimple &pass)
{
  //   return pass.serialize();
  return "";
}

template<int64_t len> std::string Manager::serialize(const PassImmutable<len> &pass)
{
  return pass.serialize();
}

}  // namespace blender::draw
