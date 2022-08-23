/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#include "GPU_debug.h"

#include "draw_pass.hh"

namespace blender::draw {

/* -------------------------------------------------------------------- */
/** \name Pass Submission
 * \{ */

void PassSimple::submit(command::RecordingState &state) const
{
  draw_commands_buf_.finalize();
  draw_commands_buf_.bind();

  GPU_debug_group_begin(debug_name);

  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::SubPass:
        sub_passes_[header.command_index].submit(state);
        break;
      default:
        commands_[header.command_index].execute(header.type, state);
        break;
    }
  }

  GPU_debug_group_end();
}

void PassSimple::Sub::submit(command::RecordingState &state) const
{
  GPU_debug_group_begin(debug_name);

  for (const command::Header &header : headers_) {
    commands_[header.command_index].execute(header.type, state);
  }

  GPU_debug_group_end();
}

void PassMain::submit(command::RecordingState &state) const
{
  GPU_debug_group_begin(debug_name);

  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::SubPass:
        sub_passes_[header.command_index].submit(state);
        break;
      case Type::MultiDraw:
        /* TODO */
        break;
      default:
        commands_[header.command_index].execute(header.type, state);
        break;
    }
  }

  GPU_debug_group_end();
}

void PassMain::Sub::submit(command::RecordingState &state) const
{
  draw_commands_buf_.finalize();
  draw_commands_buf_.bind();

  GPU_debug_group_begin(debug_name);

  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::MultiDraw:
        /* TODO */
        break;
      default:
        commands_[header.command_index].execute(header.type, state);
        break;
    }
  }

  GPU_debug_group_end();
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Pass Serialization
 * \{ */

std::string PassSimple::serialize() const
{
  std::stringstream ss;
  ss << "PassSimple(" << debug_name << ")" << std::endl;
  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::None:
        break;
      case Type::SubPass:
        ss << sub_passes_[header.command_index].serialize();
        break;
      default:
        ss << commands_[header.command_index].serialize(header.type) << std::endl;
        break;
    }
  }
  return ss.str();
}

std::string PassSimple::Sub::serialize() const
{
  std::stringstream ss;
  ss << ".sub(" << debug_name << ")" << std::endl;
  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::None:
        break;
      default:
        ss << "  " << commands_[header.command_index].serialize(header.type) << std::endl;
        break;
    }
  }
  return ss.str();
}

std::string PassMain::serialize() const
{
  std::stringstream ss;
  ss << "PassMain(" << debug_name << ")" << std::endl;
  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::None:
        break;
      case Type::SubPass:
        ss << sub_passes_[header.command_index].serialize();
        break;
      case Type::MultiDraw:
        /* TODO */
        break;
      default:
        ss << commands_[header.command_index].serialize(header.type) << std::endl;
        break;
    }
  }
  return ss.str();
}

std::string PassMain::Sub::serialize() const
{
  std::stringstream ss;
  ss << ".sub(" << debug_name << ")" << std::endl;
  for (const command::Header &header : headers_) {
    switch (header.type) {
      case Type::None:
        break;
      case Type::MultiDraw:
        /* TODO */
        break;
      default:
        ss << "  " << commands_[header.command_index].serialize(header.type) << std::endl;
        break;
    }
  }
  return ss.str();
}

/** \} */

}  // namespace blender::draw
