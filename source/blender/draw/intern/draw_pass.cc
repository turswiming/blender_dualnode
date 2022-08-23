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

template<class T> inline void detail::PassBase<T>::submit(command::RecordingState &state) const
{
  draw_commands_buf_.bind();

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

/** \} */

/* -------------------------------------------------------------------- */
/** \name Pass Serialization
 * \{ */

template<class T> inline std::string detail::PassBase<T>::serialize(std::string line_prefix) const
{
  std::stringstream ss;
  ss << line_prefix << "." << debug_name << std::endl;
  line_prefix += "  ";
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
        ss << line_prefix << commands_[header.command_index].serialize(header.type) << std::endl;
        break;
    }
  }
  return ss.str();
}

/** \} */

}  // namespace blender::draw
