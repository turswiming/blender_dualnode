/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 * A unique identifier for each object component.
 * It is used to access each component data such as matrices and object attributes.
 * It is valid only for the current draw, it is not persistent.
 *
 * The most significant bit is used to encode if the object needs to invert the front face winding
 * because of its object matrix handedness. This is handy because this means sorting inside
 * #MultiDraw command will put all inverted commands last.
 *
 * Default value of 0 points toward an non-cull-able object with unit bounding box centered at
 * the origin.
 */

#include "BKE_duplilist.h"
#include "DNA_object_types.h"

#include "draw_shader_shared.h"

namespace blender::draw {

struct ResourceHandle {
  uint raw;

  ResourceHandle() = default;
  ResourceHandle(uint raw_) : raw(raw_){};
  ResourceHandle(uint index, bool inverted_handedness)
  {
    raw = index;
    SET_FLAG_FROM_TEST(raw, inverted_handedness, 0x80000000u);
  }

  bool has_inverted_handedness() const
  {
    return (raw & 0x80000000u) != 0;
  }

  uint resource_index() const
  {
    return (raw & 0x7FFFFFFFu);
  }
};

};  // namespace blender::draw
