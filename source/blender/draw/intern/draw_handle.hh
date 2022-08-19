/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 */

namespace blender::draw {

/* -------------------------------------------------------------------- */
/** \name DRWState to GPU state conversion
 * \{ */

/**
 * A unique identifier for each object component.
 * It is used to access each component data such as matrices and object attributes.
 * It is valid only for the current draw, it is not persistent.
 *
 * The most significant bit is used to encode if the object needs to invert the front face winding
 * because of its object matrix handedness.
 *
 * Default value of 0 points toward an non-cull-able object with unit bounding box centered at
 * the origin.
 */
struct ResourceHandle {
  uint raw;

  bool has_inverted_handedness() const
  {
    return (raw & 0x80000000u) != 0;
  }

  uint resource_index() const
  {
    return (raw & 0x7FFFFFFFu);
  }
};

/** \} */

};  // namespace blender::draw
