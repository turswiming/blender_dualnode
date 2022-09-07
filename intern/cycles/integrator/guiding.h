/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

CCL_NAMESPACE_BEGIN

struct GuidingParams {
  /* The subset of path guiding parameters that can trigger a creation/rebuild
   * of the guiding field. */
  bool use = false;
  GuidingDistributionType type = GUIDING_TYPE_PARALLAX_AWARE_VMM;
  int training_iterations = 128;
  GuidingParams() = default;

  bool modified(const GuidingParams &other) const
  {
    return !((use == other.use) && (type == other.type) &&
             (training_iterations == other.training_iterations));
  }
};

CCL_NAMESPACE_END
