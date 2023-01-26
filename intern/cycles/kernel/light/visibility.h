/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#pragma once

CCL_NAMESPACE_BEGIN

ccl_device_inline void light_visibility_correction(KernelGlobals kg,
                                                   IntegratorState state,
                                                   const int shader,
                                                   Spectrum *light_visibility)
{
  Spectrum visible_components = zero_spectrum();
  const BsdfEval scatter_eval = INTEGRATOR_STATE(state, path, scatter_eval);
  Spectrum transmission = scatter_eval.sum - (scatter_eval.glossy + scatter_eval.diffuse);
  if (!(shader & SHADER_EXCLUDE_DIFFUSE)) {
    visible_components += scatter_eval.diffuse;
  }

  if (!(shader & SHADER_EXCLUDE_GLOSSY)) {
    visible_components += scatter_eval.glossy;
  }

  if (!(shader & SHADER_EXCLUDE_TRANSMIT)) {
    visible_components += transmission;
  }

  *light_visibility = safe_divide(visible_components, scatter_eval.sum);
}

CCL_NAMESPACE_END