/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2022 Blender Foundation */

#ifndef __UTIL_GUIDING_H__
#define __UTIL_GUIDING_H__

#include "util/system.h"

CCL_NAMESPACE_BEGIN

static inline bool guiding_supported()
{
#ifdef WITH_PATH_GUIDING
#  if defined(__ARM_NEON)
  return true;
#  else
  return system_cpu_support_sse41();
#  endif
#else
  return false;
#endif
}

CCL_NAMESPACE_END

#endif /* __UTIL_GUIDING_H__ */
