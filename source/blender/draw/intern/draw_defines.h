/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation.
 */

/** \file
 * \ingroup draw
 *
 * List of defines that are shared with the GPUShaderCreateInfos. We do this to avoid
 * dragging larger headers into the createInfo pipeline which would cause problems.
 */

#pragma once

#define DRW_OBJ_MAT_SLOT 7
#define DRW_OBJ_ATTR_SLOT 6

#define DRW_FINALIZE_GROUP_SIZE 64
#define DRW_VISIBILITY_GROUP_SIZE 64