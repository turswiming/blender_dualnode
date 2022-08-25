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

#define DRW_VIEW_UBO_SLOT 0

#define DRW_OBJ_MAT_SLOT 7
#define DRW_COMMAND_SLOT 6
#define DRW_OBJ_INFOS_SLOT 5
#define DRW_OBJ_ATTR_SLOT 4

#define DRW_COMMAND_GROUP_SIZE 64
#define DRW_FINALIZE_GROUP_SIZE 64
/* Must be multiple of 32. Set to 64 for optimal scheduling on GCN hardware. */
#define DRW_VISIBILITY_GROUP_SIZE 64
