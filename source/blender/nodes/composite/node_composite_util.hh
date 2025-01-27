/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2006 Blender Foundation */

/** \file
 * \ingroup nodes
 */

#pragma once

#include "DNA_ID.h"
#include "DNA_node_types.h"

#include "BLT_translation.h"

#include "node_composite_register.hh"
#include "node_util.h"

#include "NOD_composite.h"
#include "NOD_socket.h"
#include "NOD_socket_declarations.hh"

#define CMP_SCALE_MAX 12000

bool cmp_node_poll_default(const struct bNodeType *ntype,
                           const struct bNodeTree *ntree,
                           const char **r_disabled_hint);
void cmp_node_update_default(struct bNodeTree *ntree, struct bNode *node);
void cmp_node_type_base(struct bNodeType *ntype, int type, const char *name, short nclass);
