/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

#pragma once

/** \file
 * \ingroup sequencer
 */

#include "BLI_span.hh"

struct Scene;
struct Sequence;
struct SeqRetimingHandle;

blender::MutableSpan<SeqRetimingHandle> SEQ_retiming_handles_get(const Sequence *seq);
void SEQ_retiming_data_ensure(const struct Scene *scene, struct Sequence *seq);
void SEQ_retiming_data_clear(struct Sequence *seq);
bool SEQ_retiming_is_allowed(const struct Sequence *seq);
void SEQ_retiming_add_handle(struct Sequence *seq, const int timeline_frame);
void SEQ_retiming_remove_handle(struct Sequence *seq, SeqRetimingHandle *handle);
void SEQ_retiming_offset_handle(const struct Scene *scene,
                                struct Sequence *seq,
                                SeqRetimingHandle *handle,
                                const int offset);
float SEQ_retiming_handle_speed_get(const struct Scene *scene,
                                    const struct Sequence *seq,
                                    const int point_index);
