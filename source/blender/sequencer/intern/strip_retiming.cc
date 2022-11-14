/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup bke
 */

#include "MEM_guardedalloc.h"

#include "DNA_scene_types.h"
#include "DNA_sequence_types.h"

#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_span.hh"

#include "BKE_fcurve.h"
#include "BKE_movieclip.h"
#include "BKE_scene.h"
#include "BKE_sound.h"

#include "DNA_anim_types.h"
#include "DNA_sound_types.h"

#include "IMB_imbuf.h"

#include "RNA_prototypes.h"

#include "SEQ_channels.h"
#include "SEQ_iterator.h"
#include "SEQ_relations.h"
#include "SEQ_render.h"
#include "SEQ_retiming.hh"
#include "SEQ_sequencer.h"
#include "SEQ_time.h"
#include "SEQ_transform.h"

#include "sequencer.h"
#include "strip_time.h"
#include "utils.h"


using blender::MutableSpan;

MutableSpan<SeqRetimingHandle> SEQ_retiming_handles_get(const Sequence *seq)
{
  blender::MutableSpan<SeqRetimingHandle> handles(seq->retiming_handles,
                                                  seq->retiming_handle_count);
  return handles;
}

struct SeqRetimingHandle *SEQ_retiming_last_handle_get(const struct Sequence *seq)
{
  return seq->retiming_handles + seq->retiming_handle_count - 1;
}

int SEQ_retiming_handle_index_get(const Sequence *seq, const SeqRetimingHandle *handle)
{
  return handle - seq->retiming_handles;
}

static const SeqRetimingHandle *retiming_find_closest_handle_up_to_frame(const Sequence *seq,
                                                                         const int frame_index)
{
  const SeqRetimingHandle *prev = NULL;
  for (auto const &handle : SEQ_retiming_handles_get(seq)) {
    if (handle.strip_frame_index <= frame_index) {
      prev = &handle;
    }
    if (handle.strip_frame_index > frame_index) {
      break;
    }
  }

  return prev;
}

int SEQ_retiming_handles_count(const Sequence *seq)
{
  return seq->retiming_handle_count;
}

void SEQ_retiming_data_ensure(const Scene *scene, Sequence *seq)
{
  if (!SEQ_retiming_is_allowed(seq)) {
    return;
  }

  if (seq->retiming_handles != NULL) {
    return;
  }

  seq->retiming_handles = (SeqRetimingHandle *)MEM_calloc_arrayN(
      2, sizeof(SeqRetimingHandle), __func__);
  SeqRetimingHandle *handle = seq->retiming_handles + 1;
  handle->strip_frame_index = seq_time_strip_original_content_length_get(scene, seq) - 1;
  handle->retiming_factor = 1.0f;
  seq->retiming_handle_count = 2;
}

void SEQ_retiming_data_clear(Sequence *seq)
{
  seq->retiming_handles = NULL;
  seq->retiming_handle_count = 0;
}

bool SEQ_retiming_is_active(const Scene *scene, const Sequence *seq)
{
  if (seq->retiming_handle_count < 2) {
    return false;
  }
  if (seq->retiming_handle_count > 2) {
    return true;
  }

  SeqRetimingHandle *last_handle = SEQ_retiming_last_handle_get(seq);
  const int last_frame_index = seq->len / seq_time_media_playback_rate_factor_get(scene, seq);
  if (last_handle->strip_frame_index != last_frame_index - 1) {
    return true;
  }

  return false;
}

bool SEQ_retiming_is_allowed(const Sequence *seq)
{
  return ELEM(seq->type,
              SEQ_TYPE_SOUND_RAM,
              SEQ_TYPE_IMAGE,
              SEQ_TYPE_META,
              SEQ_TYPE_SCENE,
              SEQ_TYPE_MOVIE,
              SEQ_TYPE_MOVIECLIP,
              SEQ_TYPE_MASK);
}

float seq_retiming_evaluate(const Sequence *seq, const int frame_index)
{
  const SeqRetimingHandle *previous_handle = retiming_find_closest_handle_up_to_frame(seq,
                                                                                      frame_index);
  const SeqRetimingHandle *next_handle = previous_handle + 1;

  if (next_handle == NULL) {
    return 1.0f;
  }

  const int segment_length = next_handle->strip_frame_index - previous_handle->strip_frame_index;
  const int segment_index = frame_index - previous_handle->strip_frame_index;
  const float segment_fac = segment_index / (float)segment_length;

  const float target_diff = next_handle->retiming_factor - previous_handle->retiming_factor;
  return previous_handle->retiming_factor + (target_diff * segment_fac);
}

SeqRetimingHandle *SEQ_retiming_add_handle(Sequence *seq, const int timeline_frame)
{
  float frame_index = timeline_frame - SEQ_time_start_frame_get(seq);
  float value = seq_retiming_evaluate(seq, frame_index);

  const SeqRetimingHandle *closest_handle = retiming_find_closest_handle_up_to_frame(seq,
                                                                                     frame_index);
  if (closest_handle->strip_frame_index == frame_index) {
    return NULL; /* Retiming handle already exists. */
  }

  SeqRetimingHandle *handles = seq->retiming_handles;
  size_t handle_count = SEQ_retiming_handles_count(seq);
  const int new_handle_index = closest_handle - handles + 1;
  BLI_assert(new_handle_index >= 0);
  BLI_assert(new_handle_index < handle_count);

  SeqRetimingHandle *new_handles = (SeqRetimingHandle *)MEM_callocN(
      (handle_count + 1) * sizeof(SeqRetimingHandle), __func__);
  if (new_handle_index > 0) {
    memcpy(new_handles, handles, new_handle_index * sizeof(SeqRetimingHandle));
  }
  if (new_handle_index < handle_count) {
    memcpy(new_handles + new_handle_index + 1,
           handles + new_handle_index,
           (handle_count - new_handle_index) * sizeof(SeqRetimingHandle));
  }
  MEM_freeN(handles);
  seq->retiming_handles = new_handles;
  seq->retiming_handle_count++;

  SeqRetimingHandle *added_handle = (new_handles + new_handle_index);
  added_handle->strip_frame_index = frame_index;
  added_handle->retiming_factor = value;

  return added_handle;
}

void SEQ_retiming_offset_handle(const Scene *scene,
                                Sequence *seq,
                                SeqRetimingHandle *handle,
                                const int offset)
{
  if (handle->strip_frame_index == 0) {
    return; /* First handle can not be moved. */
  }

  MutableSpan handles = SEQ_retiming_handles_get(seq);
  for (; handle < handles.end(); handle++) {
    handle->strip_frame_index += offset;
  }

  SEQ_time_update_meta_strip_range(scene, seq_sequence_lookup_meta_by_seq(scene, seq));
  seq_time_update_effects_strip_range(scene, seq_sequence_lookup_effects_by_seq(scene, seq));
}

void SEQ_retiming_remove_handle(Sequence *seq, SeqRetimingHandle *handle)
{
  SeqRetimingHandle *last_handle = SEQ_retiming_last_handle_get(seq);
  if (handle->strip_frame_index == 0 || handle == last_handle) {
    return; /* First and last handle can not be removed. */
  }

  size_t handle_count = SEQ_retiming_handles_count(seq);
  SeqRetimingHandle *handles = (SeqRetimingHandle *)MEM_callocN(
      (handle_count - 1) * sizeof(SeqRetimingHandle), __func__);

  const int handle_index = handle - seq->retiming_handles;
  memcpy(handles, seq->retiming_handles, (handle_index) * sizeof(SeqRetimingHandle));
  memcpy(handles + handle_index,
         seq->retiming_handles + handle_index + 1,
         (handle_count - handle_index - 1) * sizeof(SeqRetimingHandle));
  MEM_freeN(seq->retiming_handles);
  seq->retiming_handles = handles;
  seq->retiming_handle_count--;
}

float SEQ_retiming_handle_speed_get(const Scene *scene,
                                    const Sequence *seq,
                                    const SeqRetimingHandle *handle)
{
  if (handle->strip_frame_index == 0) {
    return 1.0f;
  }

  const SeqRetimingHandle *handle_prev = handle - 1;

  const int frame_index_max = seq_time_strip_original_content_length_get(scene, seq) - 1;
  const int frame_retimed_prev = round_fl_to_int(handle_prev->retiming_factor * frame_index_max);
  const int frame_index_prev = handle_prev->strip_frame_index;
  const int frame_retimed = round_fl_to_int(handle->retiming_factor * frame_index_max);
  const int frame_index = handle->strip_frame_index;

  const int fragment_length_retimed = frame_retimed - frame_retimed_prev;
  const int fragment_length_original = frame_index - frame_index_prev;

  const float speed = (float)fragment_length_retimed / (float)fragment_length_original;
  return speed;
}

using std::vector;

class RetimingRange {
 public:
  int start, end;
  float speed;

  enum eIntersectType {
    FULL,
    PARTIAL_START,
    PARTIAL_END,
    INSIDE,
    NONE,
  };

  RetimingRange(int start_frame, int end_frame, float speed)
      : start(start_frame), end(end_frame), speed(speed)
  {
  }

  const eIntersectType intersect_type(const RetimingRange &other) const
  {
    if (other.start <= start && other.end >= end) {
      return FULL;
    }
    if (other.start > start && other.start < end && other.end > start && other.end < end) {
      return INSIDE;
    }
    if (other.start > start && other.start < end) {
      return PARTIAL_END;
    }
    if (other.end > start && other.end < end) {
      return PARTIAL_START;
    }
    return NONE;
  }
};

class RetimingRangeData {
 public:
  vector<RetimingRange> ranges;
  RetimingRangeData(const Scene *scene, const Sequence *seq)
  {
    MutableSpan handles = SEQ_retiming_handles_get(seq);
    for (const SeqRetimingHandle &handle : handles) {
      if (handle.strip_frame_index == 0) {
        continue;
      }
      const SeqRetimingHandle *handle_prev = &handle - 1;
      float speed = SEQ_retiming_handle_speed_get(scene, seq, &handle);
      int frame_start = SEQ_time_start_frame_get(seq) + handle_prev->strip_frame_index;
      int frame_end = SEQ_time_start_frame_get(seq) + handle.strip_frame_index;

      RetimingRange range = RetimingRange(frame_start, frame_end, speed);
      ranges.push_back(range);
    }
  }

  RetimingRangeData &operator*=(const RetimingRangeData rhs)
  {
    if (ranges.size() == 0) {
      for (const RetimingRange &rhs_range : rhs.ranges) {
        RetimingRange range = RetimingRange(
            rhs_range.start, rhs_range.end, rhs_range.speed);
        ranges.push_back(range);
      }
    }
    else {
      for (int i = 0; i < ranges.size(); i++) {
        RetimingRange &range = ranges[i];
        for (const RetimingRange &rhs_range : rhs.ranges) {
          if (range.intersect_type(rhs_range) == range.NONE) {
            continue;
          }
          else if (range.intersect_type(rhs_range) == range.FULL) {
            range.speed *= rhs_range.speed;
          }
          else if (range.intersect_type(rhs_range) == range.PARTIAL_START) {
            RetimingRange range_left = RetimingRange(
                range.start, rhs_range.end, range.speed * rhs_range.speed);
            range.start = rhs_range.end + 1;
            ranges.insert(ranges.begin() + i, range_left);
          }
          else if (range.intersect_type(rhs_range) == range.PARTIAL_END) {
            RetimingRange range_left = RetimingRange(
                range.start, rhs_range.start - 1, range.speed);
            range.start = rhs_range.start;
            ranges.insert(ranges.begin() + i, range_left);
          }
          else if (range.intersect_type(rhs_range) == range.INSIDE) {
            RetimingRange range_left = RetimingRange(
                range.start, rhs_range.start - 1, range.speed);
            RetimingRange range_mid = RetimingRange(
                rhs_range.start, rhs_range.start, rhs_range.speed * range.speed);
            range.start = rhs_range.end + 1;
            ranges.insert(ranges.begin() + i, range_left);
            ranges.insert(ranges.begin() + i, range_mid);
            break;
          }
        }
      }
    }
    return *this;
  }
};

static RetimingRangeData seq_retiming_range_data_get(const Scene *scene, const Sequence *seq)
{
  RetimingRangeData strip_retiming_data = RetimingRangeData(scene, seq);

  const Sequence *meta_parent = seq_sequence_lookup_meta_by_seq(scene, seq);
  if (meta_parent == nullptr) {
    return strip_retiming_data;
  }

  RetimingRangeData meta_retiming_data = RetimingRangeData(scene, meta_parent);
  strip_retiming_data *= meta_retiming_data;
  return strip_retiming_data;
}

#include <BKE_sound.h>
void SEQ_retiming_sound_animation_data_set(const Scene *scene, const Sequence *seq)
{
  RetimingRangeData retiming_data = seq_retiming_range_data_get(scene, seq);
  for (const RetimingRange &range : retiming_data.ranges) {
    BKE_sound_set_scene_sound_pitch_constant_range(
        seq->scene_sound, range.start, range.end, range.speed);
  }
}
