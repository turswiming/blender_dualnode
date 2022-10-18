/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup spseq
 */

#include "MEM_guardedalloc.h"

#include "BLI_blenlib.h"
#include "BLI_span.hh"

#include "DNA_anim_types.h"

#include "BKE_context.h"
#include "BKE_fcurve.h"
#include "BKE_scene.h"

#include "BLF_api.h"

#include "GPU_batch.h"
#include "GPU_batch_utils.h"
#include "GPU_immediate.h"
#include "GPU_immediate_util.h"
#include "GPU_matrix.h"
#include "GPU_select.h"
#include "GPU_state.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"

#include "WM_api.h"
#include "WM_types.h"

#include "ED_gizmo_library.h"
#include "ED_screen.h"
#include "ED_view3d.h"

#include "UI_interface.h"
#include "UI_interface_icons.h"
#include "UI_resources.h"
#include "UI_view2d.h"

#include "SEQ_iterator.h"
#include "SEQ_retiming.hh"
#include "SEQ_sequencer.h"
#include "SEQ_time.h"

/* Own include. */
#include "sequencer_intern.h"

#define REMOVE_GIZMO_HEIGHT 12.0f * U.dpi_fac;        /* Pixels from bottom of strip.  */
#define RETIME_HANDLE_TRIANGLE_SIZE 10.0f * U.dpi_fac /* Also used for mouseover test. */
#define RETIME_BUTTON_SIZE 0.6f                       /* Factor based on icon size. */

static float strip_y_rescale(const Sequence *seq, const float y_value)
{
  const float y_range = SEQ_STRIP_OFSTOP - SEQ_STRIP_OFSBOTTOM;
  return (y_value * y_range) + seq->machine + SEQ_STRIP_OFSBOTTOM;
}

class RetimingHandle {
 public:
  int index;
  bool is_last_handle;

  RetimingHandle(const Sequence *seq, const SeqRetimingHandle *handle, const int handle_index)
  {
    seq_ = seq;
    handle_ = handle;
    index = handle_index;

    const int handle_count = SEQ_retiming_handles_count(seq);
    is_last_handle = (handle_index == handle_count - 1);
  }

  float x()
  {
    return SEQ_time_start_frame_get(seq_) + handle_->strip_frame_index + (is_last_handle ? 1 : 0);
  }

  float y()
  {
    return strip_y_rescale(seq_, handle_->retiming_factor);
  }

  RetimingHandle next_handle_get()
  {
    const int handle_count = SEQ_retiming_handles_count(seq_);
    BLI_assert(index < handle_count - 1);
    const SeqRetimingHandle *next = handle_ + 1;
    return RetimingHandle(seq_, next, index + 1);
  }

 private:
  const Sequence *seq_;
  const SeqRetimingHandle *handle_;
};

template<typename T> class RetimingHandlesIterator {
 public:
  using value_type = SeqRetimingHandle;
  using pointer = value_type *;
  using reference = T;

  struct Iterator {
    Iterator(const Sequence *seq, pointer ptr, int handle_index)
    {
      ptr_ = ptr;
      seq_ = seq;
      handle_index_ = handle_index;
    }

    reference operator*() const
    {
      return reference(seq_, ptr_, handle_index_);
    }
    pointer operator->()
    {
      return ptr_;
    }
    Iterator &operator++()
    {
      ptr_++;
      handle_index_++;
      return *this;
    }
    Iterator operator++(int)
    {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    friend bool operator==(const Iterator &a, const Iterator &b)
    {
      return a.ptr_ == b.ptr_;
    };
    friend bool operator!=(const Iterator &a, const Iterator &b)
    {
      return a.ptr_ != b.ptr_;
    };

   private:
    pointer ptr_;
    int handle_index_;
    const Sequence *seq_;
  };

  Iterator begin()
  {
    return Iterator(seq_, handles_array_, 0);
  }
  Iterator end()
  {
    const int count = this->count();
    return Iterator(seq_, handles_array_ + count, count - 1);
  }

 public:
  RetimingHandlesIterator(const Sequence *seq)
  {
    handles_array_ = seq->retiming_handles;
    seq_ = seq;
  }

  reference *mouse_over_handle_get(const View2D *v2d, const int mval[2])
  {
    int best_distance = INT_MAX;
    int best_handle_index = -1;

    RetimingHandlesIterator handles = RetimingHandlesIterator(seq_);
    for (auto handle : handles) {
      int distance = round_fl_to_int(fabsf(UI_view2d_view_to_region_x(v2d, handle.x()) - mval[0]));

      if (distance < RETIME_HANDLE_TRIANGLE_SIZE && distance < best_distance) {
        best_distance = distance;
        best_handle_index = handle.index;
      }
    }

    if (best_handle_index == -1) {
      return nullptr;
    }

    return new reference(seq_, handles_array_ + best_handle_index, best_handle_index);
  }

  int count()
  {
    return SEQ_retiming_handles_count(seq_);
  }

 private:
  const Sequence *seq_;
  SeqRetimingHandle *handles_array_;
};

static float pixels_to_view_width(const bContext *C, const float width)
{
  const View2D *v2d = UI_view2d_fromcontext(C);
  float scale_x = UI_view2d_view_to_region_x(v2d, 1) - UI_view2d_view_to_region_x(v2d, 0.0f);
  return width / scale_x;
}

static float pixels_to_view_height(const bContext *C, const float height)
{
  const View2D *v2d = UI_view2d_fromcontext(C);
  float scale_y = UI_view2d_view_to_region_y(v2d, 1) - UI_view2d_view_to_region_y(v2d, 0.0f);
  return height / scale_y;
}

static float strip_start_screenspace_get(const bContext *C, const Sequence *seq)
{
  const View2D *v2d = UI_view2d_fromcontext(C);
  const Scene *scene = CTX_data_scene(C);
  return UI_view2d_view_to_region_x(v2d, SEQ_time_left_handle_frame_get(scene, seq));
}

static float strip_end_screenspace_get(const bContext *C, const Sequence *seq)
{
  const View2D *v2d = UI_view2d_fromcontext(C);
  const Scene *scene = CTX_data_scene(C);
  return UI_view2d_view_to_region_x(v2d, SEQ_time_right_handle_frame_get(scene, seq));
}

static Sequence *active_seq_from_context(const bContext *C)
{
  const Editing *ed = SEQ_editing_get(CTX_data_scene(C));
  return ed->act_seq;
}

static rctf strip_box_get(const bContext *C, const Sequence *seq)
{
  const View2D *v2d = UI_view2d_fromcontext(C);
  rctf rect;
  rect.xmin = strip_start_screenspace_get(C, seq);
  rect.xmax = strip_end_screenspace_get(C, seq);
  rect.ymin = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 0));
  rect.ymax = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 1));
  return rect;
}

static rctf remove_box_get(const bContext *C, const Sequence *seq)
{
  rctf rect = strip_box_get(C, seq);
  rect.ymax = rect.ymin + REMOVE_GIZMO_HEIGHT;
  return rect;
}

static bool mouse_is_inside_box(const rctf *box, const int mval[2])
{
  return mval[0] >= box->xmin && mval[0] <= box->xmax && mval[1] >= box->ymin &&
         mval[1] <= box->ymax;
}

/* -------------------------------------------------------------------- */
/** \name Retiming Add Handle Gizmo
 * \{ */

typedef struct RetimeButtonGizmo {
  wmGizmo gizmo;
  int icon_id;
  const Sequence *seq_under_mouse;
  bool is_mouse_over_gizmo;
} RetimeButtonGizmo;

typedef struct ButtonDimensions {
  float height;
  float width;
  float x;
  float y;
} ButtonDimensions;

static ButtonDimensions button_dimensions_get(const bContext *C, const RetimeButtonGizmo *gizmo)
{
  const Scene *scene = CTX_data_scene(C);
  const View2D *v2d = UI_view2d_fromcontext(C);
  const Sequence *seq = gizmo->seq_under_mouse;

  const float icon_height = UI_icon_get_height(gizmo->icon_id) * U.dpi_fac;
  const float icon_width = UI_icon_get_width(gizmo->icon_id) * U.dpi_fac;
  const float icon_x = UI_view2d_view_to_region_x(v2d, BKE_scene_frame_get(scene)) +
                       icon_width / 2;
  const float icon_y = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 0.5)) -
                       icon_height / 2;
  const ButtonDimensions dimensions = {icon_height, icon_width, icon_x, icon_y};
  return dimensions;
}

static rctf button_box_get(const bContext *C, const RetimeButtonGizmo *gizmo)
{
  ButtonDimensions button_dimensions = button_dimensions_get(C, gizmo);

  float x_range = button_dimensions.width;
  float y_range = button_dimensions.height;

  rctf rect;
  rect.xmin = button_dimensions.x;
  rect.xmax = button_dimensions.x + x_range;
  rect.ymin = button_dimensions.y;
  rect.ymax = button_dimensions.y + y_range;

  return rect;
}

static void gizmo_retime_handle_add_draw(const bContext *C, wmGizmo *gz)
{
  RetimeButtonGizmo *gizmo = (RetimeButtonGizmo *)gz;
  if (gizmo->seq_under_mouse == NULL) {
    return;
  }

  const ButtonDimensions button = button_dimensions_get(C, gizmo);
  const rctf strip_box = strip_box_get(C, gizmo->seq_under_mouse);
  if (!BLI_rctf_isect_pt(&strip_box, button.x, button.y)) {
    return;
  }

  wmOrtho2_region_pixelspace(CTX_wm_region(C));
  GPU_blend(GPU_BLEND_ALPHA);
  uint pos = GPU_vertformat_attr_add(immVertexFormat(), "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
  immBindBuiltinProgram(GPU_SHADER_3D_UNIFORM_COLOR);

  const float alpha = gizmo->is_mouse_over_gizmo ? 1.0f : 0.6f;

  immUniformColor4f(0.2f, 0.2f, 0.2f, alpha);
  imm_draw_circle_fill_2d(pos,
                          button.x + button.width / 2,
                          button.y + button.height / 2,
                          button.width * RETIME_BUTTON_SIZE,
                          32);
  immUnbindProgram();

  GPU_polygon_smooth(false);
  UI_icon_draw_alpha(button.x, button.y, gizmo->icon_id, alpha);
  GPU_polygon_smooth(true);

  GPU_blend(GPU_BLEND_NONE);
}

static int gizmo_retime_handle_add_test_select(bContext *C, wmGizmo *gz, const int mval[2])
{
  RetimeButtonGizmo *gizmo = (RetimeButtonGizmo *)gz;
  Sequence *seq = active_seq_from_context(C);

  Sequence *mouse_over_seq = NULL;
  gizmo->is_mouse_over_gizmo = false;

  /* Store strip under mouse cursor. */
  const rctf strip_box = strip_box_get(C, seq);
  if (mouse_is_inside_box(&strip_box, mval)) {
    mouse_over_seq = seq;
  }

  if (gizmo->seq_under_mouse != mouse_over_seq) {
    gizmo->seq_under_mouse = mouse_over_seq;
    WM_event_add_notifier(C, NC_SCENE | ND_SEQUENCER, CTX_data_scene(C));
  }

  if (gizmo->seq_under_mouse == NULL) {
    return -1;
  }

  const rctf button_box = button_box_get(C, gizmo);
  if (!mouse_is_inside_box(&button_box, mval)) {
    return -1;
  }

  gizmo->is_mouse_over_gizmo = true;
  const Scene *scene = CTX_data_scene(C);
  wmGizmoOpElem *op_elem = WM_gizmo_operator_get(gz, 0);
  RNA_int_set(&op_elem->ptr, "timeline_frame", BKE_scene_frame_get(scene));

  WM_event_add_notifier(C, NC_SCENE | ND_SEQUENCER, CTX_data_scene(C));
  return 0;
}

static void gizmo_retime_handle_add_setup(wmGizmo *gz)
{
  RetimeButtonGizmo *gizmo = (RetimeButtonGizmo *)gz;
  gizmo->icon_id = ICON_ADD;
}

void GIZMO_GT_retime_handle_add(wmGizmoType *gzt)
{
  /* Identifiers. */
  gzt->idname = "GIZMO_GT_retime_handle_add";

  /* Api callbacks. */
  gzt->setup = gizmo_retime_handle_add_setup;
  gzt->draw = gizmo_retime_handle_add_draw;
  gzt->test_select = gizmo_retime_handle_add_test_select;
  gzt->struct_size = sizeof(RetimeButtonGizmo);

  /* Currently only used for cursor display. */
  RNA_def_boolean(gzt->srna, "show_drag", true, "Show Drag", "");
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Retiming Move Handle Gizmo
 * \{ */

typedef struct RetimeHandleMoveGizmo {
  wmGizmo gizmo;
  const Sequence *mouse_over_seq;
  int mouse_over_handle_x;
} RetimeHandleMoveGizmo;

static void retime_handle_draw(const bContext *C,
                               const RetimeHandleMoveGizmo *gizmo,
                               uint pos,
                               const Sequence *seq,
                               RetimingHandle *handle)
{
  const Scene *scene = CTX_data_scene(C);
  if (handle->x() == SEQ_time_left_handle_frame_get(scene, seq) ||
      handle->index == 0) {  // xxx .is first_handle?
    return;                  /* Ignore first handle. */
  }

  const View2D *v2d = UI_view2d_fromcontext(C);
  const rctf strip_box = strip_box_get(C, seq);
  if (!BLI_rctf_isect_x(&strip_box,
                        UI_view2d_view_to_region_x(v2d, handle->x()))) {  // xxx .is_visible
    return; /* Handle out of strip bounds. */
  }

  const int ui_triangle_size = RETIME_HANDLE_TRIANGLE_SIZE;
  const float bottom = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 0.0f)) + 2;
  const float middle = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 0.5f));
  const float top = UI_view2d_view_to_region_y(v2d, strip_y_rescale(seq, 1.0f)) - 2;
  const float handle_position = UI_view2d_view_to_region_x(v2d, handle->x());
  const int right_handle_frame = SEQ_time_right_handle_frame_get(scene, seq);

  if (seq == gizmo->mouse_over_seq && handle->x() == gizmo->mouse_over_handle_x) {
    immUniformColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  }
  else {
    immUniformColor4f(0.65f, 0.65f, 0.65f, 1.0f);
  }

  immBegin(GPU_PRIM_TRI_FAN, 3);
  if (handle->x() == right_handle_frame) {
    /* Offset handle position by 1 frame, so it does not look out of place. */  // xxx hmmm could
                                                                                // this be part of
                                                                                // .x() function?
    immVertex2f(pos, handle_position - ui_triangle_size, middle - ui_triangle_size / 2);
    immVertex2f(pos, handle_position - ui_triangle_size, middle + ui_triangle_size / 2);
    immVertex2f(pos, handle_position, middle);
  }
  else {

    immVertex2f(pos, handle_position - ui_triangle_size / 2, bottom);
    immVertex2f(pos, handle_position + ui_triangle_size / 2, bottom);
    immVertex2f(pos, handle_position, bottom + ui_triangle_size);
  }

  immEnd();

  immBegin(GPU_PRIM_LINES, 2);
  immVertex2f(pos, handle_position, bottom);
  immVertex2f(pos, handle_position, top);
  immEnd();
}

static void retime_speed_text_draw(const bContext *C, const Sequence *seq, RetimingHandle *handle)
{
  if (handle->is_last_handle) {
    return;
  }

  const Scene *scene = CTX_data_scene(C);
  const int start_frame = SEQ_time_left_handle_frame_get(scene, seq);
  const int end_frame = SEQ_time_right_handle_frame_get(scene, seq);

  RetimingHandle next_handle = handle->next_handle_get();
  if (next_handle.x() < start_frame || handle->x() > end_frame) {
    return; /* Label out of strip bounds. */
  }

  const float speed = SEQ_retiming_handle_speed_get(scene, seq, next_handle.index);

  char label_str[20];
  const size_t label_len = BLI_snprintf_rlen(
      label_str, sizeof(label_str), "%d%%", round_fl_to_int(speed * 100.0f));

  const float width = pixels_to_view_width(C, BLF_width(BLF_default(), label_str, label_len));

  const float xmin = max_ff(SEQ_time_left_handle_frame_get(scene, seq), handle->x());
  const float xmax = min_ff(SEQ_time_right_handle_frame_get(scene, seq), next_handle.x());

  const float text_x = (xmin + xmax - width) / 2;
  const float text_y = strip_y_rescale(seq, 0) + pixels_to_view_height(C, 5);

  if (width > xmax - xmin) {
    return; /* Not enough space to draw label. */
  }

  const uchar col[4] = {255, 255, 255, 255};
  UI_view2d_text_cache_add(UI_view2d_fromcontext(C), text_x, text_y, label_str, label_len, col);
}

static void gizmo_retime_handle_draw(const bContext *C, wmGizmo *gz)
{
  const RetimeHandleMoveGizmo *gizmo = (RetimeHandleMoveGizmo *)gz;
  const View2D *v2d = UI_view2d_fromcontext(C);

  wmOrtho2_region_pixelspace(CTX_wm_region(C));
  GPU_blend(GPU_BLEND_ALPHA);
  uint pos = GPU_vertformat_attr_add(immVertexFormat(), "pos", GPU_COMP_F32, 2, GPU_FETCH_FLOAT);
  immBindBuiltinProgram(GPU_SHADER_3D_UNIFORM_COLOR);

  Sequence *seq = active_seq_from_context(C);
  SEQ_retiming_data_ensure(CTX_data_scene(C), seq);
  RetimingHandlesIterator handles = RetimingHandlesIterator<RetimingHandle>(seq);

  for (auto handle : handles) {
    retime_handle_draw(C, gizmo, pos, seq, &handle);
    retime_speed_text_draw(C, seq, &handle);
  }

  immUnbindProgram();
  GPU_blend(GPU_BLEND_NONE);

  UI_view2d_text_cache_draw(CTX_wm_region(C));
  UI_view2d_view_ortho(v2d); /* 'UI_view2d_text_cache_draw()' messes up current view. */
}

static int gizmo_retime_handle_test_select(bContext *C, wmGizmo *gz, const int mval[2])
{
  RetimeHandleMoveGizmo *gizmo = (RetimeHandleMoveGizmo *)gz;
  Scene *scene = CTX_data_scene(C);

  gizmo->mouse_over_seq = NULL;

  Sequence *seq = active_seq_from_context(C);
  SEQ_retiming_data_ensure(CTX_data_scene(C), seq);
  RetimingHandlesIterator handles = RetimingHandlesIterator<RetimingHandle>(seq);
  RetimingHandle *handle = handles.mouse_over_handle_get(UI_view2d_fromcontext(C), mval);

  if (handle == NULL) {
    return -1;
  }

  if (handle->x() == SEQ_time_left_handle_frame_get(scene, seq) || handle->index == 0) {
    return -1;
  }

  rctf strip_box = strip_box_get(C, seq);
  BLI_rctf_resize_x(&strip_box, BLI_rctf_size_x(&strip_box) + 2 * RETIME_HANDLE_TRIANGLE_SIZE);
  if (!mouse_is_inside_box(&strip_box, mval)) {
    return -1;
  }

  gizmo->mouse_over_seq = seq;
  gizmo->mouse_over_handle_x = handle->x();

  wmGizmoOpElem *op_elem = WM_gizmo_operator_get(gz, 0);
  RNA_int_set(&op_elem->ptr, "handle_index", handle->index);

  WM_event_add_notifier(C, NC_SCENE | ND_SEQUENCER, scene);
  return 0;
}

static int gizmo_retime_handle_cursor_get(wmGizmo *gz)
{
  if (RNA_boolean_get(gz->ptr, "show_drag")) {
    return WM_CURSOR_EW_SCROLL;
  }
  return WM_CURSOR_DEFAULT;
}

static void gizmo_retime_handle_setup(wmGizmo *gz)
{
  gz->flag = WM_GIZMO_DRAW_MODAL;
}

void GIZMO_GT_retime_handle(wmGizmoType *gzt)
{
  /* Identifiers. */
  gzt->idname = "GIZMO_GT_retime_handle_move";

  /* Api callbacks. */
  gzt->setup = gizmo_retime_handle_setup;
  gzt->draw = gizmo_retime_handle_draw;
  gzt->test_select = gizmo_retime_handle_test_select;
  gzt->cursor_get = gizmo_retime_handle_cursor_get;
  gzt->struct_size = sizeof(RetimeHandleMoveGizmo);

  /* Currently only used for cursor display. */
  RNA_def_boolean(gzt->srna, "show_drag", true, "Show Drag", "");
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Retiming Remove Handle Gizmo
 * \{ */

static void gizmo_retime_remove_draw(const bContext *UNUSED(C), wmGizmo *UNUSED(gz))
{
  /* Pass. */
}

static int gizmo_retime_remove_test_select(bContext *C, wmGizmo *gz, const int mval[2])
{
  Scene *scene = CTX_data_scene(C);
  Sequence *seq = active_seq_from_context(C);

  SEQ_retiming_data_ensure(CTX_data_scene(C), seq);
  RetimingHandlesIterator handles = RetimingHandlesIterator<RetimingHandle>(seq);
  RetimingHandle *handle = handles.mouse_over_handle_get(UI_view2d_fromcontext(C), mval);

  if (handle == NULL) {
    return -1;
  }

  if (handle->x() == SEQ_time_left_handle_frame_get(scene, seq) || handle->index == 0) {
    return -1; /* Ignore first handle. */
  }

  if (handle->is_last_handle) {
    return -1; /* Last handle can not be removed. */
  }

  rctf box = remove_box_get(C, seq);

  BLI_rctf_resize_x(&box, BLI_rctf_size_x(&box) + 2 * RETIME_HANDLE_TRIANGLE_SIZE);
  if (!mouse_is_inside_box(&box, mval)) {
    return -1;
  }

  wmGizmoOpElem *op_elem = WM_gizmo_operator_get(gz, 0);
  RNA_int_set(&op_elem->ptr, "handle_index", handle->index);

  WM_event_add_notifier(C, NC_SCENE | ND_SEQUENCER, scene);
  return 0;
}

static int gizmo_retime_remove_cursor_get(wmGizmo *gz)
{
  if (RNA_boolean_get(gz->ptr, "show_drag")) {
    return WM_CURSOR_ERASER;
  }
  return WM_CURSOR_DEFAULT;
}

void GIZMO_GT_retime_remove(wmGizmoType *gzt)
{
  /* Identifiers. */
  gzt->idname = "GIZMO_GT_retime_handle_remove";

  /* Api callbacks. */
  gzt->draw = gizmo_retime_remove_draw;
  gzt->test_select = gizmo_retime_remove_test_select;
  gzt->cursor_get = gizmo_retime_remove_cursor_get;
  gzt->struct_size = sizeof(wmGizmo);

  /* Currently only used for cursor display. */
  RNA_def_boolean(gzt->srna, "show_drag", true, "Show Drag", "");
}

/** \} */
