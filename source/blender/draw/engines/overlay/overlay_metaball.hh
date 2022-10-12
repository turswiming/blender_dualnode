/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#pragma once

#include "DEG_depsgraph_query.h"
#include "DNA_camera_types.h"
#include "DNA_space_types.h"
#include "ED_view3d.h"
#include "UI_resources.h"

#include "draw_cache.h"
#include "draw_pass.hh"
#include "overlay_private.hh"
#include "overlay_shader_shared.h"

namespace blender::draw::overlay {

class Metaballs {

 private:
  PassSimple metaball_ps_ = {"MetaBalls"};
  PassSimple metaball_in_front_ps_ = {"MetaBalls_In_front"};

  ArmatureSphereBuf data_buf_ = {"metaball_data_buf"};
  ArmatureSphereBuf data_in_front_buf_ = {"metaball_data_buf"};

 public:
  void begin_sync()
  {
    data_buf_.clear();
    data_in_front_buf_.clear();
  }

  void metaball_instance_data_set(BoneInstanceData *data,
                                  Object *ob,
                                  const float *pos,
                                  const float radius,
                                  const float color[4])
  {
    /* Bone point radius is 0.05. Compensate for that. */
    mul_v3_v3fl(data->mat[0], ob->obmat[0], radius / 0.05f);
    mul_v3_v3fl(data->mat[1], ob->obmat[1], radius / 0.05f);
    mul_v3_v3fl(data->mat[2], ob->obmat[2], radius / 0.05f);
    mul_v3_m4v3(data->mat[3], ob->obmat, pos);
    /* WATCH: Reminder, alpha is wire-size. */
    OVERLAY_bone_instance_data_set_color(data, color);
  }

  void edit_object_sync(const ObjectRef &ob_ref, const Resources &res)
  {
    ArmatureSphereBuf &data_buf = (ob_ref.object->dtx & OB_DRAW_IN_FRONT) != 0 ?
                                      data_in_front_buf_ :
                                      data_buf_;
    MetaBall *mb = static_cast<MetaBall *>(ob_ref.object->data);

    const float *color;
    const float *col_radius = res.theme_settings.color_mball_radius;
    const float *col_radius_select = res.theme_settings.color_mball_radius_select;
    const float *col_stiffness = res.theme_settings.color_mball_stiffness;
    const float *col_stiffness_select = res.theme_settings.color_mball_stiffness_select;

    LISTBASE_FOREACH (MetaElem *, ml, mb->editelems) {
      const bool is_selected = (ml->flag & SELECT) != 0;
      const bool is_scale_radius = (ml->flag & MB_SCALE_RAD) != 0;
      float stiffness_radius = ml->rad * atanf(ml->s) / float(M_PI_2);
      BoneInstanceData instdata;

      color = (is_selected && is_scale_radius) ? col_radius_select : col_radius;
      metaball_instance_data_set(&instdata, ob_ref.object, &ml->x, ml->rad, color);
      data_buf.append(*reinterpret_cast<float4x4 *>(&instdata));

      color = (is_selected && !is_scale_radius) ? col_stiffness_select : col_stiffness;
      metaball_instance_data_set(&instdata, ob_ref.object, &ml->x, stiffness_radius, color);
      data_buf.append(*reinterpret_cast<float4x4 *>(&instdata));
    }
  }

  void object_sync(const ObjectRef &ob_ref, const Resources &res, const State &state)
  {
    ArmatureSphereBuf &data_buf = (ob_ref.object->dtx & OB_DRAW_IN_FRONT) != 0 ?
                                      data_in_front_buf_ :
                                      data_buf_;
    MetaBall *mb = static_cast<MetaBall *>(ob_ref.object->data);

    float *color;
    /* TODO(fclem): Remove DRW global usage. */
    UNUSED_VARS(res);
    DRW_object_wire_theme_get(ob_ref.object, state.view_layer, &color);

    LISTBASE_FOREACH (MetaElem *, ml, &mb->elems) {
      /* Draw radius only. */
      BoneInstanceData instdata;
      metaball_instance_data_set(&instdata, ob_ref.object, &ml->x, ml->rad, color);
      data_buf.append(*reinterpret_cast<float4x4 *>(&instdata));
    }
  }

  void end_sync(Resources &res, const State &state)
  {
    auto init_pass = [&](PassSimple &pass, ArmatureSphereBuf &data_buf) {
      data_buf.push_update();

      pass.init();
      pass.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     state.clipping_state);
      pass.shader_set(OVERLAY_shader_armature_sphere(true));
      pass.bind_ubo("globalsBlock", &res.globals_buf);
      pass.bind_ssbo("data_buf", &data_buf);
      pass.draw(DRW_cache_bone_point_wire_outline_get(), data_buf.size());
    };
    init_pass(metaball_ps_, data_buf_);
    init_pass(metaball_in_front_ps_, data_in_front_buf_);
  }

  void draw(Resources &res, Manager &manager, View &view)
  {
    GPU_framebuffer_bind(res.overlay_line_fb);
    manager.submit(metaball_ps_, view);
  }

  void draw_in_front(Resources &res, Manager &manager, View &view)
  {
    GPU_framebuffer_bind(res.overlay_line_in_front_fb);
    manager.submit(metaball_in_front_ps_, view);
  }
};

}  // namespace blender::draw::overlay
