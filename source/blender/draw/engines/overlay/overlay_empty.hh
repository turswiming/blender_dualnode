/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#pragma once

#include "DNA_object_types.h"

#include "draw_pass.hh"
#include "draw_shader_shared.h"
#include "overlay_private.hh"
#include "overlay_shape.hh"

namespace blender::draw::overlay {

class Empties {

 private:
  PassSimple empty_ps_ = {"Empties"};
  PassSimple empty_in_front_ps_ = {"Empties_In_front"};

  using EmptyInstanceBuf = StorageVectorBuffer<ExtraInstanceData>;

  struct CallBuffers {
    EmptyInstanceBuf plain_axes_buf = {"plain_axes_buf"};
    EmptyInstanceBuf single_arrow_buf = {"single_arrow_buf"};
    EmptyInstanceBuf cube_buf = {"cube_buf"};
    EmptyInstanceBuf circle_buf = {"circle_buf"};
    EmptyInstanceBuf sphere_buf = {"sphere_buf"};
    EmptyInstanceBuf cone_buf = {"cone_buf"};
    EmptyInstanceBuf arrows_buf = {"arrows_buf"};
    EmptyInstanceBuf image_buf = {"image_buf"};
  } call_buffers_[2];

 public:
  void begin_sync()
  {
    for (int i = 0; i < 2; i++) {
      call_buffers_[i].plain_axes_buf.clear();
      call_buffers_[i].single_arrow_buf.clear();
      call_buffers_[i].cube_buf.clear();
      call_buffers_[i].circle_buf.clear();
      call_buffers_[i].sphere_buf.clear();
      call_buffers_[i].cone_buf.clear();
      call_buffers_[i].arrows_buf.clear();
      call_buffers_[i].image_buf.clear();
    }
  }

  void object_sync(const ObjectRef &ob_ref, const Resources &res, const State &state)
  {
    CallBuffers &call_bufs = call_buffers_[int((ob_ref.object->dtx & OB_DRAW_IN_FRONT) != 0)];

    float4 color = res.object_wire_color(ob_ref, state);
    ExtraInstanceData data(
        float4x4(ob_ref.object->object_to_world), color, ob_ref.object->empty_drawsize);

    switch (ob_ref.object->empty_drawtype) {
      case OB_PLAINAXES:
        call_bufs.plain_axes_buf.append(data);
        break;
      case OB_SINGLE_ARROW:
        call_bufs.single_arrow_buf.append(data);
        break;
      case OB_CUBE:
        call_bufs.cube_buf.append(data);
        break;
      case OB_CIRCLE:
        call_bufs.circle_buf.append(data);
        break;
      case OB_EMPTY_SPHERE:
        call_bufs.sphere_buf.append(data);
        break;
      case OB_EMPTY_CONE:
        call_bufs.cone_buf.append(data);
        break;
      case OB_ARROWS:
        call_bufs.arrows_buf.append(data);
        break;
      case OB_EMPTY_IMAGE:
        /* This only show the frame. See OVERLAY_image_empty_cache_populate() for the image. */
        call_bufs.image_buf.append(data);
        break;
    }
  }

  void end_sync(Resources &res, ShapeCache &shapes, const State &state)
  {
    auto init_pass = [&](PassSimple &pass, CallBuffers &call_bufs) {
      pass.init();
      pass.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     state.clipping_state);
      pass.shader_set(OVERLAY_shader_extra(false));
      pass.bind_ubo("globalsBlock", &res.globals_buf);

      call_bufs.plain_axes_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.plain_axes_buf);
      pass.draw(shapes.plain_axes, call_bufs.plain_axes_buf.size());

      call_bufs.single_arrow_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.single_arrow_buf);
      pass.draw(shapes.single_arrow, call_bufs.single_arrow_buf.size());

      call_bufs.cube_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.cube_buf);
      pass.draw(shapes.cube, call_bufs.cube_buf.size());

      call_bufs.circle_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.circle_buf);
      pass.draw(shapes.circle, call_bufs.circle_buf.size());

      call_bufs.sphere_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.sphere_buf);
      pass.draw(shapes.empty_sphere, call_bufs.sphere_buf.size());

      call_bufs.cone_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.cone_buf);
      pass.draw(shapes.empty_cone, call_bufs.cone_buf.size());

      call_bufs.arrows_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.arrows_buf);
      pass.draw(shapes.arrows, call_bufs.arrows_buf.size());

      call_bufs.image_buf.push_update();
      pass.bind_ssbo("data_buf", &call_bufs.image_buf);
      pass.draw(shapes.quad_wire, call_bufs.image_buf.size());
    };
    init_pass(empty_ps_, call_buffers_[0]);
    init_pass(empty_in_front_ps_, call_buffers_[1]);
  }

  void draw(Resources &res, Manager &manager, View &view)
  {
    GPU_framebuffer_bind(res.overlay_line_fb);
    manager.submit(empty_ps_, view);
  }

  void draw_in_front(Resources &res, Manager &manager, View &view)
  {
    GPU_framebuffer_bind(res.overlay_line_in_front_fb);
    manager.submit(empty_in_front_ps_, view);
  }
};

}  // namespace blender::draw::overlay
