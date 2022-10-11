/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#pragma once

#include "DNA_space_types.h"
#include "ED_view3d.h"

#include "draw_cache.h"
#include "draw_pass.hh"
#include "overlay_private.hh"
#include "overlay_shader_shared.h"

namespace blender::draw::overlay {

class Grid {
 private:
  UniformBuffer<OVERLAY_GridData> data_;

  bool enabled = false;

  PassSimple grid_ps_ = {"grid_ps_"};

  float3 grid_axes = float3(0.0f);
  float3 zplane_axes = float3(0.0f);
  OVERLAY_GridBits grid_flag_;
  OVERLAY_GridBits zneg_flag_;
  OVERLAY_GridBits zpos_flag_;

 public:
  void update_ubo()
  {
    float grid_steps[SI_GRID_STEPS_LEN] = {
        0.001f, 0.01f, 0.1f, 1.0f, 10.0f, 100.0f, 1000.0f, 10000.0f};
    float grid_steps_y[SI_GRID_STEPS_LEN] = {0.0f}; /* When zero, use value from grid_steps. */
    data_.line_size = max_ff(0.0f, U.pixelsize - 1.0f) * 0.5f;
    /* Default, nothing is drawn. */
    grid_flag_ = zneg_flag_ = zpos_flag_ = OVERLAY_GridBits(0);

    /* SPACE_VIEW3D */
    Scene *scene = DRW_context_state_get()->scene;
    View3D *v3d = DRW_context_state_get()->v3d;
    RegionView3D *rv3d = DRW_context_state_get()->rv3d;

    const bool show_axis_x = true;      //(pd->v3d_gridflag & V3D_SHOW_X) != 0;
    const bool show_axis_y = true;      //(pd->v3d_gridflag & V3D_SHOW_Y) != 0;
    const bool show_axis_z = true;      //(pd->v3d_gridflag & V3D_SHOW_Z) != 0;
    const bool show_floor = true;       //(pd->v3d_gridflag & V3D_SHOW_FLOOR) != 0;
    const bool show_ortho_grid = true;  //(pd->v3d_gridflag & V3D_SHOW_ORTHO_GRID) != 0;

    // if (pd->hide_overlays || !(pd->v3d_gridflag & (V3D_SHOW_X | V3D_SHOW_Y | V3D_SHOW_Z |
    //                                                V3D_SHOW_FLOOR | V3D_SHOW_ORTHO_GRID))) {
    //   return;
    // }

    float viewinv[4][4], wininv[4][4];
    float viewmat[4][4], winmat[4][4];
    DRW_view_winmat_get(nullptr, winmat, false);
    DRW_view_winmat_get(nullptr, wininv, true);
    DRW_view_viewmat_get(nullptr, viewmat, false);
    DRW_view_viewmat_get(nullptr, viewinv, true);

    /* If perspective view or non-axis aligned view. */
    if (winmat[3][3] == 0.0f || rv3d->view == RV3D_VIEW_USER) {
      if (show_axis_x) {
        grid_flag_ |= PLANE_XY | SHOW_AXIS_X;
      }
      if (show_axis_y) {
        grid_flag_ |= PLANE_XY | SHOW_AXIS_Y;
      }
      if (show_floor) {
        grid_flag_ |= PLANE_XY | SHOW_GRID;
      }
    }
    else {
      if (show_ortho_grid && ELEM(rv3d->view, RV3D_VIEW_RIGHT, RV3D_VIEW_LEFT)) {
        grid_flag_ = PLANE_YZ | SHOW_AXIS_Y | SHOW_AXIS_Z | SHOW_GRID | GRID_BACK;
      }
      else if (show_ortho_grid && ELEM(rv3d->view, RV3D_VIEW_TOP, RV3D_VIEW_BOTTOM)) {
        grid_flag_ = PLANE_XY | SHOW_AXIS_X | SHOW_AXIS_Y | SHOW_GRID | GRID_BACK;
      }
      else if (show_ortho_grid && ELEM(rv3d->view, RV3D_VIEW_FRONT, RV3D_VIEW_BACK)) {
        grid_flag_ = PLANE_XZ | SHOW_AXIS_X | SHOW_AXIS_Z | SHOW_GRID | GRID_BACK;
      }
    }

    grid_axes[0] = float((grid_flag_ & (PLANE_XZ | PLANE_XY)) != 0);
    grid_axes[1] = float((grid_flag_ & (PLANE_YZ | PLANE_XY)) != 0);
    grid_axes[2] = float((grid_flag_ & (PLANE_YZ | PLANE_XZ)) != 0);

    /* Z axis if needed */
    if (((rv3d->view == RV3D_VIEW_USER) || (rv3d->persp != RV3D_ORTHO)) && show_axis_z) {
      zpos_flag_ = SHOW_AXIS_Z;

      float zvec[3], campos[3];
      negate_v3_v3(zvec, viewinv[2]);
      copy_v3_v3(campos, viewinv[3]);

      /* z axis : chose the most facing plane */
      if (fabsf(zvec[0]) < fabsf(zvec[1])) {
        zpos_flag_ |= PLANE_XZ;
      }
      else {
        zpos_flag_ |= PLANE_YZ;
      }

      zneg_flag_ = zpos_flag_;

      /* Perspective: If camera is below floor plane, we switch clipping.
       * Orthographic: If eye vector is looking up, we switch clipping. */
      if (((winmat[3][3] == 0.0f) && (campos[2] > 0.0f)) ||
          ((winmat[3][3] != 0.0f) && (zvec[2] < 0.0f))) {
        zpos_flag_ |= CLIP_ZPOS;
        zneg_flag_ |= CLIP_ZNEG;
      }
      else {
        zpos_flag_ |= CLIP_ZNEG;
        zneg_flag_ |= CLIP_ZPOS;
      }

      zplane_axes[0] = float((zpos_flag_ & (PLANE_XZ | PLANE_XY)) != 0);
      zplane_axes[1] = float((zpos_flag_ & (PLANE_YZ | PLANE_XY)) != 0);
      zplane_axes[2] = float((zpos_flag_ & (PLANE_YZ | PLANE_XZ)) != 0);
    }
    else {
      zneg_flag_ = zpos_flag_ = CLIP_ZNEG | CLIP_ZPOS;
    }

    float dist;
    // if (rv3d->persp == RV3D_CAMOB && v3d->camera && v3d->camera->type == OB_CAMERA) {
    //   Object *camera_object = DEG_get_evaluated_object(draw_ctx->depsgraph, v3d->camera);
    //   dist = ((Camera *)(camera_object->data))->clip_end;
    //   grid_flag_ |= GRID_CAMERA;
    //   zneg_flag_ |= GRID_CAMERA;
    //   zpos_flag_ |= GRID_CAMERA;
    // }
    // else {
    dist = v3d->clip_end;
    // }

    if (winmat[3][3] == 0.0f) {
      copy_v3_fl(data_.size, dist);
    }
    else {
      float viewdist = 1.0f / min_ff(fabsf(winmat[0][0]), fabsf(winmat[1][1]));
      copy_v3_fl(data_.size, viewdist * dist);
    }

    data_.distance = dist / 2.0f;

    ED_view3d_grid_steps(scene, v3d, rv3d, grid_steps);

    if ((v3d->flag & (V3D_XR_SESSION_SURFACE | V3D_XR_SESSION_MIRROR)) != 0) {
      /* The calculations for the grid parameters assume that the view matrix has no scale
       * component, which may not be correct if the user is "shrunk" or "enlarged" by zooming in or
       * out. Therefore, we need to compensate the values here. */
      /* Assumption is uniform scaling (all column vectors are of same length). */
      float viewinvscale = len_v3(viewinv[0]);
      data_.distance *= viewinvscale;
    }

    /* Convert to UBO alignment. */
    for (int i = 0; i < SI_GRID_STEPS_LEN; i++) {
      data_.steps[i][0] = grid_steps[i];
      data_.steps[i][1] = (grid_steps_y[i] != 0.0f) ? grid_steps_y[i] : grid_steps[i];
    }

    data_.push_update();
  }

  void begin_sync()
  {
    update_ubo();

    grid_ps_.init();
    grid_ps_.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_BLEND_ALPHA);

    grid_ps_.shader_set(OVERLAY_shader_grid());
    grid_ps_.bind_ubo("grid_buf", &data_);
    grid_ps_.bind_ubo("globalsBlock", &G_draw.block_ubo);
    grid_ps_.bind_texture("depth_tx", &DRW_viewport_texture_list_get()->depth);

    grid_ps_.push_constant("grid_flag", grid_flag_);
    grid_ps_.push_constant("plane_axes", grid_axes);
    grid_ps_.clear_color(float4(0.0f, 0.0f, 0.0f, 0.0f));
    grid_ps_.draw(DRW_cache_grid_get());
  }

  void draw(Manager &manager, View &view)
  {
    manager.submit(grid_ps_, view);
  }
};

}  // namespace blender::draw::overlay
