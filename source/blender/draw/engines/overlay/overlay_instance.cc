/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#include "overlay_instance.hh"
#include "draw_debug.hh"

namespace blender::draw::overlay {

void Instance::init()
{
  resources.depth_tx.wrap(DRW_viewport_texture_list_get()->depth);
  resources.color_overlay_tx.wrap(DRW_viewport_texture_list_get()->color_overlay);
  resources.color_render_tx.wrap(DRW_viewport_texture_list_get()->color);

  /* TODO(fclem): Remove DRW global usage. */
  const DRWContextState *ctx = DRW_context_state_get();

  state.depsgraph = ctx->depsgraph;
  state.scene = ctx->scene;
  state.v3d = ctx->v3d;
  state.rv3d = ctx->rv3d;

  state.pixelsize = U.pixelsize;
  state.ctx_mode = CTX_data_mode_enum_ex(ctx->object_edit, ctx->obact, ctx->object_mode);
  state.clear_in_front = (state.v3d->shading.type != OB_SOLID);
  state.use_in_front = (state.v3d->shading.type <= OB_SOLID) ||
                       BKE_scene_uses_blender_workbench(state.scene);
  state.is_wireframe_mode = (state.v3d->shading.type == OB_WIRE);
  state.hide_overlays = (state.v3d->flag2 & V3D_HIDE_OVERLAYS) != 0;
  state.xray_enabled = XRAY_ACTIVE(state.v3d);
  state.xray_enabled_and_not_wire = state.xray_enabled && (state.v3d->shading.type > OB_WIRE);
  state.xray_opacity = XRAY_ALPHA(state.v3d);
  state.cfra = DEG_get_ctime(state.depsgraph);
  state.clipping_state = RV3D_CLIPPING_ENABLED(state.v3d, state.rv3d) ? DRW_STATE_CLIP_PLANES :
                                                                        DRWState(0);

  if (!state.hide_overlays) {
    state.overlay = state.v3d->overlay;
    state.v3d_flag = state.v3d->flag;
    state.v3d_gridflag = state.v3d->gridflag;
  }
  else {
    memset(&state.overlay, 0, sizeof(state.overlay));
    state.v3d_flag = 0;
    state.v3d_gridflag = 0;
    state.overlay.flag = V3D_OVERLAY_HIDE_TEXT | V3D_OVERLAY_HIDE_MOTION_PATHS |
                         V3D_OVERLAY_HIDE_BONES | V3D_OVERLAY_HIDE_OBJECT_XTRAS |
                         V3D_OVERLAY_HIDE_OBJECT_ORIGINS;
    state.overlay.wireframe_threshold = state.v3d->overlay.wireframe_threshold;
    state.overlay.wireframe_opacity = state.v3d->overlay.wireframe_opacity;
  }

  /* TODO(fclem): Remove DRW global usage. */
  resources.globals_buf = G_draw.block_ubo;
  resources.theme_settings = G_draw.block;
}

void Instance::begin_sync()
{
  const DRWView *view_legacy = DRW_view_default_get();
  View view("OverlayView", view_legacy);

  background.begin_sync(resources, state);
  grid.begin_sync(resources, state, view);
}

void Instance::object_sync(ObjectRef &ob_ref)
{
  UNUSED_VARS(ob_ref);
}

void Instance::end_sync()
{
}

void Instance::draw(Manager &manager)
{
  /* WORKAROUND: This is to prevent crashes when using depth picking or selection.
   * The selection engine should handle theses cases instead. */
  if (!DRW_state_is_fbo()) {
    return;
  }

  const DRWView *view_legacy = DRW_view_default_get();
  View view("OverlayView", view_legacy);

  resources.line_tx.acquire(int2(resources.depth_tx.size()), GPU_RGBA8);

  resources.overlay_color_only_fb.ensure(GPU_ATTACHMENT_NONE,
                                         GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx));
  resources.overlay_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx),
                              GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx));
  resources.overlay_line_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx),
                                   GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx),
                                   GPU_ATTACHMENT_TEXTURE(resources.line_tx));

  GPU_framebuffer_bind(resources.overlay_color_only_fb);

  float4 clear_color(0.0f);
  GPU_framebuffer_clear_color(resources.overlay_color_only_fb, clear_color);

  background.draw(resources, manager);
  grid.draw(resources, manager, view);

  // anti_aliasing.draw(resources, manager, view);

  resources.line_tx.release();
}

}  // namespace blender::draw::overlay
