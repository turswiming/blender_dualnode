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
  resources.depth_in_front_tx.wrap(DRW_viewport_texture_list_get()->depth_in_front);
  resources.color_overlay_tx.wrap(DRW_viewport_texture_list_get()->color_overlay);
  resources.color_render_tx.wrap(DRW_viewport_texture_list_get()->color);

  /* TODO(fclem): Remove DRW global usage. */
  const DRWContextState *ctx = DRW_context_state_get();

  state.depsgraph = ctx->depsgraph;
  state.view_layer = ctx->view_layer;
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
  metaballs.begin_sync();
  grid.begin_sync(resources, state, view);
}

void Instance::object_sync(ObjectRef &ob_ref)
{
  const bool in_edit_mode = object_is_edit_mode(ob_ref.object);

  if (in_edit_mode && !state.hide_overlays) {
    switch (ob_ref.object->type) {
      case OB_MESH:
        break;
      case OB_ARMATURE:
        break;
      case OB_CURVES_LEGACY:
        break;
      case OB_SURF:
        break;
      case OB_LATTICE:
        break;
      case OB_MBALL:
        metaballs.edit_object_sync(ob_ref, resources);
        break;
      case OB_FONT:
        break;
      case OB_CURVES:
        break;
    }
  }

  if (!state.hide_overlays) {
    switch (ob_ref.object->type) {
      case OB_ARMATURE:
        break;
      case OB_MBALL:
        if (!in_edit_mode) {
          metaballs.object_sync(ob_ref, resources, state);
        }
        break;
      case OB_GPENCIL:
        break;
    }
  }
}

void Instance::end_sync()
{
  metaballs.end_sync(resources, state);
}

void Instance::draw(Manager &manager)
{
  /* WORKAROUND: This is to prevent crashes when using depth picking or selection.
   * The selection engine should handle theses cases instead. */
  if (!DRW_state_is_fbo()) {
    return;
  }

  int2 render_size = int2(resources.depth_tx.size());

  const DRWView *view_legacy = DRW_view_default_get();
  View view("OverlayView", view_legacy);

  resources.line_tx.acquire(render_size, GPU_RGBA8);

  resources.overlay_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx),
                              GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx));
  resources.overlay_line_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx),
                                   GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx),
                                   GPU_ATTACHMENT_TEXTURE(resources.line_tx));
  resources.overlay_color_only_fb.ensure(GPU_ATTACHMENT_NONE,
                                         GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx));

  /* TODO(fclem): Remove mandatory allocation. */
  if (!resources.depth_in_front_tx.is_valid()) {
    resources.depth_in_front_alloc_tx.acquire(render_size, GPU_DEPTH_COMPONENT24);
    resources.depth_in_front_tx.wrap(resources.depth_in_front_alloc_tx);
  }

  resources.overlay_in_front_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_in_front_tx),
                                       GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx));
  resources.overlay_line_in_front_fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_in_front_tx),
                                            GPU_ATTACHMENT_TEXTURE(resources.color_overlay_tx),
                                            GPU_ATTACHMENT_TEXTURE(resources.line_tx));

  GPU_framebuffer_bind(resources.overlay_color_only_fb);

  float4 clear_color(0.0f);
  GPU_framebuffer_clear_color(resources.overlay_color_only_fb, clear_color);

  background.draw(resources, manager);

  metaballs.draw(resources, manager, view);

  grid.draw(resources, manager, view);

  metaballs.draw_in_front(resources, manager, view);

  // anti_aliasing.draw(resources, manager, view);

  resources.line_tx.release();
  resources.depth_in_front_alloc_tx.release();
}

}  // namespace blender::draw::overlay
