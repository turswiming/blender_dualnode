/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup overlay
 */

#include "overlay_instance.hh"

namespace blender::draw::overlay {

void Instance::init()
{
  GPUTexture *viewport_depth_tx = DRW_viewport_texture_list_get()->depth;
  GPUTexture *viewport_color_tx = DRW_viewport_texture_list_get()->color_overlay;
  overlay_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(viewport_color_tx));
}

void Instance::begin_sync()
{
  grid.begin_sync();
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
  const DRWView *view_old = DRW_view_default_get();
  View view("OverlayView", view_old);

  GPU_framebuffer_bind(overlay_fb);
  grid.draw(manager, view);
}

}  // namespace blender::draw::overlay
