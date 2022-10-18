/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

namespace blender::workbench {

MeshPass::MeshPass(const char *name) : PassMain(name){};

/* Move to draw::Pass */
bool MeshPass::is_empty() const
{
  return _is_empty;
}

void MeshPass::init_pass(SceneResources &resources, DRWState state)
{
  _is_empty = true;
  PassMain::init();
  state_set(state);
  bind_texture(WB_MATCAP_SLOT, resources.matcap_tx);
  bind_ssbo(WB_MATERIAL_SLOT, &resources.material_buf);
  bind_ubo(WB_WORLD_SLOT, resources.world_buf);
}

void MeshPass::init_subpasses(ePipelineType pipeline, eShadingType shading, ShaderCache &shaders)
{
  texture_subpass_map.clear();

  for (auto geom : IndexRange(geometry_type_len)) {
    for (auto color : IndexRange(color_type_len)) {
      eGeometryType geom_type = static_cast<eGeometryType>(geom);
      eColorType color_type = static_cast<eColorType>(color);
      std::string name = std::string(get_name(geom_type)) + std::string(get_name(color_type));
      GPUShader *sh = shaders.prepass_shader_get(pipeline, geom_type, color_type, shading);
      PassMain::Sub *pass = &sub(name.c_str());
      pass->shader_set(sh);
      passes_[geom][color] = pass;
    }
  }
}

PassMain::Sub &MeshPass::sub_pass_get(ObjectRef &ref,
                                      ::Image *image /* = nullptr */,
                                      eGPUSamplerState sampler_state /* = GPU_SAMPLER_DEFAULT */,
                                      ImageUser *iuser /* = nullptr */)
{
  /*TODO(pragma37): For now we assume retrieving a subpass means it's not empty anymore*/
  _is_empty = false;

  eGeometryType geometry_type = geometry_type_from_object(ref.object);
  if (image) {
    GPUTexture *texture = nullptr;
    GPUTexture *tilemap = nullptr;
    if (image->source == IMA_SRC_TILED) {
      texture = BKE_image_get_gpu_tiles(image, iuser, nullptr);
      tilemap = BKE_image_get_gpu_tilemap(image, iuser, nullptr);
    }
    else {
      texture = BKE_image_get_gpu_texture(image, iuser, nullptr);
    }
    /* TODO(pragma37): Should be lib.name + name ??? */
    StringRefNull name = image->id.name;
    if (texture) {
      auto add_cb = [&] {
        PassMain::Sub *sub_pass =
            passes_[static_cast<int>(geometry_type)][static_cast<int>(eColorType::TEXTURE)];
        sub_pass = &sub_pass->sub(name.c_str());
        if (tilemap) {
          sub_pass->bind_texture(WB_TILE_ARRAY_SLOT, texture, sampler_state);
          sub_pass->bind_texture(WB_TILE_DATA_SLOT, tilemap);
        }
        else {
          sub_pass->bind_texture(WB_TEXTURE_SLOT, texture, sampler_state);
        }
        sub_pass->push_constant("isImageTile", tilemap != nullptr);
        sub_pass->push_constant("imagePremult", image && image->alpha_mode == IMA_ALPHA_PREMUL);
        /*TODO(pragma37): What's the point? This could be a constant in the shader. */
        sub_pass->push_constant("imageTransparencyCutoff", 0.1f);
        return sub_pass;
      };

      return *texture_subpass_map.lookup_or_add_cb(TextureSubPassKey(texture, geometry_type),
                                                   add_cb);
    }
  }
  return *passes_[static_cast<int>(geometry_type)][static_cast<int>(eColorType::MATERIAL)];
}

OpaquePass::OpaquePass(Texture &color_tx, Texture &depth_tx)
    : color_tx(color_tx), depth_tx(depth_tx){};

void OpaquePass::sync(DRWState cull_state,
                      DRWState clip_state,
                      eShadingType shading_type,
                      SceneResources &resources)
{
  TextureFromPool &color_tx = resources.color_tx;
  ShaderCache &shaders = resources.shader_cache;
  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                   cull_state | clip_state;

  gbuffer_ps_.init_pass(resources, state);
  gbuffer_ps_.init_subpasses(ePipelineType::OPAQUE, shading_type, resources.shader_cache);

  deferred_ps_.init();
  deferred_ps_.shader_set(shaders.resolve_shader_get(ePipelineType::OPAQUE, shading_type));
  deferred_ps_.bind_ubo(WB_WORLD_SLOT, resources.world_buf);
  deferred_ps_.bind_texture(WB_MATCAP_SLOT, resources.matcap_tx);
  deferred_ps_.bind_texture("normal_tx", &gbuffer_normal_tx);
  deferred_ps_.bind_texture("material_tx", &gbuffer_material_tx);
  deferred_ps_.bind_texture("depth_tx", &depth_tx);
  deferred_ps_.bind_image("out_color_img", &color_tx);
  deferred_ps_.dispatch(math::divide_ceil(int2(depth_tx.size()), int2(WB_RESOLVE_GROUP_SIZE)));
  deferred_ps_.barrier(GPU_BARRIER_TEXTURE_FETCH);
}

void OpaquePass::draw(Manager &manager, View &view)
{
  if (is_empty()) {
    return;
  }

  gbuffer_material_tx.acquire(int2(depth_tx.size()), GPU_RGBA16F);
  gbuffer_normal_tx.acquire(int2(depth_tx.size()), GPU_RG16F);
  gbuffer_object_id_tx.acquire(int2(depth_tx.size()), GPU_R16UI);

  opaque_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_material_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_normal_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_object_id_tx));
  opaque_fb.bind();

  manager.submit(gbuffer_ps_, view);

  manager.submit(deferred_ps_, view);

  gbuffer_normal_tx.release();
  gbuffer_material_tx.release();
  gbuffer_object_id_tx.release();
}

bool OpaquePass::is_empty() const
{
  return gbuffer_ps_.is_empty();
}

TransparentPass::TransparentPass(Texture &color_tx, Texture &depth_tx)
    : color_tx(color_tx), depth_tx(depth_tx){};

void TransparentPass::sync(DRWState cull_state,
                           DRWState clip_state,
                           eShadingType shading_type,
                           SceneResources &resources)
{
  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_DEPTH_LESS_EQUAL | DRW_STATE_BLEND_OIT |
                   cull_state | clip_state;

  accumulation_ps_.init_pass(resources, state);
  accumulation_ps_.clear_color(float4(0.0f, 0.0f, 0.0f, 1.0f));
  accumulation_ps_.init_subpasses(
      ePipelineType::TRANSPARENT, shading_type, resources.shader_cache);

  resolve_ps_.init();
  /*TODO(pragma37): Use shaders.resolve_shader_get()*/
  static GPUShader *resolve_shader = GPU_shader_create_from_info_name(
      "workbench_transparent_resolve");
  resolve_ps_.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_BLEND_ALPHA);
  resolve_ps_.shader_set(resolve_shader);
  resolve_ps_.bind_texture("transparentAccum", &accumulation_tx);
  resolve_ps_.bind_texture("transparentRevealage", &reveal_tx);
  resolve_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
}

void TransparentPass::draw(Manager &manager, View &view)
{
  if (is_empty()) {
    return;
  }

  accumulation_tx.acquire(int2(depth_tx.size()), GPU_RGBA16F);
  reveal_tx.acquire(int2(depth_tx.size()), GPU_R16F);
  object_id_tx.acquire(int2(depth_tx.size()), GPU_R16UI);

  transparent_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                        GPU_ATTACHMENT_TEXTURE(accumulation_tx),
                        GPU_ATTACHMENT_TEXTURE(reveal_tx),
                        GPU_ATTACHMENT_TEXTURE(object_id_tx));
  transparent_fb.bind();

  manager.submit(accumulation_ps_, view);

  resolve_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(color_tx));
  resolve_fb.bind();

  manager.submit(resolve_ps_, view);

  accumulation_tx.release();
  reveal_tx.release();
  object_id_tx.release();
}

bool TransparentPass::is_empty() const
{
  return accumulation_ps_.is_empty();
}

}  // namespace blender::workbench
