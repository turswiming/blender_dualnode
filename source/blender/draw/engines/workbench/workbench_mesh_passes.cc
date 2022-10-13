/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

namespace blender::workbench {

MeshPass::MeshPass(const char *name) : PassMain(name){};

/* Move to draw::Pass */
bool MeshPass::is_empty() const
{
  return false; /* TODO */
}

void MeshPass::init(ePipelineType pipeline,
                    eShadingType shading,
                    SceneResources &resources,
                    DRWState state)
{
  ShaderCache &shaders = resources.shader_cache;

  this->PassMain::init();
  this->state_set(state);
  this->bind_texture(WB_MATCAP_SLOT, resources.matcap_tx);
  this->bind_ssbo(WB_MATERIAL_SLOT, &resources.material_buf);
  this->bind_ubo(WB_WORLD_SLOT, resources.world_buf);

  texture_subpass_map.clear();

  for (auto geom : IndexRange(geometry_type_len)) {
    for (auto color : IndexRange(color_type_len)) {
      eGeometryType geom_type = static_cast<eGeometryType>(geom);
      eColorType color_type = static_cast<eColorType>(color);
      std::string name = std::string(get_name(geom_type)) + std::string(get_name(color_type));
      GPUShader *sh = shaders.prepass_shader_get(pipeline, geom_type, color_type, shading);
      PassMain::Sub *pass = &this->sub(name.c_str());
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

void OpaquePass::sync(DRWState cull_state,
                      DRWState clip_state,
                      eShadingType shading_type,
                      SceneResources &resources)
{
  Texture &depth_tx = resources.depth_tx;
  Texture &depth_in_front_tx = resources.depth_in_front_tx;
  TextureFromPool &color_tx = resources.color_tx;
  ShaderCache &shaders = resources.shader_cache;
  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                   cull_state | clip_state;

  gbuffer_ps_.init(ePipelineType::OPAQUE, shading_type, resources, state);

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

void OpaquePass::draw_prepass(Manager &manager, View &view, Texture &depth_tx)
{
  gbuffer_material_tx.acquire(int2(depth_tx.size()), GPU_RGBA16F);
  gbuffer_normal_tx.acquire(int2(depth_tx.size()), GPU_RG16F);
  gbuffer_object_id_tx.acquire(int2(depth_tx.size()), GPU_R16UI);

  opaque_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_material_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_normal_tx),
                   GPU_ATTACHMENT_TEXTURE(gbuffer_object_id_tx));
  opaque_fb.bind();
  opaque_fb.clear_depth(1.0f);

  manager.submit(gbuffer_ps_, view);
}

void OpaquePass::draw_resolve(Manager &manager, View &view)
{
  manager.submit(deferred_ps_, view);

  gbuffer_normal_tx.release();
  gbuffer_material_tx.release();
  gbuffer_object_id_tx.release();
}

bool OpaquePass::is_empty() const
{
  return gbuffer_ps_.is_empty();
}

void TransparentPass::sync(DRWState cull_state,
                           DRWState clip_state,
                           eShadingType shading_type,
                           SceneResources &resources)
{
  ShaderCache &shaders = resources.shader_cache;
  Texture &depth_tx = resources.depth_tx;
  DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                   cull_state | clip_state;

  accumulation_ps_.init(ePipelineType::TRANSPARENT, shading_type, resources, state);

  resolve_ps_.init();
  resolve_ps_.shader_set(
      shaders.resolve_shader_get(ePipelineType::TRANSPARENT, eShadingType::FLAT));
  resolve_ps_.bind_texture("accumulation_tx", accumulation_tx);
  resolve_ps_.bind_texture("reveal_tx", reveal_tx);
  resolve_ps_.dispatch(math::divide_ceil(depth_tx.size(), int3(WB_RESOLVE_GROUP_SIZE)));
}

void TransparentPass::draw_prepass(Manager &manager, View &view, Texture &depth_tx)
{
  accumulation_tx.acquire(int2(depth_tx.size()), GPU_RGBA16F);
  reveal_tx.acquire(int2(depth_tx.size()), GPU_R8);

  transparent_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                        GPU_ATTACHMENT_TEXTURE(accumulation_tx),
                        GPU_ATTACHMENT_TEXTURE(reveal_tx));
  transparent_fb.bind();

  manager.submit(accumulation_ps_, view);
}

void TransparentPass::draw_resolve(Manager &manager, View &view)
{
  manager.submit(resolve_ps_, view);

  accumulation_tx.release();
  reveal_tx.release();
}

bool TransparentPass::is_empty() const
{
  return accumulation_ps_.is_empty();
}

}  // namespace blender::workbench
