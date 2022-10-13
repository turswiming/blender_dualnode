/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "DRW_render.h"
#include "draw_manager.hh"
#include "draw_pass.hh"

#include "workbench_defines.hh"
#include "workbench_enums.hh"
#include "workbench_shader_shared.h"

namespace blender::workbench {

using namespace draw;

struct Material {
  float3 base_color;
  /* Packed data into a int. Decoded in the shader. */
  uint packed_data;

  Material();
  Material(float3 color);
  Material(::Object &ob, bool random = false);
  Material(::Material &mat);

  static uint32_t pack_data(float metallic, float roughness, float alpha);
};

void get_material_image(Object *ob,
                        int material_index,
                        ::Image *&image,
                        ImageUser *&iuser,
                        eGPUSamplerState &sampler_state);

class ShaderCache {
 private:
  /* TODO(fclem): We might want to change to a Map since most shader will never be compiled. */
  GPUShader *prepass_shader_cache_[shading_type_len][color_type_len][geometry_type_len]
                                  [pipeline_type_len] = {{{{nullptr}}}};
  GPUShader *resolve_shader_cache_[shading_type_len][pipeline_type_len] = {{nullptr}};

 public:
  ~ShaderCache();

  GPUShader *prepass_shader_get(ePipelineType pipeline_type,
                                eGeometryType geometry_type,
                                eColorType color_type,
                                eShadingType shading_type);

  GPUShader *resolve_shader_get(ePipelineType pipeline_type, eShadingType shading_type);
};

struct SceneResources {
  ShaderCache shader_cache;

  Texture matcap_tx = "matcap_tx";

  TextureFromPool color_tx = "wb_color_tx";
  Texture depth_tx = "wb_depth_tx";
  Texture depth_in_front_tx = "wb_depth_in_front_tx";

  StorageVectorBuffer<Material> material_buf = {"material_buf"};
  UniformBuffer<WorldData> world_buf;
};

class MeshPass : public PassMain {
 private:
  PassMain::Sub *passes_[geometry_type_len][color_type_len];

  using TextureSubPassKey = std::pair<GPUTexture *, eGeometryType>;
  Map<TextureSubPassKey, PassMain::Sub *> texture_subpass_map;

 public:
  MeshPass(const char *name);

  /* Move to draw::Pass */
  bool is_empty() const;

  void init(ePipelineType pipeline,
            eShadingType shading,
            SceneResources &resources,
            DRWState state);

  PassMain::Sub &sub_pass_get(
      ObjectRef &ref,
      ::Image *image = nullptr,
      eGPUSamplerState sampler_state = eGPUSamplerState::GPU_SAMPLER_DEFAULT,
      ImageUser *iuser = nullptr);
};

class OpaquePass {
 public:
  TextureFromPool gbuffer_normal_tx = {"gbuffer_normal_tx"};
  TextureFromPool gbuffer_material_tx = {"gbuffer_material_tx"};
  TextureFromPool gbuffer_object_id_tx = {"gbuffer_object_id_tx"};
  Framebuffer opaque_fb;

  MeshPass gbuffer_ps_ = {"Opaque.Gbuffer"};
  PassSimple deferred_ps_ = {"Opaque.Deferred"};

  void sync(DRWState cull_state,
            DRWState clip_state,
            eShadingType shading_type,
            SceneResources &resources);

  void draw_prepass(Manager &manager, View &view, Texture &depth_tx);

  void draw_resolve(Manager &manager, View &view);

  bool is_empty() const;
};

class TransparentPass {
 public:
  TextureFromPool accumulation_tx;
  TextureFromPool reveal_tx;
  Framebuffer transparent_fb;

  MeshPass accumulation_ps_ = {"Transparent.Accumulation"};
  PassSimple resolve_ps_ = {"Transparent.Resolve"};

  void sync(DRWState cull_state,
            DRWState clip_state,
            eShadingType shading_type,
            SceneResources &resources);

  void draw_prepass(Manager &manager, View &view, Texture &depth_tx);

  void draw_resolve(Manager &manager, View &view);

  bool is_empty() const;
};

class AntiAliasingPass {
 private:
  bool is_playback;
  bool is_navigating;
  int taa_sample = 0;

 public:
  Texture smaa_search_tx = {"smaa_search_tx"};
  Texture smaa_area_tx = {"smaa_area_tx"};
  TextureFromPool smaa_edge_tx = {"smaa_edge_tx"};
  TextureFromPool smaa_weight_tx = {"smaa_weight_tx"};

  Framebuffer smaa_edge_fb = {"smaa_edge_fb"};
  Framebuffer smaa_weight_fb = {"smaa_weight_fb"};
  Framebuffer smaa_resolve_fb = {"smaa_resolve_fb"};

  float4 smaa_viewport_metrics = {0.0f, 0.0f, 0.0f, 0.0f};
  float smaa_mix_factor = 0.0f;
  float taa_weight_accum = 1.0f;

  GPUShader *smaa_edge_detect_sh = nullptr;
  GPUShader *smaa_aa_weight_sh = nullptr;
  GPUShader *smaa_resolve_sh = nullptr;

  PassSimple smaa_edge_detect_ps_ = {"SMAA.EdgeDetect"};
  PassSimple smaa_aa_weight_ps_ = {"SMAA.BlendWeights"};
  PassSimple smaa_resolve_ps_ = {"SMAA.Resolve"};

  AntiAliasingPass();

  ~AntiAliasingPass();

  void init(bool reset_taa);
  void sync(SceneResources &resources);

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx);
};

}  // namespace blender::workbench
