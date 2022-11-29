/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "DNA_camera_types.h"
#include "DRW_render.h"
#include "draw_manager.hh"
#include "draw_pass.hh"

#include "workbench_defines.hh"
#include "workbench_enums.hh"
#include "workbench_shader_shared.h"

namespace blender::workbench {

using namespace draw;

class ShaderCache {
 private:
  /* TODO(fclem): We might want to change to a Map since most shader will never be compiled. */
  GPUShader *prepass_shader_cache_[pipeline_type_len][geometry_type_len][shader_type_len]
                                  [lighting_type_len][2] = {{{{{nullptr}}}}};
  GPUShader *resolve_shader_cache_[pipeline_type_len][lighting_type_len][2][2] = {{{{nullptr}}}};

 public:
  ~ShaderCache();

  GPUShader *prepass_shader_get(ePipelineType pipeline_type,
                                eGeometryType geometry_type,
                                eShaderType shader_type,
                                eLightingType lighting_type,
                                bool clip);

  GPUShader *resolve_shader_get(ePipelineType pipeline_type,
                                eLightingType lighting_type,
                                bool cavity = false,
                                bool curvature = false);
};

struct Material {
  float3 base_color;
  /* Packed data into a int. Decoded in the shader. */
  uint packed_data;

  Material();
  Material(float3 color);
  Material(::Object &ob, bool random = false);
  Material(::Material &mat);

  bool is_transparent();

  static uint32_t pack_data(float metallic, float roughness, float alpha);
};

void get_material_image(Object *ob,
                        int material_index,
                        ::Image *&image,
                        ImageUser *&iuser,
                        eGPUSamplerState &sampler_state);

struct SceneState {
  Scene *scene;

  Object *camera_object;
  Camera *camera;
  float4x4 view_projection_matrix;
  int2 resolution;

  eContextObjectMode object_mode;

  View3DShading shading;
  eLightingType lighting_type = eLightingType::STUDIO;
  bool xray_mode;

  DRWState cull_state;
  Vector<float4> clip_planes = {};

  float4 background_color;

  bool draw_cavity;
  bool draw_curvature;
  bool draw_outline;
  bool draw_dof;
  bool draw_aa;

  bool draw_object_id;
  bool draw_transparent_depth;

  int sample;
  int samples_len;
  bool reset_taa_next_sample;
  bool render_finished;

  /* Used when material_type == eMaterialType::SINGLE */
  Material material_override = Material(float3(1.0f));
  /* When r == -1.0 the shader uses the vertex color */
  Material material_attribute_color = Material(float3(-1.0f));

  void init(Object *camera_ob = nullptr);
};

struct ObjectState {
  eV3DShadingColorType color_type;
  bool sculpt_pbvh;
  bool texture_paint_mode;
  ::Image *image_paint_override;
  eGPUSamplerState override_sampler_state;
  bool draw_shadow;
  bool use_per_material_batches;

  ObjectState(const SceneState &scene_state, Object *ob);
};

class CavityEffect {
 private:
  int sample_;
  int sample_count_;
  bool curvature_enabled_;
  bool cavity_enabled_;

  static const int max_samples_ = 512;  // This value must be kept in sync with the one declared at
                                        // workbench_composite_info.hh (cavity_samples)
  UniformArrayBuffer<float4, max_samples_> samples_buf;

  void load_samples_buf(int ssao_samples);

 public:
  void init(const SceneState &scene_state, struct SceneResources &resources);
  void setup_resolve_pass(PassSimple &pass, struct SceneResources &resources);
};

struct SceneResources {
  ShaderCache shader_cache;

  StringRefNull current_matcap;
  Texture matcap_tx = "matcap_tx";

  TextureFromPool color_tx = "wb_color_tx";
  TextureFromPool object_id_tx = "wb_object_id_tx";
  TextureFromPool depth_tx = "wb_depth_tx";
  TextureFromPool depth_in_front_tx = "wb_depth_in_front_tx";

  StorageVectorBuffer<Material> material_buf = {"material_buf"};
  UniformBuffer<WorldData> world_buf;
  UniformArrayBuffer<float4, 6> clip_planes_buf;

  static const int jitter_tx_size = 64;
  Texture jitter_tx = "wb_jitter_tx";
  void load_jitter_tx(int total_samples);

  CavityEffect cavity;

  void init(const SceneState &scene_state);
};

class MeshPass : public PassMain {
 private:
  PassMain::Sub *passes_[geometry_type_len][shader_type_len];

  using TextureSubPassKey = std::pair<GPUTexture *, eGeometryType>;
  Map<TextureSubPassKey, PassMain::Sub *> texture_subpass_map_;

  bool is_empty_;

 public:
  MeshPass(const char *name);

  /* Move to draw::Pass */
  bool is_empty() const;

  void init_pass(SceneResources &resources, DRWState state, int clip_planes);
  void init_subpasses(ePipelineType pipeline,
                      eLightingType lighting,
                      bool clip,
                      ShaderCache &shaders);

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
  Framebuffer opaque_fb;

  MeshPass gbuffer_ps_ = {"Opaque.Gbuffer"};
  MeshPass gbuffer_in_front_ps_ = {"Opaque.GbufferInFront"};
  PassSimple deferred_ps_ = {"Opaque.Deferred"};

  void sync(const SceneState &scene_state, SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
  bool is_empty() const;
};

class TransparentPass {
 private:
  GPUShader *resolve_sh_;

 public:
  TextureFromPool accumulation_tx = {"accumulation_accumulation_tx"};
  TextureFromPool reveal_tx = {"accumulation_reveal_tx"};
  Framebuffer transparent_fb;

  MeshPass accumulation_ps_ = {"Transparent.Accumulation"};
  MeshPass accumulation_in_front_ps_ = {"Transparent.AccumulationInFront"};
  PassSimple resolve_ps_ = {"Transparent.Resolve"};
  Framebuffer resolve_fb;

  void sync(const SceneState &scene_state, SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
  bool is_empty() const;
};

class TransparentDepthPass {
 private:
  GPUShader *merge_sh_;

 public:
  MeshPass main_ps_ = {"TransparentDepth.Main"};
  Framebuffer main_fb = {"TransparentDepth.Main"};
  MeshPass in_front_ps_ = {"TransparentDepth.InFront"};
  Framebuffer in_front_fb = {"TransparentDepth.InFront"};
  PassSimple merge_ps_ = {"TransparentDepth.Merge"};
  Framebuffer merge_fb = {"TransparentDepth.Merge"};

  void sync(const SceneState &scene_state, SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
  bool is_empty() const;
};

class OutlinePass {
  bool enabled_;

  PassSimple ps_ = PassSimple("Workbench.Outline");
  GPUShader *sh_;
  Framebuffer fb_ = Framebuffer("Workbench.Outline");

 public:
  void init(const SceneState &scene_state);
  void sync(SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
};

class DofPass {
  bool enabled_;

  float offset_;

  static const int kernel_radius_ = 3;
  static const int samples_len_ = (kernel_radius_ * 2 + 1) * (kernel_radius_ * 2 + 1);

  UniformArrayBuffer<float4, samples_len_> samples_buf_;

  Texture source_tx_;
  Texture coc_halfres_tx_;
  TextureFromPool blur_tx_;

  Framebuffer downsample_fb_, blur1_fb_, blur2_fb_, resolve_fb_;

  GPUShader *prepare_sh_, *downsample_sh_, *blur1_sh_, *blur2_sh_, *resolve_sh_;

  PassSimple down_ps_ = {"Workbench.DoF.DownSample"};
  PassSimple down2_ps_ = {"Workbench.DoF.DownSample2"};
  PassSimple blur_ps_ = {"Workbench.DoF.Blur"};
  PassSimple blur2_ps_ = {"Workbench.DoF.Blur2"};
  PassSimple resolve_ps_ = {"Workbench.DoF.Resolve"};

  float aperture_size_;
  float distance_;
  float invsensor_size_;
  float near_;
  float far_;
  float blades_;
  float rotation_;
  float ratio_;

  void setup_samples();

 public:
  void init(const SceneState &scene_state);
  void sync(SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
  bool is_enabled();
};

class AntiAliasingPass {
 private:
  bool enabled_;
  /* Current TAA sample index in [0..samples_len_] range. */
  int sample_;
  /* Total number of samples to after which TAA stops accumulating samples. */
  int samples_len_;
  /* Weight accumulated. */
  float weight_accum_;
  /* Samples weight for this iteration. */
  float weights_[9];
  /* Sum of weights. */
  float weights_sum_;

  Texture sample0_depth_tx_ = {"sample0_depth_tx"};

  Texture taa_accumulation_tx_ = {"taa_accumulation_tx"};
  Texture smaa_search_tx_ = {"smaa_search_tx"};
  Texture smaa_area_tx_ = {"smaa_area_tx"};
  TextureFromPool smaa_edge_tx_ = {"smaa_edge_tx"};
  TextureFromPool smaa_weight_tx_ = {"smaa_weight_tx"};

  Framebuffer taa_accumulation_fb_ = {"taa_accumulation_fb"};
  Framebuffer smaa_edge_fb_ = {"smaa_edge_fb"};
  Framebuffer smaa_weight_fb_ = {"smaa_weight_fb"};
  Framebuffer smaa_resolve_fb_ = {"smaa_resolve_fb"};

  float4 smaa_viewport_metrics_ = {0.0f, 0.0f, 0.0f, 0.0f};
  float smaa_mix_factor_ = 0.0f;

  GPUShader *taa_accumulation_sh_ = nullptr;
  GPUShader *smaa_edge_detect_sh_ = nullptr;
  GPUShader *smaa_aa_weight_sh_ = nullptr;
  GPUShader *smaa_resolve_sh_ = nullptr;

  PassSimple taa_accumulation_ps_ = {"TAA.Accumulation"};
  PassSimple smaa_edge_detect_ps_ = {"SMAA.EdgeDetect"};
  PassSimple smaa_aa_weight_ps_ = {"SMAA.BlendWeights"};
  PassSimple smaa_resolve_ps_ = {"SMAA.Resolve"};

 public:
  AntiAliasingPass();
  ~AntiAliasingPass();

  void init(const SceneState &scene_state);
  void sync(SceneResources &resources, int2 resolution);
  void setup_view(View &view, int2 resolution);
  void draw(Manager &manager,
            View &view,
            SceneResources &resources,
            int2 resolution,
            GPUTexture *depth_tx,
            GPUTexture *color_tx);
};

}  // namespace blender::workbench
