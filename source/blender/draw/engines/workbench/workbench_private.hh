/* SPDX-License-Identifier: GPL-2.0-or-later */

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
  GPUShader *prepass_shader_cache_[pipeline_type_len][geometry_type_len][color_type_len]
                                  [shading_type_len] = {{{{nullptr}}}};
  GPUShader *resolve_shader_cache_[pipeline_type_len][shading_type_len][2][2] = {{{{nullptr}}}};

 public:
  ~ShaderCache();

  GPUShader *prepass_shader_get(ePipelineType pipeline_type,
                                eGeometryType geometry_type,
                                eColorType color_type,
                                eShadingType shading_type);

  GPUShader *resolve_shader_get(ePipelineType pipeline_type,
                                eShadingType shading_type,
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
  float4x4 view_projection_matrix;
  int2 resolution;

  eContextObjectMode object_mode;

  View3DShading shading;
  eShadingType shading_type = eShadingType::STUDIO;
  bool xray_mode;

  DRWState cull_state;
  DRWState clip_state;
  Vector<float4> clip_planes = {};

  float4 background_color;

  bool draw_cavity;
  bool draw_curvature;
  bool draw_outline;
  bool draw_dof;

  bool draw_object_id;
  bool draw_transparent_depth;

  int aa_samples;
  bool reset_taa;
  bool reset_taa_next_sample;

  /** Used when material_subtype == eMaterialSubType::SINGLE */
  Material material_override = Material(float3(1.0f));
  /* When r == -1.0 the shader uses the vertex color */
  Material material_attribute_color = Material(float3(-1.0f));

  void init();
};

class ObjectState {
 private:
  void setup_material_state();

 public:
  eV3DShadingColorType color_type;
  bool sculpt_pbvh;
  bool texture_paint_mode;
  ::Image *image_paint_override;
  eGPUSamplerState override_sampler_state;
  bool draw_shadow;

  eColorType material_type;
  eMaterialSubType material_subtype;
  bool use_per_material_batches;

  ObjectState(const SceneState &scene_state, Object *ob);
};

class CavityEffect {
 private:
  int sample;
  int sample_count;
  bool curvature_enabled;
  bool cavity_enabled;

 public:
  static const int JITTER_TEX_SIZE = 64;
  static const int MAX_SAMPLES = 512;  // This value must be kept in sync with the one declared at
                                       // workbench_composite_info.hh (cavity_samples)

  UniformArrayBuffer<float4, MAX_SAMPLES> samples_buf;
  /*TODO(Miguel Pozo): Move to SceneResources (Used by DoF too)*/
  Texture jitter_tx;

  void setup_resources(int iteration_samples, int total_samples);

  void init(const SceneState &scene_state, UniformBuffer<WorldData> &world_buf);

  void setup_resolve_pass(PassSimple &pass, Texture &object_id_tx);
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

  CavityEffect cavity;

  void init(const SceneState &scene_state);
};

class MeshPass : public PassMain {
 private:
  PassMain::Sub *passes_[geometry_type_len][color_type_len];

  using TextureSubPassKey = std::pair<GPUTexture *, eGeometryType>;
  Map<TextureSubPassKey, PassMain::Sub *> texture_subpass_map;

  bool _is_empty;

 public:
  MeshPass(const char *name);

  /* Move to draw::Pass */
  bool is_empty() const;

  void init_pass(SceneResources &resources, DRWState state);
  void init_subpasses(ePipelineType pipeline, eShadingType shading, ShaderCache &shaders);

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

class DofPass {
  bool enabled = false;

  static const int KERNEL_RADIUS = 3;
  static const int SAMPLES_LEN = (KERNEL_RADIUS * 2 + 1) * (KERNEL_RADIUS * 2 + 1);

  UniformArrayBuffer<float4, SAMPLES_LEN> samples_buf;

  Texture source_tx;
  Texture coc_halfres_tx;
  TextureFromPool blur_tx;

  Framebuffer downsample_fb, blur1_fb, blur2_fb, resolve_fb;

  GPUShader *prepare_sh, *downsample_sh, *blur1_sh, *blur2_sh, *resolve_sh;

  PassSimple down_ps = {"Workbench.DoF.DownSample"};
  PassSimple down2_ps = {"Workbench.DoF.DownSample2"};
  PassSimple blur_ps = {"Workbench.DoF.Blur"};
  PassSimple blur2_ps = {"Workbench.DoF.Blur2"};
  PassSimple resolve_ps = {"Workbench.DoF.Resolve"};

  float aperture_size;
  float distance;
  float invsensor_size;
  float near;
  float far;
  float blades;
  float rotation;
  float ratio;

  void setup_samples();

 public:
  void init(const SceneState &scene_state);
  void sync(SceneResources &resources);
  void draw(Manager &manager, View &view, SceneResources &resources, int2 resolution);
  bool is_enabled();
};

class AntiAliasingPass {
 private:
  /** Total number of samples to after which TAA stops accumulating samples. */
  int sample_len;
  /** Current TAA sample index in [0..sample_len] range. */
  int sample;
  /** Weight accumulated. */
  float weight_accum;
  /** Samples weight for this iteration. */
  float weights[9];
  /** Sum of weights. */
  float weights_sum;
  /** True if the history buffer contains relevant data and false if it could contain garbage. */
  // bool valid_history;

  Texture sample0_depth_tx = {"sample0_depth_tx"};

  Texture taa_accumulation_tx = {"taa_accumulation_tx"};
  Texture smaa_search_tx = {"smaa_search_tx"};
  Texture smaa_area_tx = {"smaa_area_tx"};
  TextureFromPool smaa_edge_tx = {"smaa_edge_tx"};
  TextureFromPool smaa_weight_tx = {"smaa_weight_tx"};

  Framebuffer taa_accumulation_fb = {"taa_accumulation_fb"};
  Framebuffer smaa_edge_fb = {"smaa_edge_fb"};
  Framebuffer smaa_weight_fb = {"smaa_weight_fb"};
  Framebuffer smaa_resolve_fb = {"smaa_resolve_fb"};

  float4 smaa_viewport_metrics = {0.0f, 0.0f, 0.0f, 0.0f};
  float smaa_mix_factor = 0.0f;

  GPUShader *taa_accumulation_sh = nullptr;
  GPUShader *smaa_edge_detect_sh = nullptr;
  GPUShader *smaa_aa_weight_sh = nullptr;
  GPUShader *smaa_resolve_sh = nullptr;

  PassSimple taa_accumulation_ps_ = {"TAA.Accumulation"};
  PassSimple smaa_edge_detect_ps_ = {"SMAA.EdgeDetect"};
  PassSimple smaa_aa_weight_ps_ = {"SMAA.BlendWeights"};
  PassSimple smaa_resolve_ps_ = {"SMAA.Resolve"};

 public:
  AntiAliasingPass();

  ~AntiAliasingPass();

  void init(const SceneState &scene_state);
  void sync(SceneResources &resources, int2 resolution);
  bool setup_view(View &view, int2 resolution);
  void draw(Manager &manager,
            View &view,
            SceneResources &resources,
            int2 resolution,
            GPUTexture *depth_tx,
            GPUTexture *color_tx);
};

}  // namespace blender::workbench
