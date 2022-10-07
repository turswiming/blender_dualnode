/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_studiolight.h"
#include "DEG_depsgraph_query.h"
#include "DRW_render.h"
#include "GPU_capabilities.h"

#include "smaa_textures.h"

#include "draw_manager.hh"
#include "draw_pass.hh"

#include "workbench_defines.h"
#include "workbench_shader_shared.h"

namespace blender::workbench {

using namespace draw;

#define WB_RESOLVE_GROUP_SIZE 8

enum class eGeometryType {
  MESH = 0,
  CURVES,
  POINTCLOUD,
};
static constexpr int geometry_type_len = static_cast<int>(eGeometryType::POINTCLOUD) + 1;

static inline const char *get_name(eGeometryType type)
{
  switch (type) {
    case eGeometryType::MESH:
      return "Mesh";
    case eGeometryType::CURVES:
      return "Curves";
    case eGeometryType::POINTCLOUD:
      return "PointCloud";
    default:
      BLI_assert_unreachable();
      return "";
  }
}

static inline eGeometryType geometry_type_from_object(Object *ob)
{
  switch (ob->type) {
    case OB_CURVES:
      return eGeometryType::CURVES;
    case OB_POINTCLOUD:
      return eGeometryType::POINTCLOUD;
    default:
      return eGeometryType::MESH;
  }
}

enum class ePipelineType {
  OPAQUE = 0,
  TRANSPARENT,
  SHADOW,
};
static constexpr int pipeline_type_len = static_cast<int>(ePipelineType::SHADOW) + 1;

enum class eShadingType {
  FLAT = 0, /* V3D_LIGHTING_STUDIO */
  STUDIO,   /* V3D_LIGHTING_MATCAP */
  MATCAP,   /* V3D_LIGHTING_FLAT */
};
static constexpr int shading_type_len = static_cast<int>(eShadingType::MATCAP) + 1;

enum class eColorType {
  MATERIAL = 0,
  TEXTURE,
};
static constexpr int color_type_len = static_cast<int>(eColorType::TEXTURE) + 1;

enum class eMaterialSubType {
  NONE = 0,
  MATERIAL,
  RANDOM,
  SINGLE,
  OBJECT,
  ATTRIBUTE,
};
static constexpr int material_subtype_len = static_cast<int>(eMaterialSubType::ATTRIBUTE) + 1;

struct Material {
  float3 base_color;
  /* Packed data into a int. Decoded in the shader. */
  uint packed_data;

  Material() = default;

  Material(float3 color)
  {
    base_color = color;
    packed_data = Material::pack_data(0.0f, 0.4f, 1.0f);
  }
  Material(::Object &ob, bool random = false)
  {
    if (random) {
      uint hash = BLI_ghashutil_strhash_p_murmur(ob.id.name);
      if (ob.id.lib) {
        hash = (hash * 13) ^ BLI_ghashutil_strhash_p_murmur(ob.id.lib->filepath);
      }
      float3 hsv = float3(BLI_hash_int_01(hash), 0.5f, 0.8f);
      hsv_to_rgb_v(hsv, base_color);
    }
    else {
      base_color = ob.color;
    }
    packed_data = Material::pack_data(0.0f, 0.4f, ob.color[3]);
  }
  Material(::Material &mat)
  {
    base_color = &mat.r;
    packed_data = Material::pack_data(mat.metallic, mat.roughness, mat.a);
  }

  static uint32_t pack_data(float metallic, float roughness, float alpha)
  {
    /* Remap to Disney roughness. */
    roughness = sqrtf(roughness);
    uint32_t packed_roughness = unit_float_to_uchar_clamp(roughness);
    uint32_t packed_metallic = unit_float_to_uchar_clamp(metallic);
    uint32_t packed_alpha = unit_float_to_uchar_clamp(alpha);
    return (packed_alpha << 16u) | (packed_roughness << 8u) | packed_metallic;
  }
};

class ShaderCache {
 private:
  /* TODO(fclem): We might want to change to a Map since most shader will never be compiled. */
  GPUShader *prepass_shader_cache_[shading_type_len][color_type_len][geometry_type_len]
                                  [pipeline_type_len] = {{{{nullptr}}}};
  GPUShader *resolve_shader_cache_[shading_type_len][pipeline_type_len] = {{nullptr}};

 public:
  ~ShaderCache()
  {
    for (auto i : IndexRange(shading_type_len)) {
      for (auto j : IndexRange(color_type_len)) {
        for (auto k : IndexRange(geometry_type_len)) {
          for (auto l : IndexRange(pipeline_type_len)) {
            DRW_SHADER_FREE_SAFE(prepass_shader_cache_[i][j][k][l]);
          }
        }
      }
    }
    for (auto i : IndexRange(shading_type_len)) {
      for (auto j : IndexRange(pipeline_type_len)) {
        DRW_SHADER_FREE_SAFE(resolve_shader_cache_[i][j]);
      }
    }
  }

  GPUShader *prepass_shader_get(ePipelineType pipeline_type,
                                eGeometryType geometry_type,
                                eColorType color_type,
                                eShadingType shading_type)
  {
    GPUShader *&shader_ptr =
        prepass_shader_cache_[static_cast<int>(pipeline_type)][static_cast<int>(geometry_type)]
                             [static_cast<int>(color_type)][static_cast<int>(shading_type)];

    if (shader_ptr != nullptr) {
      return shader_ptr;
    }
    std::string info_name = "workbench_next_prepass_";
    switch (geometry_type) {
      case eGeometryType::MESH:
        info_name += "mesh_";
        break;
      case eGeometryType::CURVES:
        info_name += "curves_";
        break;
      case eGeometryType::POINTCLOUD:
        info_name += "ptcloud_";
        break;
    }
    switch (pipeline_type) {
      case ePipelineType::OPAQUE:
        info_name += "opaque_";
        break;
      case ePipelineType::TRANSPARENT:
        info_name += "transparent_";
        break;
      case ePipelineType::SHADOW:
        info_name += "shadow_";
        break;
    }
    switch (shading_type) {
      case eShadingType::FLAT:
        info_name += "flat_";
        break;
      case eShadingType::STUDIO:
        info_name += "studio_";
        break;
      case eShadingType::MATCAP:
        info_name += "matcap_";
        break;
    }
    switch (color_type) {
      case eColorType::MATERIAL:
        info_name += "material";
        break;
      case eColorType::TEXTURE:
        info_name += "texture";
        break;
    }
    /* TODO Clipping */
    info_name += "_no_clip";
    shader_ptr = GPU_shader_create_from_info_name(info_name.c_str());
    prepass_shader_cache_[static_cast<int>(pipeline_type)][static_cast<int>(
        geometry_type)][static_cast<int>(color_type)][static_cast<int>(shading_type)] = shader_ptr;
    return shader_ptr;
  }

  GPUShader *resolve_shader_get(ePipelineType pipeline_type, eShadingType shading_type)
  {
    GPUShader *&shader_ptr =
        resolve_shader_cache_[static_cast<int>(shading_type)][static_cast<int>(pipeline_type)];

    if (shader_ptr != nullptr) {
      return shader_ptr;
    }
    std::string info_name = "workbench_next_resolve_";
    switch (pipeline_type) {
      case ePipelineType::OPAQUE:
        info_name += "opaque_";
        break;
      case ePipelineType::TRANSPARENT:
        info_name += "transparent_";
        break;
      case ePipelineType::SHADOW:
        BLI_assert_unreachable();
        break;
    }
    switch (shading_type) {
      case eShadingType::FLAT:
        info_name += "flat";
        break;
      case eShadingType::STUDIO:
        info_name += "studio";
        break;
      case eShadingType::MATCAP:
        info_name += "matcap";
        break;
    }
    shader_ptr = GPU_shader_create_from_info_name(info_name.c_str());
    return shader_ptr;
  }
};

struct SceneResources {
  ShaderCache shader_cache;

  Texture matcap_tx = "matcap_tx";

  TextureFromPool color_tx = "wb_color_tx";
  Texture depth_tx = "wb_depth_tx";
  Texture depth_in_front_tx = "wb_depth_in_front_tx";

  StorageVectorBuffer<workbench::Material> material_buf = {"material_buf"};
  UniformBuffer<WorldData> world_buf;

  void init(int2 output_res, float4 background_color)
  {
    world_buf.background_color = background_color;
    matcap_tx.ensure_2d_array(GPU_RGBA16F, int2(1), 1);
    depth_tx.ensure_2d(GPU_DEPTH24_STENCIL8, output_res);
  }

  void world_sync(StudioLight *studio_light)
  {
    float4x4 rot_matrix = float4x4::identity();

    // if (USE_WORLD_ORIENTATION(wpd)) {
    /* TODO */
    // }

    if (U.edit_studio_light) {
      studio_light = BKE_studiolight_studio_edit_get();
    }

    /* Studio Lights. */
    for (int i = 0; i < 4; i++) {
      LightData &light = world_buf.lights[i];

      SolidLight *sl = (studio_light) ? &studio_light->light[i] : nullptr;
      if (sl && sl->flag) {
        float3 direction = rot_matrix * float3(sl->vec);
        light.direction = float4(UNPACK3(direction), 0.0f);
        /* We should pre-divide the power by PI but that makes the lights really dim. */
        light.specular_color = float4(UNPACK3(sl->spec), 0.0f);
        light.diffuse_color_wrap = float4(UNPACK3(sl->col), sl->smooth);
      }
      else {
        light.direction = float4(1.0f, 0.0f, 0.0f, 0.0f);
        light.specular_color = float4(0.0f);
        light.diffuse_color_wrap = float4(0.0f);
      }
    }

    if (studio_light != nullptr) {
      world_buf.ambient_color = float4(UNPACK3(studio_light->light_ambient), 0.0f);
    }
    else {
      world_buf.ambient_color = float4(1.0f, 1.0f, 1.0f, 0.0f);
    }

    /* TODO */
    world_buf.use_specular = true;

    world_buf.push_update();
  }
};

class MeshPass : public PassMain {
 private:
  std::array<PassMain::Sub *, geometry_type_len> geometry_passes_;
  eColorType color_type_;

  using TextureSubPassKey = std::pair<GPUTexture *, eGeometryType>;
  Map<TextureSubPassKey, PassMain::Sub *> texture_subpass_map;

 public:
  MeshPass(const char *name) : PassMain(name){};

  /* Move to draw::Pass */
  bool is_empty() const
  {
    return false; /* TODO */
  }

  void init(ePipelineType pipeline,
            eColorType color_type,
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

    color_type_ = color_type;
    texture_subpass_map.clear();

    for (int geom = 0; geom < geometry_type_len; geom++) {
      eGeometryType geom_type = static_cast<eGeometryType>(geom);
      GPUShader *sh = shaders.prepass_shader_get(pipeline, geom_type, color_type, shading);
      PassMain::Sub *pass = &this->sub(get_name(geom_type));
      pass->shader_set(sh);
      geometry_passes_[geom] = pass;
    }
  }

  PassMain::Sub &sub_pass_get(eGeometryType geometry_type,
                              ObjectRef & /*ref*/,
                              ::Material * /*material*/)
  {
    if (color_type_ == eColorType::TEXTURE) {
      /* TODO(fclem): Always query a layered texture so we can use only a single shader. */
      GPUTexture *texture = nullptr;  // ref.object->texture_get();
      GPUTexture *tilemap = nullptr;  // ref.object->texture_get();

      auto add_cb = [&] {
        PassMain::Sub *sub_pass = geometry_passes_[static_cast<int>(geometry_type)];
        sub_pass = &sub_pass->sub("Blender Texture Name" /* texture.name */);
        sub_pass->bind_texture(WB_TEXTURE_SLOT, texture);
        sub_pass->bind_texture(WB_TILEMAP_SLOT, tilemap);
        return sub_pass;
      };

      return *texture_subpass_map.lookup_or_add_cb(TextureSubPassKey(texture, geometry_type),
                                                   add_cb);
    }
    return *geometry_passes_[static_cast<int>(geometry_type)];
  }
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
            eColorType color_type,
            SceneResources &resources)
  {
    Texture &depth_tx = resources.depth_tx;
    Texture &depth_in_front_tx = resources.depth_in_front_tx;
    TextureFromPool &color_tx = resources.color_tx;
    ShaderCache &shaders = resources.shader_cache;
    DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     cull_state | clip_state;

    gbuffer_ps_.init(ePipelineType::OPAQUE, color_type, shading_type, resources, state);

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

  void draw_prepass(Manager &manager, View &view, Texture &depth_tx)
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

  void draw_resolve(Manager &manager, View &view)
  {
    manager.submit(deferred_ps_, view);

    gbuffer_normal_tx.release();
    gbuffer_material_tx.release();
    gbuffer_object_id_tx.release();
  }

  bool is_empty() const
  {
    return gbuffer_ps_.is_empty();
  }
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
            eColorType color_type,
            SceneResources &resources)
  {
    ShaderCache &shaders = resources.shader_cache;
    Texture &depth_tx = resources.depth_tx;
    DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     cull_state | clip_state;

    accumulation_ps_.init(ePipelineType::TRANSPARENT, color_type, shading_type, resources, state);

    resolve_ps_.init();
    resolve_ps_.shader_set(
        shaders.resolve_shader_get(ePipelineType::TRANSPARENT, eShadingType::FLAT));
    resolve_ps_.bind_texture("accumulation_tx", accumulation_tx);
    resolve_ps_.bind_texture("reveal_tx", reveal_tx);
    resolve_ps_.dispatch(math::divide_ceil(depth_tx.size(), int3(WB_RESOLVE_GROUP_SIZE)));
  }

  void draw_prepass(Manager &manager, View &view, Texture &depth_tx)
  {
    accumulation_tx.acquire(int2(depth_tx.size()), GPU_RGBA16F);
    reveal_tx.acquire(int2(depth_tx.size()), GPU_R8);

    transparent_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                          GPU_ATTACHMENT_TEXTURE(accumulation_tx),
                          GPU_ATTACHMENT_TEXTURE(reveal_tx));
    transparent_fb.bind();

    manager.submit(accumulation_ps_, view);
  }

  void draw_resolve(Manager &manager, View &view)
  {
    manager.submit(resolve_ps_, view);

    accumulation_tx.release();
    reveal_tx.release();
  }

  bool is_empty() const
  {
    return accumulation_ps_.is_empty();
  }
};

class AntiAliasingPass {
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

  AntiAliasingPass()
  {
    smaa_edge_detect_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_0");
    smaa_aa_weight_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_1");
    smaa_resolve_sh = GPU_shader_create_from_info_name("workbench_smaa_stage_2");

    smaa_search_tx.ensure_2d(GPU_R8, {SEARCHTEX_WIDTH, SEARCHTEX_HEIGHT});
    GPU_texture_update(smaa_search_tx, GPU_DATA_UBYTE, searchTexBytes);
    GPU_texture_filter_mode(smaa_search_tx, true);

    smaa_area_tx.ensure_2d(GPU_RG8, {AREATEX_WIDTH, AREATEX_HEIGHT});
    GPU_texture_update(smaa_area_tx, GPU_DATA_UBYTE, areaTexBytes);
    GPU_texture_filter_mode(smaa_area_tx, true);
  }

  ~AntiAliasingPass()
  {
    if (smaa_edge_detect_sh) {
      GPU_shader_free(smaa_edge_detect_sh);
    }
    if (smaa_aa_weight_sh) {
      GPU_shader_free(smaa_aa_weight_sh);
    }
    if (smaa_resolve_sh) {
      GPU_shader_free(smaa_resolve_sh);
    }
  }

  void sync(SceneResources &resources)
  {
    {
      smaa_edge_detect_ps_.init();
      smaa_edge_detect_ps_.state_set(DRW_STATE_WRITE_COLOR);
      smaa_edge_detect_ps_.shader_set(smaa_edge_detect_sh);
      smaa_edge_detect_ps_.bind_texture("colorTex", &resources.color_tx);
      smaa_edge_detect_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
      smaa_edge_detect_ps_.clear_color(float4(0.0f));
      smaa_edge_detect_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
    }
    {
      smaa_aa_weight_ps_.init();
      smaa_aa_weight_ps_.state_set(DRW_STATE_WRITE_COLOR);
      smaa_aa_weight_ps_.shader_set(smaa_aa_weight_sh);
      smaa_aa_weight_ps_.bind_texture("edgesTex", &smaa_edge_tx);
      smaa_aa_weight_ps_.bind_texture("areaTex", smaa_area_tx);
      smaa_aa_weight_ps_.bind_texture("searchTex", smaa_search_tx);
      smaa_aa_weight_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
      smaa_aa_weight_ps_.clear_color(float4(0.0f));
      smaa_aa_weight_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
    }
    {
      smaa_resolve_ps_.init();
      smaa_resolve_ps_.state_set(DRW_STATE_WRITE_COLOR);
      smaa_resolve_ps_.shader_set(smaa_resolve_sh);
      smaa_resolve_ps_.bind_texture("blendTex", &smaa_weight_tx);
      smaa_resolve_ps_.bind_texture("colorTex", &resources.color_tx);
      smaa_resolve_ps_.push_constant("viewportMetrics", &smaa_viewport_metrics, 1);
      smaa_resolve_ps_.push_constant("mixFactor", &smaa_mix_factor, 1);
      smaa_resolve_ps_.push_constant("taaAccumulatedWeight", &taa_weight_accum, 1);
      smaa_resolve_ps_.clear_color(float4(0.0f));
      smaa_resolve_ps_.draw_procedural(GPU_PRIM_TRIS, 1, 3);
    }
  }

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    int2 size = {GPU_texture_width(depth_tx), GPU_texture_height(depth_tx)};

    taa_weight_accum = 1.0f; /* TODO */

    smaa_viewport_metrics = float4(1.0f / size.x, 1.0f / size.y, size.x, size.y);
    smaa_mix_factor = 1.0f; /* TODO */

    smaa_edge_tx.acquire(size, GPU_RG8);
    smaa_edge_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(smaa_edge_tx));
    smaa_edge_fb.bind();
    manager.submit(smaa_edge_detect_ps_, view);

    smaa_weight_tx.acquire(size, GPU_RGBA8);
    smaa_weight_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(smaa_weight_tx));
    smaa_weight_fb.bind();
    manager.submit(smaa_aa_weight_ps_, view);
    smaa_edge_tx.release();

    smaa_resolve_fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(color_tx));
    smaa_resolve_fb.bind();
    manager.submit(smaa_resolve_ps_, view);
    smaa_weight_tx.release();
  }
};

class Instance {
 public:
  SceneResources resources;

  OpaquePass opaque_ps;
  // OpaquePass opaque_in_front_ps;

  // TransparentPass transparent_ps;
  // TransparentPass transparent_in_front_ps;

  AntiAliasingPass anti_aliasing_ps;

  DRWState clip_state;
  DRWState cull_state;

  bool use_per_material_batches = false;
  eColorType color_type = eColorType::MATERIAL;
  eMaterialSubType material_subtype = eMaterialSubType::MATERIAL;
  /** Used when material_subtype == eMaterialSubType::SINGLE */
  Material material_override = Material(float3(1.0f));

  eShadingType shading_type = eShadingType::STUDIO;
  /** Chosen studiolight or matcap. */
  StudioLight *studio_light;

  void init(const int2 &output_res,
            const Depsgraph *depsgraph,
            const Object * /*camera*/,
            const View3D *v3d,
            const RegionView3D *rv3d)
  {
    Scene *scene = DEG_get_evaluated_scene(depsgraph);
    View3DShading &shading = scene->display.shading;

    cull_state = DRW_STATE_NO_DRAW;
    if (shading.flag & V3D_SHADING_BACKFACE_CULLING) {
      cull_state |= DRW_STATE_CULL_BACK;
    }

    use_per_material_batches = ELEM(
        shading.color_type, V3D_SHADING_TEXTURE_COLOR, V3D_SHADING_MATERIAL_COLOR);

    color_type = shading.color_type == V3D_SHADING_TEXTURE_COLOR ? eColorType::TEXTURE :
                                                                   eColorType::MATERIAL;

    if (color_type == eColorType::MATERIAL) {
      switch (shading.color_type) {
        case V3D_SHADING_MATERIAL_COLOR:
          material_subtype = eMaterialSubType::MATERIAL;
          break;
        case V3D_SHADING_RANDOM_COLOR:
          material_subtype = eMaterialSubType::RANDOM;
          break;
        case V3D_SHADING_SINGLE_COLOR:
          material_subtype = eMaterialSubType::SINGLE;
          break;
        case V3D_SHADING_TEXTURE_COLOR:
          BLI_assert_msg(false, "V3D_SHADING_TEXTURE_COLOR is not an eMaterialSubType");
          break;
        case V3D_SHADING_OBJECT_COLOR:
          material_subtype = eMaterialSubType::OBJECT;
          break;
        case V3D_SHADING_VERTEX_COLOR:
          material_subtype = eMaterialSubType::ATTRIBUTE;
          break;
        default:
          BLI_assert_msg(false, "Unhandled V3D_SHADING type");
      }
    }
    else {
      material_subtype = eMaterialSubType::NONE;
    }

    if (material_subtype == eMaterialSubType::SINGLE) {
      material_override = Material(shading.single_color);
    }
    else if (material_subtype == eMaterialSubType::ATTRIBUTE) {
      /* TODO(pragma37): Don't use override, check per object if it has color attribute */
      /* When r == -1.0 the shader uses the vertex color */
      material_override = Material(float3(-1.0f, 1.0f, 1.0f));
    }

    switch (shading.light) {
      case V3D_LIGHTING_FLAT:
        shading_type = eShadingType::FLAT;
        break;
      case V3D_LIGHTING_MATCAP:
        shading_type = eShadingType::MATCAP;
        break;
      case V3D_LIGHTING_STUDIO:
        shading_type = eShadingType::STUDIO;
        break;
    }

    float4 background_color = float4(0.0f);

    /*TODO(pragma37): Replace when Workbench next is complete*/
    // bool is_workbench_render = BKE_scene_uses_blender_workbench(scene);
    bool is_workbench_render = std::string(scene->r.engine) ==
                               std::string("BLENDER_WORKBENCH_NEXT");

    /* TODO(pragma37):
     * Check why Workbench Next exposes OB_MATERIAL, and Workbench exposes OB_RENDER */
    if (!v3d || (ELEM(v3d->shading.type, OB_RENDER, OB_MATERIAL) && is_workbench_render)) {
      if (scene->r.alphamode != R_ALPHAPREMUL) {
        World *w = scene->world;
        background_color = w ? float4(w->horr, w->horg, w->horb, 1.0f) :
                               float4(0.0f, 0.0f, 0.0f, 1.0f);
      }
    }

    resources.init(output_res, background_color);

    studio_light = nullptr;
    if (shading_type == eShadingType::MATCAP) {
      studio_light = BKE_studiolight_find(shading.matcap, STUDIOLIGHT_TYPE_MATCAP);
    }
    /* If matcaps are missing, use this as fallback. */
    if (studio_light == nullptr) {
      studio_light = BKE_studiolight_find(shading.studio_light, STUDIOLIGHT_TYPE_STUDIO);
    }
  }

  void begin_sync()
  {
    resources.world_sync(studio_light);

    opaque_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // opaque_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // transparent_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // transparent_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    anti_aliasing_ps.sync(resources);
  }

  void object_sync(Manager &manager, ObjectRef &ob_ref)
  {
    if (ob_ref.object->type != OB_MESH) {
      // TODO(pragma37)
      return;
    }

    if (use_per_material_batches) {
      Vector<::Material *> materials = materials_get(ob_ref);
      Span<GPUBatch *> batches = geometry_get(ob_ref, materials.size());

      if (batches.size() == materials.size()) {
        for (auto i : materials.index_range()) {
          /* TODO(fclem): This create a cull-able instance for each sub-object. This is done for
           * simplicity to reduce complexity. But this increase the overhead per object. Instead,
           * we should use an indirection buffer to the material buffer. */
          ResourceHandle handle = manager.resource_handle(ob_ref);

          Material &mat = resources.material_buf.get_or_resize(handle.resource_index());
          mat = Material(*materials[i]);

          pipeline_get(ob_ref, materials[i]).draw(batches[i], handle);
        }
      }
    }
    else {
      ResourceHandle handle = manager.resource_handle(ob_ref);

      Material &mat = resources.material_buf.get_or_resize(handle.resource_index());

      if (material_subtype == eMaterialSubType::OBJECT) {
        mat = Material(*ob_ref.object);
      }
      else if (material_subtype == eMaterialSubType::RANDOM) {
        mat = Material(*ob_ref.object, true);
      }
      else { /* SINGLE OR ATTRIBUTE */
        mat = material_override;
      }

      GPUBatch *batch = geometry_get(ob_ref);
      if (batch) {
        pipeline_get(ob_ref).draw(batch, handle);
      }
    }
  }

  PassMain::Sub &pipeline_get(ObjectRef &ob_ref, ::Material *material = nullptr)
  {
    return opaque_ps.gbuffer_ps_.sub_pass_get(
        geometry_type_from_object(ob_ref.object), ob_ref, material);
  }

  Vector<::Material *> materials_get(ObjectRef &ob_ref)
  {
    const int material_count = DRW_cache_object_material_count_get(ob_ref.object);
    Vector<::Material *> materials(material_count);
    for (auto i : IndexRange(material_count)) {
      ::Material *mat = BKE_object_material_get_eval(ob_ref.object, i + 1);
      if (mat == nullptr) {
        mat = BKE_material_default_empty();
      }
      materials[i] = mat;
    }
    return materials;
  }

  Span<GPUBatch *> geometry_get(ObjectRef &ob_ref, int material_count)
  {
    Vector<GPUMaterial *> gpu_materials(material_count, nullptr, {});
    return {DRW_cache_object_surface_material_get(
                ob_ref.object, gpu_materials.begin(), material_count),
            material_count};
  }

  GPUBatch *geometry_get(ObjectRef &ob_ref)
  {
    if (material_subtype == eMaterialSubType::ATTRIBUTE) {
      /* TODO(pragma37): Should check for vertex paint mode as well */
      return DRW_cache_mesh_surface_vertpaint_get(ob_ref.object);
    }
    return DRW_cache_object_surface_get(ob_ref.object);
  }

  void end_sync()
  {
    resources.material_buf.push_update();
  }

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    resources.color_tx.acquire(int2(resources.depth_tx.size()), GPU_RGBA16F);

    opaque_ps.draw_prepass(manager, view, resources.depth_tx);
    // volume_ps.draw_prepass(manager, view, resources.depth_tx);
    // transparent_ps.draw_prepass(manager, view, resources.depth_tx);

    // if (opaque_in_front_ps.is_empty() == false || transparent_in_front_ps.is_empty() == false) {
    //   opaque_in_front_ps.draw_prepass(manager, view, resources.depth_in_front_tx);
    //   transparent_in_front_ps.draw_prepass(manager, view, resources.depth_in_front_tx);
    // }

    opaque_ps.draw_resolve(manager, view);
    // transparent_ps.draw_resolve(manager, view);

    anti_aliasing_ps.draw(manager, view, depth_tx, color_tx);

    resources.color_tx.release();
  }

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    this->draw(manager, view, depth_tx, color_tx);
  }
};

}  // namespace blender::workbench

/* -------------------------------------------------------------------- */
/** \name Interface with legacy C DRW manager
 * \{ */

using namespace blender;

struct WORKBENCH_Data {
  DrawEngineType *engine_type;
  DRWViewportEmptyList *fbl;
  DRWViewportEmptyList *txl;
  DRWViewportEmptyList *psl;
  DRWViewportEmptyList *stl;
  workbench::Instance *instance;

  char info[GPU_INFO_SIZE];
};

static void workbench_engine_init(void *vedata)
{
  /* TODO(fclem): Remove once it is minimum required. */
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }

  WORKBENCH_Data *ved = reinterpret_cast<WORKBENCH_Data *>(vedata);
  if (ved->instance == nullptr) {
    ved->instance = new workbench::Instance();
  }

  const DRWContextState *ctx_state = DRW_context_state_get();
  View3D *v3d = ctx_state->v3d;
  RegionView3D *rv3d = ctx_state->rv3d;

  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
  int2 size = int2(GPU_texture_width(dtxl->color), GPU_texture_height(dtxl->color));

  Object *camera = nullptr;
  if (v3d) {
    if (rv3d && (rv3d->persp == RV3D_CAMOB)) {
      camera = v3d->camera;
    }
  }

  ved->instance->init(size, ctx_state->depsgraph, camera, v3d, rv3d);
}

static void workbench_cache_init(void *vedata)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  reinterpret_cast<WORKBENCH_Data *>(vedata)->instance->begin_sync();
}

static void workbench_cache_populate(void *vedata, Object *object)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  draw::Manager *manager = DRW_manager_get();

  draw::ObjectRef ref;
  ref.object = object;
  ref.dupli_object = DRW_object_get_dupli(object);
  ref.dupli_parent = DRW_object_get_dupli_parent(object);

  reinterpret_cast<WORKBENCH_Data *>(vedata)->instance->object_sync(*manager, ref);
}

static void workbench_cache_finish(void *vedata)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  reinterpret_cast<WORKBENCH_Data *>(vedata)->instance->end_sync();
}

static void workbench_draw_scene(void *vedata)
{
  WORKBENCH_Data *ved = reinterpret_cast<WORKBENCH_Data *>(vedata);
  if (!GPU_shader_storage_buffer_objects_support()) {
    STRNCPY(ved->info, "Error: No shader storage buffer support");
    return;
  }
  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
  const DRWView *default_view = DRW_view_default_get();
  draw::Manager *manager = DRW_manager_get();
  draw::View view("DefaultView", default_view);
  ved->instance->draw_viewport(*manager, view, dtxl->depth, dtxl->color);
}

static void workbench_instance_free(void *instance)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  delete reinterpret_cast<workbench::Instance *>(instance);
}

static void workbench_view_update(void *vedata)
{
  UNUSED_VARS(vedata);
}

static void workbench_id_update(void *vedata, struct ID *id)
{
  UNUSED_VARS(vedata, id);
}

static void workbench_render_to_image(void *vedata,
                                      struct RenderEngine *engine,
                                      struct RenderLayer *layer,
                                      const struct rcti *UNUSED(rect))
{
  UNUSED_VARS(vedata, engine, layer);
}

static void workbench_render_update_passes(RenderEngine *engine,
                                           Scene *scene,
                                           ViewLayer *view_layer)
{
  RE_engine_register_pass(engine, scene, view_layer, RE_PASSNAME_COMBINED, 4, "RGBA", SOCK_RGBA);
}

extern "C" {

static const DrawEngineDataSize workbench_data_size = DRW_VIEWPORT_DATA_SIZE(WORKBENCH_Data);

DrawEngineType draw_engine_workbench_next = {
    nullptr,
    nullptr,
    N_("Workbench"),
    &workbench_data_size,
    &workbench_engine_init,
    nullptr,
    &workbench_instance_free,
    &workbench_cache_init,
    &workbench_cache_populate,
    &workbench_cache_finish,
    &workbench_draw_scene,
    &workbench_view_update,
    &workbench_id_update,
    &workbench_render_to_image,
    nullptr,
};

RenderEngineType DRW_engine_viewport_workbench_next_type = {
    nullptr,
    nullptr,
    "BLENDER_WORKBENCH_NEXT",
    N_("Workbench Next"),
    RE_INTERNAL | RE_USE_STEREO_VIEWPORT | RE_USE_GPU_CONTEXT,
    nullptr,
    &DRW_render_to_image,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    &workbench_render_update_passes,
    &draw_engine_workbench_next,
    {nullptr, nullptr, nullptr},
};
}

/** \} */
