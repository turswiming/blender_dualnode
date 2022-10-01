/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "DRW_render.h"
#include "GPU_capabilities.h"

#include "draw_manager.hh"
#include "draw_pass.hh"

namespace blender::workbench {

using namespace draw;

#define WB_MATCAP_SLOT 0
#define WB_TEXTURE_SLOT 1
#define WB_TILEMAP_SLOT 2
#define WB_MATERIAL_SLOT 0

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
      return "Curves";
    case eGeometryType::CURVES:
      return "PointCloud";
    case eGeometryType::POINTCLOUD:
      return "Mesh";
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
  Material(::Object &ob)
  {
    base_color = ob.color;
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
    std::string info_name = "workbench_next_";
    switch (geometry_type) {
      case eGeometryType::MESH:
        info_name += "mesh_";
        break;
      case eGeometryType::CURVES:
        info_name += "curves_";
        break;
      case eGeometryType::POINTCLOUD:
        info_name += "pointcloud_";
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
    shader_ptr = GPU_shader_create_from_info_name(info_name.c_str());
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
  ShaderCache *shader_cache;

  Texture matcap_tx;
  Texture depth_tx;
  Texture depth_in_front_tx;

  StorageBuffer<workbench::Material> material_buf = {"material_buf"};
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
    this->PassMain::init();
    this->state_set(state);
    this->bind_texture(WB_MATCAP_SLOT, resources.matcap_tx);
    this->bind_ssbo(WB_MATERIAL_SLOT, resources.material_buf);

    color_type_ = color_type;
    texture_subpass_map.clear();

    for (int geom = 0; geom < geometry_type_len; geom++) {
      eGeometryType geom_type = static_cast<eGeometryType>(geom);
      GPUShader *sh = resources.shader_cache->prepass_shader_get(
          pipeline, geom_type, color_type, shading);
      PassMain::Sub *pass = &this->sub(get_name(geom_type));
      pass->shader_set(sh);
      geometry_passes_[geom] = pass;
    }
  }

  PassMain::Sub &sub_pass_get(eGeometryType geometry_type, ObjectRef &ref)
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
  TextureFromPool gbuffer_normal_tx;
  TextureFromPool gbuffer_albedo_tx;
  TextureFromPool gbuffer_object_id_tx;
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
    DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     cull_state | clip_state;

    gbuffer_ps_.init(ePipelineType::OPAQUE, color_type, shading_type, resources, state);

    deferred_ps_.init();
    deferred_ps_.shader_set(
        resources.shader_cache->resolve_shader_get(ePipelineType::OPAQUE, shading_type));
    deferred_ps_.bind_texture("depth_tx", depth_tx);
    deferred_ps_.bind_texture("stencil_tx", depth_tx.stencil_view());
    deferred_ps_.bind_texture("normal_tx", gbuffer_normal_tx);
    deferred_ps_.bind_texture("albedo_tx", gbuffer_albedo_tx);
    deferred_ps_.bind_texture("object_id_tx", gbuffer_object_id_tx);
    deferred_ps_.dispatch(math::divide_ceil(depth_tx.size(), int3(WB_RESOLVE_GROUP_SIZE)));
  }

  void draw_prepass(Manager &manager, View &view, Texture &depth_tx)
  {
    gbuffer_normal_tx.acquire(int2(depth_tx.size()), GPU_RG16F);
    gbuffer_albedo_tx.acquire(int2(depth_tx.size()), GPU_RGBA8);
    gbuffer_object_id_tx.acquire(int2(depth_tx.size()), GPU_R16UI);

    opaque_fb.ensure(GPU_ATTACHMENT_TEXTURE(depth_tx),
                     GPU_ATTACHMENT_TEXTURE(gbuffer_normal_tx),
                     GPU_ATTACHMENT_TEXTURE(gbuffer_albedo_tx),
                     GPU_ATTACHMENT_TEXTURE(gbuffer_object_id_tx));
    opaque_fb.bind();
    opaque_fb.clear_depth(1.0f);

    manager.submit(gbuffer_ps_, view);
  }

  void draw_resolve(Manager &manager, View &view)
  {
    manager.submit(deferred_ps_, view);

    gbuffer_normal_tx.release();
    gbuffer_albedo_tx.release();
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
    Texture &depth_tx = resources.depth_tx;
    DRWState state = DRW_STATE_WRITE_COLOR | DRW_STATE_WRITE_DEPTH | DRW_STATE_DEPTH_LESS_EQUAL |
                     cull_state | clip_state;

    accumulation_ps_.init(ePipelineType::TRANSPARENT, color_type, shading_type, resources, state);

    resolve_ps_.init();
    resolve_ps_.shader_set(resources.shader_cache->resolve_shader_get(ePipelineType::TRANSPARENT,
                                                                      eShadingType::FLAT));
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

struct Instance {
  SceneResources resources;

  OpaquePass opaque_ps;
  OpaquePass opaque_in_front_ps;

  TransparentPass transparent_ps;
  TransparentPass transparent_in_front_ps;

  DRWState clip_state;
  DRWState cull_state;

  bool use_per_material_batches;
  eColorType color_type;
  eShadingType shading_type;

  void init(const int2 &output_res,
            const rcti *output_rect,
            RenderEngine *render_,
            Depsgraph *depsgraph_,
            Object *camera_object_,
            const RenderLayer *render_layer_,
            const View3D *v3d_,
            const RegionView3D *rv3d_)
  {
    use_per_material_batches = false;  // ELEM(color_type, V3D_COLOR_TEXTURE, V3D_COLOR_MATERIAL);
    color_type = eColorType::MATERIAL;
    shading_type = eShadingType::STUDIO;
  }

  void begin_sync(Manager &manager)
  {
    opaque_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    opaque_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    transparent_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    transparent_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
  }

  void object_sync(Manager &manager, ObjectRef &ref)
  {
    ResourceHandle handle = manager.resource_handle(ref);

    if (use_per_material_batches) {
      Span<::Material *> materials = materials_get(ref);
      Span<GPUBatch *> batches = geometry_get(ref, materials);

      for (auto i : materials.index_range()) {
        PassMain::Sub &pass = pipeline_get(ref, materials[i]);
        pass.draw(batches[i], handle);
      }
    }
    else {
      GPUBatch *geom = geometry_get(ref);

      PassMain::Sub &pass = pipeline_get(ref);
      pass.draw(geom, handle);
    }
  }

  PassMain::Sub &pipeline_get(ObjectRef &ref, ::Material *material = nullptr)
  {
    return opaque_ps.gbuffer_ps_.sub_pass_get(geometry_type_from_object(ref.object), ref);
  }

  Span<GPUBatch *> geometry_get(ObjectRef &ref, Span<::Material *> materials)
  {
    return {};
  }

  GPUBatch *geometry_get(ObjectRef &ref)
  {
    return {};
  }

  Span<::Material *> materials_get(ObjectRef &ref)
  {
    return {};
  }

  void draw(Manager &manager, View &view)
  {
    opaque_ps.draw_prepass(manager, view, resources.depth_tx);
    // volume_ps.draw_prepass(manager, view, resources.depth_tx);
    transparent_ps.draw_prepass(manager, view, resources.depth_tx);

    if (opaque_in_front_ps.is_empty() == false || transparent_in_front_ps.is_empty() == false) {
      opaque_in_front_ps.draw_prepass(manager, view, resources.depth_in_front_tx);
      transparent_in_front_ps.draw_prepass(manager, view, resources.depth_in_front_tx);
    }

    opaque_ps.draw_resolve(manager, view);
    transparent_ps.draw_resolve(manager, view);
  }

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    this->draw(manager, view);

    // anti_aliasing_ps.draw(manager, view, depth_tx, color_tx);
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
  Depsgraph *depsgraph = ctx_state->depsgraph;
  View3D *v3d = ctx_state->v3d;
  RegionView3D *rv3d = ctx_state->rv3d;

  DefaultTextureList *dtxl = DRW_viewport_texture_list_get();
  int2 size = int2(GPU_texture_width(dtxl->color), GPU_texture_height(dtxl->color));

  Object *camera = nullptr;
  /* Get render borders. */
  rcti rect;
  BLI_rcti_init(&rect, 0, size[0], 0, size[1]);
  if (v3d) {
    if (rv3d && (rv3d->persp == RV3D_CAMOB)) {
      camera = v3d->camera;
    }
  }

  ved->instance->init(size, &rect, nullptr, depsgraph, camera, nullptr, v3d, rv3d);
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

static void workbench_cache_init(void *vedata)
{
  if (!GPU_shader_storage_buffer_objects_support()) {
    return;
  }
  draw::Manager *manager = DRW_manager_get();
  reinterpret_cast<WORKBENCH_Data *>(vedata)->instance->begin_sync(*manager);
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
  UNUSED_VARS(vedata);
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

void workbench_render_update_passes(RenderEngine *engine, Scene *scene, ViewLayer *view_layer)
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
    "BLENDER_WORKBENCH",
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
