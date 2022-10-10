/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_studiolight.h"
#include "DEG_depsgraph_query.h"
#include "ED_view3d.h"
#include "GPU_capabilities.h"

#include "workbench_private.hh"

namespace blender::workbench {

using namespace draw;

class Instance {
 public:
  SceneResources resources;

  OpaquePass opaque_ps;
  // OpaquePass opaque_in_front_ps;

  // TransparentPass transparent_ps;
  // TransparentPass transparent_in_front_ps;

  AntiAliasingPass anti_aliasing_ps;

  bool use_per_material_batches = false;
  eColorType color_type = eColorType::MATERIAL;
  eMaterialSubType material_subtype = eMaterialSubType::MATERIAL;
  eShadingType shading_type = eShadingType::STUDIO;

  /** Used when material_subtype == eMaterialSubType::SINGLE */
  Material material_override = Material(float3(1.0f));
  /* When r == -1.0 the shader uses the vertex color */
  Material material_attribute_color = Material(float3(-1.0f));

  DRWState cull_state;
  DRWState clip_state;
  bool view_updated; /* TODO(pragma37): move to AntialiasingPass ? */

  eContextObjectMode ob_mode;
  eGPUShaderConfig clip_mode;

  View3DShading shading;

  void init(const int2 &output_res,
            const Depsgraph *depsgraph,
            const Object * /*camera*/,
            const View3D *v3d,
            const RegionView3D *rv3d)
  {
    Scene *scene = DEG_get_evaluated_scene(depsgraph);
    const DRWContextState *context = DRW_context_state_get();

    /* TODO(pragma37):
     * Check why Workbench Next exposes OB_MATERIAL, and Workbench exposes OB_RENDER */
    bool is_render_mode = !v3d || ELEM(v3d->shading.type, OB_RENDER, OB_MATERIAL);
    // const View3DShading &shading = is_render_mode ? scene->display.shading : v3d->shading;
    const View3DShading previous_shading = shading;
    shading = is_render_mode ? scene->display.shading : v3d->shading;

    ob_mode = CTX_data_mode_enum_ex(context->object_edit, context->obact, context->object_mode);
    clip_mode = context->sh_cfg;

    cull_state = shading.flag & V3D_SHADING_BACKFACE_CULLING ? DRW_STATE_CULL_BACK :
                                                               DRW_STATE_NO_DRAW;

    bool reset_taa = false;

    /* FIXME: This reproduce old behavior when workbench was separated in 2 engines.
     * But this is a workaround for a missing update tagging. */
    DRWState new_clip_state = RV3D_CLIPPING_ENABLED(v3d, rv3d) ? DRW_STATE_CLIP_PLANES :
                                                                 DRW_STATE_NO_DRAW;
    if (clip_state != new_clip_state) {
      reset_taa = true;
    }
    clip_state = new_clip_state;

    if (rv3d && rv3d->rflag & RV3D_GPULIGHT_UPDATE) {
      reset_taa = true;
    }

    if (SHADING_XRAY_FLAG_ENABLED(shading)) {
      /* Disable shading options that aren't supported in transparency mode. */
      shading.flag &= ~(V3D_SHADING_SHADOW | V3D_SHADING_CAVITY | V3D_SHADING_DEPTH_OF_FIELD);
    }
    if (SHADING_XRAY_ENABLED(shading) != SHADING_XRAY_ENABLED(previous_shading) ||
        shading.flag != previous_shading.flag) {
      reset_taa = true;
    }

    if (!is_render_mode) {
      if (shading.type < OB_SOLID) {
        /* TODO(pragma37): Shouldn't we just skip any rendering at all ??? */
        shading.light = V3D_LIGHTING_FLAT;
        shading.color_type = V3D_SHADING_OBJECT_COLOR;
        shading.xray_alpha = 0.0f;
      }
      else if (SHADING_XRAY_ENABLED(shading)) {
        shading.xray_alpha = SHADING_XRAY_ALPHA(shading);
      }
      else {
        shading.xray_alpha = 1.0f;
      }
    }

    material_override = Material(shading.single_color);

    use_per_material_batches = ELEM(
        shading.color_type, V3D_SHADING_TEXTURE_COLOR, V3D_SHADING_MATERIAL_COLOR);

    color_type = color_type_from_v3d_shading(shading.color_type);
    material_subtype = material_subtype_from_v3d_shading(shading.color_type);
    shading_type = shading_type_from_v3d_lighting(shading.light);

    UniformBuffer<WorldData> &world_buf = resources.world_buf;

    float4x4 rot_matrix = float4x4::identity();
    if (shading.flag & V3D_SHADING_WORLD_ORIENTATION) {
      /* TODO(pragma37) C++ API ? */
      float V[4][4], R[4][4];
      DRW_view_viewmat_get(nullptr, V, false);
      axis_angle_to_mat4_single(R, 'Z', -shading.studiolight_rot_z);
      mul_m4_m4m4(R, V, R);
      swap_v3_v3(R[2], R[1]);
      negate_v3(R[2]);
      rot_matrix = float4x4(R);
    }

    StudioLight *studio_light = nullptr;
    if (U.edit_studio_light) {
      studio_light = BKE_studiolight_studio_edit_get();
    }
    else {
      if (shading_type == eShadingType::MATCAP) {
        studio_light = BKE_studiolight_find(shading.matcap, STUDIOLIGHT_TYPE_MATCAP);
      }
      /* If matcaps are missing, use this as fallback. */
      if (studio_light == nullptr) {
        studio_light = BKE_studiolight_find(shading.studio_light, STUDIOLIGHT_TYPE_STUDIO);
      }
    }

    for (int i = 0; i < 4; i++) {
      LightData &light = world_buf.lights[i];

      SolidLight *sl = (studio_light) ? &studio_light->light[i] : nullptr;
      if (sl && sl->flag) {
        float3 direction = rot_matrix.ref_3x3() * float3(sl->vec);
        light.direction = float4(direction, 0.0f);
        /* We should pre-divide the power by PI but that makes the lights really dim. */
        light.specular_color = float4(float3(sl->spec), 0.0f);
        light.diffuse_color_wrap = float4(float3(sl->col), sl->smooth);
      }
      else {
        light.direction = float4(1.0f, 0.0f, 0.0f, 0.0f);
        light.specular_color = float4(0.0f);
        light.diffuse_color_wrap = float4(0.0f);
      }
    }

    world_buf.ambient_color = float4(1.0f, 1.0f, 1.0f, 0.0f);
    world_buf.use_specular = false;

    if (studio_light != nullptr) {
      world_buf.ambient_color = float4(float3(studio_light->light_ambient), 0.0f);
      world_buf.use_specular = shading.flag & V3D_SHADING_SPECULAR_HIGHLIGHT &&
                               studio_light->flag & STUDIOLIGHT_SPECULAR_HIGHLIGHT_PASS;
    }

    world_buf.background_color = float4(0.0f);

    if (is_render_mode && scene->r.alphamode != R_ALPHAPREMUL) {
      if (World *w = scene->world) {
        world_buf.background_color = float4(w->horr, w->horg, w->horb, 1.0f);
      }
    }

    world_buf.object_outline_color = shading.object_outline_color;
    world_buf.object_outline_color.w = 1.0f;
    world_buf.ui_scale = DRW_state_is_image_render() ? 1.0f : G_draw.block.size_pixel;
    world_buf.matcap_orientation = (shading.flag & V3D_SHADING_MATCAP_FLIP_X) != 0;

    /* TODO(pragma37) volumes_do */

    resources.matcap_tx.ensure_2d_array(GPU_RGBA16F, int2(1), 1);
    resources.depth_tx.ensure_2d(GPU_DEPTH24_STENCIL8, output_res);

    anti_aliasing_ps.init(reset_taa);
    /* TODO(pragma37) taa_sample_len */
  }

  void begin_sync()
  {
    resources.world_buf.push_update();
    opaque_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // opaque_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // transparent_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    // transparent_in_front_ps.sync(cull_state, clip_state, shading_type, color_type, resources);
    anti_aliasing_ps.sync(resources);
  }

  void end_sync()
  {
    resources.material_buf.push_update();
  }

  void object_sync(Manager &manager, ObjectRef &ob_ref)
  {
    if (ob_ref.object->type != OB_MESH) {
      // TODO(pragma37)
      return;
    }

    if (use_per_material_batches) {
      const int material_count = DRW_cache_object_material_count_get(ob_ref.object);
      Span<GPUBatch *> batches = geometry_get(ob_ref, material_count);
      /* TODO(pragma37): Could this ever be false??? */
      if (batches.size() == material_count) {
        for (auto i : IndexRange(material_count)) {
          /* TODO(fclem): This create a cull-able instance for each sub-object. This is done for
           * simplicity to reduce complexity. But this increase the overhead per object. Instead,
           * we should use an indirection buffer to the material buffer. */
          ::Material *mat = BKE_object_material_get_eval(ob_ref.object, i + 1);
          if (mat == nullptr) {
            mat = BKE_material_default_empty();
          }
          ResourceHandle handle = manager.resource_handle(ob_ref);
          resources.material_buf.get_or_resize(handle.resource_index()) = Material(*mat);
          pipeline_get(ob_ref, mat).draw(batches[i], handle);
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
      else if (material_subtype == eMaterialSubType::SINGLE) {
        mat = material_override;
      }
      else if (material_subtype == eMaterialSubType::ATTRIBUTE) {
        mat = material_attribute_color;
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

  Span<GPUBatch *> geometry_get(ObjectRef &ob_ref, int material_count)
  {
    /* This is never used, but it's required by DRW_cache_object_surface_material_get */
    static Vector<GPUMaterial *> dummy_gpu_materials(1, nullptr, {});
    if (material_count > dummy_gpu_materials.size()) {
      dummy_gpu_materials.resize(material_count, nullptr);
    }
    return {DRW_cache_object_surface_material_get(
                ob_ref.object, dummy_gpu_materials.begin(), material_count),
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

  /* FIXME: This reproduce old behavior when workbench was separated in 2 engines.
   * But this is a workaround for a missing update tagging. */
  rv3d->rflag &= ~RV3D_GPULIGHT_UPDATE;
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
