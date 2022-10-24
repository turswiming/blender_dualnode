/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_editmesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_paint.h"
#include "BKE_particle.h"
#include "BKE_pbvh.h"
#include "DEG_depsgraph_query.h"
#include "DNA_fluid_types.h"
#include "ED_paint.h"
#include "ED_view3d.h"
#include "GPU_capabilities.h"

#include "workbench_private.hh"

namespace blender::workbench {

using namespace draw;

// Utils

const CustomData *get_loop_custom_data(const Mesh *mesh);
const CustomData *get_vert_custom_data(const Mesh *mesh);
GPUMaterial **get_dummy_gpu_materials(int material_count);

class Instance {
 public:
  Scene *scene;
  eContextObjectMode object_mode;
  View3DShading shading;
  eShadingType shading_type = eShadingType::STUDIO;
  bool xray_mode;
  bool draw_outline;
  bool draw_dof;
  bool draw_object_id;
  bool draw_transparent_depth;

  int2 resolution;
  SceneResources resources;

  DRWState cull_state;
  DRWState clip_state;
  Vector<float4> clip_planes = {};

  OpaquePass opaque_ps;
  TransparentPass transparent_ps;
  TransparentDepthPass transparent_depth_ps;

  /** Used when material_subtype == eMaterialSubType::SINGLE */
  Material material_override = Material(float3(1.0f));
  /* When r == -1.0 the shader uses the vertex color */
  Material material_attribute_color = Material(float3(-1.0f));

  DofPass dof_ps;
  AntiAliasingPass anti_aliasing_ps;

  void init(const int2 &output_res,
            const Depsgraph *depsgraph,
            const Object * /*camera*/,
            const View3D *v3d,
            const RegionView3D *rv3d)
  {
    const DRWContextState *context = DRW_context_state_get();
    scene = DEG_get_evaluated_scene(depsgraph);
    object_mode = CTX_data_mode_enum_ex(
        context->object_edit, context->obact, context->object_mode);
    resolution = output_res;
    /* TODO(Miguel Pozo):
     * Check why Workbench Next exposes OB_MATERIAL, and Workbench exposes OB_RENDER */
    bool is_render_mode = !v3d || ELEM(v3d->shading.type, OB_RENDER, OB_MATERIAL);
    const View3DShading previous_shading = shading;
    shading = is_render_mode ? scene->display.shading : v3d->shading;

    bool reset_taa = false;

    cull_state = shading.flag & V3D_SHADING_BACKFACE_CULLING ? DRW_STATE_CULL_BACK :
                                                               DRW_STATE_NO_DRAW;

    /* FIXME: This reproduce old behavior when workbench was separated in 2 engines.
     * But this is a workaround for a missing update tagging. */
    DRWState new_clip_state = RV3D_CLIPPING_ENABLED(v3d, rv3d) ? DRW_STATE_CLIP_PLANES :
                                                                 DRW_STATE_NO_DRAW;
    if (clip_state != new_clip_state) {
      clip_state = new_clip_state;
      reset_taa = true;
    }
    clip_planes.clear();
    if (clip_state & DRW_STATE_CLIP_PLANES) {
      int plane_len = (RV3D_LOCK_FLAGS(rv3d) & RV3D_BOXCLIP) ? 4 : 6;
      for (auto i : IndexRange(plane_len)) {
        clip_planes.append(rv3d->clip[i]);
      }
    }

    if (!is_render_mode) {
      if (shading.type < OB_SOLID) {
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
    xray_mode = !is_render_mode && shading.xray_alpha != 1.0f;

    if (SHADING_XRAY_FLAG_ENABLED(shading)) {
      /* Disable shading options that aren't supported in transparency mode. */
      shading.flag &= ~(V3D_SHADING_SHADOW | V3D_SHADING_CAVITY | V3D_SHADING_DEPTH_OF_FIELD);
    }
    if (SHADING_XRAY_ENABLED(shading) != SHADING_XRAY_ENABLED(previous_shading) ||
        shading.flag != previous_shading.flag) {
      reset_taa = true;
    }

    if (rv3d && rv3d->rflag & RV3D_GPULIGHT_UPDATE) {
      reset_taa = true;
    }

    shading_type = shading_type_from_v3d_lighting(shading.light);
    material_override = Material(shading.single_color);

    float4 background_color = float4(0.0f);
    if (is_render_mode && scene->r.alphamode != R_ALPHAPREMUL) {
      if (World *w = scene->world) {
        background_color = float4(w->horr, w->horg, w->horb, 1.0f);
      }
    }

    /* TODO(Miguel Pozo) volumes_do */

    resources.init(shading, scene->display, resolution, background_color);
    dof_ps.init(shading, resolution);
    anti_aliasing_ps.init(reset_taa);
    /* TODO(Miguel Pozo) taa_sample_len */

    draw_outline = shading.flag & V3D_SHADING_OBJECT_OUTLINE;
    draw_dof = dof_ps.is_enabled();
    draw_transparent_depth = draw_outline || draw_dof;
    draw_object_id = draw_outline || resources.cavity.curvature_enabled;
  }

  void begin_sync()
  {
    opaque_ps.sync(cull_state, clip_state, shading_type, resources, resolution);
    transparent_ps.sync(cull_state, clip_state, shading_type, resources);
    transparent_depth_ps.sync(cull_state, clip_state, resources);

    dof_ps.sync(resources);
    anti_aliasing_ps.sync(resources);
  }

  void end_sync()
  {
    resources.material_buf.push_update();
  }

  struct DrawModeInfo {
    eV3DShadingColorType color_type;
    bool sculpt_pbvh;
    bool texture_paint_mode;
    ::Image *image_paint_override;
    eGPUSamplerState override_sampler_state;
    bool draw_shadow;

    eColorType material_type;
    eMaterialSubType material_subtype;
    bool use_per_material_batches;
    void compute_info()
    {
      material_type = color_type_from_v3d_shading(color_type);
      material_subtype = material_subtype_from_v3d_shading(color_type);
      use_per_material_batches = image_paint_override == nullptr &&
                                 ELEM(color_type,
                                      V3D_SHADING_TEXTURE_COLOR,
                                      V3D_SHADING_MATERIAL_COLOR);
    }
  };

  const DrawModeInfo get_draw_mode_info(Object *ob)
  {
    const DRWContextState *draw_ctx = DRW_context_state_get();
    const Mesh *me = (ob->type == OB_MESH) ? static_cast<Mesh *>(ob->data) : nullptr;
    const bool is_active = (ob == draw_ctx->obact);
    /* TODO(Miguel Pozo) Is the double check needed?
     * If it is, wouldn't be needed for sculpt_pbvh too?
     */
    const bool is_render = DRW_state_is_image_render() && (draw_ctx->v3d == nullptr);

    DrawModeInfo dmi = {};
    dmi.color_type = (eV3DShadingColorType)shading.color_type;
    if (!(is_active && DRW_object_use_hide_faces(ob))) {
      dmi.draw_shadow = (ob->dtx & OB_DRAW_NO_SHADOW_CAST) == 0 &&
                        shading.flag & V3D_SHADING_SHADOW;
    }
    if (me == nullptr) {
      if (dmi.color_type == V3D_SHADING_TEXTURE_COLOR) {
        dmi.color_type = V3D_SHADING_MATERIAL_COLOR;
      }
      else if (dmi.color_type == V3D_SHADING_VERTEX_COLOR) {
        dmi.color_type = V3D_SHADING_OBJECT_COLOR;
      }
      /* Early return */
      dmi.compute_info();
      return dmi;
    }

    dmi.sculpt_pbvh = BKE_sculptsession_use_pbvh_draw(ob, draw_ctx->v3d) &&
                      !DRW_state_is_image_render();

    if (dmi.sculpt_pbvh) {
      /* Shadows are unsupported in sculpt mode. We could revert to the slow
       * method in this case but I'm not sure if it's a good idea given that
       * sculpted meshes are heavy to begin with. */
      dmi.draw_shadow = false;

      if (dmi.color_type == V3D_SHADING_TEXTURE_COLOR &&
          BKE_pbvh_type(ob->sculpt->pbvh) != PBVH_FACES) {
        /* Force use of material color for sculpt. */
        dmi.color_type = V3D_SHADING_MATERIAL_COLOR;
      }

      /* Bad call C is required to access the tool system that is context aware. Cast to non-const
       * due to current API. */
      bContext *C = (bContext *)DRW_context_state_get()->evil_C;
      if (C != NULL) {
        dmi.color_type = ED_paint_shading_color_override(
            C, &scene->toolsettings->paint_mode, ob, dmi.color_type);
      }
    }
    else {
      const CustomData *cd_vdata = get_vert_custom_data(me);
      const CustomData *cd_ldata = get_loop_custom_data(me);

      bool has_color = (CustomData_has_layer(cd_vdata, CD_PROP_COLOR) ||
                        CustomData_has_layer(cd_vdata, CD_PROP_BYTE_COLOR) ||
                        CustomData_has_layer(cd_ldata, CD_PROP_COLOR) ||
                        CustomData_has_layer(cd_ldata, CD_PROP_BYTE_COLOR));

      if (dmi.color_type == V3D_SHADING_TEXTURE_COLOR) {
        if (ob->dt < OB_TEXTURE || !CustomData_has_layer(cd_ldata, CD_MLOOPUV)) {
          dmi.color_type = V3D_SHADING_MATERIAL_COLOR;
        }
      }
      else if (dmi.color_type == V3D_SHADING_VERTEX_COLOR && !has_color) {
        dmi.color_type = V3D_SHADING_OBJECT_COLOR;
      }

      if (!is_render) {
        /* Force texture or vertex mode if object is in paint mode. */
        const bool is_vertpaint_mode = is_active && (object_mode == CTX_MODE_PAINT_VERTEX);
        const bool is_texpaint_mode = is_active && (object_mode == CTX_MODE_PAINT_TEXTURE);
        if (is_vertpaint_mode && has_color) {
          dmi.color_type = V3D_SHADING_VERTEX_COLOR;
        }
        else if (is_texpaint_mode && CustomData_has_layer(cd_ldata, CD_MLOOPUV)) {
          dmi.color_type = V3D_SHADING_TEXTURE_COLOR;
          dmi.texture_paint_mode = true;

          const ImagePaintSettings *imapaint = &scene->toolsettings->imapaint;
          if (imapaint->mode == IMAGEPAINT_MODE_IMAGE) {
            dmi.image_paint_override = imapaint->canvas;
            dmi.override_sampler_state = GPU_SAMPLER_REPEAT;
            SET_FLAG_FROM_TEST(dmi.override_sampler_state,
                               imapaint->interp == IMAGEPAINT_INTERP_LINEAR,
                               GPU_SAMPLER_FILTER);
          }
        }
      }
    }

    dmi.compute_info();
    return dmi;
  }

  void object_sync(Manager &manager, ObjectRef &ob_ref)
  {
    Object *ob = ob_ref.object;
    if (!DRW_object_is_renderable(ob)) {
      return;
    }

    const DrawModeInfo dmi = get_draw_mode_info(ob);

    /* Needed for mesh cache validation, to prevent two copies of
     * of vertex color arrays from being sent to the GPU (e.g.
     * when switching from eevee to workbench).
     */
    if (ob_ref.object->sculpt && ob_ref.object->sculpt->pbvh) {
      BKE_pbvh_is_drawing_set(ob_ref.object->sculpt->pbvh, dmi.sculpt_pbvh);
    }

    if (ob->type == OB_MESH && ob->modifiers.first != nullptr) {

      LISTBASE_FOREACH (ModifierData *, md, &ob->modifiers) {
        if (md->type != eModifierType_ParticleSystem) {
          continue;
        }
        ParticleSystem *psys = ((ParticleSystemModifierData *)md)->psys;
        if (!DRW_object_is_visible_psys_in_active_context(ob, psys)) {
          continue;
        }
        ParticleSettings *part = psys->part;
        const int draw_as = (part->draw_as == PART_DRAW_REND) ? part->ren_as : part->draw_as;

        if (draw_as == PART_DRAW_PATH) {
          /* TODO(Miguel Pozo):
          workbench_cache_hair_populate(
              wpd, ob, psys, md, dmi.color_type, dmi.texture_paint_mode, part->omat);
          */
        }
      }
    }

    if (!(ob->base_flag & BASE_FROM_DUPLI)) {
      ModifierData *md = BKE_modifiers_findby_type(ob, eModifierType_Fluid);
      if (md && BKE_modifier_is_enabled(scene, md, eModifierMode_Realtime)) {
        FluidModifierData *fmd = (FluidModifierData *)md;
        if (fmd->domain) {
          /* TODO(Miguel Pozo):
          workbench_volume_cache_populate(vedata, wpd->scene, ob, md, V3D_SHADING_SINGLE_COLOR);
          */
          if (fmd->domain->type == FLUID_DOMAIN_TYPE_GAS) {
            return; /* Do not draw solid in this case. */
          }
        }
      }
    }

    if (!(DRW_object_visibility_in_active_context(ob) & OB_VISIBLE_SELF)) {
      return;
    }

    if ((ob->dt < OB_SOLID) && !DRW_state_is_scene_render()) {
      return;
    }

    if (ELEM(ob->type, OB_MESH, OB_POINTCLOUD)) {
      mesh_sync(manager, ob_ref, dmi);
    }
    else if (ob->type == OB_CURVES) {
      /* TODO(Miguel Pozo):
      DRWShadingGroup *grp = workbench_material_hair_setup(
          wpd, ob, CURVES_MATERIAL_NR, dmi.color_type);
      DRW_shgroup_curves_create_sub(ob, grp, NULL);
      */
    }
    else if (ob->type == OB_VOLUME) {
      if (shading.type != OB_WIRE) {
        /* TODO(Miguel Pozo):
        workbench_volume_cache_populate(vedata, wpd->scene, ob, NULL, dmi.color_type);
        */
      }
    }
  }

  void mesh_sync(Manager &manager, ObjectRef &ob_ref, const DrawModeInfo dmi)
  {
    if (dmi.sculpt_pbvh) {
      /* TODO(Miguel Pozo):
      workbench_cache_sculpt_populate(wpd, ob, dmi.color_type);
      */
    }
    else {
      /* workbench_cache_common_populate && workbench_cache_texpaint_populate */
      if (dmi.use_per_material_batches) {
        const int material_count = DRW_cache_object_material_count_get(ob_ref.object);

        struct GPUBatch **batches;
        if (dmi.material_type == eColorType::TEXTURE) {
          batches = DRW_cache_mesh_surface_texpaint_get(ob_ref.object);
        }
        else {
          batches = DRW_cache_object_surface_material_get(
              ob_ref.object, get_dummy_gpu_materials(material_count), material_count);
        }

        if (batches) {
          for (auto i : IndexRange(material_count)) {
            if (batches[i] == nullptr) {
              continue;
            }
            /* TODO(fclem): This create a cull-able instance for each sub-object. This is done
             * for simplicity to reduce complexity. But this increase the overhead per object.
             * Instead, we should use an indirection buffer to the material buffer. */

            ResourceHandle handle = manager.resource_handle(ob_ref);
            Material &mat = resources.material_buf.get_or_resize(handle.resource_index());

            if (::Material *_mat = BKE_object_material_get_eval(ob_ref.object, i + 1)) {
              mat = Material(*_mat);
            }
            else {
              mat = Material(*BKE_material_default_empty());
            }

            ::Image *image = nullptr;
            ImageUser *iuser = nullptr;
            eGPUSamplerState sampler_state = eGPUSamplerState::GPU_SAMPLER_DEFAULT;
            if (dmi.material_type == eColorType::TEXTURE) {
              get_material_image(ob_ref.object, i + 1, image, iuser, sampler_state);
            }

            draw_mesh(ob_ref, mat, batches[i], handle, image, sampler_state, iuser);
          }
        }
      }
      else {
        struct GPUBatch *batch;
        if (dmi.material_type == eColorType::TEXTURE) {
          batch = DRW_cache_mesh_surface_texpaint_single_get(ob_ref.object);
        }
        else if (dmi.material_subtype == eMaterialSubType::ATTRIBUTE) {
          if (ob_ref.object->mode & OB_MODE_VERTEX_PAINT) {
            batch = DRW_cache_mesh_surface_vertpaint_get(ob_ref.object);
          }
          else {
            batch = DRW_cache_mesh_surface_sculptcolors_get(ob_ref.object);
          }
        }
        else {
          batch = DRW_cache_object_surface_get(ob_ref.object);
        }

        if (batch) {
          ResourceHandle handle = manager.resource_handle(ob_ref);
          Material &mat = resources.material_buf.get_or_resize(handle.resource_index());

          if (dmi.material_subtype == eMaterialSubType::OBJECT) {
            mat = Material(*ob_ref.object);
          }
          else if (dmi.material_subtype == eMaterialSubType::RANDOM) {
            mat = Material(*ob_ref.object, true);
          }
          else if (dmi.material_subtype == eMaterialSubType::SINGLE) {
            mat = material_override;
          }
          else if (dmi.material_subtype == eMaterialSubType::ATTRIBUTE) {
            mat = material_attribute_color;
          }
          else {
            mat = Material(*BKE_material_default_empty());
          }

          draw_mesh(
              ob_ref, mat, batch, handle, dmi.image_paint_override, dmi.override_sampler_state);
        }
      }
    }

    if (dmi.draw_shadow) {
      /* TODO(Miguel Pozo):
      workbench_shadow_cache_populate(vedata, ob, has_transp_mat);
      */
    }
  }

  void draw_mesh(ObjectRef &ob_ref,
                 Material &material,
                 GPUBatch *batch,
                 ResourceHandle handle,
                 ::Image *image = nullptr,
                 eGPUSamplerState sampler_state = GPU_SAMPLER_DEFAULT,
                 ImageUser *iuser = nullptr)
  {
    const bool in_front = (ob_ref.object->dtx & OB_DRAW_IN_FRONT) != 0;

    auto draw = [&](MeshPass &pass) {
      pass.sub_pass_get(ob_ref, image, sampler_state, iuser).draw(batch, handle);
    };

    if (xray_mode || material.is_transparent()) {
      if (in_front) {
        draw(transparent_ps.accumulation_in_front_ps_);
        if (draw_transparent_depth) {
          draw(transparent_depth_ps.in_front_ps_);
        }
      }
      else {
        draw(transparent_ps.accumulation_ps_);
        if (draw_transparent_depth) {
          draw(transparent_depth_ps.main_ps_);
        }
      }
    }
    else {
      if (in_front) {
        draw(opaque_ps.gbuffer_in_front_ps_);
      }
      else {
        draw(opaque_ps.gbuffer_ps_);
      }
    }
  }

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    if (!clip_planes.is_empty()) {
      view.set_clip_planes(clip_planes);
    }

    resources.color_tx.acquire(resolution, GPU_RGBA16F);
    resources.color_tx.clear(resources.world_buf.background_color);
    if (draw_object_id) {
      resources.object_id_tx.acquire(resolution, GPU_R16UI);
      resources.object_id_tx.clear(uint4(0));
    }

    resources.depth_tx.acquire(resolution, GPU_DEPTH24_STENCIL8);
    Framebuffer fb = Framebuffer("Workbench.Clear");
    fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx));
    fb.bind();
    GPU_framebuffer_clear_depth_stencil(fb, 1.0f, 0x00);

    if (!opaque_ps.gbuffer_in_front_ps_.is_empty() ||
        !transparent_ps.accumulation_in_front_ps_.is_empty()) {
      resources.depth_in_front_tx.acquire(resolution, GPU_DEPTH24_STENCIL8);
      Framebuffer fb = Framebuffer("Workbench.Clear");
      fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_in_front_tx));
      fb.bind();
      GPU_framebuffer_clear_depth_stencil(fb, 1.0f, 0x00);
    }

    opaque_ps.draw(manager, view, resources, resolution);
    transparent_ps.draw(manager, view, resources, resolution);
    transparent_depth_ps.draw(manager, view, resources, resolution);

    if (draw_outline) {
      PassSimple outline_ps = PassSimple("Workbench.Outline");
      outline_ps.state_set(DRW_STATE_WRITE_COLOR | DRW_STATE_BLEND_ALPHA_PREMUL);
      static GPUShader *outline_shader = GPU_shader_create_from_info_name(
          "workbench_effect_outline");
      outline_ps.shader_set(outline_shader);
      outline_ps.bind_ubo("world_data", resources.world_buf);
      outline_ps.bind_texture("objectIdBuffer", &resources.object_id_tx);
      outline_ps.draw_procedural(GPU_PRIM_TRIS, 1, 3);

      Framebuffer fb = Framebuffer("Workbench.Outline");
      fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(resources.color_tx));
      fb.bind();
      manager.submit(outline_ps);
    }

    // volume_ps.draw_prepass(manager, view, resources.depth_tx);

    dof_ps.draw(manager, view, resources, resolution);
    anti_aliasing_ps.draw(manager, view, depth_tx, color_tx);

    resources.color_tx.release();
    resources.object_id_tx.release();
    resources.depth_tx.release();
    resources.depth_in_front_tx.release();
  }

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    this->draw(manager, view, depth_tx, color_tx);
  }
};

// Utils

const CustomData *get_loop_custom_data(const Mesh *mesh)
{
  if (mesh->runtime.wrapper_type == ME_WRAPPER_TYPE_BMESH) {
    BLI_assert(mesh->edit_mesh != nullptr);
    BLI_assert(mesh->edit_mesh->bm != nullptr);
    return &mesh->edit_mesh->bm->ldata;
  }
  return &mesh->ldata;
}

const CustomData *get_vert_custom_data(const Mesh *mesh)
{
  if (mesh->runtime.wrapper_type == ME_WRAPPER_TYPE_BMESH) {
    BLI_assert(mesh->edit_mesh != nullptr);
    BLI_assert(mesh->edit_mesh->bm != nullptr);
    return &mesh->edit_mesh->bm->vdata;
  }
  return &mesh->vdata;
}

/* This returns an array of nullptr GPUMaterial pointers so we can call
 * DRW_cache_object_surface_material_get. They never get actually used.
 */
GPUMaterial **get_dummy_gpu_materials(int material_count)
{
  static Vector<GPUMaterial *> dummy_gpu_materials(1, nullptr, {});
  if (material_count > dummy_gpu_materials.size()) {
    dummy_gpu_materials.resize(material_count, nullptr);
  }
  return dummy_gpu_materials.begin();
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
