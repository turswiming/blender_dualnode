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
GPUMaterial **get_dummy_gpu_materials(int material_count);

class Instance {
 public:
  SceneState scene_state;

  SceneResources resources;

  OpaquePass opaque_ps;
  TransparentPass transparent_ps;
  TransparentDepthPass transparent_depth_ps;

  OutlinePass outline_ps;
  DofPass dof_ps;
  AntiAliasingPass anti_aliasing_ps;

  void init()
  {
    scene_state.init();
    resources.init(scene_state);

    outline_ps.init(scene_state);
    dof_ps.init(scene_state);
    anti_aliasing_ps.init(scene_state);
  }

  void begin_sync()
  {
    opaque_ps.sync(scene_state, resources);
    transparent_ps.sync(scene_state, resources);
    transparent_depth_ps.sync(scene_state, resources);

    outline_ps.sync(resources);
    dof_ps.sync(resources);
    anti_aliasing_ps.sync(resources, scene_state.resolution);
  }

  void end_sync()
  {
    resources.material_buf.push_update();
  }

  void object_sync(Manager &manager, ObjectRef &ob_ref)
  {
    if (scene_state.render_finished) {
      return;
    }

    Object *ob = ob_ref.object;
    if (!DRW_object_is_renderable(ob)) {
      return;
    }

    const ObjectState object_state = ObjectState(scene_state, ob);

    /* Needed for mesh cache validation, to prevent two copies of
     * of vertex color arrays from being sent to the GPU (e.g.
     * when switching from eevee to workbench).
     */
    if (ob_ref.object->sculpt && ob_ref.object->sculpt->pbvh) {
      BKE_pbvh_is_drawing_set(ob_ref.object->sculpt->pbvh, object_state.sculpt_pbvh);
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
              wpd, ob, psys, md, object_state.color_type, object_state.texture_paint_mode,
          part->omat);
          */
        }
      }
    }

    if (!(ob->base_flag & BASE_FROM_DUPLI)) {
      ModifierData *md = BKE_modifiers_findby_type(ob, eModifierType_Fluid);
      if (md && BKE_modifier_is_enabled(scene_state.scene, md, eModifierMode_Realtime)) {
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
      mesh_sync(manager, ob_ref, object_state);
    }
    else if (ob->type == OB_CURVES) {
      /* TODO(Miguel Pozo):
      DRWShadingGroup *grp = workbench_material_hair_setup(
          wpd, ob, CURVES_MATERIAL_NR, object_state.color_type);
      DRW_shgroup_curves_create_sub(ob, grp, NULL);
      */
    }
    else if (ob->type == OB_VOLUME) {
      if (scene_state.shading.type != OB_WIRE) {
        /* TODO(Miguel Pozo):
        workbench_volume_cache_populate(vedata, wpd->scene, ob, NULL, object_state.color_type);
        */
      }
    }
  }

  void mesh_sync(Manager &manager, ObjectRef &ob_ref, const ObjectState &object_state)
  {
    if (object_state.sculpt_pbvh) {
      /* TODO(Miguel Pozo):
      workbench_cache_sculpt_populate(wpd, ob, object_state.color_type);
      */
    }
    else {
      /* workbench_cache_common_populate && workbench_cache_texpaint_populate */
      if (object_state.use_per_material_batches) {
        const int material_count = DRW_cache_object_material_count_get(ob_ref.object);

        struct GPUBatch **batches;
        if (object_state.color_type == V3D_SHADING_TEXTURE_COLOR) {
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
            if (object_state.color_type == V3D_SHADING_TEXTURE_COLOR) {
              get_material_image(ob_ref.object, i + 1, image, iuser, sampler_state);
            }

            draw_mesh(ob_ref, mat, batches[i], handle, image, sampler_state, iuser);
          }
        }
      }
      else {
        struct GPUBatch *batch;
        if (object_state.color_type == V3D_SHADING_TEXTURE_COLOR) {
          batch = DRW_cache_mesh_surface_texpaint_single_get(ob_ref.object);
        }
        else if (object_state.color_type == V3D_SHADING_VERTEX_COLOR) {
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

          if (object_state.color_type == V3D_SHADING_OBJECT_COLOR) {
            mat = Material(*ob_ref.object);
          }
          else if (object_state.color_type == V3D_SHADING_RANDOM_COLOR) {
            mat = Material(*ob_ref.object, true);
          }
          else if (object_state.color_type == V3D_SHADING_SINGLE_COLOR) {
            mat = scene_state.material_override;
          }
          else if (object_state.color_type == V3D_SHADING_VERTEX_COLOR) {
            mat = scene_state.material_attribute_color;
          }
          else {
            mat = Material(*BKE_material_default_empty());
          }

          draw_mesh(ob_ref,
                    mat,
                    batch,
                    handle,
                    object_state.image_paint_override,
                    object_state.override_sampler_state);
        }
      }
    }

    if (object_state.draw_shadow) {
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

    if (scene_state.xray_mode || material.is_transparent()) {
      if (in_front) {
        draw(transparent_ps.accumulation_in_front_ps_);
        if (scene_state.draw_transparent_depth) {
          draw(transparent_depth_ps.in_front_ps_);
        }
      }
      else {
        draw(transparent_ps.accumulation_ps_);
        if (scene_state.draw_transparent_depth) {
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
    int2 resolution = scene_state.resolution;

    if (scene_state.render_finished) {
      /* Just copy back the already rendered result */
      anti_aliasing_ps.draw(manager, view, resources, resolution, depth_tx, color_tx);
      return;
    }

    anti_aliasing_ps.setup_view(view, resolution);

    if (!scene_state.clip_planes.is_empty()) {
      view.set_clip_planes(scene_state.clip_planes);
    }

    resources.color_tx.acquire(resolution, GPU_RGBA16F);
    resources.color_tx.clear(resources.world_buf.background_color);
    if (scene_state.draw_object_id) {
      resources.object_id_tx.acquire(resolution, GPU_R16UI);
      resources.object_id_tx.clear(uint4(0));
    }

    resources.depth_tx.acquire(resolution, GPU_DEPTH24_STENCIL8);
    Framebuffer fb = Framebuffer("Workbench.Clear");
    fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx));
    fb.bind();
    GPU_framebuffer_clear_depth_stencil(fb, 1.0f, 0x00);

    if (!transparent_ps.accumulation_in_front_ps_.is_empty()) {
      resources.depth_in_front_tx.acquire(resolution, GPU_DEPTH24_STENCIL8);
      if (opaque_ps.gbuffer_in_front_ps_.is_empty()) {
        /* Clear only if it wont be overwitten by opaque_ps */
        Framebuffer fb = Framebuffer("Workbench.Clear");
        fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_in_front_tx));
        fb.bind();
        GPU_framebuffer_clear_depth_stencil(fb, 1.0f, 0x00);
      }
    }

    opaque_ps.draw(manager, view, resources, resolution);
    transparent_ps.draw(manager, view, resources, resolution);
    transparent_depth_ps.draw(manager, view, resources, resolution);

    // volume_ps.draw_prepass(manager, view, resources.depth_tx);

    outline_ps.draw(manager, view, resources, resolution);
    dof_ps.draw(manager, view, resources, resolution);
    anti_aliasing_ps.draw(manager, view, resources, resolution, depth_tx, color_tx);

    resources.color_tx.release();
    resources.object_id_tx.release();
    resources.depth_tx.release();
    resources.depth_in_front_tx.release();
  }

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    this->draw(manager, view, depth_tx, color_tx);

    if (scene_state.sample + 1 < scene_state.samples_len) {
      DRW_viewport_request_redraw();
    }
  }
};

// Utils

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

  ved->instance->init();
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
  WORKBENCH_Data *ved = reinterpret_cast<WORKBENCH_Data *>(vedata);
  if (ved->instance) {
    ved->instance->scene_state.reset_taa_next_sample = true;
  }
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
