/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_editmesh.h"
#include "BKE_modifier.h"
#include "BKE_object.h"
#include "BKE_paint.h"
#include "BKE_particle.h"
#include "BKE_pbvh.h"
#include "BKE_studiolight.h"
#include "DEG_depsgraph_query.h"
#include "DNA_fluid_types.h"
#include "ED_paint.h"
#include "ED_view3d.h"
#include "GPU_capabilities.h"
#include "IMB_imbuf_types.h"

#include "workbench_private.hh"

namespace blender::workbench {

using namespace draw;

// Utils

bool get_matcap_tx(Texture &matcap_tx, const StudioLight &studio_light);
float4x4 get_world_shading_rotation_matrix(float studiolight_rot_z);
LightData get_light_data_from_studio_solidlight(const SolidLight *sl,
                                                float4x4 world_shading_rotation);

const CustomData *get_loop_custom_data(const Mesh *mesh);
const CustomData *get_vert_custom_data(const Mesh *mesh);

GPUMaterial **get_dummy_gpu_materials(int material_count);

class Instance {
 public:
  SceneResources resources;

  OpaquePass opaque_ps = OpaquePass(resources.color_tx, resources.depth_tx);
  OpaquePass opaque_in_front_ps = OpaquePass(resources.color_tx, resources.depth_in_front_tx);
  TransparentPass transparent_ps = TransparentPass(resources.color_tx, resources.depth_tx);
  TransparentPass transparent_in_front_ps = TransparentPass(resources.color_tx,
                                                            resources.depth_in_front_tx);

  AntiAliasingPass anti_aliasing_ps;

  eShadingType shading_type = eShadingType::STUDIO;

  /** Used when material_subtype == eMaterialSubType::SINGLE */
  Material material_override = Material(float3(1.0f));
  /* When r == -1.0 the shader uses the vertex color */
  Material material_attribute_color = Material(float3(-1.0f));

  DRWState cull_state;
  DRWState clip_state;

  Vector<float4> clip_planes = {};

  eContextObjectMode ob_mode;

  View3DShading shading;
  bool xray_mode;

  StringRefNull current_matcap;
  Scene *scene;

  void init(const int2 &output_res,
            const Depsgraph *depsgraph,
            const Object * /*camera*/,
            const View3D *v3d,
            const RegionView3D *rv3d)
  {
    const DRWContextState *context = DRW_context_state_get();
    scene = DEG_get_evaluated_scene(depsgraph);

    /* TODO(pragma37):
     * Check why Workbench Next exposes OB_MATERIAL, and Workbench exposes OB_RENDER */
    bool is_render_mode = !v3d || ELEM(v3d->shading.type, OB_RENDER, OB_MATERIAL);
    const View3DShading previous_shading = shading;
    shading = is_render_mode ? scene->display.shading : v3d->shading;

    ob_mode = CTX_data_mode_enum_ex(context->object_edit, context->obact, context->object_mode);

    bool reset_taa = false;

    UniformBuffer<WorldData> &world_buf = resources.world_buf;

    world_buf.viewport_size = DRW_viewport_size_get();
    world_buf.viewport_size_inv = DRW_viewport_invert_size_get();

    cull_state = shading.flag & V3D_SHADING_BACKFACE_CULLING ? DRW_STATE_CULL_BACK :
                                                               DRW_STATE_NO_DRAW;

    /* FIXME: This reproduce old behavior when workbench was separated in 2 engines.
     * But this is a workaround for a missing update tagging. */
    DRWState new_clip_state = RV3D_CLIPPING_ENABLED(v3d, rv3d) ? DRW_STATE_CLIP_PLANES :
                                                                 DRW_STATE_NO_DRAW;
    if (clip_state != new_clip_state) {
      reset_taa = true;
    }
    clip_state = new_clip_state;

    clip_planes.clear();
    if (clip_state & DRW_STATE_CLIP_PLANES) {
      int plane_len = (RV3D_LOCK_FLAGS(rv3d) & RV3D_BOXCLIP) ? 4 : 6;
      for (auto i : IndexRange(plane_len)) {
        clip_planes.append(rv3d->clip[i]);
      }
    }

    if (rv3d && rv3d->rflag & RV3D_GPULIGHT_UPDATE) {
      reset_taa = true;
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
    world_buf.xray_alpha = shading.xray_alpha;

    xray_mode = !is_render_mode && shading.xray_alpha != 1.0f;

    if (SHADING_XRAY_FLAG_ENABLED(shading)) {
      /* Disable shading options that aren't supported in transparency mode. */
      shading.flag &= ~(V3D_SHADING_SHADOW | V3D_SHADING_CAVITY | V3D_SHADING_DEPTH_OF_FIELD);
    }
    if (SHADING_XRAY_ENABLED(shading) != SHADING_XRAY_ENABLED(previous_shading) ||
        shading.flag != previous_shading.flag) {
      reset_taa = true;
    }

    shading_type = shading_type_from_v3d_lighting(shading.light);
    material_override = Material(shading.single_color);

    StudioLight *studio_light = nullptr;
    if (U.edit_studio_light) {
      studio_light = BKE_studiolight_studio_edit_get();
    }
    else {
      if (shading_type == eShadingType::MATCAP) {
        studio_light = BKE_studiolight_find(shading.matcap, STUDIOLIGHT_TYPE_MATCAP);
        if (studio_light && studio_light->name != current_matcap) {
          if (get_matcap_tx(resources.matcap_tx, *studio_light)) {
            current_matcap = studio_light->name;
          }
        }
      }
      /* If matcaps are missing, use this as fallback. */
      if (studio_light == nullptr) {
        studio_light = BKE_studiolight_find(shading.studio_light, STUDIOLIGHT_TYPE_STUDIO);
      }
    }
    if (!resources.matcap_tx.is_valid()) {
      resources.matcap_tx.ensure_2d_array(GPU_RGBA16F, int2(1), 1);
    }
    world_buf.matcap_orientation = (shading.flag & V3D_SHADING_MATCAP_FLIP_X) != 0;

    float4x4 world_shading_rotation = float4x4::identity();
    if (shading.flag & V3D_SHADING_WORLD_ORIENTATION) {
      world_shading_rotation = get_world_shading_rotation_matrix(shading.studiolight_rot_z);
    }

    for (int i = 0; i < 4; i++) {
      SolidLight *sl = (studio_light) ? &studio_light->light[i] : nullptr;
      world_buf.lights[i] = get_light_data_from_studio_solidlight(sl, world_shading_rotation);
    }

    if (studio_light != nullptr) {
      world_buf.ambient_color = float4(float3(studio_light->light_ambient), 0.0f);
      world_buf.use_specular = shading.flag & V3D_SHADING_SPECULAR_HIGHLIGHT &&
                               studio_light->flag & STUDIOLIGHT_SPECULAR_HIGHLIGHT_PASS;
    }
    else {
      world_buf.ambient_color = float4(1.0f, 1.0f, 1.0f, 0.0f);
      world_buf.use_specular = false;
    }

    world_buf.background_color = float4(0.0f);
    if (is_render_mode && scene->r.alphamode != R_ALPHAPREMUL) {
      if (World *w = scene->world) {
        world_buf.background_color = float4(w->horr, w->horg, w->horb, 1.0f);
      }
    }

    world_buf.object_outline_color = float4(float3(shading.object_outline_color), 1.0f);
    world_buf.ui_scale = DRW_state_is_image_render() ? 1.0f : G_draw.block.size_pixel;

    /* TODO(pragma37) volumes_do */

    resources.depth_tx.ensure_2d(GPU_DEPTH24_STENCIL8, output_res);
    resources.depth_in_front_tx.ensure_2d(GPU_DEPTH24_STENCIL8, output_res);

    anti_aliasing_ps.init(reset_taa);
    /* TODO(pragma37) taa_sample_len */
  }

  void begin_sync()
  {
    resources.world_buf.push_update();

    opaque_ps.sync(cull_state, clip_state, shading_type, resources);
    opaque_in_front_ps.sync(cull_state, clip_state, shading_type, resources);

    transparent_ps.sync(cull_state, clip_state, shading_type, resources);
    transparent_in_front_ps.sync(cull_state, clip_state, shading_type, resources);

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
    /* TODO(pragma37) Is the double check needed?
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
        const bool is_vertpaint_mode = is_active && (ob_mode == CTX_MODE_PAINT_VERTEX);
        const bool is_texpaint_mode = is_active && (ob_mode == CTX_MODE_PAINT_TEXTURE);
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
          /* TODO(pragma37):
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
          /* TODO(pragma37):
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
      /* TODO(pragma37):
      DRWShadingGroup *grp = workbench_material_hair_setup(
          wpd, ob, CURVES_MATERIAL_NR, dmi.color_type);
      DRW_shgroup_curves_create_sub(ob, grp, NULL);
      */
    }
    else if (ob->type == OB_VOLUME) {
      if (shading.type != OB_WIRE) {
        /* TODO(pragma37):
        workbench_volume_cache_populate(vedata, wpd->scene, ob, NULL, dmi.color_type);
        */
      }
    }
  }

  void mesh_sync(Manager &manager, ObjectRef &ob_ref, const DrawModeInfo dmi)
  {
    if (dmi.sculpt_pbvh) {
      /* TODO(pragma37):
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

            pipeline_get(ob_ref, mat, image, sampler_state, iuser).draw(batches[i], handle);
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

          pipeline_get(ob_ref, mat, dmi.image_paint_override, dmi.override_sampler_state)
              .draw(batch, handle);
        }
      }
    }

    if (dmi.draw_shadow) {
      /* TODO(pragma37):
      workbench_shadow_cache_populate(vedata, ob, has_transp_mat);
      */
    }
  }

  PassMain::Sub &pipeline_get(ObjectRef &ob_ref,
                              Material &material,
                              ::Image *image = nullptr,
                              eGPUSamplerState sampler_state = GPU_SAMPLER_DEFAULT,
                              ImageUser *iuser = nullptr)
  {
    const bool in_front = (ob_ref.object->dtx & OB_DRAW_IN_FRONT) != 0;
    if (xray_mode || material.is_transparent()) {
      if (in_front) {
        return transparent_in_front_ps.accumulation_ps_.sub_pass_get(
            ob_ref, image, sampler_state, iuser);
      }
      else {
        return transparent_ps.accumulation_ps_.sub_pass_get(ob_ref, image, sampler_state, iuser);
      }
    }
    else {
      if (in_front) {
        return opaque_in_front_ps.gbuffer_ps_.sub_pass_get(ob_ref, image, sampler_state, iuser);
      }
      else {
        return opaque_ps.gbuffer_ps_.sub_pass_get(ob_ref, image, sampler_state, iuser);
      }
    }
  }

  void draw(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    resources.color_tx.acquire(int2(resources.depth_tx.size()), GPU_RGBA16F);

    if (!clip_planes.is_empty()) {
      view.set_clip_planes(clip_planes);
    }

    PassSimple clear_ps = PassSimple("Workbench.Clear");
    clear_ps.init();
    clear_ps.clear_color(resources.world_buf.background_color);
    Framebuffer fb = Framebuffer("Workbench.Clear");
    fb.ensure(GPU_ATTACHMENT_NONE, GPU_ATTACHMENT_TEXTURE(resources.color_tx));
    fb.bind();
    manager.submit(clear_ps);

    if (!opaque_ps.is_empty() || !transparent_ps.is_empty()) {
      fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_tx));
      fb.bind();
      fb.clear_depth(1.0f);
    }
    if (!opaque_in_front_ps.is_empty() || !transparent_in_front_ps.is_empty()) {
      fb.ensure(GPU_ATTACHMENT_TEXTURE(resources.depth_in_front_tx));
      fb.bind();
      fb.clear_depth(1.0f);
    }

    opaque_ps.draw(manager, view);
    transparent_ps.draw(manager, view);
    opaque_in_front_ps.draw(manager, view);
    transparent_in_front_ps.draw(manager, view);

    // volume_ps.draw_prepass(manager, view, resources.depth_tx);

    anti_aliasing_ps.draw(manager, view, depth_tx, color_tx);

    resources.color_tx.release();
  }

  void draw_viewport(Manager &manager, View &view, GPUTexture *depth_tx, GPUTexture *color_tx)
  {
    this->draw(manager, view, depth_tx, color_tx);
  }
};

// Utils

bool get_matcap_tx(Texture &matcap_tx, const StudioLight &studio_light)
{
  ImBuf *matcap_diffuse = studio_light.matcap_diffuse.ibuf;
  ImBuf *matcap_specular = studio_light.matcap_specular.ibuf;
  if (matcap_diffuse && matcap_diffuse->rect_float) {
    int layers = 1;
    float *buffer = matcap_diffuse->rect_float;
    Vector<float> combined_buffer = {};

    if (matcap_specular && matcap_specular->rect_float) {
      int size = matcap_diffuse->x * matcap_diffuse->y * 4;
      combined_buffer.extend(matcap_diffuse->rect_float, size);
      combined_buffer.extend(matcap_specular->rect_float, size);
      buffer = combined_buffer.begin();
      layers++;
    }

    matcap_tx = Texture(studio_light.name,
                        GPU_RGBA16F,
                        int2(matcap_diffuse->x, matcap_diffuse->y),
                        layers,
                        buffer);
    return true;
  }
  return false;
}

float4x4 get_world_shading_rotation_matrix(float studiolight_rot_z)
{
  /* TODO(pragma37) C++ API ? */
  float V[4][4], R[4][4];
  DRW_view_viewmat_get(nullptr, V, false);
  axis_angle_to_mat4_single(R, 'Z', -studiolight_rot_z);
  mul_m4_m4m4(R, V, R);
  swap_v3_v3(R[2], R[1]);
  negate_v3(R[2]);
  return float4x4(R);
}

LightData get_light_data_from_studio_solidlight(const SolidLight *sl,
                                                float4x4 world_shading_rotation)
{
  LightData light = {};
  if (sl && sl->flag) {
    float3 direction = world_shading_rotation.ref_3x3() * float3(sl->vec);
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
  return light;
}

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
