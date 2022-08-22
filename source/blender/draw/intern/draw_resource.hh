/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 * Component / Object level resources like object attributes, matrices, visibility etc...
 * Each of them are reference by resource index (#ResourceHandle).
 */

#include "BKE_curve.h"
#include "BKE_duplilist.h"
#include "BKE_mesh.h"
#include "BKE_object.h"
#include "BKE_volume.h"
#include "BLI_hash.h"
#include "DNA_curve_types.h"
#include "DNA_meta_types.h"
#include "DNA_object_types.h"

#include "draw_handle.hh"
#include "draw_manager.hh"
#include "draw_shader_shared.h"

/* -------------------------------------------------------------------- */
/** \name ObjectMatrices
 * \{ */

inline void ObjectMatrices::sync(const Object &object)
{
  model = object.obmat;
  model_inverse = object.imat;
}

inline void ObjectMatrices::sync(const float4x4 &model_matrix)
{
  model = model_matrix;
  model_inverse = model_matrix.inverted();
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name ObjectInfos
 * \{ */

ENUM_OPERATORS(eObjectInfoFlag, OBJECT_NEGATIVE_SCALE)

inline void ObjectInfos::sync(const blender::draw::ObjectRef ref, bool is_active_object)
{
  color = ref.object->color;
  index = ref.object->index;
  SET_FLAG_FROM_TEST(flag, is_active_object, eObjectInfoFlag::OBJECT_ACTIVE);
  SET_FLAG_FROM_TEST(
      flag, ref.object->base_flag & BASE_SELECTED, eObjectInfoFlag::OBJECT_SELECTED);
  SET_FLAG_FROM_TEST(
      flag, ref.object->base_flag & BASE_FROM_DUPLI, eObjectInfoFlag::OBJECT_FROM_DUPLI);
  SET_FLAG_FROM_TEST(
      flag, ref.object->base_flag & BASE_FROM_SET, eObjectInfoFlag::OBJECT_FROM_SET);
  SET_FLAG_FROM_TEST(
      flag, ref.object->transflag & OB_NEG_SCALE, eObjectInfoFlag::OBJECT_NEGATIVE_SCALE);

  if (ref.dupli_object == nullptr) {
    /* TODO(fclem): this is rather costly to do at draw time. Maybe we can
     * put it in ob->runtime and make depsgraph ensure it is up to date. */
    random = BLI_hash_int_2d(BLI_hash_string(ref.object->id.name + 2), 0) * (1.0f / 0xFFFFFFFF);
  }
  else {
    random = ref.dupli_object->random_id * (1.0f / 0xFFFFFFFF);
  }
  /* Default values. Set if needed. */
  random = 0.0f;

  if (ref.object->data == nullptr) {
    orco_add = float3(0.0f);
    orco_mul = float3(1.0f);
    return;
  }

  switch (GS(reinterpret_cast<ID *>(ref.object->data)->name)) {
    case ID_VO: {
      BoundBox &bbox = *BKE_volume_boundbox_get(ref.object);
      orco_add = (float3(bbox.vec[6]) + float3(bbox.vec[0])) * 0.5f; /* Center. */
      orco_mul = float3(bbox.vec[6]) - float3(bbox.vec[0]);          /* Size. */
      break;
    }
    case ID_ME: {
      BKE_mesh_texspace_get((Mesh *)ref.object->data, orco_add, orco_mul);
      break;
    }
    case ID_CU_LEGACY: {
      Curve &cu = *(Curve *)ref.object->data;
      BKE_curve_texspace_ensure(&cu);
      orco_add = cu.loc;
      orco_mul = cu.size;
      break;
    }
    case ID_MB: {
      MetaBall &mb = *(MetaBall *)ref.object->data;
      orco_add = mb.loc;
      orco_mul = mb.size;
      break;
    }
    default:
      orco_add = float3(0.0f);
      orco_mul = float3(1.0f);
      break;
  }
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name ObjectBounds
 * \{ */

inline void ObjectBounds::sync(Object &ob)
{
  const BoundBox *bbox = BKE_object_boundbox_get(&ob);
  if (bbox == nullptr) {
    bounding_sphere.w = -1.0f; /* Disable test. */
    return;
  }
  bounding_corners[0] = bbox->vec[0];
  bounding_corners[1] = bbox->vec[4];
  bounding_corners[2] = bbox->vec[3];
  bounding_corners[3] = bbox->vec[1];
  bounding_sphere.w = 0.0f; /* Enable test. */
}

/** \} */
