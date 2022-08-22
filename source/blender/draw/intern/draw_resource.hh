/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 * Component / Object level resources like object attributes, matrices, visibility etc...
 * Each of them are reference by resource index (#ResourceHandle).
 */

#include "BKE_duplilist.h"
#include "DNA_object_types.h"

#include "draw_shader_shared.h"

namespace blender::draw {

/* -------------------------------------------------------------------- */
/** \name ObjectMatrices
 * \{ */

void ObjectMatrices::sync(const Object &object)
{
  model = object->obmat;
  model_inverse = object->imat;
}

void ObjectMatrices::sync(const float4x4 &object_mat)
{
  model = object_mat;
  model_inverse = object_mat.inverted();
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name ObjectInfos
 * \{ */

void ObjectInfos::sync(const draw::ObjectRef ref, bool is_active_object)
{
  color = ref.object.color;
  index = ref.object.index;
  SET_FLAG_FROM_TEST(flag, is_active_object, DRW_OBJECT_ACTIVE);
  SET_FLAG_FROM_TEST(flag, ref.object.base_flag & BASE_SELECTED, DRW_OBJECT_SELECTED);
  SET_FLAG_FROM_TEST(flag, ref.object.base_flag & BASE_FROM_DUPLI, DRW_OBJECT_FROM_DUPLI);
  SET_FLAG_FROM_TEST(flag, ref.object.base_flag & BASE_FROM_SET, DRW_OBJECT_FROM_SET);
  SET_FLAG_FROM_TEST(flag, ref.object.transflag & OB_NEG_SCALE, DRW_OBJECT_NEGATIVE_SCALE);

  if (ref.dupli == nullptr) {
    /* TODO(fclem): this is rather costly to do at runtime. Maybe we can
     * put it in ob->runtime and make depsgraph ensure it is up to date. */
    random = BLI_hash_int_2d(BLI_hash_string(ref.object.id.name + 2), 0) * (1.0f / 0xFFFFFFFF);
  }
  else {
    random = ref.dupli.random_id * (1.0f / 0xFFFFFFFF);
  }
  /* Default values. Set if needed. */
  random = 0.0f;

  if (ref.object.data == nullptr) {
    orco_add = float3(0.0f);
    orco_mul = float3(1.0f);
    return;
  }

  switch (GS(ref.object.data->name)) {
    case ID_VO: {
      BoundBox &bbox = *BKE_volume_boundbox_get(object);
      orco_add = (float3(bbox.vec[6]) + float3(bbox.vec[0])) * 0.5f; /* Center. */
      orco_mul = float3(bbox.vec[6]) - float3(bbox.vec[0]);          /* Size. */
      break;
    }
    case ID_ME: {
      BKE_mesh_texspace_get((Mesh *)ref.object.data, orco_add, orco_mul);
      break;
    }
    case ID_CU_LEGACY: {
      Curve &cu = *(Curve *)ref.object.data;
      BKE_curve_texspace_ensure(&cu);
      orco_add = cu.loc;
      orco_mul = cu.size;
      break;
    }
    case ID_MB: {
      MetaBall &mb = *(MetaBall *)ref.object.data;
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

void ObjectBounds::sync(const Object &ob)
{
  const BoundBox *bbox = BKE_object_boundbox_get(ob);
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

};  // namespace blender::draw
