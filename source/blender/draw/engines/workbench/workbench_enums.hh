/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BLI_assert.h"
#include "DNA_object_types.h"
#include "DNA_view3d_enums.h"

namespace blender::workbench {

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
  FLAT = 0,
  STUDIO,
  MATCAP,
};
static constexpr int shading_type_len = static_cast<int>(eShadingType::MATCAP) + 1;

static inline eShadingType shading_type_from_v3d_lighting(char lighting)
{
  switch (lighting) {
    case V3D_LIGHTING_FLAT:
      return eShadingType::FLAT;
    case V3D_LIGHTING_MATCAP:
      return eShadingType::MATCAP;
    case V3D_LIGHTING_STUDIO:
      return eShadingType::STUDIO;
    default:
      BLI_assert_unreachable();
      return static_cast<eShadingType>(-1);
  }
}

enum class eColorType {
  MATERIAL = 0,
  TEXTURE,
};
static constexpr int color_type_len = static_cast<int>(eColorType::TEXTURE) + 1;

static inline eColorType color_type_from_v3d_shading(char shading)
{
  return shading == V3D_SHADING_TEXTURE_COLOR ? eColorType::TEXTURE : eColorType::MATERIAL;
}

static inline const char *get_name(eColorType type)
{
  switch (type) {
    case eColorType::MATERIAL:
      return "Material";
    case eColorType::TEXTURE:
      return "Texture";
    default:
      BLI_assert_unreachable();
      return "";
  }
}

enum class eMaterialSubType {
  NONE = 0,
  MATERIAL,
  RANDOM,
  SINGLE,
  OBJECT,
  ATTRIBUTE,
};
static constexpr int material_subtype_len = static_cast<int>(eMaterialSubType::ATTRIBUTE) + 1;

static inline eMaterialSubType material_subtype_from_v3d_shading(char shading)
{
  switch (shading) {
    case V3D_SHADING_MATERIAL_COLOR:
      return eMaterialSubType::MATERIAL;
    case V3D_SHADING_RANDOM_COLOR:
      return eMaterialSubType::RANDOM;
    case V3D_SHADING_SINGLE_COLOR:
      return eMaterialSubType::SINGLE;
    case V3D_SHADING_TEXTURE_COLOR:
      return eMaterialSubType::NONE;
    case V3D_SHADING_OBJECT_COLOR:
      return eMaterialSubType::OBJECT;
    case V3D_SHADING_VERTEX_COLOR:
      return eMaterialSubType::ATTRIBUTE;
    default:
      BLI_assert_unreachable();
      return static_cast<eMaterialSubType>(-1);
  }
}

}  // namespace blender::workbench
