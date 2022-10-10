/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "workbench_private.hh"

#include "BLI_hash.h"

namespace blender::workbench {

Material::Material() = default;

Material::Material(float3 color)
{
  base_color = color;
  packed_data = Material::pack_data(0.0f, 0.4f, 1.0f);
}

Material::Material(::Object &ob, bool random)
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

Material::Material(::Material &mat)
{
  base_color = &mat.r;
  packed_data = Material::pack_data(mat.metallic, mat.roughness, mat.a);
}

uint32_t Material::pack_data(float metallic, float roughness, float alpha)
{
  /* Remap to Disney roughness. */
  roughness = sqrtf(roughness);
  uint32_t packed_roughness = unit_float_to_uchar_clamp(roughness);
  uint32_t packed_metallic = unit_float_to_uchar_clamp(metallic);
  uint32_t packed_alpha = unit_float_to_uchar_clamp(alpha);
  return (packed_alpha << 16u) | (packed_roughness << 8u) | packed_metallic;
}

}  // namespace blender::workbench
