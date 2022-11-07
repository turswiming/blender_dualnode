/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "DNA_ID.h"
#include "DNA_asset_types.h"

#include "BKE_asset.h"
#include "BKE_asset_representation.hh"

namespace blender::bke {

AssetRepresentation::AssetRepresentation(std::unique_ptr<AssetMetaData> metadata)
    : metadata_(std::move(metadata))
{
}

AssetRepresentation::AssetRepresentation(const ID &id) : local_id_metadata_(id.asset_data)
{
}

AssetMetaData &AssetRepresentation::get_metadata() const
{
  return local_id_metadata_ ? *local_id_metadata_ : *metadata_;
}

bool AssetRepresentation::is_local_id() const
{
  return local_id_metadata_ != nullptr;
}

}  // namespace blender::bke

/* ---------------------------------------------------------------------- */
/** \name C-API
 * \{ */

using namespace blender;

AssetMetaData *BKE_asset_representation_metadata_get(const AssetRepresentation *asset_handle)
{
  const bke::AssetRepresentation *asset = reinterpret_cast<const bke::AssetRepresentation *>(
      asset_handle);
  return &asset->get_metadata();
}

bool BKE_asset_representation_is_local_id(const AssetRepresentation *asset_handle)
{
  const bke::AssetRepresentation *asset = reinterpret_cast<const bke::AssetRepresentation *>(
      asset_handle);
  return asset->is_local_id();
}

/** \} */
