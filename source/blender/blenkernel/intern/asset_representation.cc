/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include <stdexcept>

#include "DNA_ID.h"
#include "DNA_asset_types.h"

#include "BKE_asset.h"
#include "BKE_asset_representation.hh"

namespace blender::bke {

AssetRepresentation::AssetRepresentation(std::unique_ptr<AssetMetaData> metadata)
    : metadata_(std::move(metadata))
{
}

AssetRepresentation::AssetRepresentation(ID &id) : local_asset_id_(&id)
{
  if (!id.asset_data) {
    throw std::invalid_argument("Passed ID is not an asset");
  }
}

AssetMetaData &AssetRepresentation::get_metadata() const
{
  return local_asset_id_ ? *local_asset_id_->asset_data : *metadata_;
}

bool AssetRepresentation::is_local_id() const
{
  return local_asset_id_ != nullptr;
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
