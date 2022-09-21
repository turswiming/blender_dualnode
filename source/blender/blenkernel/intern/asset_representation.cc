/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "DNA_ID.h"
#include "DNA_asset_types.h"

#include "BKE_asset_representation.hh"

namespace blender::bke {

AssetRepresentation::AssetRepresentation(std::unique_ptr<AssetMetaData> metadata)
    : metadata_(std::move(metadata))
{
}

AssetRepresentation::AssetRepresentation(ID &id) : local_id_metadata_(id.asset_data)
{
}

AssetRepresentation::AssetRepresentation(AssetMetaData &&metadata)
    : metadata_(std::make_unique<AssetMetaData>(metadata))
{
}

AssetMetaData &AssetRepresentation::get_metadata()
{
  return local_id_metadata_ ? *local_id_metadata_ : *metadata_;
}

bool AssetRepresentation::is_local()
{
  return local_id_metadata_ != nullptr;
}

}  // namespace blender::bke
