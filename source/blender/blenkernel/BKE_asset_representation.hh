/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#pragma once

#include <memory>

struct AssetMetaData;
struct ID;

namespace blender::bke {

/**
 * \brief Abstraction to reference an asset, with necessary data for display & interaction.
 *
 * #AssetRepresentation is the data-structure to store information about a single asset. It doesn't
 * contain the asset itself, but information like the metadata and preview, as well as methods to
 * interact with them. Think of it like a view on an asset.
 */
class AssetRepresentation {
  friend class AssetLibrary;

  /** Null if the asset represents a local ID, in which case the ID owns the metadata. */
  std::unique_ptr<AssetMetaData> metadata_ = nullptr;
  /** If this asset represents an ID in the current file, this references the ID directly. This
   * means the represented asset is a fully local, and editable entity. */
  ID *local_asset_id_ = nullptr; /* Non-owning. */

 public:
  explicit AssetRepresentation(std::unique_ptr<AssetMetaData> metadata);
  /** Constructs an asset representation for an ID stored in the current file. This makes the asset
   * local and fully editable. */
  explicit AssetRepresentation(ID &id);

  AssetMetaData &get_metadata() const;
  /** Returns if this asset is stored inside this current file, and as such fully editable. */
  bool is_local_id() const;
};

}  // namespace blender::bke
