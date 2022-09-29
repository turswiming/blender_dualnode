/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "DNA_space_types.h"

#include "BKE_asset.h"
#include "BKE_asset_catalog.hh"
#include "BKE_asset_library.hh"
#include "BKE_idprop.h"
#include "BKE_screen.h"

#include "ED_asset.h"

#include "node_intern.hh"

namespace blender::ed::space_node {

bool node_add_catalog_assets_poll(const bContext *C, MenuType *mt)
{
  return CTX_wm_space_node(C);
}

static void gather_items_for_asset_library(const bContext &C,
                                           const bNodeTree &node_tree,
                                           const bke::CatalogID &catalog,
                                           const AssetLibraryReference &library_ref,
                                           Vector<AssetHandle> &assets,
                                           Vector<bke::AssetCatalogPath> &child_catalogs)
{
  AssetFilterSettings filter_settings{};
  filter_settings.id_types = FILTER_ID_NT;

  const bke::AssetLibrary &library = {};
  const bke::AssetCatalogService *service = library.catalog_service.get();
  if (!service) {
    return;
  }

  const bke::AssetCatalogFilter filter = service->create_catalog_filter(catalog);

  //   bke::AssetCatalogTree *catalog_tree = BKE_asset_library_get_catalog_tree(&library);

  /* Iterate over all the assets in the current catalog. */
  ED_assetlist_storage_fetch(&library_ref, &C);
  ED_assetlist_ensure_previews_job(&library_ref, &C);
  ED_assetlist_iterate(library_ref, [&](AssetHandle asset) {
    if (!ED_asset_filter_matches_asset(&filter_settings, &asset)) {
      return true;
    }
    const AssetMetaData *meta_data = ED_asset_handle_get_metadata(&asset);
    if (!filter.contains(meta_data->catalog_id)) {
      return true;
    }
    const AssetMetaData &asset_data = *ED_asset_handle_get_metadata(&asset);
    const IDProperty *tree_type = BKE_asset_metadata_idprop_find(&asset_data, "type");
    if (tree_type == nullptr || IDP_Int(tree_type) != node_tree.type) {
      return true;
    }
    assets.append(asset);
    return true;
  });
}

static void gather_items_for_all_assets(const bContext &C,
                                        const bNodeTree &node_tree,
                                        const bke::CatalogID &catalog,
                                        Vector<AssetHandle> &assets,
                                        Vector<bke::AssetCatalogPath> &child_catalogs)
{
  int i;
  LISTBASE_FOREACH_INDEX (const bUserAssetLibrary *, asset_library, &U.asset_libraries, i) {
    AssetLibraryReference library_ref{};
    library_ref.custom_library_index = i;
    library_ref.type = ASSET_LIBRARY_CUSTOM;
    /* Skip local assets to avoid duplicates when the asset is part of the local file library. */
    gather_items_for_asset_library(C, node_tree, catalog, library_ref, assets, child_catalogs);
  }

  AssetLibraryReference library_ref{};
  library_ref.custom_library_index = -1;
  library_ref.type = ASSET_LIBRARY_LOCAL;
  gather_items_for_asset_library(C, node_tree, catalog, library_ref, assets, child_catalogs);
}

void node_add_catalog_assets_draw(const bContext *C, Menu *menu)
{
  const SpaceNode &snode = *CTX_wm_space_node(C);
  const bNodeTree *edit_tree = snode.edittree;
  if (!edit_tree) {
    return;
  }

  const bke::AssetCatalogPath menu_catalog_path("Test Catalog");

  Vector<AssetHandle> assets_in_catalog;
  Vector<bke::AssetCatalogPath> child_catalogs;
  gather_items_for_all_assets(
      *C, *edit_tree, menu_catalog_path, assets_in_catalog, child_catalogs);
  if (assets_in_catalog.is_empty()) {
    return;
  }
  uiLayout *layout = menu->layout;
  uiItemS(layout);
  for (AssetHandle asset : assets_in_catalog) {
    uiItemO(layout, ED_asset_handle_get_name(&asset), ICON_NONE, "NODE_add_node");
  }
}

static MenuType NODE_MT_node_add_catalog_assets = []() {
  MenuType type{};
  BLI_strncpy(type.idname, "NODE_MT_node_add_catalog_assets", sizeof(type.idname));
  type.poll = node_add_catalog_assets_poll;
  type.draw = node_add_catalog_assets_draw;
  return type;
}();

}  // namespace blender::ed::space_node