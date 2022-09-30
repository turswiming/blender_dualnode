/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "DNA_space_types.h"

#include "BKE_asset.h"
#include "BKE_asset_catalog.hh"
#include "BKE_asset_library.hh"
#include "BKE_idprop.h"
#include "BKE_screen.h"

#include "RNA_access.h"
#include "RNA_prototypes.h"

#include "ED_asset.h"

#include "node_intern.hh"

namespace blender::ed::space_node {

static bool node_add_menu_poll(const bContext *C, MenuType * /*mt*/)
{
  return CTX_wm_space_node(C);
}

static void gather_items_for_asset_library(const bContext &C,
                                           const bNodeTree &node_tree,
                                           const AssetLibraryReference &library_ref,
                                           bke::AssetLibrary &library,
                                           const bke::AssetCatalogPath &path,
                                           Vector<AssetHandle> &assets)
{
  const bke::AssetCatalogFilter filter = library.catalog_service->create_catalog_filter(path);
  AssetFilterSettings type_filter{};
  type_filter.id_types = FILTER_ID_NT;

  /* Iterate over all the assets in the current catalog. */
  ED_assetlist_storage_fetch(&library_ref, &C);
  ED_assetlist_ensure_previews_job(&library_ref, &C);
  ED_assetlist_iterate(library_ref, [&](AssetHandle asset) {
    if (!ED_asset_filter_matches_asset(&type_filter, &asset)) {
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

static void gather_child_catalogs(bke::AssetLibrary &library,
                                  const bke::AssetCatalogPath &catalog_path,
                                  Vector<bke::AssetCatalogPath> &child_catalogs)
{
  bke::AssetCatalogTree *tree = library.catalog_service->get_catalog_tree();
  if (!tree) {
    return;
  }
  tree->foreach_item([&](bke::AssetCatalogTreeItem &catalog_item) {
    if (catalog_item.catalog_path().parent() == catalog_path) {
      child_catalogs.append(catalog_item.catalog_path());
    }
  });
}

static void gather_root_catalogs(const bke::AssetLibrary &library,
                                 Vector<bke::AssetCatalogPath> &root_catalogs)
{
  bke::AssetCatalogTree *tree = library.catalog_service->get_catalog_tree();
  if (!tree) {
    return;
  }
  tree->foreach_root_item([&](bke::AssetCatalogTreeItem &catalog_item) {
    root_catalogs.append(catalog_item.catalog_path());
  });
}

static void node_add_catalog_assets_draw(const bContext *C, Menu *menu)
{
  const Main &bmain = *CTX_data_main(C);
  const SpaceNode &snode = *CTX_wm_space_node(C);
  const bNodeTree *edit_tree = snode.edittree;
  if (!edit_tree) {
    return;
  }

  const bke::AssetCatalogPath menu_path("Test Catalog");

  Vector<AssetHandle> assets;
  Vector<bke::AssetCatalogPath> child_catalogs;
  for (const AssetLibraryReference &ref : bke::all_asset_library_refs()) {
    if (bke::AssetLibrary *library = BKE_asset_library_load(&bmain, ref)) {
      gather_items_for_asset_library(*C, *edit_tree, ref, *library, menu_path, assets);
      gather_child_catalogs(*library, menu_path, child_catalogs);
    }
  }
  if (assets.is_empty() && child_catalogs.is_empty()) {
    return;
  }
  uiLayout *layout = menu->layout;
  uiItemS(layout);
  for (AssetHandle asset : assets) {
    PointerRNA props;
    uiItemFullO(layout,
                "NODE_OT_add_group_asset",
                ED_asset_handle_get_name(&asset),
                ICON_NONE,
                nullptr,
                WM_OP_INVOKE_DEFAULT,
                0,
                &props);
    /* TODO: Really use null for PointerRNA.owner_id? */
    /* TODO: Really use null for PointerRNA.owner_id? */
    RNA_pointer_set(&props,
                    "asset_handle",
                    PointerRNA{
                        nullptr,
                        &RNA_AssetHandle,
                        &asset  // AGH! NO! THIS DOES NOT WORK!
                    });
  }
  for (const bke::AssetCatalogPath &child_catalog : child_catalogs) {
    uiItemM(layout, "NODE_MT_node_add_catalog_assets", child_catalog.name().c_str(), ICON_NONE);
  }
}

static void add_root_catalogs_draw(const bContext *C, Menu *menu)
{
  const Main &bmain = *CTX_data_main(C);

  Vector<bke::AssetCatalogPath> catalogs;
  for (const AssetLibraryReference &ref : bke::all_asset_library_refs()) {
    if (bke::AssetLibrary *library = BKE_asset_library_load(&bmain, ref)) {
      gather_root_catalogs(*library, catalogs);
    }
  }

  uiLayout *layout = menu->layout;
  for (const bke::AssetCatalogPath &path : catalogs) {
    uiItemM(layout, "NODE_MT_node_add_catalog_assets", path.name().c_str(), ICON_NONE);
  }
}

MenuType add_catalog_assets_menu_type()
{
  MenuType type{};
  BLI_strncpy(type.idname, "NODE_MT_node_add_catalog_assets", sizeof(type.idname));
  type.poll = node_add_menu_poll;
  type.draw = node_add_catalog_assets_draw;
  return type;
}

MenuType add_root_catalogs_menu_type()
{
  MenuType type{};
  BLI_strncpy(type.idname, "NODE_MT_node_add_root_catalogs", sizeof(type.idname));
  type.poll = node_add_menu_poll;
  type.draw = add_root_catalogs_draw;
  return type;
}

}  // namespace blender::ed::space_node