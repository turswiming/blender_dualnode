/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spoutliner
 */

#include "BLI_listbase.h"
#include "BLI_listbase_wrapper.hh"

#include "BKE_collection.h"
#include "BKE_main.h"

#include "DNA_collection_types.h"
#include "DNA_space_types.h"

#include "BLT_translation.h"

#include "../outliner_intern.hh"
#include "common.hh"
#include "tree_display.hh"
#include "tree_element.hh"

namespace blender::ed::outliner {

/* Convenience/readability. */
template<typename T> using List = ListBaseWrapper<T>;

TreeDisplayOverrideLibraryProperties::TreeDisplayOverrideLibraryProperties(
    SpaceOutliner &space_outliner)
    : AbstractTreeDisplay(space_outliner)
{
}

SubTree TreeDisplayOverrideLibraryProperties::buildTree(const TreeSourceData &source_data)
{
  SubTree tree = add_library_contents(*source_data.bmain);

  for (TreeElement &top_level_te : tree) {
    TreeStoreElem *tselem = TREESTORE(&top_level_te);
    if (!tselem->used) {
      tselem->flag &= ~TSE_CLOSED;
    }
  }

  return tree;
}

SubTree TreeDisplayOverrideLibraryProperties::add_library_contents(Main &mainvar)
{
  SubTree tree;

  const short filter_id_type = id_filter_get();

  ListBase *lbarray[INDEX_ID_MAX];
  int tot;
  if (filter_id_type) {
    lbarray[0] = which_libbase(&mainvar, space_outliner_.filter_id_type);
    tot = 1;
  }
  else {
    tot = set_listbasepointers(&mainvar, lbarray);
  }

  for (int a = 0; a < tot; a++) {
    if (!lbarray[a] || !lbarray[a]->first) {
      continue;
    }

    ID *id = nullptr;

    /* check if there's data in current id list */
    for (ID *id_iter : List<ID>(lbarray[a])) {
      if (ID_IS_OVERRIDE_LIBRARY_REAL(id_iter) && !ID_IS_LINKED(id_iter)) {
        id = id_iter;
        break;
      }
    }

    if (id == nullptr) {
      continue;
    }

    /* Create data-block list parent element on demand. */
    TreeElement *id_base_te = nullptr;
    SubTree *subtree_to_expand = &tree;

    if (!filter_id_type) {
      id_base_te = outliner_add_element(&space_outliner_, lbarray[a], tree, TSE_ID_BASE, 0);
      id_base_te->directdata = lbarray[a];
      id_base_te->name = outliner_idcode_to_plural(GS(id->name));

      subtree_to_expand = &id_base_te->child_elements;
    }

    for (ID *id : List<ID>(lbarray[a])) {
      if (ID_IS_OVERRIDE_LIBRARY_REAL(id) && !ID_IS_LINKED(id)) {
        TreeElement *override_tree_element = outliner_add_element(
            &space_outliner_, id, id_base_te, TSE_LIBRARY_OVERRIDE_BASE, 0);

        if (override_tree_element->child_elements.is_empty()) {
          subtree_to_expand->remove(*override_tree_element);
        }
      }
    }
  }

  /* Remove ID base elements that turn out to be empty. */
  for (TreeElement &te : tree) {
    if (te.child_elements.is_empty()) {
      tree.remove(te);
    }
  }

  return tree;
}

short TreeDisplayOverrideLibraryProperties::id_filter_get() const
{
  if (space_outliner_.filter & SO_FILTER_ID_TYPE) {
    return space_outliner_.filter_id_type;
  }
  return 0;
}

}  // namespace blender::ed::outliner
