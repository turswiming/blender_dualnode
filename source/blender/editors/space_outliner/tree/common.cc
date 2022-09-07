/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spoutliner
 *
 * Functions and helpers shared between tree-display types or other tree related code.
 */

#include "BLI_listbase.h"

#include "BKE_idtype.h"

#include "DNA_anim_types.h"
#include "DNA_object_types.h"
#include "DNA_outliner_types.h"

#include "RNA_access.h"
#include "RNA_prototypes.h"

#include "../outliner_intern.hh"
#include "common.hh"
#include "tree_display.hh"

namespace blender::ed::outliner {

/* -------------------------------------------------------------------- */
/** \name ID Helpers.
 * \{ */

const char *outliner_idcode_to_plural(short idcode)
{
  const char *propname = BKE_idtype_idcode_to_name_plural(idcode);
  PropertyRNA *prop = RNA_struct_type_find_property(&RNA_BlendData, propname);
  return (prop) ? RNA_property_ui_name(prop) : "UNKNOWN";
}

/** \} */

void outliner_make_object_parent_hierarchy(SubTree &subtree)
{
  /* build hierarchy */
  /* XXX also, set extents here... */
  for (TreeElement &te : subtree) {
    TreeStoreElem *tselem = TREESTORE(&te);

    if ((tselem->type != TSE_SOME_ID) || te.idcode != ID_OB) {
      continue;
    }

    Object *ob = (Object *)tselem->id;
    if (!ob->parent || !ob->parent->id.newid) {
      continue;
    }

    TreeElement *new_parent = reinterpret_cast<TreeElement *>(ob->parent->id.newid);

    std::unique_ptr floating_te = subtree.remove(te);
    new_parent->child_elements.add_back(std::move(floating_te));
  }
}

bool outliner_animdata_test(const AnimData *adt)
{
  if (adt) {
    return (adt->action || adt->drivers.first || adt->nla_tracks.first);
  }
  return false;
}

}  // namespace blender::ed::outliner
