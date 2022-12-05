/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "BKE_gpencil.h"
#include "DNA_gpencil_types.h"
#include "gpencil_new_proposal.hh"

namespace blender::bke {

GPData convert_old_to_new_gpencil_data(bGPdata *old_gpd)
{
  GPData new_gpd;

  /* Add all layers */
  Vector<std::string> layer_names;
  LISTBASE_FOREACH (bGPDlayer *, lay, &old_gpd->layers) {
    layer_names.append(std::string(lay->info));
  }
  new_gpd.add_layers(layer_names.as_span());

  /* Add all frames */
  int layer_index{-1};
  LISTBASE_FOREACH (bGPDlayer *, lay, &old_gpd->layers) {
    Vector<int> frame_indices;
    LISTBASE_FOREACH (bGPDframe *, frm, &lay->frames) {
      frame_indices.append(frm->framenum);
    }
    new_gpd.add_frames_on_layer(++layer_index, frame_indices.as_span());
  }

  return new_gpd;
}

bGPdata *convert_new_to_old_gpencil_data(const GPData &new_gpd)
{
  bGPdata *gpd = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));

  return gpd;
}

}  // namespace blender::bke