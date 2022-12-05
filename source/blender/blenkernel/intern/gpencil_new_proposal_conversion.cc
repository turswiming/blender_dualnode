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

  BLI_listbase_clear(&gpd->layers);
  gpd->totlayer = gpd->totframe = gpd->totstroke = 0;

  int frm_offset{0};
  gpd->totlayer = new_gpd.layers_size;
  for (int lay_i = 0; lay_i < new_gpd.layers_size; lay_i++) {
    bGPDlayer *gpl = reinterpret_cast<bGPDlayer *>(MEM_mallocN(sizeof(bGPDlayer), __func__));
    const ::GPLayer *lay{new_gpd.layers_array + lay_i};
    sprintf(gpl->info, "%s", lay->name);

    BLI_listbase_clear(&gpl->mask_layers);
    BLI_listbase_clear(&gpl->frames);
    BLI_addtail(&gpd->layers, gpl);

    /* Add frames of correct layer index.
       Assumes that frames in new data structure are sorted by layer index.
    */
    while ((frm_offset < new_gpd.frames_size) &&
           (new_gpd.frames_array[frm_offset].layer_index == lay_i)) {
      bGPDframe *gpf = reinterpret_cast<bGPDframe *>(MEM_mallocN(sizeof(bGPDframe), __func__));
      const ::GPFrame *frm{new_gpd.frames_array + frm_offset};
      gpf->framenum = frm->start_time;

      BLI_listbase_clear(&gpf->strokes);

      BLI_addtail(&gpl->frames, gpf);
      ++(gpd->totframe);
      ++frm_offset;
    }
  }

  return gpd;
}

}  // namespace blender::bke