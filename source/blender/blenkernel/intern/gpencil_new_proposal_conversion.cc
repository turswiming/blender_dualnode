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

  int layer_index{0};
  LISTBASE_FOREACH (bGPDlayer *, old_gpl, &old_gpd->layers) {
    new_gpd.add_layer(std::string(old_gpl->info));

    LISTBASE_FOREACH (bGPDframe *, old_gpf, &old_gpl->frames) {
      int new_gpf_index{new_gpd.add_frame_on_layer(layer_index, old_gpf->framenum)};
      GPFrame &new_gpf{new_gpd.frames_for_write(new_gpf_index)};
    }

    ++layer_index;
  }

  return new_gpd;
}

bGPdata *convert_new_to_old_gpencil_data(const GPData &new_gpd)
{
  bGPdata *old_gpd = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));

  BLI_listbase_clear(&old_gpd->layers);
  old_gpd->totlayer = old_gpd->totframe = old_gpd->totstroke = 0;

  int frm_offset{0};
  for (int lay_i = 0; lay_i < new_gpd.layers_size; lay_i++) {
    bGPDlayer *old_gpl = reinterpret_cast<bGPDlayer *>(MEM_mallocN(sizeof(bGPDlayer), __func__));
    const ::GPLayer *new_gpl{new_gpd.layers_array + lay_i};

    sprintf(old_gpl->info, "%s", new_gpl->name);

    BLI_listbase_clear(&old_gpl->mask_layers);
    BLI_listbase_clear(&old_gpl->frames);

    /* Add frames of correct layer index.
       Assumes that frames in new data structure are sorted by layer index.
    */
    while ((frm_offset < new_gpd.frames_size) &&
           (new_gpd.frames_array[frm_offset].layer_index == lay_i)) {
      bGPDframe *gpf = reinterpret_cast<bGPDframe *>(MEM_mallocN(sizeof(bGPDframe), __func__));
      const ::GPFrame *new_gpf{new_gpd.frames_array + frm_offset};
      gpf->framenum = new_gpf->start_time;

      BLI_listbase_clear(&gpf->strokes);

      ++(old_gpd->totframe);
      BLI_addtail(&old_gpl->frames, gpf);
      ++frm_offset;
    }

    ++(old_gpd->totlayer);
    BLI_addtail(&old_gpd->layers, old_gpl);
  }

  return old_gpd;
}

}  // namespace blender::bke