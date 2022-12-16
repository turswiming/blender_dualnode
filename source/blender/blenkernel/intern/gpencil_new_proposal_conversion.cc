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

  int layer_index;
  LISTBASE_FOREACH_INDEX (const bGPDlayer *, old_gpl, &old_gpd->layers, layer_index) {
    new_gpd.add_layer(std::string(old_gpl->info));

    LISTBASE_FOREACH (const bGPDframe *, old_gpf, &old_gpl->frames) {
      int new_gpf_index{new_gpd.add_frame_on_layer(layer_index, old_gpf->framenum)};
      GPFrame &new_gpf{new_gpd.frames_for_write(new_gpf_index)};

      Vector<int> offsets({0});
      LISTBASE_FOREACH (const bGPDstroke *, old_gps, &old_gpf->strokes) {
        offsets.append(old_gps->totpoints + offsets.last());
      }

      CurvesGeometry &new_gps{new_gpf.strokes_as_curves()};
      MutableAttributeAccessor attributes = new_gps.attributes_for_write();
      new_gps.resize(offsets.last(), offsets.size() - 1);
      new_gps.offsets_for_write().copy_from(offsets);
      new_gps.curve_types_for_write().fill(CURVE_TYPE_POLY);

      MutableSpan<float3> new_gps_positions{new_gps.positions_for_write()};
      SpanAttributeWriter<float> radii = attributes.lookup_or_add_for_write_only_span<float>(
          ATTR_RADIUS, ATTR_DOMAIN_POINT);
      int stroke_index;
      LISTBASE_FOREACH_INDEX (const bGPDstroke *, old_gps, &old_gpf->strokes, stroke_index) {
        const IndexRange point_index_in_curve{new_gps.points_for_curve(stroke_index)};

        for (int point_index = 0; point_index < old_gps->totpoints; point_index++) {
          const bGPDspoint *old_pt{old_gps->points + point_index};
          new_gps_positions[point_index_in_curve[point_index]] = {old_pt->x, old_pt->y, old_pt->z};
          radii.span[point_index_in_curve[point_index]] = old_pt->pressure;
        }
      }

      radii.finish();
    }
  }

  return new_gpd;
}

bGPdata *convert_new_to_old_gpencil_data(const GPData &new_gpd)
{
  bGPdata *old_gpd = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));

  BLI_listbase_clear(&old_gpd->layers);
  old_gpd->totlayer = old_gpd->totframe = old_gpd->totstroke = 0;

  int frame_index{0};
  for (const int layer_index : IndexRange(new_gpd.layers_size)) {
    bGPDlayer *old_gpl = reinterpret_cast<bGPDlayer *>(MEM_mallocN(sizeof(bGPDlayer), __func__));
    const ::GPLayer *new_gpl{new_gpd.layers_array + layer_index};

    sprintf(old_gpl->info, "%s", new_gpl->name);

    BLI_listbase_clear(&old_gpl->mask_layers);
    BLI_listbase_clear(&old_gpl->frames);

    /* Add frames of correct layer index.
     * Assumes that frames in new data structure are sorted by layer index.
     */
    while ((frame_index < new_gpd.frames_size) &&
           (new_gpd.frames_array[frame_index].layer_index == layer_index)) {
      bGPDframe *old_gpf = reinterpret_cast<bGPDframe *>(MEM_mallocN(sizeof(bGPDframe), __func__));
      const GPFrame &new_gpf{new_gpd.frames(frame_index)};
      old_gpf->framenum = new_gpf.start_time;

      BLI_listbase_clear(&old_gpf->strokes);
      const CurvesGeometry &new_gps{new_gpf.strokes_as_curves()};
      AttributeAccessor attributes = new_gps.attributes();

      Span<float3> new_gps_positions = new_gps.positions();
      VArray<float> new_gps_radii = attributes.lookup_or_default<float>(ATTR_RADIUS, ATTR_DOMAIN_POINT, 0);

      for (int stroke_index = 0; stroke_index < new_gpf.strokes_num(); stroke_index++) {
        bGPDstroke *old_gps = reinterpret_cast<bGPDstroke *>(
            MEM_mallocN(sizeof(bGPDstroke), __func__));

        int point_num{new_gps.points_num_for_curve(stroke_index)};
        old_gps->points = reinterpret_cast<bGPDspoint *>(
            MEM_calloc_arrayN(point_num, sizeof(bGPDspoint), __func__));
        old_gps->totpoints = point_num;
        old_gps->triangles = nullptr;
        old_gps->editcurve = nullptr;
        old_gps->dvert = nullptr;

        int point_index{0};
        for (int new_gps_point_index : new_gps.points_for_curve(stroke_index)) {
          bGPDspoint *pt = &old_gps->points[point_index];
          copy_v3_v3(&pt->x, new_gps_positions[new_gps_point_index]);
          pt->pressure = new_gps_radii[new_gps_point_index];
          ++point_index;
        }

        BLI_addtail(&old_gpf->strokes, old_gps);
      }

      ++(old_gpd->totframe);
      BLI_addtail(&old_gpl->frames, old_gpf);
      ++frame_index;
    }

    ++(old_gpd->totlayer);
    BLI_addtail(&old_gpd->layers, old_gpl);
  }

  return old_gpd;
}

}  // namespace blender::bke
