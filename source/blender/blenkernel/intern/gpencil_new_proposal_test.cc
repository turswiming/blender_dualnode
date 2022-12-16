/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "testing/testing.h"

#include "BLI_index_mask_ops.hh"
#include "BLI_math_vec_types.hh"

#include "BKE_curves.hh"
#include "BKE_curves_utils.hh"
#include "BKE_gpencil.h"

#include "DNA_gpencil_types.h"

#include "gpencil_new_proposal.hh"

namespace blender::bke::gpencil::tests {

static GPData build_gpencil_data(int num_layers,
                                 int frames_per_layer,
                                 int strokes_per_frame,
                                 int points_per_stroke)
{
  GPData gpd;

  Vector<std::string> test_names;
  for (const int i : IndexRange(num_layers)) {
    test_names.append(std::string("GPLayer") + std::to_string(i));
  }
  gpd.add_layers(test_names.as_span());

  Array<int> test_start_frames(IndexRange(frames_per_layer).as_span());
  for (const int i : gpd.layers().index_range()) {
    gpd.add_frames_on_layer(i, test_start_frames);
  }

  for (const int frame_i : gpd.frames().index_range()) {
    CurvesGeometry &strokes = gpd.frames_for_write(frame_i).strokes_as_curves();
    strokes.resize(points_per_stroke * strokes_per_frame, strokes_per_frame);
    strokes.offsets_for_write().drop_back(1).fill(points_per_stroke);
    curves::accumulate_counts_to_offsets(strokes.offsets_for_write());
    MutableSpan<float3> positions = strokes.positions_for_write();
    for (const int curve_i : strokes.curves_range()) {
      const IndexRange points = strokes.points_for_curve(curve_i);
      for (const int i : IndexRange(points.size())) {
        positions[points[i]] = {float(i), float(curve_i * i % points.size()), float(i + curve_i)};
      }
    }
  }

  return gpd;
}

static bGPdata *build_old_gpencil_data(int num_layers,
                                       int frames_per_layer,
                                       int strokes_per_frame,
                                       int points_per_stroke)
{
  bGPdata *gpd = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));
  BLI_listbase_clear(&gpd->layers);
  gpd->totlayer = gpd->totframe = gpd->totstroke = 0;
  for (int i = 0; i < num_layers; i++) {
    bGPDlayer *gpl = reinterpret_cast<bGPDlayer *>(MEM_mallocN(sizeof(bGPDlayer), __func__));
    sprintf(gpl->info, "%s%d", "GPLayer", i);
    gpl->flag = 0;

    BLI_listbase_clear(&gpl->mask_layers);
    BLI_listbase_clear(&gpl->frames);
    for (int j = 0; j < frames_per_layer; j++) {
      bGPDframe *gpf = reinterpret_cast<bGPDframe *>(MEM_mallocN(sizeof(bGPDframe), __func__));
      gpf->framenum = j;

      BLI_listbase_clear(&gpf->strokes);
      for (int k = 0; k < strokes_per_frame; k++) {
        bGPDstroke *gps = reinterpret_cast<bGPDstroke *>(
            MEM_mallocN(sizeof(bGPDstroke), __func__));
        gps->totpoints = points_per_stroke;
        gps->points = reinterpret_cast<bGPDspoint *>(
            MEM_calloc_arrayN(points_per_stroke, sizeof(bGPDspoint), __func__));
        gps->triangles = nullptr;
        gps->editcurve = nullptr;
        gps->dvert = nullptr;

        for (int l = 0; l < points_per_stroke; l++) {
          float pos[3] = {(float)l, (float)((l * k) % points_per_stroke), (float)(l + k)};
          bGPDspoint *pt = &gps->points[l];
          copy_v3_v3(&pt->x, pos);
        }

        BLI_addtail(&gpf->strokes, gps);
      }
      gpd->totstroke += strokes_per_frame;
      BLI_addtail(&gpl->frames, gpf);
    }
    gpd->totframe += frames_per_layer;
    BLI_addtail(&gpd->layers, gpl);
  }
  gpd->totlayer += num_layers;

  return gpd;
}

static bGPdata *copy_old_gpencil_data(bGPdata *gpd_src)
{
  bGPdata *gpd_dst = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));
  BLI_listbase_clear(&gpd_dst->layers);
  LISTBASE_FOREACH (bGPDlayer *, gpl_src, &gpd_src->layers) {
    bGPDlayer *gpl_dst = BKE_gpencil_layer_duplicate(gpl_src, true, true);
    BLI_addtail(&gpd_dst->layers, gpl_dst);
  }

  return gpd_dst;
}

static void insert_new_frame_old_gpencil_data(bGPdata *gpd, int frame_num)
{
  bGPDlayer *gpl_active = BKE_gpencil_layer_active_get(gpd);
  BKE_gpencil_frame_addnew(gpl_active, frame_num);
}

static void free_old_gpencil_data(bGPdata *gpd)
{
  BKE_gpencil_free_layers(&gpd->layers);
  MEM_SAFE_FREE(gpd);
}

static void insert_new_stroke_old_gpencil_data(bGPDframe *gpf,
                                               int point_num,
                                               float *position,
                                               float *pressure)
{

  bGPDstroke *gps = reinterpret_cast<bGPDstroke *>(MEM_mallocN(sizeof(bGPDstroke), __func__));
  gps->totpoints = point_num;
  gps->points = reinterpret_cast<bGPDspoint *>(
      MEM_calloc_arrayN(point_num, sizeof(bGPDspoint), __func__));
  gps->triangles = nullptr;
  gps->editcurve = nullptr;
  gps->dvert = nullptr;

  for (int pt_i = 0; pt_i < point_num; pt_i++) {
    bGPDspoint *pt = &gps->points[pt_i];
    copy_v3_v3(&pt->x, position + 3 * pt_i);
    pt->pressure = pressure[pt_i];
  }

  BLI_addtail(&gpf->strokes, gps);
}

static void insert_new_stroke_new_gpencil_data(GPFrame &gpf,
                                               int point_num,
                                               float *position,
                                               float *pressure)
{
  int stroke_index{gpf.strokes_num()};
  gpf.add_new_stroke(point_num);
  CurvesGeometry &curves = gpf.strokes_as_curves();

  int point_index{0};
  MutableSpan<float3> pos = curves.positions_for_write();
  for (int point_index_in_curve : curves.points_for_curve(stroke_index)) {
    pos[point_index_in_curve].x = position[3 * point_index];
    pos[point_index_in_curve].y = position[3 * point_index + 1];
    pos[point_index_in_curve].z = position[3 * point_index + 2];
    ++point_index;
  }
}

static void compare_gpencil_stroke_data(const CurvesGeometry &curves,
                                        int curve_index,
                                        const bGPDstroke *stk)
{
  // Stroke/Curve length
  int curve_point_num{curves.points_num_for_curve(curve_index)};
  EXPECT_EQ(curve_point_num, stk->totpoints);
  if (curve_point_num != stk->totpoints) {
    return;
  }

  // Get curve attributes
  Span<float3> curve_positions{curves.positions()};

  IndexRange curve_id{curves.points_for_curve(curve_index)};
  int stk_point_id{0};
  for (int curve_point_id : curve_id) {
    const bGPDspoint *stk_point{stk->points + stk_point_id};
    const float3 &curve_pos{curve_positions[curve_point_id]};

    // Point positions
    EXPECT_EQ(curve_pos.x, stk_point->x);
    EXPECT_EQ(curve_pos.y, stk_point->y);
    EXPECT_EQ(curve_pos.z, stk_point->z);

    ++stk_point_id;
  }
}

static void compare_gpencil_data_structures(const GPData &new_gpd, const bGPdata *old_gpd)
{
  EXPECT_EQ(new_gpd.layers_size, old_gpd->totlayer);
  EXPECT_EQ(new_gpd.frames_size, old_gpd->totframe);

  int layer_id{0};
  int frame_id{0};
  LISTBASE_FOREACH (bGPDlayer *, old_gpl, &old_gpd->layers) {
    if (layer_id >= new_gpd.layers_size) {
      break;
    }

    // Compare layer data
    const ::GPLayer *new_gpl = &(new_gpd.layers_array[layer_id]);
    EXPECT_EQ(std::strcmp(new_gpl->name, old_gpl->info), 0);

    LISTBASE_FOREACH (bGPDframe *, old_gpf, &old_gpl->frames) {
      if (frame_id >= new_gpd.frames_size) {
        break;
      }
      const GPFrame &new_gpf{new_gpd.frames(frame_id)};

      // Compare frame data
      EXPECT_EQ(new_gpf.layer_index, layer_id);
      EXPECT_EQ(new_gpf.start_time, old_gpf->framenum);
      EXPECT_EQ(new_gpf.strokes_num(), BLI_listbase_count(&old_gpf->strokes));

      int curve_index{0};
      LISTBASE_FOREACH (const bGPDstroke *, stk, &old_gpf->strokes) {
        if (curve_index >= new_gpf.strokes_num()) {
          break;
        }

        // Compare stroke data
        compare_gpencil_stroke_data(new_gpf.strokes_as_curves(), curve_index, stk);

        ++curve_index;
      }

      ++frame_id;
    }
    ++layer_id;
  }
}

TEST(gpencil_proposal, EmptyGPData)
{
  GPData data;
  EXPECT_EQ(data.layers_size, 0);
  EXPECT_EQ(data.frames_size, 0);
}

TEST(gpencil_proposal, OneLayer)
{
  GPData data(1, 0);
  EXPECT_EQ(data.layers_size, 1);
  EXPECT_EQ(data.frames_size, 0);
}

TEST(gpencil_proposal, LayerName)
{
  GPLayer layer1;
  EXPECT_STREQ(layer1.name, "GP_Layer");

  GPLayer layer2("FooLayer");
  EXPECT_STREQ(layer2.name, "FooLayer");
}

TEST(gpencil_proposal, AddOneLayer)
{
  GPData data;

  const int layer_index = data.add_layer("FooLayer");
  EXPECT_EQ(data.layers_size, 1);
  EXPECT_STREQ(data.layers(layer_index).name, "FooLayer");
}

TEST(gpencil_proposal, AddLayers)
{
  GPData data;
  StringRefNull layer_names[3] = {"TestLayer1", "TestLayer2", "TestLayer3"};

  for (int i : IndexRange(3)) {
    data.add_layer(layer_names[i]);
  }
  EXPECT_EQ(data.layers_size, 3);

  for (int i : IndexRange(3)) {
    EXPECT_STREQ(data.layers(i).name, layer_names[i].c_str());
  }
}

TEST(gpencil_proposal, ChangeLayerName)
{
  GPData data;

  const int layer_index = data.add_layer("FooLayer");
  EXPECT_EQ(data.layers_size, 1);
  EXPECT_STREQ(data.layers(layer_index).name, "FooLayer");

  strcpy(data.layers_for_write(layer_index).name, "BarLayer");

  EXPECT_EQ(data.layers_size, 1);
  EXPECT_STREQ(data.layers(layer_index).name, "BarLayer");
}

TEST(gpencil_proposal, AddFrameToLayer)
{
  GPData data;

  data.add_layer("TestLayer1");
  const int layer2_index = data.add_layer("TestLayer2");

  const int frame_index = data.add_frame_on_layer(layer2_index, 0);
  EXPECT_NE(frame_index, -1);

  EXPECT_EQ(data.frames_size, 1);
  EXPECT_EQ(data.frames().last().layer_index, 1);
  EXPECT_EQ(data.frames(frame_index).layer_index, 1);

  data.frames_for_write(frame_index).start_time = 20;
  EXPECT_EQ(data.frames(frame_index).start_time, 20);
}

TEST(gpencil_proposal, CheckFramesSorted1)
{
  GPData data;

  const int frame_numbers1[5] = {10, 5, 6, 1, 3};
  const int frame_numbers_sorted1[5] = {1, 3, 5, 6, 10};

  int layer1_index = data.add_layer("TestLayer1");
  for (int i : IndexRange(5)) {
    const int frame_index = data.add_frame_on_layer(layer1_index, frame_numbers1[i]);
    EXPECT_NE(frame_index, -1);
    EXPECT_EQ(data.frames(frame_index).start_time, frame_numbers1[i]);
  }

  for (const int i : data.frames().index_range()) {
    EXPECT_EQ(data.frames(i).start_time, frame_numbers_sorted1[i]);
  }
}

TEST(gpencil_proposal, CheckFramesSorted2)
{
  GPData data;

  const int frame_numbers_layer1[5] = {10, 5, 6, 1, 3};
  const int frame_numbers_layer2[5] = {8, 5, 7, 1, 4};
  const int frame_numbers_sorted2[10][2] = {
      {0, 1}, {0, 3}, {0, 5}, {0, 6}, {0, 10}, {1, 1}, {1, 4}, {1, 5}, {1, 7}, {1, 8}};

  const int layer1_index = data.add_layer("TestLayer1");
  const int layer2_index = data.add_layer("TestLayer2");
  for (int i : IndexRange(5)) {
    data.add_frame_on_layer(layer1_index, frame_numbers_layer1[i]);
    data.add_frame_on_layer(layer2_index, frame_numbers_layer2[i]);
  }

  for (const int i : data.frames().index_range()) {
    EXPECT_EQ(data.frames(i).layer_index, frame_numbers_sorted2[i][0]);
    EXPECT_EQ(data.frames(i).start_time, frame_numbers_sorted2[i][1]);
  }
}

TEST(gpencil_proposal, IterateOverFramesOnLayer)
{
  GPData data;

  const int frame_numbers_layer1[5] = {10, 5, 6, 1, 3};
  const int frame_numbers_layer2[5] = {8, 5, 7, 1, 4};

  const int frame_numbers_sorted[10] = {1, 3, 5, 6, 10, 1, 4, 5, 7, 8};

  const int layer1_index = data.add_layer("TestLayer1");
  const int layer2_index = data.add_layer("TestLayer2");
  for (int i : IndexRange(5)) {
    data.add_frame_on_layer(layer1_index, frame_numbers_layer1[i]);
    data.add_frame_on_layer(layer2_index, frame_numbers_layer2[i]);
  }

  data.frames_on_layer(layer1_index);
  EXPECT_TRUE(data.runtime->frames_index_range_cache.contains(layer1_index));
  for (const int i : data.frames_on_layer(layer1_index)) {
    EXPECT_EQ(data.frames(i).start_time, frame_numbers_sorted[i]);
  }

  data.frames_on_layer(layer2_index);
  EXPECT_TRUE(data.runtime->frames_index_range_cache.contains(layer2_index));
  for (const int i : data.frames_on_layer(layer2_index)) {
    EXPECT_EQ(data.frames(i).start_time, frame_numbers_sorted[i]);
  }
}

TEST(gpencil_proposal, AddSingleStroke)
{
  GPData data;
  const int layer1_index = data.add_layer("TestLayer1");

  const int frame_index = data.add_frame_on_layer(layer1_index, 0);
  EXPECT_NE(frame_index, -1);

  IndexRange point_indices = data.frames_for_write(frame_index).add_new_stroke(100);

  EXPECT_EQ(data.strokes_num(), 1);
  EXPECT_EQ(data.frames(frame_index).strokes_num(), 1);
  EXPECT_EQ(point_indices.size(), 100);
}

TEST(gpencil_proposal, ChangeStrokePoints)
{
  GPData data;
  const int layer1_index = data.add_layer("TestLayer1");

  static const Array<float3> test_positions{{
      {1.0f, 2.0f, 3.0f},
      {4.0f, 5.0f, 6.0f},
      {7.0f, 8.0f, 9.0f},
  }};

  const int frame_index = data.add_frame_on_layer(layer1_index, 0);
  EXPECT_NE(frame_index, -1);

  CurvesGeometry &curves = data.frames_for_write(frame_index).strokes_as_curves();
  curves.resize(curves.points_num() + test_positions.size(), curves.curves_num() + 1);
  curves.offsets_for_write().last() = curves.offsets().last(1) + test_positions.size();

  const IndexRange new_points = curves.points_for_curve(curves.curves_range().last());
  MutableSpan<float3> new_positions = curves.positions_for_write().slice(new_points);
  new_positions.copy_from(test_positions);

  for (const int i : curves.positions().index_range()) {
    EXPECT_V3_NEAR(new_positions[i], test_positions[i], 1e-5f);
  }
}

TEST(gpencil_proposal, BigGPData)
{
  GPData data = build_gpencil_data(5, 500, 100, 100);

  EXPECT_EQ(data.strokes_num(), 250e3);
  EXPECT_EQ(data.points_num(), 25e6);
}

TEST(gpencil_proposal, TimeBigGPDataCopy)
{
  int layers_num = 10, frames_num = 500, strokes_num = 100, points_num = 100;

  GPData data = build_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  GPData data_copy;

  {
    SCOPED_TIMER("BigGPDataCopy");
    data_copy = data;
  }

  bGPdata *old_data = build_old_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  bGPdata *old_data_copy;

  {
    SCOPED_TIMER("BigGPDataCopyOld");
    old_data_copy = copy_old_gpencil_data(old_data);
  }

  free_old_gpencil_data(old_data);
  free_old_gpencil_data(old_data_copy);
}

TEST(gpencil_proposal, TimeInsertFrame)
{
  int layers_num = 100, frames_num = 1000, strokes_num = 10, points_num = 10;
  GPData data = build_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  data.set_active_layer(7);

  {
    SCOPED_TIMER("TimeInsertFrame");
    data.add_frame_on_active_layer(347);
  }

  EXPECT_EQ(data.frames_on_active_layer().size(), 1001);

  bGPdata *old_data = build_old_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  int i = 0;
  bGPDlayer *gpl_active = NULL;
  LISTBASE_FOREACH_INDEX (bGPDlayer *, gpl, &old_data->layers, i) {
    if (i == 7) {
      BKE_gpencil_layer_active_set(old_data, gpl);
      gpl_active = gpl;
      break;
    }
  }
  /* Remove the frame so we can insert it again. */
  LISTBASE_FOREACH (bGPDframe *, gpf, &gpl_active->frames) {
    if (gpf->framenum == 347) {
      BKE_gpencil_layer_frame_delete(gpl_active, gpf);
      break;
    }
  }

  {
    SCOPED_TIMER("TimeOldInsertFrame");
    insert_new_frame_old_gpencil_data(old_data, 347);
  }

  EXPECT_EQ(BLI_listbase_count(&gpl_active->frames), 1000);

  free_old_gpencil_data(old_data);
}

TEST(gpencil_proposal, TimeMultiFrameTransformStrokes)
{
  int layers_num = 1, frames_num = 1000, strokes_num = 100, points_num = 100;
  GPData data = build_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  data.set_active_layer(0);

  float4x4 translate_mat = float4x4::from_location({1.0f, 2.0f, 3.0f});
  {
    SCOPED_TIMER("TimeMultiFrameTransformStrokes");
    for (const int i : data.frames_on_active_layer()) {
      GPFrame &gpf = data.frames_for_write(i);
      CurvesGeometry &curves = gpf.strokes_as_curves();
      curves.transform(translate_mat);
    }
  }

  bGPdata *old_data = build_old_gpencil_data(layers_num, frames_num, strokes_num, points_num);
  BKE_gpencil_layer_active_set(old_data, reinterpret_cast<bGPDlayer *>(old_data->layers.first));

  float matrix[4][4], loc[3] = {1.0f, 2.0f, 3.0f};
  unit_m4(matrix);
  copy_v3_v3(matrix[3], loc);
  {
    SCOPED_TIMER("TimeOldMultiFrameTransformStrokes");
    bGPDlayer *gpl_active = BKE_gpencil_layer_active_get(old_data);
    LISTBASE_FOREACH (bGPDframe *, gpf, &gpl_active->frames) {
      LISTBASE_FOREACH (bGPDstroke *, gps, &gpf->strokes) {
        for (int i = 0; i < gps->totpoints; i++) {
          bGPDspoint *pt = &gps->points[i];
          mul_m4_v3(matrix, &pt->x);
        }
      }
    }
  }

  free_old_gpencil_data(old_data);
}

TEST(gpencil_proposal, OldToNewConversion)
{
  int layers_num = 2, frames_num = 3, strokes_num = 2, points_num = 2;

  bGPdata *old_data = build_old_gpencil_data(layers_num, frames_num, strokes_num, points_num);

  GPData data = convert_old_to_new_gpencil_data(old_data);

  compare_gpencil_data_structures(data, old_data);

  free_old_gpencil_data(old_data);
}

TEST(gpencil_proposal, NewToOldConversion)
{
  int layers_num = 2, frames_num = 3, strokes_num = 2, points_num = 2;

  GPData data = build_gpencil_data(layers_num, frames_num, strokes_num, points_num);

  bGPdata *old_data = convert_new_to_old_gpencil_data(data);

  compare_gpencil_data_structures(data, old_data);

  free_old_gpencil_data(old_data);
}

TEST(gpencil_proposal, OldToNewStrokeConversion)
{
  int layers_num = 1, frames_num = 1, strokes_num = 0, points_num = 0;

  bGPdata *old_data = build_old_gpencil_data(layers_num, frames_num, strokes_num, points_num);

  const int point_num{5};
  float pos[point_num * 3] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7};
  float prs[point_num] = {0.5, 0.75, 1, 1.5, 1.75};

  LISTBASE_FOREACH (bGPDlayer *, lay, &old_data->layers) {
    LISTBASE_FOREACH (bGPDframe *, frm, &lay->frames) {
      insert_new_stroke_old_gpencil_data(frm, point_num, pos, prs);
    }
  }

  GPData data = convert_old_to_new_gpencil_data(old_data);

  compare_gpencil_data_structures(data, old_data);

  free_old_gpencil_data(old_data);
}

TEST(gpencil_proposal, NewToOldStrokeConversion)
{
  int layers_num = 1, frames_num = 1, strokes_num = 0, points_num = 0;

  GPData data = build_gpencil_data(layers_num, frames_num, strokes_num, points_num);

  const int point_num{5};
  float pos[point_num * 3] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7};
  float prs[point_num] = {0.5, 0.75, 1, 1.5, 1.75};
  insert_new_stroke_new_gpencil_data(data.frames_for_write(0), point_num, pos, prs);

  bGPdata *old_data = convert_new_to_old_gpencil_data(data);

  compare_gpencil_data_structures(data, old_data);

  free_old_gpencil_data(old_data);
}
}  // namespace blender::bke::gpencil::tests
