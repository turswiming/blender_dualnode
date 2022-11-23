/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include <algorithm>

#include "BKE_curves.hh"
#include "BKE_gpencil.h"

#include "BLI_index_mask_ops.hh"
#include "BLI_math_vec_types.hh"

#include "DNA_gpencil_types.h"

#include "gpencil_new_proposal.hh"

namespace blender::bke {

Span<float3> GPStroke::points_positions() const
{
  return {geometry_->positions().begin() + offset_, points_num_};
}

MutableSpan<float3> GPStroke::points_positions_for_write() const
{
  return {geometry_->positions_for_write().begin() + offset_, points_num_};
}

void GPStroke::transform(float4x4 matrix)
{
  threading::parallel_for(
      points_positions_for_write().index_range(), 512, [&](const IndexRange range) {
        for (float3 &position : points_positions_for_write().slice(range)) {
          position = matrix * position;
        }
      });
}

GPFrame::GPFrame(int start_frame, int end_frame)
{
  this->start_time = start_frame;
  this->end_time = end_frame;
  this->strokes = nullptr;
}

GPFrame::GPFrame(const GPFrame &other) : GPFrame(other.start_time, other.end_time)
{
  if (other.strokes != nullptr) {
    /* Make sure old strokes are freed before copying. */
    MEM_SAFE_FREE(this->strokes);
    this->strokes = MEM_new<CurvesGeometry>(__func__);

    *reinterpret_cast<CurvesGeometry *>(this->strokes) = CurvesGeometry::wrap(*other.strokes);
  }
  this->layer_index = other.layer_index;
}

GPFrame &GPFrame::operator=(const GPFrame &other)
{
  if (this != &other && other.strokes != nullptr) {
    /* Make sure old strokes are freed before copying. */
    MEM_SAFE_FREE(this->strokes);
    this->strokes = MEM_new<CurvesGeometry>(__func__);

    *reinterpret_cast<CurvesGeometry *>(this->strokes) = CurvesGeometry::wrap(*other.strokes);
  }
  this->layer_index = other.layer_index;
  this->start_time = other.start_time;
  this->end_time = other.end_time;
  return *this;
}

GPFrame::GPFrame(GPFrame &&other) : GPFrame(other.start_time, other.end_time)
{
  if (this != &other) {
    std::swap(this->strokes, other.strokes);
    other.strokes = nullptr;
  }
  this->layer_index = other.layer_index;
}

GPFrame &GPFrame::operator=(GPFrame &&other)
{
  if (this != &other) {
    std::swap(this->strokes, other.strokes);
    other.strokes = nullptr;
  }
  this->layer_index = other.layer_index;
  this->start_time = other.start_time;
  this->end_time = other.end_time;
  return *this;
}

GPFrame::~GPFrame()
{
  MEM_delete(reinterpret_cast<CurvesGeometry *>(this->strokes));
  this->strokes = nullptr;
}

bool GPFrame::operator<(const GPFrame &other) const
{
  if (this->start_time == other.start_time) {
    return this->layer_index < other.layer_index;
  }
  return this->start_time < other.start_time;
}

bool GPFrame::operator<(const std::pair<int, int> elem) const
{
  if (this->start_time == elem.second) {
    return this->layer_index < elem.first;
  }
  return this->start_time < elem.second;
}

bool GPFrame::operator==(const GPFrame &other) const
{
  return this->layer_index == other.layer_index && this->start_time == other.start_time;
}

CurvesGeometry &GPFrame::strokes_as_curves()
{
  return CurvesGeometry::wrap(*this->strokes);
}

int GPFrame::strokes_num() const
{
  if (this->strokes == nullptr) {
    return 0;
  }
  return this->strokes->curve_num;
}

int GPFrame::points_num() const
{
  if (this->strokes == nullptr) {
    return 0;
  }
  return this->strokes->point_num;
}

Vector<GPStroke> GPFrame::strokes_for_write()
{
  Vector<GPStroke> strokes;
  for (const int i : this->strokes_as_curves().offsets().drop_back(1).index_range()) {
    int offset = this->strokes_as_curves().offsets()[i];
    int length = this->strokes_as_curves().offsets()[i + 1] - offset;
    strokes.append({reinterpret_cast<CurvesGeometry *>(this->strokes), length, offset});
  }
  return strokes;
}

GPStroke GPFrame::add_new_stroke(int new_points_num)
{
  if (this->strokes == nullptr) {
    this->strokes = MEM_new<CurvesGeometry>(__func__);
  }
  CurvesGeometry &strokes = this->strokes_as_curves();
  int orig_last_offset = strokes.offsets().last();

  strokes.resize(strokes.points_num() + new_points_num, strokes.curves_num() + 1);
  strokes.offsets_for_write().last() = strokes.points_num();

  /* Use poly type by default. */
  strokes.curve_types_for_write().last() = CURVE_TYPE_POLY;

  strokes.tag_topology_changed();
  return {reinterpret_cast<CurvesGeometry *>(this->strokes), new_points_num, orig_last_offset};
}

}  // namespace blender::bke