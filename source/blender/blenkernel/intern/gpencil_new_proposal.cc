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

/* GPLayerGroup */
GPLayerGroup::GPLayerGroup()
{
  this->children = nullptr;
  this->children_size = 0;
  this->layer_indices = nullptr;
  this->layer_indices_size = 0;
}

GPLayerGroup::GPLayerGroup(const StringRefNull name) : GPLayerGroup()
{
  BLI_assert(name.size() < 128);
  strcpy(this->name, name.c_str());
}

GPLayerGroup::~GPLayerGroup()
{
  /* Recursivly free the children of this layer group first. */
  for (int i = 0; i < this->children_size; i++) {
    MEM_delete(&this->children[i]);
  }
  /* Then free its data. */
  MEM_SAFE_FREE(this->children);
  MEM_SAFE_FREE(this->layer_indices);
}

IndexMask GPLayerGroup::layers_index_mask()
{
  return {reinterpret_cast<int64_t>(this->layer_indices), this->layer_indices_size};
}

/* GPStroke */
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

/* GPFrame */
GPFrame::GPFrame(int layer_index)
{
  this->layer_index = layer_index;
  this->strokes = nullptr;
}

GPFrame::GPFrame(int layer_index, int start_time, int end_time)
{
  this->layer_index = layer_index;
  this->start_time = start_time;
  this->end_time = end_time;
  this->strokes = nullptr;
}

GPFrame::GPFrame(const GPFrame &other) : GPFrame(other.layer_index)
{
  if (other.strokes != nullptr) {
    /* Make sure old strokes are freed before copying. */
    MEM_SAFE_FREE(this->strokes);
    this->strokes = MEM_new<CurvesGeometry>(__func__);

    *reinterpret_cast<CurvesGeometry *>(this->strokes) = CurvesGeometry::wrap(*other.strokes);
  }
  this->start_time = other.start_time;
  this->end_time = other.end_time;
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

GPFrame::GPFrame(GPFrame &&other) : GPFrame(other.layer_index)
{
  if (this != &other) {
    std::swap(this->strokes, other.strokes);
    other.strokes = nullptr;
  }
  this->start_time = other.start_time;
  this->end_time = other.end_time;
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

bool GPFrame::operator<(const GPFrameKey key) const
{
  if (this->layer_index == key.layer_index) {
    return this->start_time < key.start_time;
  }
  return this->layer_index < key.layer_index;
}

bool GPFrame::operator<(const GPFrame &other) const
{
  if (this->layer_index == other.layer_index) {
    return this->start_time < other.start_time;
  }
  return this->layer_index < other.layer_index;
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

/* GPLayer */
GPLayer::GPLayer(const StringRefNull name)
{
  strcpy(this->name, name.c_str());
}

bool GPLayer::operator==(const GPLayer &other) const
{
  return STREQ(this->name, other.name);
}

/* GPData */
GPData::GPData(const int layers_size, const int frame_size)
{
  BLI_assert(layers_size >= 0);
  BLI_assert(frame_size >= 0);

  this->frames_size = frame_size;
  this->layers_size = layers_size;

  if (this->frames_size > 0) {
    this->frames_array = reinterpret_cast<::GPFrame *>(
        MEM_malloc_arrayN(this->frames_size, sizeof(::GPFrame), __func__));
    default_construct_n(reinterpret_cast<GPFrame *>(this->frames_array), this->frames_size);
  }
  else {
    this->frames_array = nullptr;
  }
  CustomData_reset(&this->frame_data);

  if (this->layers_size > 0) {
    this->layers_array = reinterpret_cast<::GPLayer *>(
        MEM_malloc_arrayN(this->layers_size, sizeof(::GPLayer), __func__));
    default_construct_n(reinterpret_cast<GPLayer *>(this->layers_array), this->layers_size);
    this->active_layer_index = 0;
  }
  else {
    this->layers_array = nullptr;
    this->active_layer_index = -1;
  }

  this->default_group = MEM_new<::GPLayerGroup>(__func__);

  this->runtime = MEM_new<GPDataRuntime>(__func__);
}

GPData::GPData(const GPData &other) : GPData(other.layers_size, other.frames_size)
{
  copy_gpdata(*this, other);
}

GPData &GPData::operator=(const GPData &other)
{
  if (this != &other) {
    copy_gpdata(*this, other);
  }
  return *this;
}

GPData::GPData(GPData &&other) : GPData(other.layers_size, other.frames_size)
{
  move_gpdata(*this, other);
}

GPData &GPData::operator=(GPData &&other)
{
  if (this != &other) {
    move_gpdata(*this, other);
  }
  return *this;
}

GPData::~GPData()
{
  /* Free frames and frame custom data. */
  destruct_n(reinterpret_cast<GPFrame *>(this->frames_array), this->frames_size);
  MEM_SAFE_FREE(this->frames_array);
  CustomData_free(&this->frame_data, this->frames_size);

  /* Free layer and layer groups. */
  destruct_n(reinterpret_cast<GPLayer *>(this->layers_array), this->layers_size);
  MEM_SAFE_FREE(this->layers_array);
  MEM_delete(reinterpret_cast<GPLayerGroup *>(this->default_group));
  this->default_group = nullptr;

  /* Free the runtime structure. */
  MEM_delete(this->runtime);
  this->runtime = nullptr;
}

Span<GPFrame> GPData::frames() const
{
  return {reinterpret_cast<const GPFrame *>(this->frames_array), this->frames_size};
}

const GPFrame &GPData::frames(int index) const
{
  return this->frames()[index];
}

MutableSpan<GPFrame> GPData::frames_for_write()
{
  return {reinterpret_cast<GPFrame *>(this->frames_array), this->frames_size};
}

GPFrame &GPData::frames_for_write(int index)
{
  return this->frames_for_write()[index];
}

IndexRange GPData::frames_on_layer(int layer_index) const
{
  if (layer_index < 0 || layer_index > this->layers_size) {
    return {};
  }

  /* If the indices are cached for this layer, use the cache. */
  if (this->runtime->frames_index_range_cache.contains(layer_index)) {
    return this->runtime->frames_index_range_cache_for_layer(layer_index);
  }

  /* A double checked lock. */
  std::scoped_lock{this->runtime->frames_index_range_cache_mutex};
  if (this->runtime->frames_index_range_cache.contains(layer_index)) {
    return this->runtime->frames_index_range_cache_for_layer(layer_index);
  }

  GPFrame search_val{layer_index};

  auto it_lower = std::lower_bound(this->frames().begin(),
                                   this->frames().end(),
                                   search_val,
                                   [](const GPFrame &frame_A, const GPFrame &frame_B) {
                                     return frame_A.layer_index < frame_B.layer_index;
                                   });
  auto it_upper = std::upper_bound(this->frames().begin(),
                                   this->frames().end(),
                                   search_val,
                                   [](const GPFrame &frame_A, const GPFrame &frame_B) {
                                     return frame_A.layer_index < frame_B.layer_index;
                                   });

  /* Could not find this layer. */
  if (it_lower == this->frames().end()) {
    return {};
  }

  /* Get the index of the first frame. */
  int start_idx = std::distance(this->frames().begin(), it_lower);
  /* Calculate size of the layer. */
  int frames_size = std::distance(it_lower, it_upper);

  /* Cache the resulting index range. */
  this->runtime->frames_index_range_cache.add(layer_index, {start_idx, frames_size});
  return {start_idx, frames_size};
}

IndexRange GPData::frames_on_layer(GPLayer &layer) const
{
  int index = this->layers().first_index_try(layer);
  if (index == -1) {
    return {};
  }
  return frames_on_layer(index);
}

IndexRange GPData::frames_on_active_layer() const
{
  return frames_on_layer(this->active_layer_index);
}

Span<GPLayer> GPData::layers() const
{
  return {reinterpret_cast<const GPLayer *>(this->layers_array), this->layers_size};
}

const GPLayer &GPData::layers(int index) const
{
  return layers()[index];
}

MutableSpan<GPLayer> GPData::layers_for_write()
{
  return {reinterpret_cast<GPLayer *>(this->layers_array), this->layers_size};
}

GPLayer &GPData::layers_for_write(int index)
{
  return layers_for_write()[index];
}

const GPLayer &GPData::active_layer()
{
  return this->layers()[this->active_layer_index];
}

GPLayer &GPData::active_layer_for_write()
{
  return this->layers_for_write()[this->active_layer_index];
}

int GPData::add_layer(StringRefNull name)
{
  /* Ensure that the layer array has enough space. */
  if (!ensure_layers_array_has_size_at_least(this->layers_size + 1)) {
    return -1;
  }

  GPLayer new_layer(name);
  /* Move new_layer to the end in the array. */
  this->layers_for_write().last() = std::move(new_layer);
  return this->layers_size - 1;
}

void GPData::add_layers(Array<StringRefNull> names)
{
  for (StringRefNull name : names) {
    this->add_layer(name);
  }
}

int GPData::find_layer_by_name(StringRefNull name)
{
  for (const int i : this->layers().index_range()) {
    if (STREQ(this->layers(i).name, name.c_str())) {
      return i;
    }
  }
  return -1;
}

int GPData::add_frame_on_layer(int layer_index, int frame_start)
{
  /* TODO: Check for collisions before resizing the array. */
  if (!ensure_frames_array_has_size_at_least(this->frames_size + 1)) {
    return -1;
  }

  return add_frame_on_layer_initialized(layer_index, frame_start, 1);
}

int GPData::add_frame_on_layer(GPLayer &layer, int frame_start)
{
  int index = this->layers().first_index_try(layer);
  if (index == -1) {
    return -1;
  }
  return add_frame_on_layer(index, frame_start);
}

int GPData::add_frame_on_active_layer(int frame_start)
{
  return add_frame_on_layer(this->active_layer_index, frame_start);
}

void GPData::add_frames_on_layer(int layer_index, Array<int> start_frames)
{
  int new_frames_size = start_frames.size();
  /* TODO: Check for collisions before resizing the array. */
  if (!ensure_frames_array_has_size_at_least(this->frames_size + new_frames_size)) {
    return;
  }

  int reserved = new_frames_size;
  for (int start_frame : start_frames) {
    add_frame_on_layer_initialized(layer_index, start_frame, reserved);
    reserved--;
  }
}

int GPData::strokes_num() const
{
  /* TODO: could be done with parallel_for */
  int count = 0;
  for (const GPFrame &gpf : this->frames()) {
    count += gpf.strokes_num();
  }
  return count;
}

int GPData::points_num() const
{
  /* TODO: could be done with parallel_for */
  int count = 0;
  for (const GPFrame &gpf : this->frames()) {
    count += gpf.points_num();
  }
  return count;
}

void GPData::set_active_layer(int layer_index)
{
  if (layer_index < 0 || layer_index >= this->layers_size) {
    return;
  }
  this->active_layer_index = layer_index;
}

const void GPData::copy_gpdata(GPData &dst, const GPData &src)
{
  /* Make sure previous frame data is freed. */
  MEM_SAFE_FREE(dst.frames_array);
  CustomData_free(&dst.frame_data, dst.frames_size);

  /* Copy frame data. */
  dst.frames_size = src.frames_size;
  dst.frames_array = reinterpret_cast<::GPFrame *>(
      MEM_malloc_arrayN(dst.frames_size, sizeof(::GPFrame), __func__));
  uninitialized_copy_n(reinterpret_cast<GPFrame *>(src.frames_array),
                       src.frames_size,
                       reinterpret_cast<GPFrame *>(dst.frames_array));
  CustomData_copy(&src.frame_data, &dst.frame_data, CD_MASK_ALL, CD_DUPLICATE, dst.frames_size);

  /* Make sure layer data is freed then copy it over. */
  MEM_SAFE_FREE(dst.layers_array);
  dst.layers_size = src.layers_size;
  dst.layers_array = reinterpret_cast<::GPLayer *>(
      MEM_malloc_arrayN(dst.layers_size, sizeof(::GPLayer), __func__));
  uninitialized_copy_n(reinterpret_cast<GPLayer *>(src.layers_array),
                       src.layers_size,
                       reinterpret_cast<GPLayer *>(dst.layers_array));
  dst.active_layer_index = src.active_layer_index;

  /* Copy layer default group. */
  *dst.default_group = *src.default_group;
}

const void GPData::move_gpdata(GPData &dst, GPData &src)
{
  /* Move frame data. */
  dst.frames_size = src.frames_size;
  std::swap(dst.frames_array, src.frames_array);
  std::swap(dst.frame_data, src.frame_data);
  MEM_SAFE_FREE(src.frames_array);
  CustomData_free(&src.frame_data, src.frames_size);
  src.frames_size = 0;

  /* Move layer data. */
  dst.layers_size = src.layers_size;
  std::swap(dst.layers_array, src.layers_array);
  dst.active_layer_index = src.active_layer_index;
  MEM_SAFE_FREE(src.layers_array);
  src.layers_size = 0;
  src.active_layer_index = -1;

  /* Move layer group and runtime pointers. */
  std::swap(dst.default_group, src.default_group);
  std::swap(dst.runtime, src.runtime);
}

const bool GPData::ensure_layers_array_has_size_at_least(int64_t size)
{
  if (this->layers_size > size) {
    return true;
  }

  int old_size = this->layers_size;
  this->layers_size = size;

  ::GPLayer *new_array = reinterpret_cast<::GPLayer *>(
      MEM_calloc_arrayN(this->layers_size, sizeof(::GPLayer), __func__));
  if (new_array == nullptr) {
    return false;
  }

  if (this->layers_array != nullptr) {
    /* Since the layers have default move constructors, we just use memcpy here. */
    memcpy(new_array, this->layers_array, old_size * sizeof(::GPLayer));
    MEM_SAFE_FREE(this->layers_array);
  }
  this->layers_array = new_array;

  return true;
}

const bool GPData::ensure_frames_array_has_size_at_least(int64_t size)
{
  if (this->frames_size > size) {
    return true;
  }

  int old_size = this->frames_size;
  this->frames_size = size;

  ::GPFrame *new_array = reinterpret_cast<::GPFrame *>(
      MEM_malloc_arrayN(this->frames_size, sizeof(::GPFrame), __func__));
  if (new_array == nullptr) {
    return false;
  }

  if (this->frames_array != nullptr) {
    uninitialized_relocate_n(reinterpret_cast<GPFrame *>(this->frames_array),
                             old_size,
                             reinterpret_cast<GPFrame *>(new_array));
    default_construct_n(reinterpret_cast<GPFrame *>(new_array + old_size),
                        this->frames_size - old_size);
    MEM_SAFE_FREE(this->frames_array);
    this->frames_array = new_array;
  }
  else {
    this->frames_array = new_array;
    default_construct_n(reinterpret_cast<GPFrame *>(this->frames_array), this->frames_size);
  }
  return true;
}

int GPData::add_frame_on_layer_initialized(int layer_index, int start_time, int reserved)
{
  /* Create the new frame. */
  GPFrame new_frame(layer_index);
  new_frame.start_time = start_time;

  int last_index = this->frames_size - reserved - 1;

  /* Check if the frame can be appended at the end. */
  if (this->frames_size == 0 || this->frames_size == reserved ||
      this->frames(last_index) < new_frame.get_frame_key()) {
    this->frames_for_write(last_index + 1) = std::move(new_frame);
    return last_index + 1;
  }

  /* Look for the first frame that is equal or greater than the new frame. */
  auto it = std::lower_bound(
      this->frames().begin(), this->frames().drop_back(reserved).end(), new_frame.get_frame_key());
  /* Get the index of the frame. */
  int index = std::distance(this->frames().begin(), it);
  /* Move all the frames and make space at index. */
  initialized_reversed_move_n(reinterpret_cast<GPFrame *>(this->frames_array + index),
                              this->frames_size - index - 1,
                              reinterpret_cast<GPFrame *>(this->frames_array + index + 1));
  /* Move the new frame into the space at index. */
  this->frames_for_write(index) = std::move(new_frame);

  return index;
}

void GPData::update_frames_array()
{
  /* Make sure frames are ordered by layers and chronologically. */
  std::sort(this->frames_for_write().begin(), this->frames_for_write().end());

  /* Clear the cached reanges since they are (probably) no longer valid. */
  this->runtime->frames_index_range_cache.clear();
}

}  // namespace blender::bke