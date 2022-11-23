/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup DNA
 */

#pragma once

#include "BKE_curves.hh"

#include "BLI_index_mask_ops.hh"

#include "DNA_ID.h"
#include "DNA_curves_types.h"
#include "DNA_customdata_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
namespace blender::bke {
class GPDataRuntime;
}  // namespace blender::bke
using GPDataRuntimeHandle = blender::bke::GPDataRuntime;
#else
typedef struct GPDataRuntimeHandle GPDataRuntimeHandle;
#endif

typedef struct GPLayerGroup {
  /**
   * An array of GPLayerGroup's. A layer group can have N >= 0 number of layer group children.
   */
  struct GPLayerGroup *children;
  int children_size;

  /**
   * An array of indices to the layers in GPData.layers_array. These are the layers contained in
   * the group.
   */
  int *layer_indices;
  int layer_indices_size;

  /**
   * The name of the layer group.
   */
  char name[128];

  /* ... */
} GPLayerGroup;

typedef struct GPLayer {
  /**
   * The name of the layer.
   */
  char name[128];

  /**
   * The layer flag.
   */
  int flag;

  /* ... */
} GPLayer;

typedef struct GPFrame {
  /**
   * The curves in this frame. Each individual curve is a single stroke. The CurvesGeometry
   * structure also stores attributes on the strokes and points.
   */
  CurvesGeometry *strokes;

  /**
   * The frame flag.
   */
  int flag;

  /**
   * The index of the layer in GPData.layers_array that this frame is in.
   */
  int layer_index;

  /**
   * The start frame in the scene that the grease pencil frame is displayed.
   */
  int start_time;
  int end_time; /* UNUSED for now. */

  /* ... */
} GPFrame;

typedef struct GPData {
  /**
   * The array of grease pencil frames. This is kept in chronological order (tiebreaks for two
   * frames on different layers are resloved by the order of the layers).
   */
  GPFrame *frames_array;
  int frames_size;

  /**
   * All attributes stored on the frames.
   */
  CustomData frame_data;

  /**
   * The array of grease pencil layers.
   */
  GPLayer *layers_array;
  int layers_size;

  /**
   * The index of the active layer in the GPData.layers_array.
   */
  int active_layer_index;

  /**
   * The root layer group. This must not be nullptr.
   */
  GPLayerGroup *default_group;

  /**
   * The runtime data.
   */
  GPDataRuntimeHandle *runtime;
} GPData;

/**
 * This would be the new Grease Pencil ID structure. This is where the animation data, materials,
 * etc. are stored. Layers, Frames, Groups and RuntimeData would be stored in GPData.
 */
typedef struct GreasePencil {
  ID id;
  /* Animation data (must be immediately after id). */
  struct AnimData *adt;

  /* Pointer to the actual data-block containing the frames, layers and layer groups. Note: This is
   * stored as a pointer to easily wrap it in a class. */
  GPData *grease_pencil_data;

  /* GreasePencil flag. */
  int flag;

  /** Materials array. */
  struct Material **mat;
  /** Total materials. */
  short totcol;

  /* ... */
} GreasePencil;

namespace blender::bke {

class GPLayerGroup : ::GPLayerGroup { /* Unused for now. Placeholder class. */
 public:
  GPLayerGroup()
  {
    this->children = nullptr;
    this->children_size = 0;
    this->layer_indices = nullptr;
    this->layer_indices_size = 0;
  }

  GPLayerGroup(const StringRefNull name) : GPLayerGroup()
  {
    BLI_assert(name.size() < 128);
    strcpy(this->name, name.c_str());
  }

  ~GPLayerGroup()
  {
    /* Recursivly free the children of this layer group first. */
    for (int i = 0; i < this->children_size; i++) {
      MEM_delete(&this->children[i]);
    }
    /* Then free its data. */
    MEM_SAFE_FREE(this->children);
    MEM_SAFE_FREE(this->layer_indices);
  }

  IndexMask layers_index_mask()
  {
    return {reinterpret_cast<int64_t>(this->layer_indices), this->layer_indices_size};
  }
};

class GPDataRuntime {
 public:
  /* mutable void *sbuffer */

  /**
   * Cache that maps the index of a layer to the index mask of the frames in that layer.
   */
  mutable Map<int, Vector<int64_t>> frame_index_masks_cache;
  mutable std::mutex frame_index_masks_cache_mutex;

  IndexMask frame_index_masks_cache_for_layer(int layer_index)
  {
    return frame_index_masks_cache.lookup(layer_index).as_span();
  }
};

/**
 * A wrapper class around a single curve in GPFrame.strokes (CurvesGeometry). It holds the offset
 * of where to find the stroke in the frame and it's size.
 * This class is only meant to facilitate the handling of individual strokes.
 */
class GPStroke {
 public:
  GPStroke(CurvesGeometry *geometry, int num_points, int offset)
      : geometry_(geometry), points_num_(num_points), offset_(offset){};

  ~GPStroke() = default;

  int points_num() const
  {
    return points_num_;
  }

  /**
   * Start index of this stroke in the points array of geometry_.
   */
  int points_offset() const
  {
    return offset_;
  }

  Span<float3> points_positions() const;
  MutableSpan<float3> points_positions_for_write() const;
  void transform(float4x4 matrix);

 private:
  CurvesGeometry *geometry_ = nullptr;
  int points_num_ = 0;
  int offset_;
};

class GPFrame : public ::GPFrame {
 public:
  GPFrame() : GPFrame(-1, -1)
  {
  }

  GPFrame(int start_frame) : GPFrame(start_frame, -1)
  {
  }

  GPFrame(int start_frame, int end_frame);

  GPFrame(const GPFrame &other);
  GPFrame &operator=(const GPFrame &other);
  GPFrame(GPFrame &&other);
  GPFrame &operator=(GPFrame &&other);

  ~GPFrame();

  bool operator<(const GPFrame &other) const;
  /* Assumes that elem.first is the layer index and elem.second is the start time. */
  bool operator<(const std::pair<int, int> elem) const;

  bool operator==(const GPFrame &other) const;

  CurvesGeometry &strokes_as_curves();

  int strokes_num() const;
  int points_num() const;

  Vector<GPStroke> strokes_for_write();
  GPStroke add_new_stroke(int new_points_num);
};

class GPLayer : public ::GPLayer {
 public:
  GPLayer() : GPLayer("GP_Layer")
  {
  }

  GPLayer(const StringRefNull name)
  {
    strcpy(this->name, name.c_str());
  }

  ~GPLayer() = default;

  bool operator==(const GPLayer &other) const
  {
    return STREQ(this->name, other.name);
  }
};

class GPData : public ::GPData {
 public:
  GPData() : GPData(0, 0)
  {
  }

  GPData(const int layers_size, const int frame_size)
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

  GPData(const GPData &other) : GPData(other.layers_size, other.frames_size)
  {
    copy_gpdata(*this, other);
  }

  GPData &operator=(const GPData &other)
  {
    if (this != &other) {
      copy_gpdata(*this, other);
    }
    return *this;
  }

  GPData(GPData &&other) : GPData(other.layers_size, other.frames_size)
  {
    move_gpdata(*this, other);
  }

  GPData &operator=(GPData &&other)
  {
    if (this != &other) {
      move_gpdata(*this, other);
    }
    return *this;
  }

  ~GPData()
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

  Span<GPFrame> frames() const
  {
    return {reinterpret_cast<const GPFrame *>(this->frames_array), this->frames_size};
  }

  const GPFrame &frames(int index) const
  {
    return this->frames()[index];
  }

  MutableSpan<GPFrame> frames_for_write()
  {
    return {reinterpret_cast<GPFrame *>(this->frames_array), this->frames_size};
  }

  GPFrame &frames_for_write(int index)
  {
    return this->frames_for_write()[index];
  }

  IndexMask frames_on_layer(int layer_index) const
  {
    if (layer_index < 0 || layer_index > this->layers_size) {
      return IndexMask();
    }

    /* If the indices are cached for this layer, use the cache. */
    if (this->runtime->frame_index_masks_cache.contains(layer_index)) {
      return this->runtime->frame_index_masks_cache_for_layer(layer_index);
    }

    /* A double checked lock. */
    std::scoped_lock{this->runtime->frame_index_masks_cache_mutex};
    if (this->runtime->frame_index_masks_cache.contains(layer_index)) {
      return this->runtime->frame_index_masks_cache_for_layer(layer_index);
    }

    Vector<int64_t> indices;
    const IndexMask mask = index_mask_ops::find_indices_based_on_predicate(
        IndexMask(this->frames_size), 1024, indices, [&](const int index) {
          return this->frames()[index].layer_index == layer_index;
        });

    /* Cache the resulting index mask. */
    this->runtime->frame_index_masks_cache.add(layer_index, std::move(indices));
    return mask;
  }

  IndexMask frames_on_layer(GPLayer &layer) const
  {
    int index = this->layers().first_index_try(layer);
    if (index == -1) {
      return IndexMask();
    }
    return frames_on_layer(index);
  }

  IndexMask frames_on_active_layer() const
  {
    return frames_on_layer(this->active_layer_index);
  }

  Span<GPLayer> layers() const
  {
    return {reinterpret_cast<const GPLayer *>(this->layers_array), this->layers_size};
  }

  const GPLayer &layers(int index) const
  {
    return layers()[index];
  }

  MutableSpan<GPLayer> layers_for_write()
  {
    return {reinterpret_cast<GPLayer *>(this->layers_array), this->layers_size};
  }

  GPLayer &layers_for_write(int index)
  {
    return layers_for_write()[index];
  }

  const GPLayer &active_layer()
  {
    return this->layers()[this->active_layer_index];
  }

  GPLayer &active_layer_for_write()
  {
    return this->layers_for_write()[this->active_layer_index];
  }

  int add_layer(StringRefNull name)
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

  void add_layers(Array<StringRefNull> names)
  {
    for (StringRefNull name : names) {
      this->add_layer(name);
    }
  }

  int find_layer_by_name(StringRefNull name)
  {
    for (const int i : this->layers().index_range()) {
      if (STREQ(this->layers(i).name, name.c_str())) {
        return i;
      }
    }
    return -1;
  }

  int add_frame_on_layer(int layer_index, int frame_start)
  {
    /* TODO: Check for collisions before resizing the array. */
    if (!ensure_frames_array_has_size_at_least(this->frames_size + 1)) {
      return -1;
    }

    return add_frame_on_layer_initialized(layer_index, frame_start, 1);
  }

  int add_frame_on_layer(GPLayer &layer, int frame_start)
  {
    int index = this->layers().first_index_try(layer);
    if (index == -1) {
      return -1;
    }
    return add_frame_on_layer(index, frame_start);
  }

  int add_frame_on_active_layer(int frame_start)
  {
    return add_frame_on_layer(active_layer_index, frame_start);
  }

  void add_frames_on_layer(int layer_index, Array<int> start_frames)
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

  int strokes_num() const
  {
    /* TODO: could be done with parallel_for */
    int count = 0;
    for (const GPFrame &gpf : this->frames()) {
      count += gpf.strokes_num();
    }
    return count;
  }

  int points_num() const
  {
    /* TODO: could be done with parallel_for */
    int count = 0;
    for (const GPFrame &gpf : this->frames()) {
      count += gpf.points_num();
    }
    return count;
  }

  void set_active_layer(int layer_index)
  {
    if (layer_index < 0 || layer_index >= this->layers_size) {
      return;
    }
    this->active_layer_index = layer_index;
  }

 private:
  const void copy_gpdata(GPData &dst, const GPData &src)
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

  const void move_gpdata(GPData &dst, GPData &src)
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

  const bool ensure_layers_array_has_size_at_least(int64_t size)
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

  const bool ensure_frames_array_has_size_at_least(int64_t size)
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

  /**
   * Creates a new frame and inserts it into the  \a frames_array so that the ordering is kept.
   * Assumes that \a frames_array is sorted and that the array has been reallocated + expaned by \a
   * reserved.
   */
  int add_frame_on_layer_initialized(int layer_index, int frame_start, int reserved)
  {
    /* Create the new frame. */
    GPFrame new_frame(frame_start);
    new_frame.layer_index = layer_index;

    int last_index = this->frames_size - reserved - 1;

    /* Check if the frame can be appended at the end. */
    if (this->frames_size == 0 || this->frames_size == reserved ||
        this->frames(last_index) < std::pair<int, int>(layer_index, frame_start)) {
      this->frames_for_write(last_index + 1) = std::move(new_frame);
      return last_index + 1;
    }

    /* Look for the first frame that is equal or greater than the new frame. */
    auto it = std::lower_bound(this->frames().begin(),
                               this->frames().drop_back(reserved).end(),
                               std::pair<int, int>(layer_index, frame_start));
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

  void update_frames_array()
  {
    /* Make sure frames are ordered chronologically and by layer order. */
    std::sort(this->frames_for_write().begin(), this->frames_for_write().end());

    /* Clear the cached indices since they are (probably) no longer valid. */
    this->runtime->frame_index_masks_cache.clear();
  }
};

}  // namespace blender::bke

#ifdef __cplusplus
}
#endif