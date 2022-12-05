/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup DNA
 */

#pragma once

#include "BKE_curves.hh"

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

/** GPLayer.flag */
typedef enum eGPLayerFlag {
  LAYER_SELECT = (1 << 0),
  LAYER_HIDE = (1 << 1),
  LAYER_LOCK = (1 << 2),
} eGPLayerFlag;

typedef struct GPLayer {
  /**
   * The name of the layer.
   */
  char name[128];

  /**
   * The layer flag (see `eGPLayerFlag`).
   */
  int flag;

  /* ... */
} GPLayer;

/** GPFrame.flag */
typedef enum eGPFrameFlag {
  FRAME_SELECT = (1 << 0),
} eGPFrameFlag;

typedef struct GPFrame {
  /**
   * The curves in this frame. Each individual curve is a single stroke. The CurvesGeometry
   * structure also stores attributes on the strokes and points.
   */
  CurvesGeometry *strokes;

  /**
   * The frame flag (see `eGPFrameFlag`).
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
   * The array of grease pencil frames.
   *
   * The array is always assumed to be sorted. First by the layer index of each frame, then in
   * chronological order, like so:
   *
   *  `[ frame 1 on layer 1, frame 2 on layer 1, frame 3 on layer 1, ...
   *     frame 1 on layer 2, frame 2 on layer 2, frame 3 on layer 2, ...]`
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
  GPLayerGroup();
  GPLayerGroup(const StringRefNull name);

  ~GPLayerGroup();

  IndexMask layers_index_mask();
};

class GPDataRuntime {
 public:
  /* mutable void *sbuffer */

  /**
   * Cache that maps the index of a layer to the index range of the frames in that layer.
   */
  mutable Map<int, IndexRange> frames_index_range_cache;
  mutable std::mutex frames_index_range_cache_mutex;

  IndexRange frames_index_range_cache_for_layer(int layer_index)
  {
    return frames_index_range_cache.lookup(layer_index);
  };
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

  /** Structure that defines a key for a frame. Used for ordering/sorting. */
  struct GPFrameKey {
    int layer_index;
    int start_time;
  };

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

  bool operator<(const GPFrameKey key) const;
  bool operator<(const GPFrame &other) const;

  bool operator==(const GPFrame &other) const;

  CurvesGeometry &strokes_as_curves();

  int strokes_num() const;
  int points_num() const;

  Vector<GPStroke> strokes_for_write();
  GPStroke add_new_stroke(int new_points_num);

  constexpr GPFrameKey get_frame_key() const
  {
    return {layer_index, start_time};
  };
};

class GPLayer : public ::GPLayer {
 public:
  GPLayer() : GPLayer("GP_Layer")
  {
  }

  GPLayer(const StringRefNull name);

  ~GPLayer() = default;

  bool operator==(const GPLayer &other) const;
};

class GPData : public ::GPData {
 public:
  GPData() : GPData(0, 0)
  {
  }

  GPData(const int layers_size, const int frame_size);

  GPData(const GPData &other);
  GPData &operator=(const GPData &other);
  GPData(GPData &&other);
  GPData &operator=(GPData &&other);

  ~GPData();

  Span<GPFrame> frames() const;
  const GPFrame &frames(int index) const;

  MutableSpan<GPFrame> frames_for_write();
  GPFrame &frames_for_write(int index);

  IndexRange frames_on_layer(int layer_index) const;
  IndexRange frames_on_layer(GPLayer &layer) const;
  IndexRange frames_on_active_layer() const;

  Span<GPLayer> layers() const;
  const GPLayer &layers(int index) const;
  MutableSpan<GPLayer> layers_for_write();
  GPLayer &layers_for_write(int index);

  const GPLayer &active_layer();
  GPLayer &active_layer_for_write();

  int add_layer(StringRefNull name);
  void add_layers(Array<StringRefNull> names);

  int find_layer_by_name(StringRefNull name);

  int add_frame_on_layer(int layer_index, int frame_start);
  int add_frame_on_layer(GPLayer &layer, int frame_start);
  int add_frame_on_active_layer(int frame_start);
  void add_frames_on_layer(int layer_index, Array<int> start_frames);

  int strokes_num() const;
  int points_num() const;

  void set_active_layer(int layer_index);

 private:
  const void copy_gpdata(GPData &dst, const GPData &src);
  const void move_gpdata(GPData &dst, GPData &src);

  const bool ensure_layers_array_has_size_at_least(int64_t size);
  const bool ensure_frames_array_has_size_at_least(int64_t size);

  /**
   * Creates a new frame and inserts it into the  \a frames_array so that the ordering is kept.
   * Assumes that \a frames_array is sorted and that the array has been reallocated + expaned by \a
   * reserved.
   */
  int add_frame_on_layer_initialized(int layer_index, int frame_start, int reserved);

  void update_frames_array();
};

GPData convert_old_to_new_gpencil_data(bGPdata *old_gpd);
bGPdata *convert_new_to_old_gpencil_data(const GPData &new_gpd);

}  // namespace blender::bke

#ifdef __cplusplus
}
#endif