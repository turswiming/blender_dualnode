/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 * Passes record draw commands.
 */

#include "DRW_gpu_wrapper.hh"

#include "draw_command_shared.hh"
#include "draw_handle.hh"
#include "draw_state.h"

namespace blender::draw::command {

class DrawCommandBuf;

/* -------------------------------------------------------------------- */
/** \name Recording State
 * \{ */

/**
 * Command recording state.
 * Keep track of several states and avoid redundant state changes.
 */
struct RecordingState {
  GPUShader *shader = nullptr;
  bool front_facing = true;
  bool inverted_view = false;
  DRWState pipeline_state = DRW_STATE_NO_DRAW;
  int view_clip_plane_count = 0;

  void front_facing_set(bool front_facing)
  {
    /* Facing is inverted if view is not in expected handedness. */
    front_facing = this->inverted_view == front_facing;
    /* Remove redundant changes. */
    if (assign_if_different(this->front_facing, front_facing)) {
      GPU_front_facing(!front_facing);
    }
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Regular Commands
 * \{ */

enum class Type : uint8_t {
  /**
   * None Type commands are either uninitialized or are repurposed as data storage.
   * They are skipped during submission.
   */
  None = 0,

  /** Commands stored as Undetermined in regular command buffer. */
  Barrier,
  Clear,
  Dispatch,
  DispatchIndirect,
  Draw,
  DrawIndirect,
  PushConstant,
  ResourceBind,
  ShaderBind,
  StateSet,
  StencilSet,

  /** Special commands stored in separate buffers. */
  SubPass,
  MultiDraw,
};

/**
 * The index of the group is implicit since it is known by the one who want to
 * access it. This also allows to have an indexed object to split the command
 * stream.
 */
struct Header {
  /** Command type. */
  Type type;
  /** Command index in command heap of this type. */
  uint command_index;
};

struct ShaderBind {
  GPUShader *shader;

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct ResourceBind {
  eGPUSamplerState sampler;
  int slot;
  bool is_reference;

  enum class Type : uint8_t {
    Sampler = 0,
    Image,
    UniformBuf,
    StorageBuf,
  } type;

  union {
    /** TODO: Use draw::Texture|StorageBuffer|UniformBuffer as resources as they will give more
     * debug info. */
    GPUUniformBuf *uniform_buf;
    GPUUniformBuf **uniform_buf_ref;
    GPUStorageBuf *storage_buf;
    GPUStorageBuf **storage_buf_ref;
    /** NOTE: Texture is used for both Sampler and Image binds. */
    GPUTexture *texture;
    GPUTexture **texture_ref;
  };

  ResourceBind() = default;

  ResourceBind(int slot_, GPUUniformBuf *res)
      : slot(slot_), is_reference(false), type(Type::UniformBuf), uniform_buf(res){};
  ResourceBind(int slot_, GPUUniformBuf **res)
      : slot(slot_), is_reference(true), type(Type::UniformBuf), uniform_buf_ref(res){};
  ResourceBind(int slot_, GPUStorageBuf *res)
      : slot(slot_), is_reference(false), type(Type::StorageBuf), storage_buf(res){};
  ResourceBind(int slot_, GPUStorageBuf **res)
      : slot(slot_), is_reference(true), type(Type::StorageBuf), storage_buf_ref(res){};
  ResourceBind(int slot_, draw::Image *res)
      : slot(slot_), is_reference(false), type(Type::Image), texture(draw::as_texture(res)){};
  ResourceBind(int slot_, draw::Image **res)
      : slot(slot_), is_reference(true), type(Type::Image), texture_ref(draw::as_texture(res)){};
  ResourceBind(int slot_, GPUTexture *res, eGPUSamplerState state)
      : sampler(state), slot(slot_), is_reference(false), type(Type::Sampler), texture(res){};
  ResourceBind(int slot_, GPUTexture **res, eGPUSamplerState state)
      : sampler(state), slot(slot_), is_reference(true), type(Type::Sampler), texture_ref(res){};

  void execute() const;
  std::string serialize() const;
};

struct PushConstant {
  int location;
  uint8_t array_len;
  uint8_t comp_len;
  enum class Type : uint8_t {
    IntValue = 0,
    FloatValue,
    IntReference,
    FloatReference,
  } type;
  /**
   * IMPORTANT: Data is at the end of the struct as it can span over the next commands.
   * These next commands are not real commands but just memory to hold the data and are not
   * referenced by any Command::Header.
   * This is a hack to support float4x4 copy.
   */
  union {
    int int1_value;
    int2 int2_value;
    int3 int3_value;
    int4 int4_value;
    float float1_value;
    float2 float2_value;
    float3 float3_value;
    float4 float4_value;
    const int *int_ref;
    const int2 *int2_ref;
    const int3 *int3_ref;
    const int4 *int4_ref;
    const float *float_ref;
    const float2 *float2_ref;
    const float3 *float3_ref;
    const float4 *float4_ref;
    const float4x4 *float4x4_ref;
  };

  PushConstant() = default;

  PushConstant(int loc, const float &val)
      : location(loc), array_len(1), comp_len(1), type(Type::FloatValue), float1_value(val){};
  PushConstant(int loc, const float2 &val)
      : location(loc), array_len(1), comp_len(2), type(Type::FloatValue), float2_value(val){};
  PushConstant(int loc, const float3 &val)
      : location(loc), array_len(1), comp_len(3), type(Type::FloatValue), float3_value(val){};
  PushConstant(int loc, const float4 &val)
      : location(loc), array_len(1), comp_len(4), type(Type::FloatValue), float4_value(val){};

  PushConstant(int loc, const int &val)
      : location(loc), array_len(1), comp_len(1), type(Type::IntValue), int1_value(val){};
  PushConstant(int loc, const int2 &val)
      : location(loc), array_len(1), comp_len(2), type(Type::IntValue), int2_value(val){};
  PushConstant(int loc, const int3 &val)
      : location(loc), array_len(1), comp_len(3), type(Type::IntValue), int3_value(val){};
  PushConstant(int loc, const int4 &val)
      : location(loc), array_len(1), comp_len(4), type(Type::IntValue), int4_value(val){};

  PushConstant(int loc, const float *val, int arr)
      : location(loc), array_len(arr), comp_len(1), type(Type::FloatReference), float_ref(val){};
  PushConstant(int loc, const float2 *val, int arr)
      : location(loc), array_len(arr), comp_len(2), type(Type::FloatReference), float2_ref(val){};
  PushConstant(int loc, const float3 *val, int arr)
      : location(loc), array_len(arr), comp_len(3), type(Type::FloatReference), float3_ref(val){};
  PushConstant(int loc, const float4 *val, int arr)
      : location(loc), array_len(arr), comp_len(4), type(Type::FloatReference), float4_ref(val){};
  PushConstant(int loc, const float4x4 *val)
      : location(loc), array_len(1), comp_len(16), type(Type::FloatReference), float4x4_ref(val){};

  PushConstant(int loc, const int *val, int arr)
      : location(loc), array_len(arr), comp_len(1), type(Type::IntReference), int_ref(val){};
  PushConstant(int loc, const int2 *val, int arr)
      : location(loc), array_len(arr), comp_len(2), type(Type::IntReference), int2_ref(val){};
  PushConstant(int loc, const int3 *val, int arr)
      : location(loc), array_len(arr), comp_len(3), type(Type::IntReference), int3_ref(val){};
  PushConstant(int loc, const int4 *val, int arr)
      : location(loc), array_len(arr), comp_len(4), type(Type::IntReference), int4_ref(val){};

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct Draw {
  GPUBatch *batch;
  uint vertex_len;
  uint vertex_first;
  uint instance_len;
  ResourceHandle handle;

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct DrawIndirect {
  GPUBatch *batch;
  GPUStorageBuf **indirect_buf;
  ResourceHandle handle;

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct Dispatch {
  bool is_reference;
  union {
    int3 size;
    int3 *size_ref;
  };

  Dispatch() = default;

  Dispatch(int3 group_len) : is_reference(false), size(group_len){};
  Dispatch(int3 *group_len) : is_reference(true), size_ref(group_len){};

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct DispatchIndirect {
  GPUStorageBuf **indirect_buf;

  void execute(RecordingState &state) const;
  std::string serialize() const;
};

struct Barrier {
  eGPUBarrier type;

  void execute() const;
  std::string serialize() const;
};

struct Clear {
  uint8_t clear_channels; /* #eGPUFrameBufferBits. But want to save some bits. */
  uint8_t stencil;
  float depth;
  float4 color;

  void execute() const;
  std::string serialize() const;
};

struct StateSet {
  DRWState state;

  void execute(RecordingState &recording_state) const;
  std::string serialize() const;
};

struct StencilSet {
  uint write_mask;
  uint compare_mask;
  uint reference;

  void execute() const;
  std::string serialize() const;
};

struct Undetermined {
  union {
    ShaderBind shader_bind;
    ResourceBind resource_bind;
    PushConstant push_constant;
    Draw draw;
    DrawIndirect draw_indirect;
    Dispatch dispatch;
    DispatchIndirect dispatch_indirect;
    Barrier barrier;
    Clear clear;
    StateSet state_set;
    StencilSet stencil_set;
  };

  void execute(const command::Type &type, command::RecordingState &state) const
  {
    switch (type) {
      case command::Type::ShaderBind:
        shader_bind.execute(state);
        break;
      case command::Type::ResourceBind:
        resource_bind.execute();
        break;
      case command::Type::PushConstant:
        push_constant.execute(state);
        break;
      case command::Type::Draw:
        draw.execute(state);
        break;
      case command::Type::DrawIndirect:
        draw_indirect.execute(state);
        break;
      case command::Type::Dispatch:
        dispatch.execute(state);
        break;
      case command::Type::DispatchIndirect:
        dispatch_indirect.execute(state);
        break;
      case command::Type::Barrier:
        barrier.execute();
        break;
      case command::Type::Clear:
        clear.execute();
        break;
      case command::Type::StateSet:
        state_set.execute(state);
        break;
      case command::Type::StencilSet:
        stencil_set.execute();
        break;
      default:
        /* Skip Type::None. */
        break;
    }
  }

  std::string serialize(const command::Type &type) const
  {
    switch (type) {
      case command::Type::ShaderBind:
        return shader_bind.serialize();
      case command::Type::ResourceBind:
        return resource_bind.serialize();
      case command::Type::PushConstant:
        return push_constant.serialize();
      case command::Type::Draw:
        return draw.serialize();
      case command::Type::DrawIndirect:
        return draw_indirect.serialize();
      case command::Type::Dispatch:
        return dispatch.serialize();
      case command::Type::DispatchIndirect:
        return dispatch_indirect.serialize();
      case command::Type::Barrier:
        return barrier.serialize();
      case command::Type::Clear:
        return clear.serialize();
      case command::Type::StateSet:
        return state_set.serialize();
      case command::Type::StencilSet:
        return stencil_set.serialize();
      default:
        /* Skip Type::None. */
        return "";
    }
  }
};

/** Try to keep the command size as low as possible for performance. */
BLI_STATIC_ASSERT(sizeof(Undetermined) <= 24, "One of the command type is too large.")

/** \} */

/* -------------------------------------------------------------------- */
/** \name Draw Commands
 *
 * A draw command buffer used to issue single draw commands without instance merging or any
 * other optimizations.
 * \{ */

class DrawCommandBuf {
 private:
  using ResourceIdBuf = StorageArrayBuffer<uint, 128>;

 public:
  void clear(){};

  void append_draw(Vector<command::Header> &headers,
                   Vector<command::Undetermined> &commands,
                   GPUBatch *batch,
                   uint instance_len,
                   uint vertex_len,
                   uint vertex_first,
                   ResourceHandle handle)
  {
    vertex_first = vertex_first != -1 ? vertex_first : 0;
    instance_len = instance_len != -1 ? instance_len : 1;

    int64_t index = commands.append_and_get_index({});
    headers_.append({type, static_cast<uint>(index)});
    commands[index].draw = {batch, instance_len, vertex_len, vertex_first, handle};
  }

  void bind(ResourceIdBuf &resource_id_buf)
  {
    uint total_instance = 0;
    for (DrawCommand &cmd : command_buf_) {
      int batch_vert_len, batch_inst_len;
      /* Now that GPUBatches are guaranteed to be finished, extract their parameters. */
      GPU_batch_draw_parameter_get(cmd.batch, &batch_vert_len, &batch_inst_len);
      /* Instancing attributes are not supported using the new pipeline since we use the base
       * instance to set the correct resource_id. Workaround is a storage_buf + gl_InstanceID. */
      BLI_assert(batch_inst_len == 1);

      cmd.v_count = max_ii(cmd.v_count, batch_vert_len);

      if (cmd.resource_id > 0) {
        /* Save correct offset to start of resource_id buffer region for this draw. */
        cmd.i_first = total_instance;
        total_instance += cmd.i_count;
        /* Ensure the buffer is big enough. */
        resource_id_buf.get_or_resize(total_instance - 1);

        /* Copy the resource id for all instances. */
        for (int i = cmd.i_first; i < (cmd.i_first + cmd.i_count); i++) {
          resource_id_buf[i] = cmd.resource_id;
        }
      }
    }

    if (total_instance > 0) {
      resource_id_buf.push_update();
    }
  }
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Multi Draw Commands
 *
 * For efficient rendering of large scene we strive to minimize the number of draw call and state
 * changes. This reduces the amount of work the CPU has to do. To this end, we group many rendering
 * commands and sort them per render state using Command::MultiDraw as a container.
 *
 * We sort by Command::MultiDraw index using a prefix sum on CPU.
 * Then we sort the MultiDrawUnit inside each MultiDraw by their drw_resource_id on GPU.
 *
 * For the sake of the example consider that:
 * - Command1/2 are rendering with different shaders.
 * - GPUBatch1/2 are two different mesh data block.
 * - Each column is a MultiDrawUnit.
 *
 * +---------------------------------------------------------------+------------------------------+
 * |                       CPU Timeline                            |          Granularity         |
 * +---------------------------------------------------------------+------------------------------+
 * |                    Command1                   |    Command2   |  < Command::Header           |
 * |       GPUBatch1       |       GPUBatch2       |   GPUBatch1   |  < Command::MultiDraw        |
 * |   1   |   0       0   |   0       0       0   |   1       1   |  < Front facing inverted     |
 * |  MDI  |      MDI      |          MDI          |      MDI      |  < MultiDrawIndirect emitted |
 * +---------------------------------------------------------------+------------------------------+
 * |                       GPU Timeline                            |          Granularity         |
 * +---------------------------------------------------------------+------------------------------+
 * |   4   |   2       5   |   1   |   3       4   |   6       7   |  < Resource_id (sorted)      |
 * |   1   |   1       1   |   1   |   0   |   1   |   1       1   |  < Visibility test result    |
 * |   4   |     2 + 5     |         1 + 4         |      6+7      |  < DrawCommand (compacted)   |
 * +---------------------------------------------------------------+------------------------------+
 *
 * In the example above, we will issue 4 multi draw indirect calls.
 *
 * \{ */

struct MultiDrawBuf {
  void clear(){};

  void bind(){};

  uint append(uint instance_len, uint vertex_len, uint vertex_first, ResourceHandle handle){};
};

/** \} */

};  // namespace blender::draw::command