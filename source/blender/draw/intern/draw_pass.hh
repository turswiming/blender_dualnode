/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

#pragma once

/** \file
 * \ingroup draw
 *
 * Passes record draw commands. There exists different pass types for different purpose but they
 * only change in resource load (memory & CPU usage). They can be swapped without any functional
 * change.
 *
 * `Pass`:
 * Should be used on heavy load passes such as ones that may contain scene objects. Draw call
 * submission is optimized for large number of draw calls. But has a significant overhead per
 * #Pass. Use many #PassSub along with a main #Pass to reduce the overhead and allow groupings of
 * commands.
 *
 * `PassSub`:
 * A lightweight #Pass that lives inside a main #Pass. It can only be created from #Pass.sub()
 * and is auto managed. This mean it can be created, filled and thrown away. A #PassSub reference
 * is valid until the next #Pass.init() of the parent pass. Commands recorded inside a #PassSub are
 * inserted inside the parent #Pass where the sub have been created durring submission.
 *
 *
 * NOTE: A pass can be recorded once and resubmitted any number of time. This can be a good
 * optimization for passes that are always the same for each frame. The only thing to be aware of
 * is the life time of external resources. If a pass contains draw-calls with non default
 * ResourceHandle (not 0) or a reference to any non static resources (GPUBatch, PushConstant ref,
 * ResourceBind ref) it will have to be re-recorded if any of these reference becomes invalid.
 */

#include "BLI_vector.hh"
#include "DRW_gpu_wrapper.hh"

#include "draw_command.hh"
#include "draw_handle.hh"
#include "draw_shader_shared.h"
#include "draw_state.h"

namespace blender::draw {

using namespace blender::draw;
using namespace blender::draw::command;

class Manager;

/* -------------------------------------------------------------------- */
/** \name Pass API
 * \{ */

/**
 * Public API of a draw pass.
 */
class Pass {
  friend Manager;

 private:
  /** Highest level of the command stream. Split command stream in different command types. */
  Vector<command::Header> headers_;
  /** Commands referenced by headers (which contains their types). */
  Vector<command::Undetermined> commands_;
  /** Sub-passes referenced by headers. */
  // Vector<Pass> sub_passes_;
  /** Currently bound shader. Used for interface queries. */
  GPUShader *shader_;

 public:
  const char *debug_name;

  Pass() = delete;
  Pass(const char *name = "Unamed") : debug_name(name){};

  /**
   * Reset the pass command pool.
   */
  void init();

  /**
   * Create a sub-pass inside this pass. The sub-pass memory is auto managed.
   */
  // Sub &sub();

  /**
   * Changes the fixed function pipeline state.
   * Starts as DRW_STATE_NO_DRAW at the start of a Pass submission.
   * SubPass inherit previous pass state.
   *
   * IMPORTANT: This does not set the stencil mask/reference values. Add a call to state_stencil()
   * to ensure correct behavior of stencil aware draws.
   */
  void state_set(DRWState state);

  /**
   * Clear the current frame-buffer.
   */
  void clear_color(float4 color);
  void clear_depth(float depth);
  void clear_stencil(uint8_t stencil);
  void clear_depth_stencil(float depth, uint8_t stencil);
  void clear_color_depth_stencil(float4 color, float depth, uint8_t stencil);

  /**
   * Reminders:
   * - (compare_mask & reference) is what is tested against (compare_mask & stencil_value)
   *   stencil_value being the value stored in the stencil buffer.
   * - (write-mask & reference) is what gets written if the test condition is fulfilled.
   */
  void state_stencil(uint8_t write_mask, uint8_t reference, uint8_t compare_mask);

  /**
   * Bind a shader. Any following bind() or push_constant() call will use its interface.
   */
  void shader_set(GPUShader *shader);

  /**
   * Record a draw call.
   * NOTE: Setting the count or first to -1 will use the values from the batch.
   * NOTE: An instance or vertex count of 0 will discard the draw call. It will not be recorded.
   */
  void draw(GPUBatch *batch,
            uint instance_len = -1,
            uint vertex_len = -1,
            uint vertex_first = -1,
            ResourceHandle handle = {0});
  void draw(GPUBatch *batch,
            StorageBuffer<DrawCommand> &indirect_buffer,
            ResourceHandle handle = {0});

  /** Shorter version for the common case. */
  void draw(GPUBatch *batch, ResourceHandle handle);

  /**
   * Record a procedural draw call. Geometry is **NOT** source from a GPUBatch.
   * NOTE: An instance or vertex count of 0 will discard the draw call. It will not be recorded.
   */
  void draw(GPUPrimType primitive,
            uint instance_len,
            uint vertex_len,
            uint vertex_first = -1,
            ResourceHandle handle = {0});
  void draw(GPUPrimType primitive,
            StorageBuffer<DrawCommand> &indirect_buffer,
            ResourceHandle handle = {0});

  /**
   * Record a compute dispatch call.
   */
  void dispatch(int3 group_len);
  void dispatch(int3 *group_len);
  void dispatch(StorageBuffer<DispatchCommand> &indirect_buffer);

  /**
   * Record a barrier call to synchronize arbitrary load/store operation between draw calls.
   */
  void barrier(eGPUBarrier type);

  /**
   * Bind a shader resource.
   *
   * Reference versions are to be used when the resource might be resize / realloc or even change
   * between the time it is referenced and the time it is dereferenced for drawing.
   *
   * IMPORTANT: Will keep a reference to the data and dereference it upon drawing. Make sure data
   * still alive until pass submission.
   */
  void bind(const char *name, GPUStorageBuf *buffer);
  void bind(const char *name, GPUUniformBuf *buffer);
  void bind(const char *name, draw::Image *image);
  void bind(const char *name, GPUTexture *texture, eGPUSamplerState state = GPU_SAMPLER_MAX);
  void bind(const char *name, GPUStorageBuf **buffer);
  void bind(const char *name, GPUUniformBuf **buffer);
  void bind(const char *name, draw::Image **image);
  void bind(const char *name, GPUTexture **texture, eGPUSamplerState state = GPU_SAMPLER_MAX);
  void bind(int slot, GPUStorageBuf *buffer);
  void bind(int slot, GPUUniformBuf *buffer);
  void bind(int slot, draw::Image *image);
  void bind(int slot, GPUTexture *texture, eGPUSamplerState state = GPU_SAMPLER_MAX);
  void bind(int slot, GPUStorageBuf **buffer);
  void bind(int slot, GPUUniformBuf **buffer);
  void bind(int slot, draw::Image **image);
  void bind(int slot, GPUTexture **texture, eGPUSamplerState state = GPU_SAMPLER_MAX);

  /**
   * Update a shader constant.
   *
   * Reference versions are to be used when the resource might change between the time it is
   * referenced and the time it is dereferenced for drawing.
   *
   * IMPORTANT: Will keep a reference to  the data and dereference it upon drawing. Make sure data
   * still alive until pass submission.
   *
   * NOTE: bool reference version is expected to take bool1 reference which is aliased to int.
   */
  void push_constant(const char *name, const float &data);
  void push_constant(const char *name, const float2 &data);
  void push_constant(const char *name, const float3 &data);
  void push_constant(const char *name, const float4 &data);
  void push_constant(const char *name, const int &data);
  void push_constant(const char *name, const int2 &data);
  void push_constant(const char *name, const int3 &data);
  void push_constant(const char *name, const int4 &data);
  void push_constant(const char *name, const bool &data);
  void push_constant(const char *name, const float4x4 &data);
  void push_constant(const char *name, const float *data, int array_len = 1);
  void push_constant(const char *name, const float2 *data, int array_len = 1);
  void push_constant(const char *name, const float3 *data, int array_len = 1);
  void push_constant(const char *name, const float4 *data, int array_len = 1);
  void push_constant(const char *name, const int *data, int array_len = 1);
  void push_constant(const char *name, const int2 *data, int array_len = 1);
  void push_constant(const char *name, const int3 *data, int array_len = 1);
  void push_constant(const char *name, const int4 *data, int array_len = 1);
  void push_constant(const char *name, const float4x4 *data);

  /**
   * Turn the pass into a string for inspection.
   */
  std::string serialize() const;

  friend std::ostream &operator<<(std::ostream &stream, const Pass &pass)
  {
    return stream << pass.serialize();
  }

 private:
  /**
   * Internal Helpers
   */

  int push_constant_offset(const char *name);

  void clear(eGPUFrameBufferBits planes, float4 color, float depth, uint8_t stencil);

  GPUBatch *procedural_batch_get(GPUPrimType primitive);

  static void submit_command(command::RecordingState &state,
                             const command::Type &type,
                             const command::Undetermined &command);

  static std::string serialize_command(const command::Type &type,
                                       const command::Undetermined &command);

  /**
   * Return a new command recorded with the given type.
   */
  command::Undetermined &create_command(command::Type type);

  /**
   * Submit a Pass as GPU commands. Should only be called by draw::Manager.
   */
  void submit(command::RecordingState &state) const;
};

/**
 * Main pass type.
 * Optimized for many draw calls and sub-pass.
 *
 * IMPORTANT: To be used only for passes containing lots of draw calls since it has a potentially
 * high overhead due to batching and culling optimizations.
 */
class PassMain : public Pass {
 private:
  /** Multi draw indirect rendering for many draw calls efficient rendering. */
  command::MultiDrawBuffer multi_draws_;

 public:
  PassMain() = delete;
  PassMain(const char *name) : Pass(name){};

#if 0
  void draw(GPUBatch *batch,
            uint instance_len = -1,
            uint vertex_len = -1,
            uint vertex_first = -1,
            ResourceHandle handle = {0});
  void draw(GPUBatch *batch,
            StorageBuffer<DrawCommand> &indirect_buffer,
            ResourceHandle handle = {0});
#endif
};

/** \} */

/* -------------------------------------------------------------------- */
/** \name Pass
 * \{ */

/**
 * Reset the pass command pool.
 */
inline void Pass::init()
{
  headers_.clear();
  commands_.clear();
}

inline command::Undetermined &Pass::create_command(command::Type type)
{
  int64_t index = commands_.append_and_get_index({});
  headers_.append({static_cast<uint>(index), type, false});
  return commands_[index];
}

inline void Pass::submit(command::RecordingState &state) const
{
#if 0
  if (command_last_ != -1) {
    uint cmd_index = 0;
    while (cmd_index != -1) {
      command::Header &header = headers_[cmd_index];
      cmd_index = header.next;

      if (header.is_subpass) {
        /** WARNING: Recursive. */
        PassSub &sub = *reinterpret_cast<PassSub *>(&commands_[header.command_index]);
        ss << submit();
      }
      else if (false) {
        /** TODO(fclem): MultiDraw */
      }
      else {
        ss << detail::PassCommon::submit_command(
                  state, header.type, commands_[header.command_index])
           << std::endl;
      }
    }
  }
#endif
}

inline std::string Pass::serialize() const
{
#if 0
  if (command_last_ != -1) {
    uint cmd_index = 0;
    while (cmd_index != -1) {
      command::Header &header = headers_[cmd_index];
      cmd_index = header.next;

      if (header.is_subpass) {
        /** WARNING: Recursive. */
        PassSub &sub = *reinterpret_cast<PassSub *>(&commands_[header.command_index]);
        ss << sub.serialize();
      }
      else if (false) {
        /** TODO(fclem): MultiDraw */
      }
      else {
        ss << detail::PassCommon::serialize_command(header.type, commands_[header.command_index])
           << std::endl;
      }
    }
  }
#endif
  return "";
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Draw Calls
 * \{ */

inline void Pass::clear(eGPUFrameBufferBits planes, float4 color, float depth, uint8_t stencil)
{
  create_command(command::Type::Clear).clear = {(uint8_t)planes, stencil, depth, color};
}

inline GPUBatch *Pass::procedural_batch_get(GPUPrimType primitive)
{
  switch (primitive) {
    case GPU_PRIM_POINTS:
      return drw_cache_procedural_points_get();
    case GPU_PRIM_LINES:
      return drw_cache_procedural_lines_get();
    case GPU_PRIM_TRIS:
      return drw_cache_procedural_triangles_get();
    case GPU_PRIM_TRI_STRIP:
      return drw_cache_procedural_triangle_strips_get();
    default:
      /* Add new one as needed. */
      BLI_assert_unreachable();
      return nullptr;
  }
}

inline void Pass::draw(
    GPUBatch *batch, uint instance_len, uint vertex_len, uint vertex_first, ResourceHandle handle)
{
  if (instance_len == 0 || vertex_len == 0) {
    return;
  }
  BLI_assert(shader_);
  create_command(Type::Draw).draw = {batch, instance_len, vertex_len, vertex_first, handle};
}

inline void Pass::draw(GPUBatch *batch,
                       StorageBuffer<DrawCommand> &indirect_buffer,
                       ResourceHandle handle)
{
  BLI_assert(shader_);
  create_command(Type::DrawIndirect).draw_indirect = {batch, &indirect_buffer, handle};
}

inline void Pass::draw(GPUPrimType primitive,
                       uint instance_len,
                       uint vertex_len,
                       uint vertex_first,
                       ResourceHandle handle)
{
  draw(procedural_batch_get(primitive), instance_len, vertex_len, vertex_first, handle);
}

inline void Pass::draw(GPUPrimType primitive,
                       StorageBuffer<DrawCommand> &indirect_buffer,
                       ResourceHandle handle)
{
  draw(procedural_batch_get(primitive), indirect_buffer, handle);
}

inline void Pass::draw(GPUBatch *batch, ResourceHandle handle)
{
  draw(batch, -1, -1, -1, handle);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Compute Dispatch
 * \{ */

inline void Pass::dispatch(int3 group_len)
{
  BLI_assert(shader_);
  create_command(Type::Dispatch).dispatch = {group_len};
}

inline void Pass::dispatch(int3 *group_len)
{
  BLI_assert(shader_);
  create_command(Type::Dispatch).dispatch = {group_len};
}

inline void Pass::dispatch(StorageBuffer<DispatchCommand> &indirect_buffer)
{
  BLI_assert(shader_);
  create_command(Type::DispatchIndirect).dispatch_indirect = {&indirect_buffer};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Clear
 * \{ */

inline void Pass::clear_color(float4 color)
{
  clear(GPU_COLOR_BIT, color, 0.0f, 0);
}

inline void Pass::clear_depth(float depth)
{
  clear(GPU_DEPTH_BIT, float4(0.0f), depth, 0);
}

inline void Pass::clear_stencil(uint8_t stencil)
{
  clear(GPU_STENCIL_BIT, float4(0.0f), 0.0f, stencil);
}

inline void Pass::clear_depth_stencil(float depth, uint8_t stencil)
{
  clear(GPU_DEPTH_BIT | GPU_STENCIL_BIT, float4(0.0f), depth, stencil);
}

inline void Pass::clear_color_depth_stencil(float4 color, float depth, uint8_t stencil)
{
  clear(GPU_DEPTH_BIT | GPU_STENCIL_BIT | GPU_COLOR_BIT, color, depth, stencil);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Barrier
 * \{ */

inline void Pass::barrier(eGPUBarrier type)
{
  create_command(Type::Barrier).barrier = {type};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name State
 * \{ */

inline void Pass::state_set(DRWState state)
{
  create_command(Type::StateSet).state_set = {state};
}

inline void Pass::state_stencil(uint8_t write_mask, uint8_t reference, uint8_t compare_mask)
{
  create_command(Type::StencilSet).stencil_set = {write_mask, reference, compare_mask};
}

inline void Pass::shader_set(GPUShader *shader)
{
  shader_ = shader;
  create_command(Type::ShaderBind).shader_bind = {shader};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Resource bind
 * \{ */

inline int Pass::push_constant_offset(const char *name)
{
  return GPU_shader_get_uniform(shader_, name);
}

inline void Pass::bind(const char *name, GPUStorageBuf *buffer)
{
  bind(GPU_shader_get_uniform_block_binding(shader_, name), buffer);
}

inline void Pass::bind(const char *name, GPUUniformBuf *buffer)
{
  bind(GPU_shader_get_ssbo(shader_, name), buffer);
}

inline void Pass::bind(const char *name, GPUTexture *texture, eGPUSamplerState state)
{
  bind(GPU_shader_get_texture_binding(shader_, name), texture, state);
}

inline void Pass::bind(const char *name, draw::Image *image)
{
  bind(GPU_shader_get_texture_binding(shader_, name), image);
}

inline void Pass::bind(int slot, GPUStorageBuf *buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

inline void Pass::bind(int slot, GPUUniformBuf *buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

inline void Pass::bind(int slot, GPUTexture *texture, eGPUSamplerState state)
{
  create_command(Type::ResourceBind).resource_bind = {slot, texture, state};
}

inline void Pass::bind(int slot, draw::Image *image)
{
  create_command(Type::ResourceBind).resource_bind = {slot, image};
}

inline void Pass::bind(const char *name, GPUStorageBuf **buffer)
{
  bind(GPU_shader_get_uniform_block_binding(shader_, name), buffer);
}

inline void Pass::bind(const char *name, GPUUniformBuf **buffer)
{
  bind(GPU_shader_get_ssbo(shader_, name), buffer);
}

inline void Pass::bind(const char *name, GPUTexture **texture, eGPUSamplerState state)
{
  bind(GPU_shader_get_texture_binding(shader_, name), texture, state);
}

inline void Pass::bind(const char *name, draw::Image **image)
{
  bind(GPU_shader_get_texture_binding(shader_, name), image);
}

inline void Pass::bind(int slot, GPUStorageBuf **buffer)
{

  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

inline void Pass::bind(int slot, GPUUniformBuf **buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

inline void Pass::bind(int slot, GPUTexture **texture, eGPUSamplerState state)
{
  create_command(Type::ResourceBind).resource_bind = {slot, texture, state};
}

inline void Pass::bind(int slot, draw::Image **image)
{
  create_command(Type::ResourceBind).resource_bind = {slot, image};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Push Constant
 * \{ */

inline void Pass::push_constant(const char *name, const float &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const float2 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const float3 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const float4 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const int &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const int2 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const int3 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const int4 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const bool &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const float *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const float2 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const float3 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const float4 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const int *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const int2 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const int3 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const int4 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

inline void Pass::push_constant(const char *name, const float4x4 *data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

inline void Pass::push_constant(const char *name, const float4x4 &data)
{
  /* WORKAROUND: Push 3 consecutive commands to hold the 64 bytes of the float4x4.
   * This assumes that all commands are always stored in flat array of memory. */
  Undetermined commands[3];

  PushConstant &cmd = commands[0].push_constant;
  cmd.location = push_constant_offset(name);
  cmd.array_len = 1;
  cmd.comp_len = 16;
  cmd.type = PushConstant::Type::FloatValue;
  /* Copy overrides the next 2 commands. We append them as Type::None to not evaluate them. */
  *reinterpret_cast<float4x4 *>(&cmd.float4_value) = data;

  create_command(Type::PushConstant) = commands[0];
  create_command(Type::None) = commands[1];
  create_command(Type::None) = commands[2];
}

inline void Pass::submit_command(command::RecordingState &state,
                                 const command::Type &type,
                                 const command::Undetermined &command)
{
  switch (type) {
    case command::Type::ShaderBind:
      command.shader_bind.execute(state);
      break;
    case command::Type::ResourceBind:
      command.resource_bind.execute();
      break;
    case command::Type::PushConstant:
      command.push_constant.execute(state);
      break;
    case command::Type::Draw:
      command.draw.execute(state);
      break;
    case command::Type::DrawIndirect:
      command.draw_indirect.execute(state);
      break;
    case command::Type::Dispatch:
      command.dispatch.execute(state);
      break;
    case command::Type::DispatchIndirect:
      command.dispatch_indirect.execute(state);
      break;
    case command::Type::Barrier:
      command.barrier.execute();
      break;
    case command::Type::Clear:
      command.clear.execute();
      break;
    case command::Type::StateSet:
      command.state_set.execute(state);
      break;
    case command::Type::StencilSet:
      command.stencil_set.execute();
      break;
    default:
      /* Skip Type::None. */
      break;
  }
}

inline std::string Pass::serialize_command(const command::Type &type,
                                           const command::Undetermined &command)
{
  switch (type) {
    case command::Type::ShaderBind:
      return command.shader_bind.serialize();
    case command::Type::ResourceBind:
      return command.resource_bind.serialize();
    case command::Type::PushConstant:
      return command.push_constant.serialize();
    case command::Type::Draw:
      return command.draw.serialize();
    case command::Type::DrawIndirect:
      return command.draw_indirect.serialize();
    case command::Type::Dispatch:
      return command.dispatch.serialize();
    case command::Type::DispatchIndirect:
      return command.dispatch_indirect.serialize();
    case command::Type::Barrier:
      return command.barrier.serialize();
    case command::Type::Clear:
      return command.clear.serialize();
    case command::Type::StateSet:
      return command.state_set.serialize();
    case command::Type::StencilSet:
      return command.stencil_set.serialize();
    default:
      /* Skip Type::None. */
      return "";
  }
}

/** \} */

}  // namespace blender::draw
