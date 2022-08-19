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
 * Example (P for Parent, s for Sub):
 * This recording order [P, s1, P, s2, s1, P, s2]
 * will yield
 * this execution order [P, s1, s1, P, s2, s2, P]
 *
 * `PassSimple`:
 * Drawing is a bit easier for the CPU but there is no batching or culling optimization. Has lower
 * memory footprint and CPU usage than #Pass for simple workload. Cannot contain sub passes.
 *
 * `PassImmutable`:
 * Commands are packed into a static array. The size must be declared upfront. To be used with
 * small static passes. Has the lowest memory footprint and CPU usage. Cannot contain sub passes.
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

class Manager;

/* -------------------------------------------------------------------- */
/** \name Pass API
 * \{ */

namespace detail {

using namespace blender::draw;

/**
 * Public API of a draw pass.
 */
class PassInterface {
 public:
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

 private:
  /** Currently bound shader. Used for interface queries. */
  GPUShader *shader_;

  /**
   * Helpers
   */

  int push_constant_offset(const char *name)
  {
    return GPU_shader_get_uniform(shader_, name);
  }

  void clear(eGPUFrameBufferBits planes, float4 color, float depth, uint8_t stencil)
  {
    create_command(command::Type::Clear).clear = {(uint8_t)planes, stencil, depth, color};
  }

  GPUBatch *procedural_batch_get(GPUPrimType primitive)
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

  virtual command::Undetermined &create_command(command::Type type) = 0;
};

class PassCommon : public PassInterface {
  friend Manager;

 public:
  const char *debug_name = "Unamed";

  PassCommon(const char *name) : debug_name(name){};

 protected:
  static void submit_command(command::RecordingState &state,
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

  static std::string serialize_command(const command::Type &type,
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
};

}  // namespace detail

/** \} */

/* -------------------------------------------------------------------- */
/** \name Passes
 * \{ */

class PassSub;

/**
 * Main pass type.
 * Optimized for many draw calls and sub-pass.
 *
 * IMPORTANT: To be used only for passes containing lots of draw calls since it has a potentially
 * high overhead due to batching and culling optimizations.
 */
class Pass : public detail::PassCommon {
 private:
  /** Highest level of the command stream. Split command stream in different command types. */
  Vector<command::Header> headers_;
  /** Commands referenced by groups (which contains their types). */
  Vector<command::Undetermined> commands_;
  /** Multi draw indirect rendering for many draw calls efficient rendering. */
  command::MultiDrawBuffer multi_draws_;
  /** Last command index for easy appending of commands. */
  uint command_last_;

 public:
  Pass() = delete;
  Pass(const char *name) : detail::PassCommon(name){};

  void init();

  /**
   * Create a sub-pass inside this pass. The sub-pass memory is auto managed.
   */
  PassSub &sub();

  /**
   * Turn the pass into a string for inspection.
   */
  std::string serialize() const;

  friend std::ostream &operator<<(std::ostream &stream, const Pass &pass)
  {
    return stream << pass.serialize();
  }

 private:
  command::Undetermined &create_command(command::Type type) override final;
  void submit(command::RecordingState &state) const;
};

/**
 * A pass that is linked to its parent pass and gets its storage from it.
 * Trivially constructible and destructible since it does not contain actual data.
 *
 * IMPORTANT: Sub pass commands will be inserted between parent pass commands where
 * the sub-pass was created. This can change the main pass state (like bound shader)
 * for the subsequent commands.
 */
class PassSub : public detail::PassCommon {
  friend Pass;

 private:
  Pass &parent_;
  uint command_first_ = -1;
  uint command_last_ = -1;

 public:
  /**
   * Create a sub-pass inside this sub-pass. The sub-pass memory is auto managed.
   */
  PassSub &sub();

  /**
   * Turn the pass into a string for inspection.
   */
  std::string serialize() const;

  friend std::ostream &operator<<(std::ostream &stream, const PassSub &pass)
  {
    return stream << pass.serialize();
  }

 private:
  /** Private constructor since a sub pass can only be created using the sub() function. */
  PassSub(const char *name, Pass &parent) : detail::PassCommon(name), parent_(parent){};

  /**
   * Return a new command recorded with the given type.
   */
  command::Undetermined &create_command(command::Type type) override final;
};

/**
 * A pass that cannot have sub passes. All commands are in sequential order.
 * Drawing is a bit easier for the CPU but there is no batching or culling optimization.
 */
class PassSimple : public detail::PassCommon {
 private:
  /** Highest level of the command stream. Split command stream in different command types. */
  Vector<command::Type> command_types_;
  /** Commands referenced by groups (which contains their types). */
  Vector<command::Undetermined> commands_;

 public:
  PassSimple() = delete;
  PassSimple(const char *name) : detail::PassCommon(name){};

  /**
   * Reset the pass command pool.
   */
  void init()
  {
    command_types_.clear();
    commands_.clear();
  }

  /**
   * Turn the pass into a string for inspection.
   */
  std::string serialize() const
  {
    std::stringstream ss;
    ss << "PassSimple(" << debug_name << ", len=" << commands_.size() << ")" << std::endl;
    for (auto i : command_types_.index_range()) {
      ss << detail::PassCommon::serialize_command(command_types_[i], commands_[i]) << std::endl;
    }
    return ss.str();
  }

  friend std::ostream &operator<<(std::ostream &stream, const PassSimple &pass)
  {
    return stream << pass.serialize();
  }

 private:
  /**
   * Return a new command recorded with the given type.
   */
  command::Undetermined &create_command(command::Type type) override final
  {
    command_types_.append(type);
    int64_t index = commands_.append_and_get_index({});
    return commands_[index];
  }

  /**
   * Submit a Pass as GPU commands. Should only be called by draw::Manager.
   */
  void submit(command::RecordingState &state) const
  {
    for (auto i : IndexRange(command_types_.size())) {
      detail::PassCommon::submit_command(state, command_types_[i], commands_[i]);
    }
  }
};

/**
 * A pass that needs to have its size specified upfront and is immutable afterwards.
 * Drawing is a bit easier for the CPU but there is no batching or culling optimization.
 * Has the lowest memory footprint.
 */
template<int64_t len> class PassImmutable : public detail::PassCommon {
 private:
  std::array<command::Type, len> command_types_;
  std::array<command::Undetermined, len> commands_;

 public:
  PassImmutable() = delete;
  PassImmutable(const char *name) : detail::PassCommon(name){};

  /**
   * Reset the pass command pool.
   */
  void init()
  {
    memset(command_types_, 0, sizeof(command_types_));
  }

  /**
   * Turn the pass into a string for inspection.
   */
  std::string serialize() const
  {
    std::stringstream ss;
    ss << "PassImmutable(" << debug_name << ", len=" << len << ")" << std::endl;
    for (auto i : IndexRange(command_types_.size())) {
      ss << detail::PassCommon::serialize_command(command_types_[i], commands_[i]) << std::endl;
    }
    return ss.str();
  }

  friend std::ostream &operator<<(std::ostream &stream, const PassImmutable &pass)
  {
    return stream << pass.serialize();
  }

 private:
  /**
   * Return a new command recorded with the given type.
   */
  command::Undetermined &create_command(command::Type type) override final
  {
    for (auto i : IndexRange(command_types_.size())) {
      if (command_types_[i] == command::Type::None) {
        command_types_[i] = type;
        return commands_[i];
      }
    }
    BLI_assert_unreachable();
    /* Return any command to avoid compiler error. But this will override another command and lead
     * to undefined behavior. */
    return commands_.back();
  }

  /**
   * Submit a Pass as GPU commands. Should only be called by draw::Manager.
   */
  void submit(command::RecordingState &state) const
  {
    for (auto i : IndexRange(command_types_.size())) {
      detail::PassCommon::submit_command(state, command_types_[i], commands_[i++]);
    }
  }
};

/** \} */

}  // namespace blender::draw
