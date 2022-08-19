/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 Blender Foundation. */

/** \file
 * \ingroup draw
 */

#include "draw_pass.hh"

namespace blender::draw {

/* -------------------------------------------------------------------- */
/** \name Pass
 * \{ */

#if 0 /* TODO */
void Pass::init()
{
  command_last_ = -1;
}

PassSub &Pass::sub()
{
  /* WORKAROUND: Push 2 consecutive commands to hold the PassSub in the command array.
   * This assumes that all commands are always stored in flat array of memory. */
  BLI_STATIC_ASSERT(sizeof(PassSub) <= sizeof(command::Undetermined) * 2, "PassSub is too large");
}

void Pass::submit(command::RecordingState &state) const
{
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
}

std::string Pass::serialize() const
{
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
}

command::Undetermined &Pass::create_command(command::Type type)
{
  command_types_.append(type);
  int64_t index = commands_.append_and_get_index({});
  return commands_[index];
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Sub Pass
 * \{ */

PassSub &PassSub::sub()
{
  return PassSub();
}

void PassSub::submit(command::RecordingState &state) const
{
}

command::Undetermined &PassSub::create_command(command::Type type)
{
}
#endif

/** \} */

}  // namespace blender::draw

namespace blender::draw::detail {

using namespace blender::draw::command;

/* -------------------------------------------------------------------- */
/** \name Draw Calls
 * \{ */

void PassInterface::draw(
    GPUBatch *batch, uint instance_len, uint vertex_len, uint vertex_first, ResourceHandle handle)
{
  if (instance_len == 0 || vertex_len == 0) {
    return;
  }
  BLI_assert(shader_);
  create_command(Type::Draw).draw = {batch, instance_len, vertex_len, vertex_first, handle};
}

void PassInterface::draw(GPUBatch *batch,
                         StorageBuffer<DrawCommand> &indirect_buffer,
                         ResourceHandle handle)
{
  BLI_assert(shader_);
  create_command(Type::DrawIndirect).draw_indirect = {batch, &indirect_buffer, handle};
}

void PassInterface::draw(GPUPrimType primitive,
                         uint instance_len,
                         uint vertex_len,
                         uint vertex_first,
                         ResourceHandle handle)
{
  draw(procedural_batch_get(primitive), instance_len, vertex_len, vertex_first, handle);
}

void PassInterface::draw(GPUPrimType primitive,
                         StorageBuffer<DrawCommand> &indirect_buffer,
                         ResourceHandle handle)
{
  draw(procedural_batch_get(primitive), indirect_buffer, handle);
}

void PassInterface::draw(GPUBatch *batch, ResourceHandle handle)
{
  draw(batch, -1, -1, -1, handle);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Compute Dispatch
 * \{ */

void PassInterface::dispatch(int3 group_len)
{
  BLI_assert(shader_);
  create_command(Type::Dispatch).dispatch = {group_len};
}

void PassInterface::dispatch(int3 *group_len)
{
  BLI_assert(shader_);
  create_command(Type::Dispatch).dispatch = {group_len};
}

void PassInterface::dispatch(StorageBuffer<DispatchCommand> &indirect_buffer)
{
  BLI_assert(shader_);
  create_command(Type::DispatchIndirect).dispatch_indirect = {&indirect_buffer};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Clear
 * \{ */

void PassInterface::clear_color(float4 color)
{
  clear(GPU_COLOR_BIT, color, 0.0f, 0);
}

void PassInterface::clear_depth(float depth)
{
  clear(GPU_DEPTH_BIT, float4(0.0f), depth, 0);
}

void PassInterface::clear_stencil(uint8_t stencil)
{
  clear(GPU_STENCIL_BIT, float4(0.0f), 0.0f, stencil);
}

void PassInterface::clear_depth_stencil(float depth, uint8_t stencil)
{
  clear(GPU_DEPTH_BIT | GPU_STENCIL_BIT, float4(0.0f), depth, stencil);
}

void PassInterface::clear_color_depth_stencil(float4 color, float depth, uint8_t stencil)
{
  clear(GPU_DEPTH_BIT | GPU_STENCIL_BIT | GPU_COLOR_BIT, color, depth, stencil);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Barrier
 * \{ */

void PassInterface::barrier(eGPUBarrier type)
{
  create_command(Type::Barrier).barrier = {type};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name State
 * \{ */

void PassInterface::state_set(DRWState state)
{
  create_command(Type::StateSet).state_set = {state};
}

void PassInterface::state_stencil(uint8_t write_mask, uint8_t reference, uint8_t compare_mask)
{
  create_command(Type::StencilSet).stencil_set = {write_mask, reference, compare_mask};
}

void PassInterface::shader_set(GPUShader *shader)
{
  shader_ = shader;
  create_command(Type::ShaderBind).shader_bind = {shader};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Resource bind
 * \{ */

void PassInterface::bind(const char *name, GPUStorageBuf *buffer)
{
  bind(GPU_shader_get_uniform_block_binding(shader_, name), buffer);
}

void PassInterface::bind(const char *name, GPUUniformBuf *buffer)
{
  bind(GPU_shader_get_ssbo(shader_, name), buffer);
}

void PassInterface::bind(const char *name, GPUTexture *texture, eGPUSamplerState state)
{
  bind(GPU_shader_get_texture_binding(shader_, name), texture, state);
}

void PassInterface::bind(const char *name, draw::Image *image)
{
  bind(GPU_shader_get_texture_binding(shader_, name), image);
}

void PassInterface::bind(int slot, GPUStorageBuf *buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

void PassInterface::bind(int slot, GPUUniformBuf *buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

void PassInterface::bind(int slot, GPUTexture *texture, eGPUSamplerState state)
{
  create_command(Type::ResourceBind).resource_bind = {slot, texture, state};
}

void PassInterface::bind(int slot, draw::Image *image)
{
  create_command(Type::ResourceBind).resource_bind = {slot, image};
}

void PassInterface::bind(const char *name, GPUStorageBuf **buffer)
{
  bind(GPU_shader_get_uniform_block_binding(shader_, name), buffer);
}

void PassInterface::bind(const char *name, GPUUniformBuf **buffer)
{
  bind(GPU_shader_get_ssbo(shader_, name), buffer);
}

void PassInterface::bind(const char *name, GPUTexture **texture, eGPUSamplerState state)
{
  bind(GPU_shader_get_texture_binding(shader_, name), texture, state);
}

void PassInterface::bind(const char *name, draw::Image **image)
{
  bind(GPU_shader_get_texture_binding(shader_, name), image);
}

void PassInterface::bind(int slot, GPUStorageBuf **buffer)
{

  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

void PassInterface::bind(int slot, GPUUniformBuf **buffer)
{
  create_command(Type::ResourceBind).resource_bind = {slot, buffer};
}

void PassInterface::bind(int slot, GPUTexture **texture, eGPUSamplerState state)
{
  create_command(Type::ResourceBind).resource_bind = {slot, texture, state};
}

void PassInterface::bind(int slot, draw::Image **image)
{
  create_command(Type::ResourceBind).resource_bind = {slot, image};
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Push Constant
 * \{ */

void PassInterface::push_constant(const char *name, const float &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const float2 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const float3 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const float4 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const int &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const int2 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const int3 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const int4 &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const bool &data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const float *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const float2 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const float3 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const float4 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const int *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const int2 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const int3 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const int4 *data, int array_len)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data, array_len};
}

void PassInterface::push_constant(const char *name, const float4x4 *data)
{
  create_command(Type::PushConstant).push_constant = {push_constant_offset(name), data};
}

void PassInterface::push_constant(const char *name, const float4x4 &data)
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

/** \} */

}  // namespace blender::draw::detail
