/* SPDX-License-Identifier: Apache-2.0
 * Copyright 2011-2022 Blender Foundation */

#include "integrator/path_trace_work_cpu.h"

#include "device/cpu/kernel.h"
#include "device/device.h"

#include "kernel/film/write_passes.h"
#include "kernel/integrator/path_state.h"

#include "integrator/pass_accessor_cpu.h"
#include "integrator/path_trace_display.h"

#include "scene/scene.h"
#include "session/buffers.h"

#include "util/atomic.h"
#include "util/log.h"
#include "util/tbb.h"

CCL_NAMESPACE_BEGIN

/* Create TBB arena for execution of path tracing and rendering tasks. */
static inline tbb::task_arena local_tbb_arena_create(const Device *device)
{
  /* TODO: limit this to number of threads of CPU device, it may be smaller than
   * the system number of threads when we reduce the number of CPU threads in
   * CPU + GPU rendering to dedicate some cores to handling the GPU device. */
  return tbb::task_arena(device->info.cpu_threads);
}

/* Get CPUKernelThreadGlobals for the current thread. */
static inline CPUKernelThreadGlobals *kernel_thread_globals_get(
    vector<CPUKernelThreadGlobals> &kernel_thread_globals)
{
  const int thread_index = tbb::this_task_arena::current_thread_index();
  DCHECK_GE(thread_index, 0);
  DCHECK_LE(thread_index, kernel_thread_globals.size());

  return &kernel_thread_globals[thread_index];
}

PathTraceWorkCPU::PathTraceWorkCPU(Device *device,
                                   Film *film,
                                   DeviceScene *device_scene,
                                   bool *cancel_requested_flag)
    : PathTraceWork(device, film, device_scene, cancel_requested_flag),
      kernels_(Device::get_cpu_kernels())
{
  DCHECK_EQ(device->info.type, DEVICE_CPU);
}

void PathTraceWorkCPU::init_execution()
{
  /* Cache per-thread kernel globals. */
  device_->get_cpu_kernel_thread_globals(kernel_thread_globals_);
}

#if defined(WITH_PATH_GUIDING)
/* Note: It seems that this is called before every rendering iteration/progression and not once per
 * rendering. May be we find a way to call it only once per renderering. */
void PathTraceWorkCPU::guiding_init_kernel_globals(void *guiding_field, void *sample_data_storage)
{
  /* Cache per-thread kernel globals. */
  // device_->get_cpu_kernel_thread_globals(kernel_thread_globals_);
#  if PATH_GUIDING_LEVEL >= 4
  openpgl::cpp::Field *field = (openpgl::cpp::Field *)guiding_field;
#  endif
  for (int thread_index = 0; thread_index < kernel_thread_globals_.size(); thread_index++) {
#  if PATH_GUIDING_LEVEL >= 3

    kernel_thread_globals_[thread_index].opgl_sample_data_storage = (openpgl::cpp::SampleStorage *)
        sample_data_storage;
#  endif

#  if PATH_GUIDING_LEVEL >= 4
    kernel_thread_globals_[thread_index].opgl_guiding_field = field;
    /* We might be able to just create a new Sampling distribution if the pointer is NULL */
    if (kernel_thread_globals_[thread_index].opgl_surface_sampling_distribution)
      delete kernel_thread_globals_[thread_index].opgl_surface_sampling_distribution;
    kernel_thread_globals_[thread_index].opgl_surface_sampling_distribution =
        new openpgl::cpp::SurfaceSamplingDistribution(field);
    if (kernel_thread_globals_[thread_index].opgl_volume_sampling_distribution)
      delete kernel_thread_globals_[thread_index].opgl_volume_sampling_distribution;
    kernel_thread_globals_[thread_index].opgl_volume_sampling_distribution =
        new openpgl::cpp::VolumeSamplingDistribution(field);
#  endif
  }
}
#endif

void PathTraceWorkCPU::render_samples(RenderStatistics &statistics,
                                      int start_sample,
                                      int samples_num,
                                      int sample_offset)
{
#if defined(WITH_PATH_GUIDING) && defined(WITH_PATH_GUIDING_DEBUG_PRINT)
  VLOG_WORK << "render_samples: start_sample = " << start_sample
            << "\t samples_num = " << samples_num;
#endif
  const int64_t image_width = effective_buffer_params_.width;
  const int64_t image_height = effective_buffer_params_.height;
  const int64_t total_pixels_num = image_width * image_height;

  if (device_->profiler.active()) {
    for (CPUKernelThreadGlobals &kernel_globals : kernel_thread_globals_) {
      kernel_globals.start_profiling();
    }
  }

  tbb::task_arena local_arena = local_tbb_arena_create(device_);
  local_arena.execute([&]() {
    parallel_for(int64_t(0), total_pixels_num, [&](int64_t work_index) {
      if (is_cancel_requested()) {
        return;
      }

      const int y = work_index / image_width;
      const int x = work_index - y * image_width;

      KernelWorkTile work_tile;
      work_tile.x = effective_buffer_params_.full_x + x;
      work_tile.y = effective_buffer_params_.full_y + y;
      work_tile.w = 1;
      work_tile.h = 1;
      work_tile.start_sample = start_sample;
      work_tile.sample_offset = sample_offset;
      work_tile.num_samples = 1;
      work_tile.offset = effective_buffer_params_.offset;
      work_tile.stride = effective_buffer_params_.stride;

      CPUKernelThreadGlobals *kernel_globals = kernel_thread_globals_get(kernel_thread_globals_);

      render_samples_full_pipeline(kernel_globals, work_tile, samples_num);
    });
  });
  if (device_->profiler.active()) {
    for (CPUKernelThreadGlobals &kernel_globals : kernel_thread_globals_) {
      kernel_globals.stop_profiling();
    }
  }

  statistics.occupancy = 1.0f;
}

void PathTraceWorkCPU::render_samples_full_pipeline(KernelGlobalsCPU *kg,
                                                    const KernelWorkTile &work_tile,
                                                    const int samples_num)
{
  const bool has_bake = device_scene_->data.bake.use;

  IntegratorStateCPU integrator_states[2];

  IntegratorStateCPU *state = &integrator_states[0];
  IntegratorStateCPU *shadow_catcher_state = nullptr;

  if (device_scene_->data.integrator.has_shadow_catcher) {
    shadow_catcher_state = &integrator_states[1];
    path_state_init_queues(shadow_catcher_state);
  }

  KernelWorkTile sample_work_tile = work_tile;
  float *render_buffer = buffers_->buffer.data();

  for (int sample = 0; sample < samples_num; ++sample) {
    if (is_cancel_requested()) {
      break;
    }

    if (has_bake) {
      if (!kernels_.integrator_init_from_bake(kg, state, &sample_work_tile, render_buffer)) {
        break;
      }
    }
    else {
      if (!kernels_.integrator_init_from_camera(kg, state, &sample_work_tile, render_buffer)) {
        break;
      }
    }

#ifdef WITH_PATH_GUIDING
    const bool use_guiding = kernel_data.integrator.use_guiding;
    if (use_guiding) {
      /* Clear path segment storage. */
      guiding_prepare_integrator_state(kg, state);
    }
#endif

    kernels_.integrator_megakernel(kg, state, render_buffer);

#if defined(WITH_PATH_GUIDING) && PATH_GUIDING_LEVEL >= 1
    if (use_guiding) {
      /* Push the generated sample data to the global sample data storage. */
      guiding_push_sample_data_to_global_storage(kg, state, render_buffer);
    }
#endif
    if (shadow_catcher_state) {
      kernels_.integrator_megakernel(kg, shadow_catcher_state, render_buffer);
    }

    ++sample_work_tile.start_sample;
  }
}

void PathTraceWorkCPU::copy_to_display(PathTraceDisplay *display,
                                       PassMode pass_mode,
                                       int num_samples)
{
  half4 *rgba_half = display->map_texture_buffer();
  if (!rgba_half) {
    /* TODO(sergey): Look into using copy_to_display() if mapping failed. Might be needed for
     * some implementations of PathTraceDisplay which can not map memory? */
    return;
  }

  const KernelFilm &kfilm = device_scene_->data.film;

  const PassAccessor::PassAccessInfo pass_access_info = get_display_pass_access_info(pass_mode);

  const PassAccessorCPU pass_accessor(pass_access_info, kfilm.exposure, num_samples);

  PassAccessor::Destination destination = get_display_destination_template(display);
  destination.pixels_half_rgba = rgba_half;

  tbb::task_arena local_arena = local_tbb_arena_create(device_);
  local_arena.execute([&]() {
    pass_accessor.get_render_tile_pixels(buffers_.get(), effective_buffer_params_, destination);
  });

  display->unmap_texture_buffer();
}

void PathTraceWorkCPU::destroy_gpu_resources(PathTraceDisplay * /*display*/)
{
}

bool PathTraceWorkCPU::copy_render_buffers_from_device()
{
  return buffers_->copy_from_device();
}

bool PathTraceWorkCPU::copy_render_buffers_to_device()
{
  buffers_->buffer.copy_to_device();
  return true;
}

bool PathTraceWorkCPU::zero_render_buffers()
{
  buffers_->zero();
  return true;
}

int PathTraceWorkCPU::adaptive_sampling_converge_filter_count_active(float threshold, bool reset)
{
  const int full_x = effective_buffer_params_.full_x;
  const int full_y = effective_buffer_params_.full_y;
  const int width = effective_buffer_params_.width;
  const int height = effective_buffer_params_.height;
  const int offset = effective_buffer_params_.offset;
  const int stride = effective_buffer_params_.stride;

  float *render_buffer = buffers_->buffer.data();

  uint num_active_pixels = 0;

  tbb::task_arena local_arena = local_tbb_arena_create(device_);

  /* Check convergency and do x-filter in a single `parallel_for`, to reduce threading overhead. */
  local_arena.execute([&]() {
    parallel_for(full_y, full_y + height, [&](int y) {
      CPUKernelThreadGlobals *kernel_globals = &kernel_thread_globals_[0];

      bool row_converged = true;
      uint num_row_pixels_active = 0;
      for (int x = 0; x < width; ++x) {
        if (!kernels_.adaptive_sampling_convergence_check(
                kernel_globals, render_buffer, full_x + x, y, threshold, reset, offset, stride)) {
          ++num_row_pixels_active;
          row_converged = false;
        }
      }

      atomic_fetch_and_add_uint32(&num_active_pixels, num_row_pixels_active);

      if (!row_converged) {
        kernels_.adaptive_sampling_filter_x(
            kernel_globals, render_buffer, y, full_x, width, offset, stride);
      }
    });
  });

  if (num_active_pixels) {
    local_arena.execute([&]() {
      parallel_for(full_x, full_x + width, [&](int x) {
        CPUKernelThreadGlobals *kernel_globals = &kernel_thread_globals_[0];
        kernels_.adaptive_sampling_filter_y(
            kernel_globals, render_buffer, x, full_y, height, offset, stride);
      });
    });
  }

  return num_active_pixels;
}

void PathTraceWorkCPU::cryptomatte_postproces()
{
  const int width = effective_buffer_params_.width;
  const int height = effective_buffer_params_.height;

  float *render_buffer = buffers_->buffer.data();

  tbb::task_arena local_arena = local_tbb_arena_create(device_);

  /* Check convergency and do x-filter in a single `parallel_for`, to reduce threading overhead. */
  local_arena.execute([&]() {
    parallel_for(0, height, [&](int y) {
      CPUKernelThreadGlobals *kernel_globals = &kernel_thread_globals_[0];
      int pixel_index = y * width;

      for (int x = 0; x < width; ++x, ++pixel_index) {
        kernels_.cryptomatte_postprocess(kernel_globals, render_buffer, pixel_index);
      }
    });
  });
}

#if defined(WITH_PATH_GUIDING)
void PathTraceWorkCPU::guiding_prepare_integrator_state(KernelGlobalsCPU *kg,
                                                        IntegratorStateCPU *state)
{
#  if PATH_GUIDING_LEVEL >= 1
  /* Linking to the thread local instances of the path segment storage and the guided sampling
   * distrubtion */
  state->guiding.path_segment_storage = kg->opgl_path_segment_storage;
  // TODO: find a better place to do this (once per rendering)
  state->guiding.path_segment_storage->Reserve(kernel_data.integrator.transparent_max_bounce +
                                               kernel_data.integrator.max_bounce + 3);
  state->guiding.path_segment_storage->Clear();
  /* Resetting the pointers to current path segment */
  state->guiding.path_segment = nullptr;
  state->ao.shadow_path.path_segment = nullptr;
  state->shadow.shadow_path.path_segment = nullptr;
#  endif

#  if PATH_GUIDING_LEVEL >= 4
  state->guiding.surface_sampling_distribution = kg->opgl_surface_sampling_distribution;
  state->guiding.volume_sampling_distribution = kg->opgl_volume_sampling_distribution;
#  endif
}

#  if PATH_GUIDING_LEVEL >= 1
void PathTraceWorkCPU::guiding_push_sample_data_to_global_storage(
    KernelGlobalsCPU *kg, IntegratorStateCPU *state, ccl_global float *ccl_restrict render_buffer)
{

#    if defined(PATH_GUIDING_DEBUG_VALIDATE) && PATH_GUIDING_LEVEL >= 5
  /* we can guard this check with a debug flag/define */
  bool validSegments = state->guiding.path_segment_storage->ValidateSegments();
  if (!validSegments)
    std::cout << "!!! validSegments = " << validSegments << " !!!" << std::endl;
#    endif

#    if defined(WITH_CYCLES_DEBUG) && PATH_GUIDING_LEVEL >= 5
  // guiding_write_pixel_estimate_buffer(kg, state, render_buffer);

  pgl_vec3f pgl_final_color = state->guiding.path_segment_storage->CalculatePixelEstimate(false);
  const uint32_t render_pixel_index = INTEGRATOR_STATE(state, path, render_pixel_index);
  const uint64_t render_buffer_offset = (uint64_t)render_pixel_index *
                                        kernel_data.film.pass_stride;
  ccl_global float *buffer = render_buffer + render_buffer_offset;
  float3 final_color = make_float3(pgl_final_color.x, pgl_final_color.y, pgl_final_color.z);
  if (kernel_data.film.pass_opgl_color != PASS_UNUSED) {
    kernel_write_pass_float3(buffer + kernel_data.film.pass_opgl_color, final_color);
  }
#    else
  (void)render_buffer;
#    endif
#    if PATH_GUIDING_LEVEL >= 2
  /* Converting the path segment representation of the random walk into radiance samples. */
  const bool use_guide_direct_light = kernel_data.integrator.use_guide_direct_light;
  const bool use_mis_weights = kernel_data.integrator.use_mis_weights;
  state->guiding.path_segment_storage->PrepareSamples(
      false, nullptr, use_mis_weights, use_guide_direct_light, false);
#    endif

#    if defined(PATH_GUIDING_DEBUG_VALIDATE) && PATH_GUIDING_LEVEL >= 5
  bool validSamples = state->guiding.path_segment_storage->ValidateSamples();
  if (!validSamples)
    std::cout << "!!! validSamples = " << validSamples << " !!!" << std::endl;
#    endif

#    if PATH_GUIDING_LEVEL >= 3
  /* Pushing the radiance data/samples from the current random walk/path to the global sample
   * stoarge. */
  size_t nSamples = 0;
  const openpgl::cpp::SampleData *samples = state->guiding.path_segment_storage->GetSamples(
      nSamples);
  kg->opgl_sample_data_storage->AddSamples(samples, nSamples);
#    endif

#    if PATH_GUIDING_LEVEL >= 1
  /* Clearing the storage and the pointer to the current segment to be ready for the next random
   * walk/path. */
  state->guiding.path_segment_storage->Clear();
  state->guiding.path_segment = nullptr;
#    endif
}
#  endif
#endif

CCL_NAMESPACE_END
