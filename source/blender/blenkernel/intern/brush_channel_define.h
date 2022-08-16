#if defined(BRUSH_CHANNEL_DEFINE_EXTERNAL) || defined(BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES) || \
    defined(BRUSH_CHANNEL_DEFINE_INTERNAL)
#  ifdef MAKE_PROP
#    undef MAKE_PROP
#  endif
#  ifdef MAKE_PROP_EX
#    undef MAKE_PROP_EX
#  endif
#endif

#ifdef BRUSH_CHANNEL_DEFINE_EXTERNAL
#  define MAKE_PROP_EX(idname, category, flag, ...) MAKE_PROP(idname, category, uiflags)
#  define MAKE_PROP(idname, category, ...) extern const char *BRUSH_BUILTIN_##idname;
#elif defined(BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES)
#  define MAKE_PROP_EX(idname, category, flag, ...) MAKE_PROP(idname, category)
#  define MAKE_PROP(idname, category, ...) const char *BRUSH_BUILTIN_##idname = #  idname;
#elif defined(BRUSH_CHANNEL_DEFINE_INTERNAL)
#  define MAKE_PROP(idname, category, ...) {#  idname, category, {__VA_ARGS__}, 0, {}},
#  define MAKE_PROP_EX(idname, category, flag, ...) {#  idname, category, {__VA_ARGS__}, flag, {}},

#endif

#ifdef SHOW_WORKSPACE
#  undef SHOW_WORKSPACE
#endif
#ifdef SHOW_CONTEXT
#  undef SHOW_CONTEXT
#endif
#ifdef SHOW_CONTEXT
#  undef SHOW_CONTEXT
#endif
#ifdef SHOW_ALL
#  undef SHOW_ALL
#endif

/*
  UI visibility flags.  Note that some brush types
  may override these in their own channels, see BKE_brush_channelset_ensure_channels
*/

#define SHOW_WORKSPACE BRUSH_CHANNEL_SHOW_IN_WORKSPACE
#define SHOW_CONTEXT BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU
#define SHOW_HEADER BRUSH_CHANNEL_SHOW_IN_HEADER
#define SHOW_ALL (SHOW_WORKSPACE | SHOW_CONTEXT | SHOW_HEADER)

/* Note that channel sets for individual brush types are built in
 * BKE_brush_channelset_ensure_channels
 */

#define SCULPT_GEO_TOOLS \
  SCULPT_TOOL_DRAW, SCULPT_TOOL_CLAY, SCULPT_TOOL_CLAY_STRIPS, SCULPT_TOOL_FLATTEN, \
      SCULPT_TOOL_PINCH, SCULPT_TOOL_FILL, SCULPT_TOOL_INFLATE, SCULPT_TOOL_THUMB, \
      SCULPT_TOOL_SNAKE_HOOK, SCULPT_TOOL_CREASE, SCULPT_TOOL_BLOB, SCULPT_TOOL_CLAY_STRIPS, \
      SCULPT_TOOL_LAYER

#ifdef GEOMETRY
#undef GEOMETRY
#endif

#define GEOMETRY(flag) UI(PAINT_MODE_SCULPT, flag, {SCULPT_GEO_TOOLS})

#ifdef AUTOMASKING
#  undef AUTOMASKING
#endif
#define AUTOMASKING UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE)

#ifdef SCULPT_PAINT
#  undef SCULPT_PAINT
#endif
#ifdef SCULPT_CLOTHABLE_TOOLS
#  undef SCULPT_CLOTHABLE_TOOLS
#endif

#define SCULPT_CLOTHABLE_TOOLS \
  SCULPT_TOOL_CLOTH, SCULPT_TOOL_BOUNDARY, SCULPT_TOOL_ELASTIC_DEFORM, SCULPT_TOOL_POSE

#define SCULPT_PAINT UI(PAINT_MODE_SCULPT, SHOW_CONTEXT | SHOW_WORKSPACE, {SCULPT_TOOL_PAINT})

MAKE_PROP(size, "Basic", UI(SHOW_ALL))
MAKE_PROP(unprojected_radius, "Basic", UI(SHOW_ALL))
MAKE_PROP(strength, "Basic", UI(SHOW_ALL))
MAKE_PROP(automasking_boundary_edges_propagation_steps, "Automasking", AUTOMASKING)
MAKE_PROP(use_automasking_boundary_edges, "Automasking", AUTOMASKING)
MAKE_PROP(use_automasking_boundary_face_sets, "Automasking", AUTOMASKING)
MAKE_PROP(use_automasking_face_sets, "Automasking", AUTOMASKING)
MAKE_PROP(use_automasking_topology, "Automasking", AUTOMASKING)
MAKE_PROP(area_radius_factor, "Basic", GEOMETRY(SHOW_WORKSPACE|SHOW_CONTEXT))
MAKE_PROP(blur_kernel_radius, "Basic")
MAKE_PROP(boundary_offset, "Basic", UI(PAINT_MODE_SCULPT, SHOW_ALL, {SCULPT_TOOL_BOUNDARY}))
MAKE_PROP(clone_alpha,
          "Basic",
          UI(PAINT_MODE_TEXTURE_2D, SHOW_WORKSPACE),
          UI(PAINT_MODE_TEXTURE_3D, SHOW_WORKSPACE))
MAKE_PROP(clone_offset,
          "Basic",
          UI(PAINT_MODE_TEXTURE_2D, SHOW_WORKSPACE),
          UI(PAINT_MODE_TEXTURE_3D, SHOW_WORKSPACE))
MAKE_PROP(crease_pinch_factor, "Basic")
MAKE_PROP(cursor_color_add, "Basic")
MAKE_PROP(cursor_color_subtract, "Basic")
MAKE_PROP(cursor_overlay_alpha, "Basic")
MAKE_PROP(disconnected_distance_max, "Basic")
MAKE_PROP(elastic_deform_volume_preservation,
          "Basic",
          UI(PAINT_MODE_SCULPT,
             SHOW_WORKSPACE,
             {SCULPT_TOOL_ELASTIC_DEFORM | SCULPT_TOOL_BOUNDARY | SCULPT_TOOL_CLOTH}))
MAKE_PROP(falloff_angle, "Basic")
MAKE_PROP(fill_threshold, "Basic")
MAKE_PROP(flow, "Basic")
MAKE_PROP(grad_spacing, "Basic")
MAKE_PROP(hardness, "Basic")
MAKE_PROP(height, "Basic", UI(PAINT_MODE_SCULPT, SHOW_ALL, {SCULPT_TOOL_LAYER}))
MAKE_PROP(invert_to_scrape_fill, "Basic")
MAKE_PROP(mask_overlay_alpha, "Basic")
MAKE_PROP(mask_stencil_dimension, "Basic")
MAKE_PROP(mask_stencil_pos, "Basic")
MAKE_PROP(multiplane_scrape_angle, "Basic")
MAKE_PROP(normal_radius_factor, "Basic")
MAKE_PROP(normal_weight, "Basic")
MAKE_PROP(plane_offset, "Basic")
MAKE_PROP(plane_trim, "Basic")
MAKE_PROP(rake_factor, "Basic")
MAKE_PROP(sharp_threshold, "Basic")
MAKE_PROP(show_multiplane_scrape_planes_preview, "Basic")
MAKE_PROP(stencil_dimension, "Basic")
MAKE_PROP(stencil_pos, "Basic")
MAKE_PROP(surface_smooth_current_vertex, "Basic")
MAKE_PROP(surface_smooth_iterations, "Basic")
MAKE_PROP(surface_smooth_shape_preservation, "Basic")
MAKE_PROP(texture_overlay_alpha, "Basic")
MAKE_PROP(texture_sample_bias, "Basic")
MAKE_PROP(tilt_strength_factor, "Basic")
MAKE_PROP(tip_roundness, "Basic")
MAKE_PROP(tip_scale_x, "Basic")
MAKE_PROP(use_accumulate, "Basic", UI(SHOW_ALL))
MAKE_PROP(use_adaptive_space, "Basic")
MAKE_PROP(use_airbrush, "Basic")
MAKE_PROP(use_alpha, "Basic")
MAKE_PROP(use_anchor, "Basic")
MAKE_PROP(use_connected_only, "Basic")
MAKE_PROP(use_cursor_overlay, "Basic")
MAKE_PROP(use_cursor_overlay_override, "Basic")
MAKE_PROP(use_curve, "Basic")
MAKE_PROP(use_custom_icon, "Basic")
MAKE_PROP(use_edge_to_edge, "Basic")
MAKE_PROP(use_frontface, "Basic")
MAKE_PROP(use_frontface_falloff, "Basic")
MAKE_PROP(use_grab_active_vertex, "Basic")
MAKE_PROP(use_grab_silhouette, "Basic")
MAKE_PROP(use_line, "Basic")
MAKE_PROP(use_multiplane_scrape_dynamic, "Basic")
MAKE_PROP(use_original_normal, "Basic")
MAKE_PROP(use_original_plane, "Basic")
MAKE_PROP(use_persistent, "Basic")
MAKE_PROP(use_plane_trim, "Basic")
MAKE_PROP(use_primary_overlay, "Basic")
MAKE_PROP(use_primary_overlay_override, "Basic")
MAKE_PROP(use_restore_mesh, "Basic")
MAKE_PROP(use_secondary_overlay, "Basic")
MAKE_PROP(use_secondary_overlay_override, "Basic")
MAKE_PROP(use_smooth_stroke, "Basic", UI(SHOW_WORKSPACE))
MAKE_PROP(use_space, "Basic", UI(SHOW_WORKSPACE))
MAKE_PROP(use_space_attenuation, "Basic", UI(SHOW_WORKSPACE))
MAKE_PROP(use_vertex_grease_pencil, "Basic")
MAKE_PROP(weight, "Basic")
MAKE_PROP(wet_paint_radius_factor, "Basic")
MAKE_PROP(cloth_constraint_softbody_strength,
          "Cloth",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(cloth_damping, "Cloth", UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(cloth_mass, "Cloth", UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(cloth_sim_falloff,
          "Cloth",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(cloth_sim_limit,
          "Cloth",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(use_cloth_collision,
          "Cloth",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(use_cloth_pin_simulation_boundary,
          "Cloth",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_CLOTHABLE_TOOLS}))
MAKE_PROP(color,
          "Paint",
          UI(PAINT_MODE_SCULPT, SHOW_ALL, {SCULPT_TOOL_PAINT}),
          UI(PAINT_MODE_TEXTURE_2D, SHOW_ALL),
          UI(PAINT_MODE_TEXTURE_3D, SHOW_ALL),
          UI(PAINT_MODE_VERTEX, SHOW_ALL))
MAKE_PROP(density, "Paint", SCULPT_PAINT)
MAKE_PROP(rate, "Paint")
MAKE_PROP(secondary_color,
          "Paint",
          UI(PAINT_MODE_SCULPT, SHOW_ALL, {SCULPT_TOOL_PAINT}),
          UI(PAINT_MODE_TEXTURE_2D, SHOW_ALL),
          UI(PAINT_MODE_TEXTURE_3D, SHOW_ALL),
          UI(PAINT_MODE_VERTEX, SHOW_ALL))
MAKE_PROP(use_paint_antialiasing, "Paint", UI(PAINT_MODE_TEXTURE_2D, SHOW_WORKSPACE))
MAKE_PROP(use_paint_weight, "Paint")
MAKE_PROP(wet_mix, "Paint", SCULPT_PAINT)
MAKE_PROP(wet_persistence, "Paint", SCULPT_PAINT)
MAKE_PROP(pose_ik_segments, "Pose", UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_TOOL_POSE}))
MAKE_PROP(pose_offset, "Pose", UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_TOOL_POSE}))
MAKE_PROP(pose_smooth_iterations,
          "Pose",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_TOOL_POSE}))
MAKE_PROP(use_pose_ik_anchored, "Pose", UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_TOOL_POSE}))
MAKE_PROP(use_pose_lock_rotation,
          "Pose",
          UI(PAINT_MODE_SCULPT, SHOW_WORKSPACE, {SCULPT_TOOL_POSE}))
MAKE_PROP(auto_smooth_factor, "Smooth", GEOMETRY(SHOW_WORKSPACE | SHOW_CONTEXT))
MAKE_PROP(topology_rake_factor, "Smooth", GEOMETRY(SHOW_WORKSPACE))
MAKE_PROP(dash_ratio, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(dash_samples, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(jitter, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(jitter_absolute, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(smooth_stroke_factor, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(smooth_stroke_radius, "Stroke", UI(SHOW_WORKSPACE))
MAKE_PROP(spacing, "Stroke", UI(SHOW_WORKSPACE))



// tst
#undef MAKE_PROP
#undef MAKE_PROP_EX
#undef SHOW_WORKSPACE
#undef SHOW_HEADER
#undef SHOW_CONTEXT
#undef SHOW_ALL
#undef AUTOMASKING
#undef SCULPT_GEO_TOOLS
#undef SCULPT_PAINT
#undef SCULPT_CLOTHABLE_TOOLS
#undef GEOMETRY
