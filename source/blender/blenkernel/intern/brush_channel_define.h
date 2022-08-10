#ifdef MAKE_PROP
#  undef MAKE_PROP
#endif
#ifdef MAKE_PROP_EX
#  undef MAKE_PROP_EX
#endif

#ifdef BRUSH_CHANNEL_DEFINE_EXTERNAL
#  define MAKE_PROP(idname, category, uiflag) extern const char *BRUSH_BUILTIN_##idname;
#  define MAKE_PROP_EX(idname, category, uiflag, flag) MAKE_PROP(idname)
#elif defined(BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES)
#  define MAKE_PROP(idname, category, uiflag) const char *BRUSH_BUILTIN_##idname = #  idname;
#  define MAKE_PROP_EX(idname, category, uiflag, flag) MAKE_PROP(idname)
#elif defined(BRUSH_CHANNEL_DEFINE_INTERNAL)
#  define MAKE_PROP_EX(idname, category, uiflag, flag) {#  idname, category, uiflag, flag, {}},
#  define MAKE_PROP(idname, category, uiflag) MAKE_PROP_EX(idname, category, uiflag, 0)

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
#undef SHOW_ALL
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

MAKE_PROP(size, "Basic", SHOW_ALL)
MAKE_PROP(unprojected_radius, "Basic", SHOW_ALL)
MAKE_PROP(strength, "Basic", SHOW_ALL)
MAKE_PROP(auto_smooth_factor, "Smooth", SHOW_WORKSPACE | SHOW_CONTEXT)
MAKE_PROP(density, "Basic", 0)
MAKE_PROP(tip_scale_x, "Basic", 0)

MAKE_PROP(color, "Color", 0)
MAKE_PROP(secondary_color, "Color", 0)
MAKE_PROP(wet_mix, "Color", 0)
MAKE_PROP(flow, "Color", 0)
MAKE_PROP(wet_persistence, "Color", 0)

#undef SHOW_WORKSPACE
#undef SHOW_HEADER
#undef SHOW_CONTEXT
#undef SHOW_ALL
