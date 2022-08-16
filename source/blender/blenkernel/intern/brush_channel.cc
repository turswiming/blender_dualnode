#include "MEM_guardedalloc.h"

#include "BLI_assert.h"
#include "BLI_compiler_attrs.h"
#include "BLI_compiler_compat.h"
#include "BLI_ghash.h"
#include "BLI_index_range.hh"
#include "BLI_listbase.h"
#include "BLI_map.hh"
#include "BLI_math.h"
#include "BLI_math_vec_types.hh"
#include "BLI_string_ref.hh"
#include "BLI_string_utils.h"
#include "BLI_vector.hh"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_path.h"
#include "RNA_prototypes.h"

#include "DNA_brush_channel_types.h"
#include "DNA_brush_enums.h"
#include "DNA_brush_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "DNA_color_types.h"
#include "DNA_material_types.h"

#include "BKE_brush.h"
#include "BKE_brush_channel.h"
#include "BKE_colortools.h"
#include "BKE_idprop.h"
#include "BKE_idprop.hh"
#include "BKE_lib_id.h"
#include "BKE_paint.h"
#include "BKE_pbvh.h"
#include "BLO_read_write.h"

#include <string>
#include <vector>

const char builtin_brush_categories[][128] = {"Basic", "Smooth", "Color"};

static BrushChannelType empty_brush_type = {"error", "error", "error", "error"};

#define BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES
#include "brush_channel_define.h"
#undef BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES

using string = std::string;

using KeyString = const char *;
using blender::float3;
using blender::float4;
using blender::IndexRange;
using blender::Map;
using blender::StringRef;
using blender::Vector;

static Map<StringRef, BrushChannelType> builtin_channels;

using ChannelNameMap = Map<StringRef, BrushChannel *>;

ATTR_NO_OPT static void init_builtin_brush_channels()
{
  struct UI {
    int mode;
    Vector<int> tools;
    int uiflag;
    bool all_tool_modes = false;
    bool all_tools = false;
    bool exists = true;

    UI(const UI &b)
    {
      tools = b.tools;
      uiflag = b.uiflag;
      all_tool_modes = b.all_tool_modes;
      all_tools = b.all_tools;
      exists = b.exists;
      mode = b.mode;
    }

    UI(int _mode, int _uiflag, int _tool) : mode(_mode), uiflag(_uiflag)
    {
      tools.append(_tool);
    }

    UI(int _uiflag) : uiflag(_uiflag)
    {
      all_tool_modes = all_tools = true;
    }

    UI(ePaintMode _mode, int _uiflag) : mode((int)_mode), uiflag(_uiflag)
    {
      all_tools = true;
    }

    UI(ePaintMode mode, int uiflag, int tool)
    {
      UI::UI((int)mode, uiflag, tool);
    }

    UI(ePaintMode _mode, int _uiflag, Vector<int> _tools)
        : mode((int)_mode), tools(_tools), uiflag(_uiflag)
    {
    }
  };

  struct ChannelProp {
    char path[512];
    char category[64];
    Vector<UI> extra_uiflags;
    int flag;
    BrushMappingPreset mappings;
  };

#ifdef BRUSH_CHANNEL_DEFINE_EXTERNAL
#  undef BRUSH_CHANNEL_DEFINE_EXTERNAL
#endif

#ifdef BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES
#  undef BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES
#endif
#ifdef MAKE_PROP
#  undef MAKE_PROP
#endif
#ifdef MAKE_PROP_EX
#  undef MAKE_PROP_EX
#endif

  //#define MAKE_PROP(idname, category, uiflags) MAKE_PROP_EX(idname, category, uiflags, 0)

  //#define SHOW_CONTEXT BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU
  //#define SHOW_WORKSPACE BRUSH_CHANNEL_SHOW_IN_WORKSPACE

#define BRUSH_CHANNEL_DEFINE_INTERNAL
  ChannelProp channel_props[] = {
#include "brush_channel_define.h"
  };

#undef SHOW_CONTEXT
#undef SHOW_WORKSPACE

  StructRNA *srna = RNA_struct_find("Brush");

  Brush dummy = {0};
  PointerRNA _ptr = {&dummy.id}, *ptr = &_ptr;

  printf("total properties: %i\n", (int)ARRAY_SIZE(channel_props));

  for (int i : IndexRange(ARRAY_SIZE(channel_props))) {
    ChannelProp &def = channel_props[i];
    BrushChannelType type = {0};

    BLI_strncpy(type.idname, def.path, sizeof(type.idname));
    BLI_strncpy(type.category, def.category, sizeof(type.category));
    type.flag = def.flag;

    for (auto &ui : def.extra_uiflags) {
      for (int mode = 0; mode < (int)PAINT_MODE_INVALID; mode++) {
        if (mode != ui.mode && !ui.all_tool_modes) {
          continue;
        }

        int uiflag = ui.uiflag ? ui.uiflag : -1;

        if (ui.all_tools) {
          for (int j = 0; j < ARRAY_SIZE(type.paint_mode_uiflags); j++) {
            type.paint_mode_uiflags[mode].tools[j] = uiflag;
          }
        }
        else {
          for (auto tool : ui.tools) {
            type.paint_mode_uiflags[mode].tools[tool] = uiflag;
          }
        }
      }
    }

    PropertyRNA *prop = RNA_struct_type_find_property(srna, def.path);
    BLI_assert(prop);

    if (!prop) {
      printf("%s: Missing property %s\n", __func__, def.path);
      continue;
    }

    PropertyType prop_type = RNA_property_type(prop);
    PropertySubType prop_subtype = RNA_property_subtype(prop);

    const char *uiname = RNA_property_ui_name(prop);
    BLI_strncpy(type.uiname, uiname, sizeof(type.uiname));

    type.min = 0.0;
    type.max = 1.0;
    type.soft_min = 0.0;
    type.soft_max = 1.0;

    switch (prop_type) {
      case PROP_BOOLEAN:
        type.type = BRUSH_CHANNEL_TYPE_BOOL;
        break;
      case PROP_INT: {
        int min, max, soft_min, soft_max;
        int step;

        RNA_property_int_range(nullptr, prop, &min, &max);
        RNA_property_int_ui_range(nullptr, prop, &soft_min, &soft_max, &step);

        type.min = (float)min;
        type.max = (float)max;
        type.soft_min = (float)soft_min;
        type.soft_max = (float)soft_max;
        type.type = BRUSH_CHANNEL_TYPE_INT;
        break;
      }
      case PROP_FLOAT: {
        float precision, step;

        RNA_property_float_range(ptr, prop, &type.min, &type.max);
        RNA_property_float_ui_range(ptr, prop, &type.soft_min, &type.soft_max, &step, &precision);

        if (!RNA_property_array_check(prop)) {
          type.type = BRUSH_CHANNEL_TYPE_FLOAT;
        }
        else {
          int dimen = RNA_property_array_length(ptr, prop);

          switch (dimen) {
            case 3:
              type.type = BRUSH_CHANNEL_TYPE_VEC3;
              break;
            case 4:
              type.type = BRUSH_CHANNEL_TYPE_VEC4;
              break;
            default:
              BLI_assert_unreachable();
          }
        }
        break;
      }
      default:
        break;
    }

    switch (prop_subtype) {
      case PROP_COLOR:
      case PROP_COLOR_GAMMA:
        type.subtype = BRUSH_CHANNEL_COLOR;
        break;
      case PROP_FACTOR:
        type.subtype = BRUSH_CHANNEL_FACTOR;
        break;
      case PROP_PERCENTAGE:
        type.subtype = BRUSH_CHANNEL_PERCENT;
        break;
      case PROP_ANGLE:
        type.subtype = BRUSH_CHANNEL_ANGLE;
        break;
      case PROP_PIXEL:
        type.subtype = BRUSH_CHANNEL_PIXEL;
        break;
      default:
        break;
    }

    builtin_channels.add(strdup(type.idname), type);
  }
#undef BRUSH_CHANNEL_DEFINE_INTERNAL
}
static void check_builtin_brush_channels()
{
  if (builtin_channels.size() == 0) {
    init_builtin_brush_channels();
  }
}

ATTR_NO_OPT ChannelNameMap *get_namemap(BrushChannelSet *chset)
{
  return reinterpret_cast<ChannelNameMap *>(chset->channelmap);
}

BrushChannelSet *BKE_brush_channelset_create()
{
  BrushChannelSet *chset = MEM_cnew<BrushChannelSet>("BrushChannelSet");

  chset->channels.first = chset->channels.last = nullptr;
  chset->channels_num = 0;

  ChannelNameMap *map = MEM_new<ChannelNameMap>("ChannelNameMap");
  chset->channelmap = static_cast<void *>(map);

  return chset;
}

void BKE_brush_channel_free_data(BrushChannel *ch)
{
  MEM_SAFE_FREE(ch->category);

  for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
    BrushMapping *mp = ch->mappings + i;

    if (mp->curve.curve) {
      BKE_curvemapping_free(mp->curve.curve);
    }
  }

  if (ch->curve.curve) {
    BKE_curvemapping_free(ch->curve.curve);
  }
}

void BKE_brush_channelset_free(BrushChannelSet *chset)
{
  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BKE_brush_channel_free_data(ch);
  }

  BLI_freelistN(&chset->channels);

  MEM_delete<ChannelNameMap>(get_namemap(chset));
  MEM_freeN(static_cast<void *>(chset));
}

static void brush_mapping_reset(BrushMapping *mp, int type)
{
  mp->type = type;
  mp->premultiply_factor = 1.0f;
  mp->min = 0.0f;
  mp->max = 1.0f;
  mp->factor = 1.0f;
  mp->blendmode = MA_RAMP_MULT;
  mp->curve.preset = BRUSH_CURVE_LIN;
}

BrushChannel *BKE_brush_channelset_add(BrushChannelSet *chset, BrushChannelType *type)
{
  BrushChannel *ch = MEM_cnew<BrushChannel>("BrushChannel");

  BLI_strncpy(ch->idname, type->idname, sizeof(ch->idname));
  BLI_strncpy(ch->uiname, type->uiname, sizeof(ch->uiname));

  ch->def = type;
  ch->type = type->type;
  ch->flag = type->flag;

  ch->ui_order = chset->channels_num;

  BLI_addtail(&chset->channels, static_cast<void *>(ch));
  get_namemap(chset)->add(ch->idname, ch);
  chset->channels_num++;

  for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
    BrushMapping *mp = ch->mappings + i;
    brush_mapping_reset(mp, i);
  }

  BKE_brush_channelset_ui_order_check(chset);

  return ch;
}

void BKE_brush_channelset_begin(BrushChannelSet *chset, BrushChannelType *type)
{
  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    ch->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;
  }
}

ATTR_NO_OPT BrushChannel *_BKE_brush_channelset_lookup(BrushChannelSet *chset, const char *idname)
{
  ChannelNameMap *namemap = get_namemap(chset);

  return namemap->lookup(StringRef(idname));
}

bool _BKE_brush_channelset_has(BrushChannelSet *chset, const char *idname)
{
  return get_namemap(chset)->contains(idname);
}

const char *BKE_brush_channel_category_get(BrushChannel *ch)
{
  return ch->category ? ch->category : ch->def->category;
}

const void BKE_brush_channel_category_set(BrushChannel *ch, const char *category)
{
  if (STREQ(category, ch->def->category)) {
    MEM_SAFE_FREE(ch->category);
    ch->category = nullptr;

    return;
  }

  ch->category = BLI_strdup(category);
}

ATTR_NO_OPT BrushChannel *_BKE_brush_channelset_ensure(BrushChannelSet *chset, const char *idname)
{
  if (!builtin_channels.contains(idname)) {
    printf("channel types:\n");
    for (StringRef key : builtin_channels.keys()) {
      printf("  %s\n", key.data());
    }

    printf("Unknown brush channel %s\n", idname);
    return nullptr;
  }

  ChannelNameMap *namemap = get_namemap(chset);
  BrushChannelType &type = builtin_channels.lookup(idname);

  if (!namemap->contains(type.idname)) {
    BKE_brush_channelset_add(chset, &type);
  }

  return _BKE_brush_channelset_lookup(chset, idname);
}

void BKE_brush_channelset_ensure_all_modes(Brush *brush)
{
  if (!brush->channels) {
    brush->channels = BKE_brush_channelset_create();
    BKE_brush_channelset_ensure_channels(brush->channels, PAINT_MODE_INVALID, 0);
  }

  if (brush->sculpt_tool) {
    BKE_brush_channelset_ensure_channels(brush->channels, PAINT_MODE_SCULPT, brush->sculpt_tool);
  }

  if (brush->vertexpaint_tool) {
    BKE_brush_channelset_ensure_channels(
        brush->channels, PAINT_MODE_VERTEX, brush->vertexpaint_tool);
  }

  if (brush->imagepaint_tool) {
    BKE_brush_channelset_ensure_channels(
        brush->channels, PAINT_MODE_TEXTURE_3D, brush->imagepaint_tool);
  }

  if (brush->curves_sculpt_tool) {
    BKE_brush_channelset_ensure_channels(
        brush->channels, PAINT_MODE_SCULPT_CURVES, brush->curves_sculpt_tool);
  }

  if (brush->weightpaint_tool) {
    BKE_brush_channelset_ensure_channels(
        brush->channels, PAINT_MODE_WEIGHT, brush->weightpaint_tool);
  }
}
void BKE_brush_channelset_ensure_channels(BrushChannelSet *chset, ePaintMode mode, int tool)
{
  check_builtin_brush_channels();

  BKE_brush_channelset_ensure(chset, size);
  BKE_brush_channelset_ensure(chset, unprojected_radius);
  BKE_brush_channelset_ensure(chset, strength);
  BKE_brush_channelset_ensure(chset, spacing);

  if (mode != PAINT_MODE_INVALID) {
    for (BrushChannelType &type : builtin_channels.values()) {
      int uiflag = type.paint_mode_uiflags[(int)mode].tools[tool];

      if (uiflag) {
        BrushChannel *ch = _BKE_brush_channelset_ensure(chset, type.idname);

        if (uiflag != -1) {
          if (!(ch->ui_flag & BRUSH_CHANNEL_SHOW_IN_HEADER_USER_SET)) {
            ch->ui_flag &= ~BRUSH_CHANNEL_SHOW_IN_HEADER;
            ch->ui_flag |= uiflag & BRUSH_CHANNEL_SHOW_IN_HEADER;
          }
          if (!(ch->ui_flag & BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU_USER_SET)) {
            ch->ui_flag &= ~BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU;
            ch->ui_flag |= uiflag & BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU;
          }
          if (!(ch->ui_flag & BRUSH_CHANNEL_SHOW_IN_WORKSPACE_USER_SET)) {
            ch->ui_flag &= ~BRUSH_CHANNEL_SHOW_IN_WORKSPACE;
            ch->ui_flag |= uiflag & BRUSH_CHANNEL_SHOW_IN_WORKSPACE;
          }
        }
      }
    }
  }
  /* Some helper lambdas */

  auto _ensure = [&](const char *idname) { return _BKE_brush_channelset_ensure(chset, idname); };

#ifdef ensure
#  undef ensure
#endif

#define ensure _ensure(make_builtin_ch_name(idname), ui_flag);

  auto _ensure_ui = [&](const char *idname, int ui_flag) {
    _ensure(idname);

    BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);
    ch->ui_flag |= ui_flag;

    return ch;
  };

#ifdef ensure_ui
#  undef ensure_ui
#endif

#define ensure_ui(idname, ui_flag) _ensure_ui(make_builtin_ch_name(idname), ui_flag)

  const int SHOW_WORKSPACE = BRUSH_CHANNEL_SHOW_IN_WORKSPACE;
  const int SHOW_CONTEXT = BRUSH_CHANNEL_SHOW_IN_CONTEXT_MENU;
  const int SHOW_HEADER = BRUSH_CHANNEL_SHOW_IN_HEADER;
  const int SHOW_ALL = (SHOW_WORKSPACE | SHOW_CONTEXT | SHOW_HEADER);

  if (mode == PAINT_MODE_SCULPT) {
    switch (tool) {
      case SCULPT_TOOL_PAINT:
        ensure_ui(wet_mix, SHOW_WORKSPACE);
        ensure_ui(wet_persistence, SHOW_WORKSPACE);
        ensure_ui(color, SHOW_ALL);
        ensure_ui(secondary_color, SHOW_ALL);
        ensure_ui(flow, SHOW_WORKSPACE | SHOW_CONTEXT);
        ensure_ui(density, SHOW_WORKSPACE | SHOW_CONTEXT);
        ensure_ui(tip_scale_x, SHOW_WORKSPACE | SHOW_CONTEXT);
        break;
    }
  }

#undef ensure
#undef ensure_ui
}

char *BKE_brush_channel_rna_path(const ID *owner, const BrushChannel *ch)
{
  switch (GS(owner->name)) {
    case ID_BR:
      return BLI_strdup(ch->idname);
    case ID_SCE: {
      string path = "tool_settings.unified_properties[\"";

      path += string(ch->idname) + "\"]";

      return BLI_strdup(path.c_str());
    }
    default:
      return BLI_strdup("");
  }
}

void brush_channel_ensure_value(const ID *id, BrushChannel *ch)
{
  if (!(ch->flag & BRUSH_CHANNEL_NEEDS_EVALUATE)) {
    return;
  }

  string path;
  StructRNA *srna;

  string rnaname = ch->idname;

  if (GS(id->name) == ID_BR) {
    path = rnaname;
    srna = &RNA_Brush;
  }
  else {
    path = "tool_settings.unified_properties[\"";
    path += rnaname;
    path += "\"]";

    rnaname = string("[\"") + rnaname + string("\"]");
    srna = &RNA_Scene;
  }

  ch->flag &= ~BRUSH_CHANNEL_NEEDS_EVALUATE;
  PointerRNA ptr, ptr2;
  PropertyRNA *prop = nullptr;

  RNA_pointer_create(const_cast<ID *>(id), srna, const_cast<ID *>(id), &ptr);
  if (RNA_path_resolve(&ptr, path.c_str(), &ptr2, &prop)) {
    PropertyType prop_type = RNA_property_type(prop);

    switch (prop_type) {
      case PROP_FLOAT:
        if (RNA_property_array_check(prop)) {
          RNA_float_get_array(&ptr2, rnaname.c_str(), ch->vector);
        }
        else {
          ch->fvalue = RNA_float_get(&ptr2, rnaname.c_str());
        }
        break;
      case PROP_INT:
        ch->ivalue = RNA_int_get(&ptr2, rnaname.c_str());
        break;
      default:
        printf("Unknown prop type %d\n", (int)prop_type);
        break;
    }
  }
  else {
    printf("Error looking up path %s", path.c_str());
    return;
  }
}

static bool channel_has_mappings(BrushChannel *ch)
{
  for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
    if (ch->mappings[i].flag & BRUSH_MAPPING_ENABLED) {
      return true;
    }
  }

  return false;
}

/* idx is used by vector channels */
double BKE_brush_channel_eval_mappings(BrushChannel *ch,
                                       BrushMappingData *mapdata,
                                       double f,
                                       int idx)
{

  if (idx == 3 && !(ch->flag & BRUSH_CHANNEL_APPLY_MAPPING_TO_ALPHA)) {
    return f;
  }

  if (mapdata) {
    double factor = f;

    for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
      BrushMapping *mp = ch->mappings + i;

      if (!(mp->flag & BRUSH_MAPPING_ENABLED)) {
        continue;
      }

      float inputf = ((float *)mapdata)[i] * mp->premultiply_factor;

      switch ((eBrushMappingFunc)mp->mapfunc) {
        case BRUSH_MAPFUNC_NONE:
          break;
        case BRUSH_MAPFUNC_SAW:
          inputf -= floorf(inputf);
          break;
        case BRUSH_MAPFUNC_TENT:
          inputf -= floorf(inputf);
          inputf = 1.0f - fabs(inputf - 0.5f) * 2.0f;
          break;
        case BRUSH_MAPFUNC_COS:
          inputf = 1.0f - (cos(inputf * (float)M_PI * 2.0f) * 0.5f + 0.5f);
          break;
        case BRUSH_MAPFUNC_CUTOFF:
          /*Cutoff is meant to create a fadeout effect,
            which requires inverting the input.  To avoid
            user confusion we just do it here instead of making
            them check the inverse checkbox.*/
          inputf = 1.0f - inputf;
          CLAMP(inputf, 0.0f, mp->func_cutoff * 2.0f);
          break;
        case BRUSH_MAPFUNC_SQUARE:
          inputf -= floorf(inputf);
          inputf = inputf > mp->func_cutoff ? 1.0f : 0.0f;
          break;
        default:
          break;
      }

      if (mp->flag & BRUSH_MAPPING_INVERT) {
        inputf = 1.0f - inputf;
      }

      double f2 = BKE_brush_curve_strength_ex(
          mp->curve.preset, mp->curve.curve, inputf, 1.0f, false);
      f2 = mp->min + (mp->max - mp->min) * f2;

      /* make sure to update blend_items in rna_brush_engine.c
        when adding new mode implementations */
      switch (mp->blendmode) {
        case MA_RAMP_BLEND:
          break;
        case MA_RAMP_MULT:
          f2 *= factor;
          break;
        case MA_RAMP_DIV:
          f2 = factor / (f2 == 0.0f ? 0.0001f : f2);
          break;
        case MA_RAMP_ADD:
          f2 += factor;
          break;
        case MA_RAMP_SUB:
          f2 = factor - f2;
          break;
        case MA_RAMP_DIFF:
          f2 = fabsf(factor - f2);
          break;
        default:
          printf("Unsupported brush mapping blend mode for %s (%s); will mix instead\n",
                 ch->uiname,
                 ch->idname);
          break;
      }

      factor += (f2 - factor) * mp->factor;
    }

    f = factor;
    CLAMP(f, ch->def->min, ch->def->max);
  }

  return f;
}

static BrushChannel *brush_channel_final(const Brush *brush,
                                         const Scene *scene,
                                         const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(brush->channels, idname);

  if (BKE_brush_channel_inherits(brush, scene->toolsettings, ch)) {
    ch = _BKE_brush_channelset_lookup(scene->toolsettings->unified_channels, idname);
    brush_channel_ensure_value(&scene->id, ch);
  }
  else {
    brush_channel_ensure_value(&brush->id, ch);
  }

  return ch;
}

float _BKE_brush_eval_float(const Brush *brush,
                            const Scene *scene,
                            const char *idname,
                            BrushMappingData *mapping)
{
  BrushChannel *ch = brush_channel_final(brush, scene, idname);

  return mapping ? BKE_brush_channel_eval_mappings(ch, mapping, ch->fvalue, 0) : ch->fvalue;
}

int _BKE_brush_eval_int(const Brush *brush,
                        const Scene *scene,
                        const char *idname,
                        BrushMappingData *mapping)
{
  BrushChannel *ch = brush_channel_final(brush, scene, idname);

  return mapping ? (int)BKE_brush_channel_eval_mappings(ch, mapping, (double)ch->ivalue, 0) :
                   ch->ivalue;
}

static void brush_channel_evaluate(const Brush *br,
                                   const Scene *scene,
                                   BrushChannelSet *chset,
                                   const char *idname,
                                   BrushMappingData *mapping)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);

  if (!mapping) {
    return;
  }

  switch (ch->type) {
    case BRUSH_CHANNEL_TYPE_FLOAT:
      ch->fvalue = BKE_brush_channel_eval_mappings(ch, mapping, ch->fvalue, 0);
      break;
    case BRUSH_CHANNEL_TYPE_INT:
    case BRUSH_CHANNEL_TYPE_BOOL:
    case BRUSH_CHANNEL_TYPE_ENUM:
      ch->ivalue = (int)BKE_brush_channel_eval_mappings(ch, mapping, ch->ivalue, 0);
      break;
    case BRUSH_CHANNEL_TYPE_VEC3:
    case BRUSH_CHANNEL_TYPE_VEC4: {
      int size = ch->type == BRUSH_CHANNEL_TYPE_VEC3 ? 3 : 4;

      for (int i = 0; i < size; i++) {
        ch->vector[i] = (int)BKE_brush_channel_eval_mappings(ch, mapping, ch->vector[i], i);
        break;
      }
      break;
    }
  }
}

float _BKE_brush_channelset_float_get(BrushChannelSet *chset, const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);
  return ch->fvalue;
}

void _BKE_brush_channelset_float_set(BrushChannelSet *chset, const char *idname, float f)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);
  ch->fvalue = f;
}

int _BKE_brush_channelset_int_get(BrushChannelSet *chset, const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);
  return ch->ivalue;
}

void _BKE_brush_channelset_int_set(BrushChannelSet *chset, const char *idname, int i)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);
  ch->ivalue = i;
}

bool brush_channelset_rna_path_resolve(const ID *id,
                                       BrushChannelSet *chset,
                                       const char *idname,
                                       PointerRNA *r_ptr,
                                       PropertyRNA **r_prop,
                                       char **final_idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);

  if (!ch) {
    return false;
  }

  bool ok = true;

  StructRNA *srna;

  switch (GS(id->name)) {
    case ID_SCE:
      srna = &RNA_Scene;
      break;
    case ID_BR:
      srna = &RNA_Brush;
      break;
    default:
      return false;
  }

  PointerRNA id_ptr;
  RNA_pointer_create(
      const_cast<ID *>(id), srna, static_cast<void *>(const_cast<ID *>(id)), &id_ptr);

  char *path = BKE_brush_channel_rna_path(id, ch);
  ok = RNA_path_resolve(&id_ptr, path, r_ptr, r_prop);
  MEM_SAFE_FREE(path);

  *final_idname = static_cast<char *>(
      MEM_mallocN(strlen(idname) + 6, "brush_channelset_rna_path_resolve.final_idname"));

  if (RNA_property_is_idprop(*r_prop)) {
    sprintf(*final_idname, "[\"%s\"]", idname);
  }
  else {
    sprintf(*final_idname, "%s", idname);
  }

  return ok;
}

int _BKE_brush_int_get(const ID *id, BrushChannelSet *chset, const char *idname)
{
  PointerRNA ptr;
  PropertyRNA *prop;
  char *final_idname;
  int ret = 0;

  if (brush_channelset_rna_path_resolve(id, chset, idname, &ptr, &prop, &final_idname)) {
    ret = RNA_int_get(&ptr, final_idname);
  }

  MEM_SAFE_FREE(final_idname);
  return ret;
}

void _BKE_brush_int_set(ID *id, BrushChannelSet *chset, const char *idname, int i)
{
  PointerRNA ptr;
  PropertyRNA *prop;
  char *final_idname;

  if (brush_channelset_rna_path_resolve(id, chset, idname, &ptr, &prop, &final_idname)) {
    RNA_int_set(&ptr, final_idname, i);
    _BKE_brush_channelset_lookup(chset, idname)->ivalue = i;
  }

  MEM_SAFE_FREE(final_idname);
}

float _BKE_brush_float_get(const ID *id, BrushChannelSet *chset, const char *idname)
{
  PointerRNA ptr;
  PropertyRNA *prop;
  char *final_idname;
  float ret = 0;

  if (brush_channelset_rna_path_resolve(id, chset, idname, &ptr, &prop, &final_idname)) {
    ret = RNA_float_get(&ptr, final_idname);
  }

  MEM_SAFE_FREE(final_idname);
  return ret;
}

void _BKE_brush_float_set(ID *id, BrushChannelSet *chset, const char *idname, float f)
{
  PointerRNA ptr;
  PropertyRNA *prop;
  char *final_idname;

  if (brush_channelset_rna_path_resolve(id, chset, idname, &ptr, &prop, &final_idname)) {
    RNA_float_set(&ptr, final_idname, f);
    _BKE_brush_channelset_lookup(chset, idname)->fvalue = f;
  }

  MEM_SAFE_FREE(final_idname);
}

void BKE_brush_mapping_copy_data(BrushMapping *dest, BrushMapping *src)
{
  *dest = *src;

  if (dest->curve.curve) {
    dest->curve.curve = BKE_curvemapping_copy(dest->curve.curve);
    BKE_curvemapping_init(dest->curve.curve);
  }
}

void BKE_brush_mapping_free_data(BrushMapping *mp)
{
  if (mp->curve.curve) {
    BKE_curvemapping_free(mp->curve.curve);
  }
}

void BKE_brush_channel_copy_data(BrushChannel *dest, BrushChannel *src)
{
  BrushChannel *prev = dest->prev, *next = dest->next;

  *dest = *src;
  dest->prev = prev;
  dest->next = next;

  if (dest->category) {
    dest->category = static_cast<char *>(MEM_dupallocN(static_cast<void *>(dest->category)));
  }

  if (dest->curve.curve) {
    dest->curve.curve = BKE_curvemapping_copy(dest->curve.curve);
    BKE_curvemapping_init(dest->curve.curve);
  }

  for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
    BKE_brush_mapping_copy_data(dest->mappings + i, src->mappings + i);
  }
}

BrushChannel *BKE_brush_channel_copy(BrushChannel *ch)
{
  BrushChannel *ch2 = MEM_cnew<BrushChannel>("BrushChannel");
  BKE_brush_channel_copy_data(ch2, ch);
  return ch2;
}

BrushChannelSet *BKE_brush_channelset_copy(BrushChannelSet *chset)
{
  BrushChannelSet *chset2 = BKE_brush_channelset_create();
  ChannelNameMap *namemap2 = get_namemap(chset2);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BrushChannel *ch2 = BKE_brush_channel_copy(ch);

    chset2->channels_num++;
    namemap2->add(ch2->idname, ch2);
    BLI_addtail(&chset2->channels, static_cast<void *>(ch2));
  }

  return chset2;
}

void BKE_brush_channelset_blend_read(BrushChannelSet *chset, BlendDataReader *reader)
{
  check_builtin_brush_channels();

  BLO_read_list(reader, &chset->channels);

  ChannelNameMap *namemap = MEM_new<ChannelNameMap>("ChannelNameMap");
  chset->channelmap = static_cast<void *>(namemap);
  chset->channels_num = 0;

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    chset->channels_num++;

    namemap->add(StringRef(ch->idname), ch);

    if (builtin_channels.contains(ch->idname)) {
      ch->def = &builtin_channels.lookup(ch->idname);
    }
    else {
      ch->def = &empty_brush_type;
    }
    /* Read user-defined category if it exists. */
    BLO_read_data_address(reader, &ch->category);
    BLO_read_data_address(reader, &ch->curve.curve);

    if (ch->curve.curve) {
      BKE_curvemapping_blend_read(reader, ch->curve.curve);
      BKE_curvemapping_init(ch->curve.curve);
    }

    for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
      BrushMapping *mp = ch->mappings + i;

      BLO_read_data_address(reader, &mp->curve.curve);

      if (mp->curve.curve) {
        BKE_curvemapping_blend_read(reader, mp->curve.curve);
        BKE_curvemapping_init(mp->curve.curve);
      }
    }
  }
}

void BKE_brush_channelset_blend_write(BrushChannelSet *chset, BlendWriter *writer)
{
  BLO_write_struct(writer, BrushChannelSet, chset);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BLO_write_struct(writer, BrushChannel, ch);

    if (ch->category) {
      BLO_write_string(writer, ch->category);
    }

    if (ch->curve.curve) {
      BKE_curvemapping_blend_write(writer, ch->curve.curve);
    }

    for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
      BrushMapping *mp = ch->mappings + i;

      if (mp->curve.curve) {
        BKE_curvemapping_blend_write(writer, mp->curve.curve);
      }
    }
  }
}

bool BKE_brush_channel_inherits(const Brush *brush,
                                const ToolSettings *tool_settings,
                                BrushChannel *ch)
{
  BrushChannel *scene_ch = _BKE_brush_channelset_lookup(tool_settings->unified_channels,
                                                        ch->idname);

  if (!scene_ch) {
    return false;
  }

  if (ch->flag & BRUSH_CHANNEL_INHERIT) {
    return true;
  }

  if (scene_ch->flag & BRUSH_CHANNEL_FORCE_INHERIT) {
    return !(ch->flag & BRUSH_CHANNEL_IGNORE_FORCE_INHERIT);
  }

  return false;
}

BrushChannelSet *BKE_brush_channelset_create_final(const Brush *brush,
                                                   const Scene *scene,
                                                   BrushMappingData *mapdata)
{
  BrushChannelSet *chset = BKE_brush_channelset_copy(brush->channels);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BrushChannel *ch1 = _BKE_brush_channelset_lookup(brush->channels, ch->idname);
    BrushChannel *ch2 = _BKE_brush_channelset_lookup(scene->toolsettings->unified_channels,
                                                     ch->idname);

    ch1->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;
    brush_channel_ensure_value(&brush->id, ch1);
    ch2->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;
    brush_channel_ensure_value(&scene->id, ch2);

    bool inherit = BKE_brush_channel_inherits(brush, scene->toolsettings, ch);

    if (inherit) {
      BKE_brush_channel_free_data(ch);
      BKE_brush_channel_copy_data(ch, ch2);
    }

    for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
      BrushMapping *mp = ch->mappings + i;
      BrushMapping *mp1 = ch1->mappings + i;
      BrushMapping *mp2 = ch2->mappings + i;

      if ((mp1->flag & BRUSH_MAPPING_INHERIT_NEVER) || !inherit) {
        BKE_brush_mapping_free_data(mp);
        BKE_brush_mapping_copy_data(mp, mp1);
      }
      else if ((mp1->flag & BRUSH_MAPPING_INHERIT_ALWAYS) || inherit) {
        BKE_brush_mapping_free_data(mp);
        BKE_brush_mapping_copy_data(mp, mp2);
      }
    }

    ch2->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;

    brush_channel_evaluate(brush, scene, chset, ch->idname, mapdata);
  }

  return chset;
}

void BKE_brush_channelset_toolsettings_init(ToolSettings *ts)
{
  check_builtin_brush_channels();

  if (!ts->unified_properties) {
    ts->unified_properties = IDP_New(IDP_GROUP, nullptr, "group");
  }

  if (!ts->unified_channels) {
    ts->unified_channels = BKE_brush_channelset_create();
  }

  for (const BrushChannelType &type : builtin_channels.values()) {
    _BKE_brush_channelset_ensure(ts->unified_channels, type.idname);
  }

  StructRNA *srna = &RNA_Brush;
  Brush defaults = {0};

  BLI_strncpy(defaults.id.name, "BRDefaults", sizeof(defaults.id.name));

  defaults.sculpt_tool = SCULPT_TOOL_DRAW;
  BKE_brush_sculpt_reset(&defaults);

  PointerRNA ptr;
  ptr.owner_id = &defaults.id;
  ptr.data = &defaults;
  ptr.type = &RNA_Brush;

  LISTBASE_FOREACH (BrushChannel *, ch, &ts->unified_channels->channels) {
    IDProperty *idprop = IDP_GetPropertyFromGroup(ts->unified_properties, ch->idname);

    if (!idprop) {
      PointerRNA prop_ptr;
      PropertyRNA *prop;
      double default_value = 0.0;
      float vector4[4] = {0.0f, 0.0f, 0.0f, 0.0f};

      if (RNA_path_resolve(&ptr, ch->idname, &prop_ptr, &prop)) {
        switch (RNA_property_type(prop)) {
          case PROP_BOOLEAN:
            default_value = RNA_boolean_get(&prop_ptr, ch->idname) ? 1.0 : 0.0;
            break;
          case PROP_INT:
            default_value = RNA_int_get(&prop_ptr, ch->idname);
            break;
          case PROP_FLOAT:
            if (RNA_property_array_check(prop)) {
              RNA_float_get_array(&prop_ptr, ch->idname, vector4);
            }
            else {
              default_value = RNA_float_get(&prop_ptr, ch->idname);
            }
            break;
          case PROP_ENUM:
            default_value = (double)RNA_enum_get(&prop_ptr, ch->idname);
            break;
          default:
            break;
        }
        printf("found property!\n");
      }

      IDPropertyTemplate tmpl;
      char type;

      switch (ch->type) {
        case BRUSH_CHANNEL_TYPE_FLOAT:
          tmpl.f = (float)default_value;
          type = IDP_FLOAT;
          break;
        case BRUSH_CHANNEL_TYPE_INT:
        case BRUSH_CHANNEL_TYPE_ENUM:
        case BRUSH_CHANNEL_TYPE_BITMASK:
        case BRUSH_CHANNEL_TYPE_BOOL:
          tmpl.i = (int)default_value;
          type = IDP_INT;
          break;
        case BRUSH_CHANNEL_TYPE_VEC4:
          tmpl.array.type = IDP_FLOAT;
          tmpl.array.len = 4;
          type = IDP_ARRAY;
          break;
        default:
          printf("%s: unsupported brush channel type for unified channel %s: %d\n",
                 __func__,
                 ch->idname,
                 ch->type);
          continue;
      }

      idprop = IDP_New(type, &tmpl, ch->idname);
      IDP_AddToGroup(ts->unified_properties, idprop);

      if (ch->type == BRUSH_CHANNEL_TYPE_VEC4) {
        memcpy(idprop->data.pointer, static_cast<void *>(vector4), sizeof(vector4));
      }
    }

    IDPropertyUIData *uidata = IDP_ui_data_ensure(idprop);

    MEM_SAFE_FREE(uidata->description);

    PropertyRNA *prop = RNA_struct_type_find_property(srna, ch->def->idname);
    BLI_assert(prop);

    uidata->description = BLI_strdup(RNA_property_description(prop));

    PropertySubType prop_subtype = RNA_property_subtype(prop);

    switch (ch->type) {
      case BRUSH_CHANNEL_TYPE_FLOAT: {
        IDPropertyUIDataFloat *uidataf = reinterpret_cast<IDPropertyUIDataFloat *>(uidata);

        float min, max, soft_min, soft_max, step, precision;

        RNA_property_float_range(nullptr, prop, &min, &max);
        RNA_property_float_ui_range(nullptr, prop, &soft_min, &soft_max, &step, &precision);

        uidataf->min = (float)min;
        uidataf->max = (float)max;
        uidataf->soft_min = (float)soft_min;
        uidataf->soft_max = (float)soft_max;
        uidataf->step = (float)step;
        uidataf->precision = (int)precision;
        break;
      }
      case BRUSH_CHANNEL_TYPE_INT: {
        IDPropertyUIDataInt *uidatai = reinterpret_cast<IDPropertyUIDataInt *>(uidata);

        RNA_property_int_range(nullptr, prop, &uidatai->min, &uidatai->max);
        RNA_property_int_ui_range(
            nullptr, prop, &uidatai->soft_min, &uidatai->soft_max, &uidatai->step);
        break;
      }
      case BRUSH_CHANNEL_TYPE_ENUM:
      case BRUSH_CHANNEL_TYPE_BITMASK:
        break;
      case BRUSH_CHANNEL_TYPE_BOOL: {
        IDPropertyUIDataInt *uidatai = reinterpret_cast<IDPropertyUIDataInt *>(uidata);

        uidatai->min = uidatai->soft_min = 0;
        uidatai->max = uidatai->soft_max = 1;
        uidatai->step = 1;

        break;
      }
    }

    uidata->rna_subtype = prop_subtype;
  }

  BKE_libblock_free_data(&defaults.id, false);
}

void _BKE_brush_channelset_mark_update(BrushChannelSet *chset, const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);

  if (ch) {
    ch->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;
  }
}

void BKE_brush_channelset_ui_order_check(BrushChannelSet *chset)
{
  Vector<BrushChannel *> channels;

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    channels.append(ch);
  }

  auto cmp = [](const BrushChannel *ch1, const BrushChannel *ch2) {
    return ch1->ui_order < ch2->ui_order;
  };

  std::sort(channels.begin(), channels.end(), cmp);

  for (int i = 0; i < channels.size(); i++) {
    channels[i]->ui_order = i;
  }
}

void BKE_brush_channelset_ui_order_move(BrushChannelSet *chset,
                                        BrushChannel *ch,
                                        int uiflag,
                                        int dir)
{
  Vector<BrushChannel *> channels;

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    channels.append(ch);
  }

  auto cmp = [](const BrushChannel *ch1, const BrushChannel *ch2) {
    return ch1->ui_order < ch2->ui_order;
  };

  std::sort(channels.begin(), channels.end(), cmp);

  const char *cat1 = BKE_brush_channel_category_get(ch);

  for (int i = 0; i < channels.size(); i++) {
    if (channels[i] == ch) {
      int j = i;
      BrushChannel *ch2 = NULL;
      const char *cat2;

      do {
        j += dir < 0 ? -1 : 1;

        if (j < 0 || j >= channels.size()) {
          break;
        }

        ch2 = channels[j];
        cat2 = BKE_brush_channel_category_get(ch2);
      } while ((!(ch2->ui_flag & uiflag) || !STREQ(cat1, cat2)));

      int neworder;

      if (ch2) {
        neworder = ch2->ui_order;
        ch2->ui_order = ch->ui_order;
      }
      else {
        neworder = dir < 0 ? 0 : channels.size();
      }

      ch->ui_order = neworder;
    }
  }

  BKE_brush_channelset_ui_order_check(chset);
}

#if 0
void BKE_brush_channels_update(Brush *active_brush, Scene *scene)
{
  if (!scene->toolsettings) {
    return;
  }

  /* Sync evaluated inheritance flags */
  LISTBASE_FOREACH (BrushChannel *, ch1, &active_brush->channels->channels) {
    BrushChannel *ch2 = _BKE_brush_channelset_lookup(scene->toolsettings->unified_channels,
                                                     ch1->idname);

    ch1->evaluated_flag = ch1->flag & ~(BRUSH_CHANNEL_INHERIT | BRUSH_CHANNEL_FORCE_INHERIT);
    ch2->evaluated_flag = ch2->flag & ~(BRUSH_CHANNEL_INHERIT | BRUSH_CHANNEL_FORCE_INHERIT);

    if (ch1->flag & BRUSH_CHANNEL_INHERIT) {
      ch1->evaluated_flag |= BRUSH_CHANNEL_INHERIT;
      ch2->evaluated_flag |= BRUSH_CHANNEL_FORCE_INHERIT;
    }

    if (ch2->flag & BRUSH_CHANNEL_FORCE_INHERIT) {
      ch1->evaluated_flag |= BRUSH_CHANNEL_INHERIT;
      ch2->evaluated_flag |= BRUSH_CHANNEL_FORCE_INHERIT;
    }
  }
}
#endif

int _BKE_brush_int_get_unified(const struct Scene *scene,
                               const struct Brush *brush,
                               const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(brush->channels, idname);

  if (!ch) {
    printf("Unknown channel %s\n", idname);
    return 0;
  }

  if (BKE_brush_channel_inherits(brush, scene->toolsettings, ch)) {
    return _BKE_brush_int_get(&scene->id, scene->toolsettings->unified_channels, idname);
  }
  else {
    return _BKE_brush_int_get(&brush->id, brush->channels, idname);
  }
}

void _BKE_brush_int_set_unified(struct Scene *scene,
                                struct Brush *brush,
                                const char *idname,
                                int i)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(brush->channels, idname);

  if (!ch) {
    printf("Unknown channel %s\n", idname);
    return;
  }

  if (BKE_brush_channel_inherits(brush, scene->toolsettings, ch)) {
    _BKE_brush_int_set(&scene->id, scene->toolsettings->unified_channels, idname, i);
  }
  else {
    _BKE_brush_int_set(&brush->id, brush->channels, idname, i);
  }
}

float _BKE_brush_float_get_unified(const struct Scene *scene,
                                   const struct Brush *brush,
                                   const char *idname)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(brush->channels, idname);

  if (!ch) {
    printf("Unknown channel %s\n", idname);
    return 0.0f;
  }

  if (BKE_brush_channel_inherits(brush, scene->toolsettings, ch)) {
    return _BKE_brush_float_get(&scene->id, scene->toolsettings->unified_channels, idname);
  }
  else {
    return _BKE_brush_float_get(&brush->id, brush->channels, idname);
  }
}

void _BKE_brush_float_set_unified(struct Scene *scene,
                                  struct Brush *brush,
                                  const char *idname,
                                  float f)
{
  BrushChannel *ch = _BKE_brush_channelset_lookup(brush->channels, idname);

  if (!ch) {
    printf("Unknown channel %s\n", idname);
    return;
  }

  if (BKE_brush_channel_inherits(brush, scene->toolsettings, ch)) {
    _BKE_brush_float_set(&scene->id, scene->toolsettings->unified_channels, idname, f);
  }
  else {
    _BKE_brush_float_set(&brush->id, brush->channels, idname, f);
  }
}

static bool brush_mapping_inherits(BrushChannel *owner,
                                   BrushChannel *unified,
                                   BrushMapping *mapping)
{
  bool inherit_ch = owner->flag & BRUSH_CHANNEL_INHERIT;

  if ((unified->flag & BRUSH_CHANNEL_FORCE_INHERIT) &&
      !(owner->flag & BRUSH_CHANNEL_IGNORE_FORCE_INHERIT)) {
    inherit_ch = true;
  }

  bool inherit = mapping->inherit_mode == BRUSH_MAPPING_INHERIT_ALWAYS;
  inherit = inherit || (mapping->inherit_mode == BRUSH_MAPPING_INHERIT_CHANNEL && inherit_ch);

  return inherit;
}
bool _BKE_brush_mapping_enabled(const struct Scene *scene,
                                const struct Brush *brush,
                                const char *idname,
                                eBrushMappingType mapping_type)
{
  BrushChannel *ch1 = _BKE_brush_channelset_lookup(brush->channels, idname);
  BrushChannel *ch2 = _BKE_brush_channelset_lookup(scene->toolsettings->unified_channels, idname);

  if (!ch1) {
    return false;
  }

  if (brush_mapping_inherits(ch1, ch2, ch1->mappings + mapping_type)) {
    return ch2->mappings[mapping_type].flag & BRUSH_MAPPING_ENABLED;
  }
  else {
    return ch1->mappings[mapping_type].flag & BRUSH_MAPPING_ENABLED;
  }
}
