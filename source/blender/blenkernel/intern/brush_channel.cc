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
#include "BLO_read_write.h"

#include <string>

#define BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES
#include "brush_channel_define.h"
#undef BRUSH_CHANNEL_DEFINE_INTERNAL_NAMES

using string = std::string;

namespace blender {

using KeyString = const char *;
using StringRef = blender::StringRef;

static Map<StringRef, BrushChannelType> builtin_channels;

using ChannelNameMap = Map<StringRef, BrushChannel *>;

static void init_builtin_brush_channels()
{
  struct ChannelProp {
    char path[512];
    char category[64];
    int flag;
    BrushMappingPreset mappings;
  };

#define BRUSH_CHANNEL_DEFINE_INTERNAL
  ChannelProp channel_props[] = {
#include "brush_channel_define.h"
  };

  StructRNA *srna = RNA_struct_find("Brush");

  Brush dummy = {0};
  PointerRNA _ptr = {&dummy.id}, *ptr = &_ptr;

  printf("total properties: %i\n", (int)ARRAY_SIZE(channel_props));

  for (int i : IndexRange(ARRAY_SIZE(channel_props))) {
    ChannelProp &def = channel_props[i];
    BrushChannelType type;

    BLI_strncpy(type.idname, def.path, sizeof(type.idname));
    BLI_strncpy(type.category, def.category, sizeof(type.category));
    type.flag = def.flag;

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

    printf("  idname: %s\n", type.idname);
    builtin_channels.add(strdup(type.idname), type);
    printf("    success: %d\n", builtin_channels.contains(type.idname) ? 1 : 0);
  }
#undef BRUSH_CHANNEL_DEFINE_INTERNAL
}
static void check_builtin_brush_channels()
{
  printf("size: %d %d\n", (int)builtin_channels.size(), (int)(builtin_channels.size() == 0));
  for (auto &key : builtin_channels.keys()) {
    printf("key: %s\n", key.data());
  }

  if (builtin_channels.size() == 0) {
    init_builtin_brush_channels();
  }
}

ChannelNameMap *get_namemap(BrushChannelSet *chset)
{
  return reinterpret_cast<ChannelNameMap *>(chset->channelmap);
}

extern "C" BrushChannelSet *BKE_brush_channelset_create()
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

extern "C" void BKE_brush_channelset_free(BrushChannelSet *chset)
{
  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BKE_brush_channel_free_data(ch);
  }

  BLI_freelistN(&chset->channels);

  MEM_delete<ChannelNameMap>(get_namemap(chset));
  MEM_freeN(chset->channelmap);
  MEM_freeN(static_cast<void *>(chset));
}

static void brush_mapping_reset(BrushMapping *mp, int type)
{
  mp->type = type;
  mp->premultiply_factor = 1.0f;
  mp->min = 0.0f;
  mp->max = 1.0f;
  mp->blendmode = MA_RAMP_MULT;
  mp->curve.preset = BRUSH_CURVE_LIN;
}

extern "C" BrushChannel *BKE_brush_channelset_add(BrushChannelSet *chset, BrushChannelType *type)
{
  BrushChannel *ch = MEM_cnew<BrushChannel>("BrushChannel");

  BLI_strncpy(ch->idname, type->idname, sizeof(ch->idname));
  BLI_strncpy(ch->uiname, type->uiname, sizeof(ch->uiname));

  ch->def = type;
  ch->type = type->type;
  ch->flag = type->flag;

  BLI_addtail(&chset->channels, static_cast<void *>(ch));
  get_namemap(chset)->add(ch->idname, ch);

  for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
    BrushMapping *mp = ch->mappings + i;
    brush_mapping_reset(mp, i);
  }

  return ch;
}

extern "C" void BKE_brush_channelset_begin(BrushChannelSet *chset, BrushChannelType *type)
{
  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    ch->flag |= BRUSH_CHANNEL_NEEDS_EVALUATE;
  }
}

extern "C" BrushChannel *_BKE_brush_channelset_lookup(BrushChannelSet *chset, const char *idname)
{
  return get_namemap(chset)->lookup(idname);
}

extern "C" bool _BKE_brush_channelset_has(BrushChannelSet *chset, const char *idname)
{
  return get_namemap(chset)->contains(idname);
}

extern "C" const char *BKE_brush_channel_category_get(BrushChannel *ch)
{
  return ch->category ? ch->category : ch->def->category;
}

extern "C" const void BKE_brush_channel_category_set(BrushChannel *ch, const char *category)
{
  if (STREQ(category, ch->def->category)) {
    MEM_SAFE_FREE(ch->category);
    ch->category = nullptr;

    return;
  }

  ch->category = BLI_strdup(category);
}

ATTR_NO_OPT extern "C" BrushChannel *_BKE_brush_channelset_ensure(BrushChannelSet *chset,
                                                                  const char *idname)
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

extern "C" void BKE_brush_channelset_ensure_channels(BrushChannelSet *chset, int sculpt_tool)
{
  check_builtin_brush_channels();

  BKE_brush_channelset_ensure(chset, size);
  BKE_brush_channelset_ensure(chset, unprojected_radius);
  BKE_brush_channelset_ensure(chset, strength);

  switch (sculpt_tool) {
  }
}

extern "C" char *BKE_brush_channel_rna_path(ID *owner, BrushChannel *ch)
{
  switch (GS(owner->name)) {
    case ID_BR:
      return BLI_strdup(ch->idname);
    case ID_SCE: {
      string path = "tool_settings.channel_properties[\"";

      path += string(ch->idname) + "\"]";

      return BLI_strdup(path.c_str());
    }
    default:
      return BLI_strdup("");
  }
}

void _brush_channel_ensure_value(Brush *brush,
                                 Scene *scene,
                                 BrushChannelSet *chset,
                                 const char *idname)
{
  ToolSettings *tool_settings = scene->toolsettings;

  BrushChannel *dest_ch = _BKE_brush_channelset_lookup(chset, idname);
  BrushChannel *brush_ch = _BKE_brush_channelset_lookup(brush->channels, idname);
  BrushChannel *scene_ch = _BKE_brush_channelset_lookup(tool_settings->unified_channels, idname);

  if (!(dest_ch->flag & BRUSH_CHANNEL_NEEDS_EVALUATE)) {
    return;
  }

  string path;
  StructRNA *srna;
  ID *id;

  const char *rnaname = dest_ch->idname;
  bool inherit = BKE_brush_channel_inherits(brush, scene->toolsettings, dest_ch);

  if (!inherit) {
    path = "";
    id = &brush->id;
    srna = RNA_struct_find("Brush");
  }
  else {
    path = "tool_settings.unified_properties[\"";
    id = &scene->id;
    srna = RNA_struct_find("Scene");
  }

  path += string(rnaname);

  if (inherit) {
    path += string("\"]");
  }

  dest_ch->flag &= ~BRUSH_CHANNEL_NEEDS_EVALUATE;
  PointerRNA ptr, ptr2;
  PropertyRNA *prop = nullptr;

  RNA_pointer_create(id, srna, nullptr, &ptr);
  if (RNA_path_resolve(&ptr, path.c_str(), &ptr2, &prop)) {
    PropertyType prop_type = RNA_property_type(prop);
    PropertySubType prop_subtype = RNA_property_subtype(prop);

    switch (prop_type) {
      case PROP_FLOAT:
        if (RNA_property_array_check(prop)) {
          RNA_float_get_array(&ptr2, rnaname, dest_ch->vector);
        }
        else {
          dest_ch->fvalue = RNA_float_get(&ptr2, rnaname);
        }
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

#define brush_channel_ensure_value(br, sd, chset, idname) \
  _brush_channel_ensure_value(br, sd, chset, make_builtin_ch_name(idname));

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

float _BKE_brush_channelset_get_float(
    Brush *br, Scene *scene, BrushChannelSet *chset, const char *idname, BrushMappingData *mapping)
{
  _brush_channel_ensure_value(br, scene, chset, idname);
  BrushChannel *ch = _BKE_brush_channelset_lookup(chset, idname);

  return BKE_brush_channel_eval_mappings(ch, mapping, ch->fvalue, 0);
}

void _BKE_brush_channelset_float_set(struct Brush *br,
                                     struct Sculpt *sd,
                                     BrushChannelSet *chset,
                                     const char *idname,
                                     float f,
                                     bool set_rna)
{
}

void BKE_brush_mapping_copy_data(BrushMapping *dest, BrushMapping *src)
{
  *dest = *src;

  if (dest->curve.curve) {
    dest->curve.curve = BKE_curvemapping_copy(dest->curve.curve);
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
  *dest = *src;

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

extern "C" BrushChannelSet *BKE_brush_channelset_copy(BrushChannelSet *chset)
{
  ChannelNameMap *namemap = get_namemap(chset);

  BrushChannelSet *chset2 = BKE_brush_channelset_create();
  ChannelNameMap *namemap2 = get_namemap(chset2);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BrushChannel *ch2 = BKE_brush_channel_copy(ch);

    namemap2->add(ch2->idname, ch2);
    BLI_addtail(&chset->channels, static_cast<void *>(ch2));
  }

  return chset2;
}

extern "C" void BKE_brush_channelset_blend_read(BrushChannelSet *chset, BlendDataReader *reader)
{
  BLO_read_list(reader, &chset->channels);

  ChannelNameMap *namemap = MEM_new<ChannelNameMap>("ChannelNameMap");
  chset->channelmap = static_cast<void *>(namemap);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
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

extern "C" bool BKE_brush_channel_inherits(Brush *brush,
                                           ToolSettings *tool_settings,
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

BrushChannelSet *BKE_brush_channelset_create_final(Brush *brush,
                                                   ToolSettings *tool_settings,
                                                   BrushMappingData *mapdata)
{
  BrushChannelSet *chset = BKE_brush_channelset_copy(brush->channels);

  LISTBASE_FOREACH (BrushChannel *, ch, &chset->channels) {
    BrushChannel *ch1 = _BKE_brush_channelset_lookup(brush->channels, ch->idname);
    BrushChannel *ch2 = _BKE_brush_channelset_lookup(tool_settings->unified_channels, ch->idname);
    bool inherit = BKE_brush_channel_inherits(brush, tool_settings, ch);

    if (inherit) {
      BKE_brush_channel_copy_data(ch, ch2);
    }

    for (int i = 0; i < BRUSH_MAPPING_MAX; i++) {
      BrushMapping *mp1 = ch1->mappings + i;
      BrushMapping *mp2 = ch2->mappings + i;

      if ((mp1->flag & BRUSH_MAPPING_INHERIT_NEVER) && inherit) {
        BKE_brush_mapping_free_data(ch->mappings + i);
        BKE_brush_mapping_copy_data(ch->mappings + i, mp1);
      }
      else if ((mp1->flag & BRUSH_MAPPING_INHERIT_ALWAYS) && !inherit) {
        BKE_brush_mapping_free_data(ch->mappings + i);
        BKE_brush_mapping_copy_data(ch->mappings + i, mp2);
      }
    }
  }

  return chset;
}

extern "C" void BKE_brush_channelset_toolsettings_init(ToolSettings *ts)
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

      if (RNA_path_resolve(&ptr, ch->idname, &prop_ptr, &prop)) {
        switch (RNA_property_type(prop)) {
          case PROP_INT:
            default_value = RNA_int_get(&prop_ptr, ch->idname);
            break;
          case PROP_FLOAT:
            default_value = RNA_float_get(&prop_ptr, ch->idname);
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
        default:
          printf(
              "%s: unsupported brush channel type for unified channel: %d\n", __func__, ch->type);
          continue;
      }

      idprop = IDP_New(type, &tmpl, ch->idname);
      IDP_AddToGroup(ts->unified_properties, idprop);
    }

    IDPropertyUIData *uidata = IDP_ui_data_ensure(idprop);

    MEM_SAFE_FREE(uidata->description);
    uidata->description = BLI_strdup(ch->def->tooltip);

    PropertyRNA *prop = RNA_struct_type_find_property(srna, ch->def->idname);
    BLI_assert(prop);

    PropertyType prop_type = RNA_property_type(prop);
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
      case BRUSH_CHANNEL_TYPE_BOOL: {
        break;
      }
    }

    uidata->rna_subtype = prop_subtype;
  }

  BKE_libblock_free_data(&defaults.id, false);
}

}  // namespace blender
