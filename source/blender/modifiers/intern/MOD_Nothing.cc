//
// Created by terrapene on 2023/4/23.
//
/** \file
 * \ingroup modifiers
 */

#include <intern/RNA_blender.h>

#include "BKE_lib_query.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_volume.h"



#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_object_types.h"
#include "DNA_screen_types.h"
#include "DNA_volume_types.h"

#include "UI_interface.h"
#include "UI_resources.h"
#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"
#include "RNA_access.h"

#include "BLI_math_vector.h"
#include "BLI_span.hh"
#include "BLI_timeit.hh"



static void initData(ModifierData *md)
{
  NothingModifierData *vmmd = (NothingModifierData*)(md);
}

static void updateDepsgraph(struct ModifierData *md,
                            const ModifierUpdateDepsgraphContext *pContext)
{
  NothingModifierData *nmd = (NothingModifierData*)(md);

}

static void foreachIDLink(void *md)
{
  NothingModifierData *nmd = (NothingModifierData*)(md);
}

static void panel_draw(const bContext *UNUSED, Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, NULL);

  uiLayoutSetPropSep(layout, true);
  uiLayoutSetPropDecorate(layout, false);

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  PanelType *panel_type = modifier_panel_register(region_type, eModifierType_Nothing, panel_draw);
}



static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *ctx, Mesh *input_mesh)
{



return input_mesh;

}
static void deformVerts(
    struct ModifierData *pData, const struct ModifierEvalContext *pContext,struct Mesh *pMesh, float (*pDouble)[3], int i){
  int vsize =pMesh->totvert;
  for (int i = 0; i < vsize; ++i) {
    pDouble[i][0]*=2;
    pDouble[i][1]*=2;

    pDouble[i][2]*=2;
  }
}
ModifierTypeInfo modifierType_Nothing = {
    /* name */ "Nothing",
    /* structName */ "NothingModifierData",
    /* structSize */ sizeof(NothingModifierData),
    /* srna */ &RNA_NothingModifier,
    /* type */ eModifierTypeType_OnlyDeform,
    /* flags */ eModifierTypeFlag_AcceptsMesh,
    /* icon */ ICON_VOLUME_DATA, /* TODO: Use correct icon. */

    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ deformVerts,
    /* deformMatrices */ NULL,
    /* deformVertsEM */ NULL,
    /* deformMatricesEM */ NULL,
    /* modifyMesh */ modifyMesh,
    /* modifyGeometrySet= */ NULL,
    /* initData= */ NULL,
    /* requiredDataMask= */ NULL,

    /* initData= */ initData,
    /* isDisabled= */ NULL,
    /* updateDepsgraph= */ updateDepsgraph,
    /* dependsOnTime= */ NULL,
    /* dependsOnNormals= */ NULL,
    /* foreachIDLink= */ NULL,
    /* foreachTexLink= */ NULL,
    /* freeRuntimeData= */ foreachIDLink,
    /* panelRegister= */ panelRegister,
    /* blendWrite= */ NULL,
    /* blendRead= */ NULL,
};
