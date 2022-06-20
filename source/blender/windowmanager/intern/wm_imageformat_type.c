/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup wm
 *
 * Operator Registry.
 */

#include "MEM_guardedalloc.h"

#include "CLG_log.h"

#include "DNA_ID.h"
#include "DNA_scene_types.h"
#include "DNA_screen_types.h"
#include "DNA_userdef_types.h"
#include "DNA_windowmanager_types.h"

#include "BLT_translation.h"

#include "BLI_blenlib.h"
#include "BLI_ghash.h"
#include "BLI_utildefines.h"

#include "BKE_context.h"
#include "BKE_idprop.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_enum_types.h"
#include "RNA_prototypes.h"

#include "WM_api.h"
#include "WM_types.h"

#include "wm.h"
#include "wm_event_system.h"

#define UNDOCUMENTED_IMAGEFORMAT_TIP N_("(undocumented image format)")

/* -------------------------------------------------------------------- */
/** \name Image Format Type Registry
 * \{ */

static GHash *global_imageformat_hash = NULL;
// TODO: no idea what this is for
/** Counter for operator-properties that should not be tagged with #OP_PROP_TAG_ADVANCED. */
static int ot_prop_basic_count = -1;

struct imfImageFormatType *WM_imageformattype_find(const char *idname, bool quiet)
{
  if (idname[0]) {
    imfImageFormatType *ift;

    /* needed to support python style names without the _OT_ syntax */
    char idname_bl[OP_MAX_TYPENAME];
    WM_operator_bl_idname(idname_bl, idname); // TODO: remove this or reimplement?? not sure if needed yet.

    ift = BLI_ghash_lookup(global_imageformat_hash, idname_bl);
    if (ift) {
      return ift;
    }

    if (!quiet) {
      CLOG_INFO(
          WM_LOG_OPERATORS, 0, "search for unknown image format '%s', '%s'\n", idname_bl, idname);
    }
  }
  else {
    if (!quiet) {
      CLOG_INFO(WM_LOG_OPERATORS, 0, "search for empty image format");
    }
  }

  return NULL;
}

void wm_imageformattype_iter(GHashIterator *ghi)
{
  BLI_ghashIterator_init(ghi, global_imageformat_hash);
}

/* -------------------------------------------------------------------- */
/** \name Image Format Type Append
 * \{ */

static imfImageFormatType *wm_imageformattype_append__begin(void)
{
  imfImageFormatType *ift = MEM_callocN(sizeof(imfImageFormatType), "imageformattype");

  ift->srna = RNA_def_struct_ptr(&BLENDER_RNA, "", &RNA_ImageFormat); // TODO: wrong type?? use other function
  /* Set the default i18n context now, so that opfunc can redefine it if needed! */
  RNA_def_struct_translation_context(ift->srna, BLT_I18CONTEXT_IMAGEFORMAT_DEFAULT);
  ift->translation_context = BLT_I18CONTEXT_IMAGEFORMAT_DEFAULT;

  return ift;
}

static void wm_imageformattype_append__end(imfImageFormatType *ift)
{
  if (ift->name == NULL) {
    CLOG_ERROR(WM_LOG_OPERATORS, "ImageFormat '%s' has no name property", ift->idname); // TODO: change error clgref
  }
  BLI_assert((ift->description == NULL) || (ift->description[0]));

  /* XXX All ops should have a description but for now allow them not to. */
  RNA_def_struct_ui_text(
      ift->srna, ift->name, ift->description ? ift->description : UNDOCUMENTED_IMAGEFORMAT_TIP);
  RNA_def_struct_identifier(&BLENDER_RNA, ift->srna, ift->idname);

  BLI_ghash_insert(global_imageformat_hash, (void *)ift->idname, ift);
}

/* All ops in 1 list (for time being... needs evaluation later). */

void wm_imageformattype_append(void (*opfunc)(imfImageFormatType *))
{
  imfImageFormatType *ift = wm_imageformattype_append__begin();
  opfunc(ift);
  wm_imageformattype_append__end(ift);
}

void wm_imageformattype_append_ptr(void (*opfunc)(imfImageFormatType *, void *), void *userdata)
{
  imfImageFormatType *ift = wm_imageformattype_append__begin();
  opfunc(ift, userdata);
  wm_imageformattype_append__end(ift);
}

/** \} */

void wm_imageformattype_remove_ptr(imfImageFormatType *ift)
{
  BLI_assert(ift == WM_imageformattype_find(otiftidname, false));

  RNA_struct_free(&BLENDER_RNA, ift->srna);

  BLI_ghash_remove(global_imageformat_hash, ift->idname, NULL, NULL);

  MEM_freeN(ift);
}

bool wm_imageformattype_remove(const char *idname)
{
  imfImageFormatType *ift = WM_imageformattype_find(idname, 0);

  if (ift == NULL) {
    return false;
  }

  wm_imageformattype_remove_ptr(ift);

  return true;
}

void wm_imageformattype_init(void)
{
  /* reserve size is set based on blender default setup */
  global_imageformat_hash = BLI_ghash_str_new_ex("wm_imageformattype_init gh", 64);
}

static void imageformattype_ghash_free_cb(imfImageFormatType *ift)
{
  if (ift->rna_ext.srna) {
    /* python image format, allocs own string */
    // TODO: delete more of these things that alloc their own strings
    MEM_freeN((void *)ift->idname);
  }

  MEM_freeN(ift);
}

void wm_imageformattype_free(void)
{
  BLI_ghash_free(global_imageformat_hash, NULL, (GHashValFreeFP)imageformattype_ghash_free_cb);
  global_imageformat_hash = NULL;
}

// TODO: probably remove
void wm_imageformattype_idname_visit_for_search(const bContext *UNUSED(C),
                                             PointerRNA *UNUSED(ptr),
                                             PropertyRNA *UNUSED(prop),
                                             const char *UNUSED(edit_text),
                                             StringPropertySearchVisitFunc visit_fn,
                                             void *visit_user_data)
{
  GHashIterator gh_iter;
  GHASH_ITER (gh_iter, global_imageformat_hash) {
    imfImageFormatType *ift = BLI_ghashIterator_getValue(&gh_iter);

    char idname_py[OP_MAX_TYPENAME];
    WM_operator_py_idname(idname_py, ift->idname);

    StringPropertySearchVisitParams visit_params = {NULL};
    visit_params.text = idname_py;
    visit_params.info = ift->name;
    visit_fn(visit_user_data, &visit_params);
  }
}

/** \} */
