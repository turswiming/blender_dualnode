/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup pythonintern
 *
 * This file exposes functionality for defining to define operators that C can call into.
 * The generic callback functions for python operators are defines in
 * 'rna_wm.c', some calling into functions here to do python specific
 * functionality.
 */

#include <Python.h>

#include "BLI_utildefines.h"

#include "WM_api.h"
#include "WM_types.h"

#include "RNA_access.h"
#include "RNA_define.h"
#include "RNA_prototypes.h"

#include "bpy_intern_string.h"
#include "bpy_imageformat_wrap.h" /* own include */
#include "bpy_rna.h"

void BPY_RNA_imageformat_wrapper(imfImageFormatType *ift, void *userdata)
{
  /* take care not to overwrite anything set in
   * WM_imageformattype_append_ptr before opfunc() is called */
  StructRNA *srna = ift->srna;
  *ift = *((imfImageFormatType *)userdata);
  ift->srna = srna; /* restore */

  /* Use i18n context from rna_ext.srna if possible (py image formats). */
  if (ift->rna_ext.srna) {
    RNA_def_struct_translation_context(ift->srna,
                                       RNA_struct_translation_context(ift->rna_ext.srna));
  }
}
