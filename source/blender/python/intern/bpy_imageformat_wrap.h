/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup pythonintern
 */

#pragma once

struct imfImageFormatType;

#ifdef __cplusplus
extern "C" {
#endif

/* exposed to rna/wm api */
/**
 * Generic function used by all Python defined image formats
 * it's passed as an argument to #WM_imageformattype_append_ptr in for image format registration.
 */
void BPY_RNA_imageformat_wrapper(struct imfImageFormatType *ift, void *userdata);

#ifdef __cplusplus
}
#endif
