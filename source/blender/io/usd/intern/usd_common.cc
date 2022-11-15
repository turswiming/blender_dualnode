/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2021 Blender Foundation. All rights reserved. */

#include "usd_common.h"
#include "usd.h"

#include <pxr/usd/ar/resolver.h>
#include <pxr/base/plug/registry.h>

#include "BKE_appdir.h"
#include "BLI_path_util.h"
#include "BLI_string.h"

#include "WM_api.h"
#include "WM_types.h"

namespace blender::io::usd {

void ensure_usd_plugin_path_registered()
{
  /* if PXR_PYTHON_SUPPORT_ENABLED is defined, we *must* be dynamic and
     the plugins are placed relative to the USD shared library hence no
     hinting is required. */
#ifndef PXR_PYTHON_SUPPORT_ENABLED
  static bool plugin_path_registered = false;
  if (plugin_path_registered) {
    return;
  }
  plugin_path_registered = true;

  /* Tell USD which directory to search for its JSON files. If 'datafiles/usd'
   * does not exist, the USD library will not be able to read or write any files. */
  const std::string blender_usd_datafiles = BKE_appdir_folder_id(BLENDER_DATAFILES, "usd");
  /* The trailing slash indicates to the USD library that the path is a directory. */
  pxr::PlugRegistry::GetInstance().RegisterPlugins(blender_usd_datafiles + "/");
#endif
}

}  // namespace blender::io::usd


void USD_path_abs(char *path, const char *basepath, bool for_import)
{
  if (BLI_path_is_rel(path)) {
    /* The path is relative to the .blend file, so use Blender's
     * standard absolute path resolution. */
    BLI_path_abs(path, basepath);
    /* Use forward slash separators, which is standard for USD. */
    BLI_str_replace_char(path, SEP, ALTSEP);
  }

  pxr::ArResolvedPath resolved_path = for_import ? pxr::ArGetResolver().Resolve(path) :
                                                   pxr::ArGetResolver().ResolveForNewAsset(path);

  std::string path_str = resolved_path.GetPathString();

  if (!path_str.empty()) {
    if (path_str.length() < FILE_MAX) {
      BLI_strncpy(path, path_str.c_str(), FILE_MAX);
      /* Use forward slash separators, which is standard for USD. */
      BLI_str_replace_char(path, SEP, ALTSEP);
    }
    else {
      WM_reportf(RPT_ERROR,
                 "USD_path_abs: resolved path %s exceeds path buffer length.",
                 path_str.c_str());
    }
  }
}
