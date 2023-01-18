/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup asset_system
 */

#include "BKE_appdir.h"

#include "AS_asset_bundled.hh"

namespace blender::asset_system {

StringRefNull bundled_assets_directory_path()
{
  static std::string path = []() {
    const char *datafiles_path = BKE_appdir_folder_id(BLENDER_DATAFILES, "assets");
    return datafiles_path;
  }();
  return path;
}

}  // namespace blender::asset_system
