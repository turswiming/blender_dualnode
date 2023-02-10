/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup asset_system
 */

#include "BLI_path_util.h"

#include "BKE_appdir.h"

#include "AS_essentials_library.hh"

namespace blender::asset_system {

StringRefNull essentials_directory_path()
{
  static std::string path = []() {
    const char *datafiles_path = BKE_appdir_folder_id(BLENDER_DATAFILES, "assets");
    return datafiles_path;
  }();
  return path;
}

}  // namespace blender::asset_system

bool ED_asset_essentials_contains_path(const char *path)
{
  const blender::StringRefNull bundled_path = blender::asset_system::essentials_directory_path();
  return BLI_path_contains(bundled_path.c_str(), path);
}
