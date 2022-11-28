/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 NVIDIA Corportation. All rights reserved. */

#include "usd_asset_utils.h"

#include <pxr/usd/ar/asset.h>
#include <pxr/usd/ar/packageUtils.h>
#include <pxr/usd/ar/resolver.h>
#include <pxr/usd/ar/writableAsset.h>

#include "BKE_main.h"

#include "BLI_fileops.h"
#include "BLI_path_util.h"
#include "BLI_string.h"

#include "WM_api.h"
#include "WM_types.h"

#include <iostream>

static const char UDIM_PATTERN[] = "<UDIM>";
static const char UDIM_PATTERN2[] = "%3CUDIM%3E";
static const int UDIM_START_TILE = 1001;
static const int UDIM_END_TILE = 1100;

namespace blender::io::usd {

static std::string get_asset_file_name(const char *src_path)
{
  char file_name[FILE_MAX];

  if (pxr::ArIsPackageRelativePath(src_path)) {
    std::pair<std::string, std::string> split = pxr::ArSplitPackageRelativePathInner(src_path);
    if (split.second.empty()) {
      WM_reportf(
          RPT_WARNING,
          "usd_import_texture(): Couldn't determine package-relative file name from path %s.",
          src_path);
      return src_path;
    }
    BLI_split_file_part(split.second.c_str(), file_name, FILE_MAX);
  }
  else {
    BLI_split_file_part(src_path, file_name, FILE_MAX);
  }

  return file_name;
}

/* Return true if the given path is an existing director
 * on the standard file system. */
static bool fs_parent_dir_exists(const char *path)
{
  char dir_path[FILE_MAX];
  BLI_split_dir_part(path, dir_path, FILE_MAX);
  bool is_dir = BLI_is_dir(dir_path);
  return is_dir;
}

static bool is_udim_path(const std::string &path)
{
  return path.find(UDIM_PATTERN) != std::string::npos ||
         path.find(UDIM_PATTERN2) != std::string::npos;
}

/* The following is copied from _SplitUdimPattern() in
 * USD library source file materialParamsUtils.cpp.
 * Split a udim file path such as /someDir/myFile.<UDIM>.exr into a
 * prefix (/someDir/myFile.) and suffix (.exr). */
static std::pair<std::string, std::string> split_udim_pattern(const std::string &path)
{
  static const std::vector<std::string> patterns = {UDIM_PATTERN, UDIM_PATTERN2};

  for (const std::string &pattern : patterns) {
    const std::string::size_type pos = path.find(pattern);
    if (pos != std::string::npos) {
      return {path.substr(0, pos), path.substr(pos + pattern.size())};
    }
  }

  return {std::string(), std::string()};
}

static std::string copy_asset_to_directory(const char *src_path, const char *dest_dir_path, bool overwrite)
{
  std::string file_name = get_asset_file_name(src_path);

  char dest_file_path[FILE_MAX];
  BLI_path_join(dest_file_path, FILE_MAX, dest_dir_path, file_name.c_str());
  BLI_str_replace_char(dest_file_path, SEP, ALTSEP);

  if (!overwrite && BLI_is_file(dest_file_path)) {
    return dest_file_path;
  }

  if (!copy_usd_asset(src_path, dest_file_path, overwrite)) {
    WM_reportf(
        RPT_WARNING, "usd_import_texture(): Couldn't copy file %s to %s.", src_path, dest_file_path);
    return src_path;
  }

  WM_reportf(RPT_INFO, "usd_import_texture(): Copied file %s to %s.", src_path, dest_file_path);

  return dest_file_path;
}

static std::string copy_udim_asset_to_directory(const char *src_path,
                                                const char *dest_dir_path,
                                                bool overwrite)
{
  /* Get prefix and suffix from udim pattern. */
  std::pair<std::string, std::string> splitPath = split_udim_pattern(src_path);
  if (splitPath.first.empty() || splitPath.second.empty()) {
    WM_reportf(RPT_ERROR, "copy_udim_asset_to_directory(): Couldn't split UDIM pattern %s.", src_path);
    return src_path;
  }

  for (int i = UDIM_START_TILE; i < UDIM_END_TILE; ++i) {
    const std::string src_udim = splitPath.first + std::to_string(i) + splitPath.second;
    if (usd_path_exists(src_udim.c_str())) {
      copy_asset_to_directory(src_udim.c_str(), dest_dir_path, overwrite);
    }
  }

  const std::string src_file_name = get_asset_file_name(src_path);
  char ret_udim_path[FILE_MAX];
  BLI_path_join(ret_udim_path, FILE_MAX, dest_dir_path, src_file_name.c_str());
  BLI_str_replace_char(ret_udim_path, SEP, ALTSEP);

  /* Blender only recognizes the <UDIM> pattern, not the
   * alternative UDIM_PATTERN2, so we make sure the returned
   * path has the former. */
  splitPath = split_udim_pattern(ret_udim_path);
  if (splitPath.first.empty() || splitPath.second.empty()) {
    WM_reportf(RPT_ERROR, "copy_udim_asset_to_directory(): Couldn't split UDIM pattern %s.", ret_udim_path);
    return ret_udim_path;
  }

  return splitPath.first + UDIM_PATTERN + splitPath.second;
}


bool copy_usd_asset(const char *src, const char *dst, bool overwrite)
{
  if (!(src && dst)) {
    return false;
  }

  pxr::ArResolver &ar = pxr::ArGetResolver();

  if (!overwrite) {
    if (!ar.Resolve(dst).IsEmpty()) {
      /* The asset exists, so this is a no-op. */
      WM_reportf(RPT_INFO, "copy_usd_asset(): Will not overwrite existing asset %s", dst);
      return true;
    }
  }

  pxr::ArResolvedPath src_path = ar.Resolve(src);

  if (src_path.IsEmpty()) {
    WM_reportf(RPT_ERROR, "copy_usd_asset(): Can't resolve path %s.", src);
    return false;
  }

  pxr::ArResolvedPath dst_path = ar.ResolveForNewAsset(dst);

  if (dst_path.IsEmpty()) {
    WM_reportf(RPT_ERROR, "copy_usd_asset(): Can't resolve path %s for writing.", dst);
    return false;
  }

  std::string why_not;
  if (!ar.CanWriteAssetToPath(dst_path, &why_not)) {
    WM_reportf(RPT_ERROR,
               "copy_usd_asset(): Can't write to asset %s.  %s.",
               dst_path.GetPathString().c_str(),
               why_not.c_str());
    return false;
  }

  std::shared_ptr<pxr::ArAsset> src_asset = ar.OpenAsset(src_path);
  if (!src_asset) {
    WM_reportf(
        RPT_ERROR, "copy_usd_asset(): Can't open source asset %s.", src_path.GetPathString().c_str());
    return false;
  }

  const size_t size = src_asset->GetSize();

  if (size == 0) {
    WM_reportf(RPT_WARNING,
               "copy_usd_asset(): Will not copy zero size source asset %s.",
               src_path.GetPathString().c_str());
    return false;
  }

  std::shared_ptr<const char> buf = src_asset->GetBuffer();

  if (!buf) {
    WM_reportf(RPT_ERROR,
               "copy_usd_asset(): Null buffer for source asset %s.",
               src_path.GetPathString().c_str());
    return false;
  }

  std::shared_ptr<pxr::ArWritableAsset> dst_asset = ar.OpenAssetForWrite(
      dst_path, pxr::ArResolver::WriteMode::Replace);
  if (!dst_asset) {
    WM_reportf(RPT_ERROR,
               "copy_usd_asset(): Can't open destination asset %s for writing.",
               src_path.GetPathString().c_str());
    return false;
  }

  size_t bytes_written = dst_asset->Write(src_asset->GetBuffer().get(), src_asset->GetSize(), 0);

  if (bytes_written == 0) {
    WM_reportf(RPT_ERROR,
               "copy_usd_asset(): Error writing to destination asset %s.",
               dst_path.GetPathString().c_str());
  }
  else {
    WM_reportf(RPT_INFO,
               "copy_usd_asset(): Copied %s to %s.",
               src_path.GetPathString().c_str(), dst_path.GetPathString().c_str());
  }

  if (!dst_asset->Close()) {
    WM_reportf(RPT_ERROR,
               "copy_usd_asset(): Couldn't close destination asset %s.",
               dst_path.GetPathString().c_str());
    return false;
  }

  return bytes_written > 0;
}

 std::string get_textures_dir(const pxr::UsdStageRefPtr stage)
{
   pxr::SdfLayerHandle layer = stage->GetRootLayer();
   std::string stage_path = layer->GetRealPath();
   if (stage_path.empty()) {
     return "";
   }

   char usd_dir_path[FILE_MAX];
   BLI_split_dir_part(stage_path.c_str(), usd_dir_path, FILE_MAX);

   return std::string(usd_dir_path) + std::string("textures/");
 }

 bool usd_path_exists(const char *src)
 {
   return src && !pxr::ArGetResolver().Resolve(src).IsEmpty();
 }

 bool usd_paths_equal(const char *p1, const char *p2)
 {
   BLI_assert_msg(!BLI_path_is_rel(p1) && !BLI_path_is_rel(p2),
                  "Paths arguments must be absolute");

   pxr::ArResolver &ar = pxr::ArGetResolver();

   std::string resolved_p1 = ar.ResolveForNewAsset(p1).GetPathString();
   std::string resolved_p2 = ar.ResolveForNewAsset(p2).GetPathString();

   return resolved_p1 == resolved_p2;
 }

 std::string usd_import_texture(const char *src, const char *import_dir, const bool overwrite)
 {
   const bool udim_src = is_udim_path(src);

   if (udim_src) {
     if (fs_parent_dir_exists(src)) {
       return src;
     }
   }
   else if (BLI_is_file(src)) {
     /* File exists in filesystem, no need to import. */
     return src;
   }

   if (strlen(import_dir) == 0) {
     WM_reportf(
         RPT_ERROR,
         "usd_import_texture(): Texture import directory path empty, couldn't import %s.",
         src);
     return src;
   }

   char dest_dir_path[FILE_MAX];
   STRNCPY(dest_dir_path, import_dir);

   if (BLI_path_is_rel(import_dir)) {
     const char *basepath = BKE_main_blendfile_path_from_global();

     if (!basepath || strlen(basepath) == 0) {
       WM_reportf(RPT_ERROR,
                  "usd_import_texture(): Texture import directory is relative "
                  "but the blend file path is empty.  "
                  "Please save the blend file before importing the USD. "
                  "Can't import %s.", src);
       return src;
     }

     BLI_path_abs(dest_dir_path, basepath);
   }

   BLI_str_replace_char(dest_dir_path, SEP, ALTSEP);

   if (!BLI_dir_create_recursive(dest_dir_path)) {
     WM_reportf(RPT_ERROR,
                "usd_import_texture(): Couldn't create texture import directory %s.",
                dest_dir_path);
     return src;
   }

   if (udim_src) {
     return copy_udim_asset_to_directory(src, dest_dir_path, overwrite);
   }

   return copy_asset_to_directory(src, dest_dir_path, overwrite);
 }

}  // namespace blender::io::usd
