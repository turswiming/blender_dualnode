/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2022 NVIDIA Corportation. All rights reserved. */

#include "usd_asset_utils.h"

#include <pxr/usd/ar/asset.h>
#include <pxr/usd/ar/resolver.h>
#include <pxr/usd/ar/writableAsset.h>

#include "BLI_path_util.h"
#include "BLI_string.h"

#include "WM_api.h"
#include "WM_types.h"

#include <iostream>

namespace blender::io::usd {

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

}  // namespace blender::io::usd
