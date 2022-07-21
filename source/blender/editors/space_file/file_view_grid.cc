/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup spfile
 */

#include "DNA_ID_enums.h"
#include "DNA_screen_types.h"
#include "DNA_space_types.h"

#include "BKE_context.h"
#include "BKE_icons.h"

#include "BLI_fileops.h"
#include "BLI_math_vector.h"
#include "BLI_rect.h"

/* TODO temp for static ImBuf -> icon_id map. */
#include "BLI_map.hh"
#include "IMB_imbuf.h"

#include "ED_fileselect.h"

#include "UI_grid_view.hh"
#include "UI_interface.h"
#include "UI_interface.hh"

#include "file_intern.h"
#include "filelist.h"

namespace blender::ed::file_browser {

class FilePreviewGridView : public ui::AbstractGridView {
  friend class FilePreviewGridItem;

  FileList &files_;

 public:
  FilePreviewGridView(FileList &filelist);
  void build_items() override;

  BIFIconID file_preview_icon_id_get(const FileDirEntry &file);
};

class FilePreviewGridItem : public ui::PreviewGridItem {
  const FileDirEntry &file_;
  std::string file_identifier_;
  /* Index in the file list. */
  const int file_idx_;

 public:
  FilePreviewGridItem(const FileDirEntry &file, int file_idx);

  void build_grid_tile(uiLayout &layout) const override;

  FileList &get_file_list() const;
  void icon_mono_color_get(uchar r_mono_color[4]);
  void add_big_combined_file_icon(uiLayout &overlap) const;
};

FilePreviewGridView::FilePreviewGridView(FileList &files) : files_(files)
{
}

void FilePreviewGridView::build_items()
{

  const int numfiles = filelist_files_ensure(&files_);

  for (int file_idx = 0; file_idx < numfiles; file_idx++) {
    const FileDirEntry *file = filelist_file(&files_, file_idx);

    add_item<FilePreviewGridItem>(*file, file_idx);
  }
}

FilePreviewGridItem::FilePreviewGridItem(const FileDirEntry &file, const int file_idx)
    : ui::PreviewGridItem(file.relpath, file.name),
      file_(file),
      /* Get a copy so the identifier is always available (the file data may be freed). */
      file_identifier_(identifier_),
      file_idx_(file_idx)
{
  /* Update reference so we don't point into the possibly freed file data. */
  /* TODO always store the identifier as std::string in the item base class? Avoids these issues.
   */
  identifier_ = file_identifier_;
}

FileList &FilePreviewGridItem::get_file_list() const
{
  const FilePreviewGridView &view = dynamic_cast<const FilePreviewGridView &>(get_view());
  return view.files_;
}

static BIFIconID file_big_file_icon_get(const FileDirEntry &file)
{
  /* TODO temp!! Needs proper lifetime/memory management */
  static Map<const ImBuf *, BIFIconID> imbuf_icon_map;

  ImBuf *imb = filelist_file_geticon_image(&file);
  return imbuf_icon_map.lookup_or_add_cb(imb,
                                         [&]() { return (BIFIconID)BKE_icon_imbuf_create(imb); });
}

static void file_icon_mono_color_get(const FileDirEntry &file, uchar r_col[4])
{
  uchar col[4] = {255, 255, 255, 255};
  if (file.typeflag & FILE_TYPE_DIR) {
    UI_GetThemeColor4ubv(TH_ICON_FOLDER, col);
  }
  else {
    UI_GetThemeColor4ubv(TH_TEXT, col);
  }

  const bool is_hidden = (file.attributes & FILE_ATTR_HIDDEN);
  if (is_hidden) {
    col[3] *= 0.3f;
  }

  copy_v4_v4_uchar(r_col, col);
}

void FilePreviewGridItem::add_big_combined_file_icon(uiLayout &layout) const
{
  uiLayout *overlap = uiLayoutOverlap(&layout);

  BIFIconID file_icon_id = file_big_file_icon_get(file_);
  uiBut *preview_but = nullptr;
  if (file_icon_id != ICON_NONE) {
    uchar mono_col[4];
    file_icon_mono_color_get(file_, mono_col);

    preview_but = add_preview_button(*overlap, file_icon_id, mono_col);
  }

  /* Smaller file type icon in the middle of image, scaled to fit container and UI scale */
  {
    float icon_opacity = 0.3f;
    uchar icon_color[4] = {0, 0, 0, 255};
    float bgcolor[4];
    UI_GetThemeColor4fv(TH_ICON_FOLDER, bgcolor);
    if (rgb_to_grayscale(bgcolor) < 0.5f) {
      icon_color[0] = 255;
      icon_color[1] = 255;
      icon_color[2] = 255;
    }
    icon_color[3] *= icon_opacity;

    FileList &files = get_file_list();
    const int icon_id = filelist_geticon(&files, file_idx_, false);

    //    uiLayoutSetAlignment(col, UI_LAYOUT_ALIGN_CENTER);
    const ui::GridViewStyle &style = get_view().get_style();
    int icon_size = style.tile_width / 3.5f;

    uiLayout *col = uiLayoutColumn(overlap, false);
    uiBlock *block = uiLayoutGetBlock(col);

    /* Add padding to vertically center the icon. */
    const rcti preview_img_rect = UI_preview_tile_but_preview_rect_get(preview_but);
    const int icon_ofs_y = style.tile_height - BLI_rcti_cent_y(&preview_img_rect) -
                           (icon_size * ((file_.typeflag & FILE_TYPE_DIR) ? 0.78f : 0.75f)) / 2;
    uiDefButPadding(block, 0, 0, 0, icon_ofs_y);
    uiDefButPreviewTile(block, icon_id, "", 0, 0, icon_size, icon_size, icon_color);
  }
}

void FilePreviewGridItem::build_grid_tile(uiLayout &layout) const
{
  uiLayout &overlap = *uiLayoutOverlap(&layout);
  add_big_combined_file_icon(overlap);
}

}  // namespace blender::ed::file_browser

namespace ui = blender::ui;
using namespace blender::ed::file_browser;

void file_grid_view_create_in_layout(FileList *files, const View2D *v2d, uiLayout *layout)
{
  uiBlock *block = uiLayoutGetBlock(layout);
  UI_block_layout_set_current(block, layout);

  ui::AbstractGridView *grid_view = UI_block_add_view(
      *block, "file preview grid view", std::make_unique<FilePreviewGridView>(*files));

  ui::GridViewBuilder builder(*block);
  builder.build_grid_view(*grid_view, *v2d);
}
