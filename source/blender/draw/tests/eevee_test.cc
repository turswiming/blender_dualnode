/* SPDX-License-Identifier: Apache-2.0 */

#include "testing/testing.h"

#include "BKE_context.h"
#include "BKE_idtype.h"
#include "BKE_main.h"
#include "BKE_node.h"
#include "BKE_object.h"
#include "DEG_depsgraph.h"
#include "RNA_define.h"

#include "draw_shader.h"
#include "draw_testing.hh"
#include "engines/eevee_next/eevee_instance.hh"

namespace blender::draw {

using namespace blender::eevee;

/* Replace with template version that is not GPU only. */
using ShadowPageCacheBuf = draw::StorageArrayBuffer<uint2, SHADOW_MAX_PAGE, false>;
using ShadowTileDataBuf = draw::StorageArrayBuffer<ShadowTileDataPacked, SHADOW_MAX_TILE, false>;

static void test_eevee_shadow_shift_clear()
{
  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  ShadowTileDataBuf tiles_data = {"tiles_data"};

  int tiles_index = 1;
  int tile_lod0 = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 5;
  int tile_lod1 = tile_lod0 + square_i(SHADOW_TILEMAP_RES);

  {
    ShadowTileMapData tilemap = {};
    tilemap.tiles_index = tiles_index;
    tilemap.grid_shift = int2(SHADOW_TILEMAP_RES);
    tilemap.is_cubeface = true;

    tilemaps_data.append(tilemap);

    tilemaps_data.push_update();
  }
  {
    ShadowTileData tile;

    tile.is_used = true;
    tile.do_update = true;
    tiles_data[tile_lod0] = shadow_tile_pack(tile);

    tile.is_used = false;
    tile.do_update = false;
    tiles_data[tile_lod1] = shadow_tile_pack(tile);

    tiles_data.push_update();
  }

  GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_tilemap_init");

  PassSimple pass("Test");
  pass.shader_set(sh);
  pass.bind_ssbo("tilemaps_buf", tilemaps_data);
  pass.bind_ssbo("tiles_buf", tiles_data);
  pass.dispatch(int3(1, 1, tilemaps_data.size()));

  Manager manager;
  manager.submit(pass);
  GPU_finish();

  tilemaps_data.read();
  tiles_data.read();

  EXPECT_EQ(tilemaps_data[0].grid_offset, int2(0));
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_lod0]).is_used, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_lod0]).do_update, true);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_lod1]).is_used, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_lod1]).do_update, true);

  GPU_shader_free(sh);
  DRW_shaders_free();
}
DRAW_TEST(eevee_shadow_shift_clear)

static void test_eevee_shadow_tag_update()
{
  StorageVectorBuffer<uint, 128> past_casters_updated = {"PastCastersUpdated"};
  StorageVectorBuffer<uint, 128> curr_casters_updated = {"CurrCastersUpdated"};

  Manager manager;
  {
    /* Simulate 1 object moving and 1 object static with changing resource index. */
    float4x4 obmat = float4x4::identity();
    float4x4 obmat2 = float4x4::from_loc_eul_scale(float3(10.0f), float3(0.0f), float3(0.5f));
    float3 half_extent = float3(0.249f, 0.249f, 0.001f);

    {
      manager.begin_sync();
      ResourceHandle hdl = manager.resource_handle(obmat, float3(0.5f, 0.5f, -1.0f), half_extent);
      manager.resource_handle(obmat2);
      manager.end_sync();
      past_casters_updated.append(hdl.resource_index());
      past_casters_updated.push_update();
    }
    {
      manager.begin_sync();
      manager.resource_handle(obmat2);
      ResourceHandle hdl = manager.resource_handle(obmat, float3(-1.0f, 0.5f, -1.0f), half_extent);
      manager.end_sync();
      curr_casters_updated.append(hdl.resource_index());
      curr_casters_updated.push_update();
    }
  }

  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  tiles_data.clear_to_zero();

  {
    ShadowTileMap tilemap(0);
    tilemap.sync_cubeface(float4x4::identity(), 0.01f, 1.0f, DEG2RADF(75.0f), Z_NEG);
    tilemaps_data.append(tilemap);
  }
  {
    ShadowTileMap tilemap(1);
    tilemap.sync_clipmap(float3(0.0f), float4x4::identity(), -1.0f, 1.0f, int2(0), 1);
    tilemaps_data.append(tilemap);
  }

  tilemaps_data.push_update();

  GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_tag_update");

  PassSimple pass("Test");
  pass.shader_set(sh);
  pass.bind_ssbo("tilemaps_buf", tilemaps_data);
  pass.bind_ssbo("tiles_buf", tiles_data);
  pass.bind_ssbo("bounds_buf", &manager.bounds_buf.previous());
  pass.bind_ssbo("resource_ids_buf", past_casters_updated);
  pass.dispatch(int3(past_casters_updated.size(), 1, tilemaps_data.size()));
  pass.bind_ssbo("bounds_buf", &manager.bounds_buf.current());
  pass.bind_ssbo("resource_ids_buf", curr_casters_updated);
  pass.dispatch(int3(curr_casters_updated.size(), 1, tilemaps_data.size()));

  manager.submit(pass);
  GPU_finish();

  tiles_data.read();

  /** The layout of these expected strings is Y down. */
  StringRefNull expected_lod0 =
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "----------------"
      "xx--------xxxx--"
      "xx--------xxxx--"
      "xx--------xxxx--"
      "xx--------xxxx--"
      "----------------"
      "----------------";
  StringRefNull expected_lod1 =
      "--------"
      "--------"
      "--------"
      "--------"
      "--------"
      "x----xx-"
      "x----xx-"
      "--------";
  StringRefNull expected_lod2 =
      "----"
      "----"
      "x-xx"
      "x-xx";
  StringRefNull expected_lod3 =
      "--"
      "xx";
  StringRefNull expected_lod4 = "x";
  const uint lod0_len = square_i(SHADOW_TILEMAP_RES);
  const uint lod1_len = square_i(SHADOW_TILEMAP_RES / 2);
  const uint lod2_len = square_i(SHADOW_TILEMAP_RES / 4);
  const uint lod3_len = square_i(SHADOW_TILEMAP_RES / 8);
  const uint lod4_len = square_i(SHADOW_TILEMAP_RES / 16);

  auto stringify_result = [&](uint start, uint len) -> std::string {
    std::string result = "";
    for (auto i : IndexRange(start, len)) {
      result += (shadow_tile_unpack(tiles_data[i]).do_update) ? "x" : "-";
    }
    return result;
  };

  EXPECT_EQ(stringify_result(0, lod0_len), expected_lod0);
  EXPECT_EQ(stringify_result(lod0_len, lod1_len), expected_lod1);
  EXPECT_EQ(stringify_result(lod0_len + lod1_len, lod2_len), expected_lod2);
  EXPECT_EQ(stringify_result(lod0_len + lod1_len + lod2_len, lod3_len), expected_lod3);
  EXPECT_EQ(stringify_result(lod0_len + lod1_len + lod2_len + lod3_len, lod4_len), expected_lod4);

  GPU_shader_free(sh);
  DRW_shaders_free();
}
DRAW_TEST(eevee_shadow_tag_update)

static void test_eevee_shadow_free()
{
  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  ShadowPageHeapBuf pages_free_data = {"PagesFreeBuf"};
  ShadowPageCacheBuf pages_cached_data = {"PagesCachedBuf"};
  ShadowPagesInfoDataBuf pages_infos_data = {"PagesInfosBuf"};

  int tiles_index = 1;
  int tile_orphaned_cached = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 5;
  int tile_orphaned_allocated = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 6;
  int tile_used_cached = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 260;
  int tile_used_allocated = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 32;
  int tile_used_unallocated = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 64;
  int tile_unused_cached = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 9;
  int tile_unused_allocated = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 8;
  int page_free_count = SHADOW_MAX_PAGE - 6;

  for (uint i : IndexRange(2, page_free_count)) {
    uint2 page = {i % SHADOW_PAGE_PER_ROW, i / SHADOW_PAGE_PER_ROW};
    pages_free_data[i] = page.x | (page.y << 16u);
  }
  pages_free_data.push_update();

  pages_infos_data.page_free_count = page_free_count;
  pages_infos_data.page_alloc_count = 0;
  pages_infos_data.page_cached_next = 2u;
  pages_infos_data.page_cached_start = 0u;
  pages_infos_data.page_cached_end = 2u;
  pages_infos_data.push_update();

  for (uint i : IndexRange(pages_cached_data.size())) {
    pages_cached_data[i] = uint2(-1, -1);
  }
  pages_cached_data[0] = uint2(0, tile_orphaned_cached);
  pages_cached_data[1] = uint2(1, tile_used_cached);
  pages_cached_data.push_update();

  {
    ShadowTileData tile;

    /* is_orphaned = true */
    tile.is_used = false;
    tile.do_update = true;

    tile.is_cached = true;
    tile.is_allocated = false;
    tiles_data[tile_orphaned_cached] = shadow_tile_pack(tile);

    tile.is_cached = false;
    tile.is_allocated = true;
    tiles_data[tile_orphaned_allocated] = shadow_tile_pack(tile);

    /* is_orphaned = false */
    tile.do_update = false;
    tile.is_used = true;

    tile.is_cached = true;
    tile.is_allocated = false;
    tiles_data[tile_used_cached] = shadow_tile_pack(tile);

    tile.is_cached = false;
    tile.is_allocated = true;
    tiles_data[tile_used_allocated] = shadow_tile_pack(tile);

    tile.is_cached = false;
    tile.is_allocated = false;
    tiles_data[tile_used_unallocated] = shadow_tile_pack(tile);

    tile.is_used = false;
    tile.is_cached = true;
    tile.is_allocated = false;
    tiles_data[tile_unused_cached] = shadow_tile_pack(tile);

    tile.is_cached = false;
    tile.is_allocated = true;
    tiles_data[tile_unused_allocated] = shadow_tile_pack(tile);

    tiles_data.push_update();
  }
  {
    ShadowTileMapData tilemap = {};
    tilemap.tiles_index = tiles_index;
    tilemaps_data.append(tilemap);
    tilemaps_data.push_update();
  }

  GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_page_free");

  PassSimple pass("Test");
  pass.shader_set(sh);
  pass.bind_ssbo("tilemaps_buf", tilemaps_data);
  pass.bind_ssbo("tiles_buf", tiles_data);
  pass.bind_ssbo("pages_infos_buf", pages_infos_data);
  pass.bind_ssbo("pages_free_buf", pages_free_data);
  pass.bind_ssbo("pages_cached_buf", pages_cached_data);
  pass.dispatch(int3(1, 1, tilemaps_data.size()));

  Manager manager;
  manager.submit(pass);
  GPU_finish();

  tiles_data.read();
  pages_infos_data.read();

  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_orphaned_cached]).is_cached, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_orphaned_cached]).is_allocated, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_orphaned_allocated]).is_cached, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_orphaned_allocated]).is_allocated, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_cached]).is_cached, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_cached]).is_allocated, true);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_allocated]).is_cached, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_allocated]).is_allocated, true);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_unallocated]).is_cached, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_used_unallocated]).is_allocated, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_unused_cached]).is_cached, true);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_unused_cached]).is_allocated, false);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_unused_allocated]).is_cached, true);
  EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_unused_allocated]).is_allocated, false);
  EXPECT_EQ(pages_infos_data.page_alloc_count, 1);
  EXPECT_EQ(pages_infos_data.page_free_count, page_free_count + 2);
  EXPECT_EQ(pages_infos_data.page_cached_next, 3);
  EXPECT_EQ(pages_infos_data.page_cached_end, 2);

  GPU_shader_free(sh);
  DRW_shaders_free();
}
DRAW_TEST(eevee_shadow_free)

class TestDefrag {
 private:
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  ShadowPageHeapBuf pages_free_data = {"PagesFreeBuf"};
  ShadowPageCacheBuf pages_cached_data = {"PagesCachedBuf"};
  ShadowPagesInfoDataBuf pages_infos_data = {"PagesInfosBuf"};

 public:
  TestDefrag(int allocation_count,
             int descriptor_offset,
             StringRefNull descriptor,
             StringRefNull expect)
  {
    for (uint i : IndexRange(SHADOW_MAX_PAGE)) {
      uint2 page = {i % SHADOW_PAGE_PER_ROW, i / SHADOW_PAGE_PER_ROW};
      pages_free_data[i] = page.x | (page.y << 16u);
    }

    for (uint i : IndexRange(tiles_data.size())) {
      tiles_data[i] = 0;
    }

    int free_count = SHADOW_MAX_PAGE;
    int tile_index = 0;

    for (uint i : IndexRange(pages_cached_data.size())) {
      pages_cached_data[i] = uint2(-1, -1);
    }

    int cached_index = descriptor_offset;
    int hole_count = 0;
    int inserted_count = 0;
    ShadowTileData tile = {};
    tile.is_cached = true;
    for (char c : descriptor) {
      switch (c) {
        case 'c':
          tile.cache_index = cached_index++ % SHADOW_MAX_PAGE;
          pages_cached_data[tile.cache_index] = uint2(pages_free_data[--free_count], tile_index);
          tiles_data[tile_index++] = shadow_tile_pack(tile);
          break;
        case 'f':
          pages_cached_data[cached_index++ % SHADOW_MAX_PAGE] = uint2(-1, -1);
          hole_count++;
          break;
        case 'i':
          tile.cache_index = (cached_index + inserted_count++) % SHADOW_MAX_PAGE;
          pages_cached_data[tile.cache_index] = uint2(pages_free_data[--free_count], tile_index);
          tiles_data[tile_index++] = shadow_tile_pack(tile);
          break;
        default:
          break;
      }
    }

    pages_infos_data.page_alloc_count = allocation_count;
    pages_infos_data.page_cached_next = cached_index + inserted_count;
    pages_infos_data.page_free_count = free_count;
    pages_infos_data.page_cached_start = descriptor_offset;
    pages_infos_data.page_cached_end = cached_index;

    tiles_data.push_update();
    pages_infos_data.push_update();
    pages_free_data.push_update();
    pages_cached_data.push_update();

    GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_page_defrag");

    PassSimple pass("Test");
    pass.shader_set(sh);
    pass.bind_ssbo("tiles_buf", tiles_data);
    pass.bind_ssbo("pages_infos_buf", pages_infos_data);
    pass.bind_ssbo("pages_free_buf", pages_free_data);
    pass.bind_ssbo("pages_cached_buf", pages_cached_data);
    pass.dispatch(int3(1, 1, 1));

    Manager manager;
    manager.submit(pass);
    GPU_finish();

    tiles_data.read();
    pages_cached_data.read();
    pages_infos_data.read();

    std::string result = "";
    int expect_cached_len = 0;
    for (auto i : IndexRange(descriptor_offset, descriptor.size())) {
      if (pages_cached_data[i % SHADOW_MAX_PAGE].y != -1) {
        result += 'c';
        expect_cached_len++;
      }
      else {
        result += 'f';
      }
    }
    EXPECT_EQ(expect, result);

    allocation_count = min_ii(allocation_count, SHADOW_MAX_PAGE);

    int additional_pages = max_ii(0, allocation_count - free_count);
    int expected_free_count = max_ii(free_count, allocation_count);
    int expected_start = descriptor_offset + hole_count + additional_pages;
    int result_cached_len = pages_infos_data.page_cached_end - pages_infos_data.page_cached_start;

    if (expected_start > SHADOW_MAX_PAGE) {
      expected_start -= SHADOW_MAX_PAGE;
    }

    EXPECT_EQ(expected_free_count, pages_infos_data.page_free_count);
    EXPECT_EQ(expected_start, pages_infos_data.page_cached_start);
    EXPECT_EQ(expect_cached_len, result_cached_len);
    EXPECT_EQ(pages_infos_data.page_cached_end, pages_infos_data.page_cached_next);

    GPU_shader_free(sh);
    DRW_shaders_free();
  }
};

static void test_eevee_shadow_defrag()
{
  TestDefrag(0, 0, "cfi", "fcc");
  TestDefrag(0, 0, "fci", "fcc");
  TestDefrag(0, 47, "ccfcffccfcfciiiii", "fffffcccccccccccc");
  TestDefrag(10, SHADOW_MAX_PAGE - 5, "ccfcffccfcfciiiii", "fffffcccccccccccc");
  TestDefrag(SHADOW_MAX_PAGE - 8, 30, "ccfcffccfcfciiiii", "fffffffffcccccccc");
  TestDefrag(SHADOW_MAX_PAGE - 4, 30, "ccfcffccfcfciiiii", "fffffffffffffcccc");
  /* Over allocation but should not crash. */
  TestDefrag(SHADOW_MAX_PAGE + 4, 30, "ccfcffccfcfciiiii", "fffffffffffffffff");
}
DRAW_TEST(eevee_shadow_defrag)

class TestAlloc {
 private:
  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  ShadowPageHeapBuf pages_free_data = {"PagesFreeBuf"};
  ShadowPageCacheBuf pages_cached_data = {"PagesCachedBuf"};
  ShadowPagesInfoDataBuf pages_infos_data = {"PagesInfosBuf"};

 public:
  TestAlloc(int page_free_count)
  {
    int tiles_index = 1;

    for (uint i : IndexRange(0, page_free_count)) {
      uint2 page = {i % SHADOW_PAGE_PER_ROW, i / SHADOW_PAGE_PER_ROW};
      pages_free_data[i] = page.x | (page.y << 16u);
    }
    pages_free_data.push_update();
    pages_cached_data.push_update();

    pages_infos_data.page_free_count = page_free_count;
    pages_infos_data.page_alloc_count = 1;
    pages_infos_data.page_cached_next = 0u;
    pages_infos_data.page_cached_start = 0u;
    pages_infos_data.page_cached_end = 0u;
    pages_infos_data.view_count = 0u;
    pages_infos_data.page_size = 256u;
    pages_infos_data.push_update();

    int tile_allocated = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 5;
    int tile_free = tiles_index * SHADOW_TILEDATA_PER_TILEMAP + 6;

    {
      ShadowTileData tile;

      tile.is_used = true;
      tile.do_update = false;

      tile.is_cached = false;
      tile.is_allocated = false;
      tiles_data[tile_free] = shadow_tile_pack(tile);

      tile.is_cached = false;
      tile.is_allocated = true;
      tiles_data[tile_allocated] = shadow_tile_pack(tile);

      tiles_data.push_update();
    }
    {
      ShadowTileMapData tilemap = {};
      tilemap.tiles_index = tiles_index;
      tilemaps_data.append(tilemap);
      tilemaps_data.push_update();
    }

    GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_page_allocate");

    PassSimple pass("Test");
    pass.shader_set(sh);
    pass.bind_ssbo("tilemaps_buf", tilemaps_data);
    pass.bind_ssbo("tiles_buf", tiles_data);
    pass.bind_ssbo("pages_infos_buf", pages_infos_data);
    pass.bind_ssbo("pages_free_buf", pages_free_data);
    pass.bind_ssbo("pages_cached_buf", pages_cached_data);
    pass.dispatch(int3(1, 1, tilemaps_data.size()));

    Manager manager;
    manager.submit(pass);
    GPU_finish();

    tiles_data.read();
    pages_infos_data.read();

    bool alloc_success = page_free_count >= 1;

    EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_free]).do_update, alloc_success);
    EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_free]).is_allocated, alloc_success);
    EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_allocated]).do_update, false);
    EXPECT_EQ(shadow_tile_unpack(tiles_data[tile_allocated]).is_allocated, true);
    EXPECT_EQ(pages_infos_data.page_free_count, page_free_count - 1);

    GPU_shader_free(sh);
    DRW_shaders_free();
  }
};

static void test_eevee_shadow_alloc()
{
  TestAlloc(SHADOW_MAX_PAGE);
  TestAlloc(1);
  TestAlloc(0);
}
DRAW_TEST(eevee_shadow_alloc)

static void test_eevee_shadow_finalize()
{
  ShadowTileMapDataBuf tilemaps_data = {"tilemaps_data"};
  ShadowTileDataBuf tiles_data = {"tiles_data"};
  ShadowPageHeapBuf pages_free_data = {"PagesFreeBuf"};
  ShadowPageCacheBuf pages_cached_data = {"PagesCachedBuf"};
  ShadowPagesInfoDataBuf pages_infos_data = {"PagesInfosBuf"};

  const uint lod0_len = square_i(SHADOW_TILEMAP_RES);
  const uint lod1_len = square_i(SHADOW_TILEMAP_RES / 2);
  const uint lod2_len = square_i(SHADOW_TILEMAP_RES / 4);
  const uint lod3_len = square_i(SHADOW_TILEMAP_RES / 8);

  const uint lod0_ofs = 0;
  const uint lod1_ofs = lod0_len;
  const uint lod2_ofs = lod1_ofs + lod1_len;
  const uint lod3_ofs = lod2_ofs + lod2_len;
  const uint lod4_ofs = lod3_ofs + lod3_len;

  for (auto i : IndexRange(SHADOW_TILEDATA_PER_TILEMAP)) {
    tiles_data[i] = 0;
  }

  {
    ShadowTileData tile;
    tile.is_used = true;
    tile.is_allocated = true;

    tile.page = uint2(1, 0);
    tile.do_update = false;
    tiles_data[lod0_ofs] = shadow_tile_pack(tile);

    tile.page = uint2(2, 0);
    tile.do_update = false;
    tiles_data[lod1_ofs] = shadow_tile_pack(tile);

    tile.page = uint2(3, 0);
    tile.do_update = true;
    tiles_data[lod2_ofs] = shadow_tile_pack(tile);

    tile.page = uint2(4, 0);
    tile.do_update = false;
    tiles_data[lod3_ofs] = shadow_tile_pack(tile);

    tile.page = uint2(5, 0);
    tile.do_update = true;
    tiles_data[lod4_ofs] = shadow_tile_pack(tile);

    tile.page = uint2(6, 0);
    tile.do_update = true;
    tiles_data[lod0_ofs + 8] = shadow_tile_pack(tile);

    tiles_data.push_update();
  }
  {
    ShadowTileMapData tilemap = {};
    tilemap.tiles_index = 0;
    tilemap.is_cubeface = true;
    tilemaps_data.append(tilemap);

    tilemaps_data.push_update();
  }
  {
    pages_infos_data.page_free_count = -5;
    pages_infos_data.page_alloc_count = 0;
    pages_infos_data.page_cached_next = 0u;
    pages_infos_data.page_cached_start = 0u;
    pages_infos_data.page_cached_end = 0u;
    pages_infos_data.view_count = 0u;
    pages_infos_data.page_size = 256u;
    pages_infos_data.push_update();
  }

  Texture tilemap_tx = {"tilemap_tx"};
  tilemap_tx.ensure_2d(GPU_R32UI,
                       int2(SHADOW_TILEMAP_RES),
                       GPU_TEXTURE_USAGE_HOST_READ | GPU_TEXTURE_USAGE_SHADER_READ |
                           GPU_TEXTURE_USAGE_SHADER_WRITE);
  tilemap_tx.clear(uint4(0));

  Texture render_map_tx = {"ShadowRenderMap",
                           GPU_R32UI,
                           GPU_TEXTURE_USAGE_HOST_READ | GPU_TEXTURE_USAGE_SHADER_READ |
                               GPU_TEXTURE_USAGE_SHADER_WRITE,
                           int2(SHADOW_TILEMAP_RES),
                           1, /* Only one layer for the test. */
                           nullptr,
                           SHADOW_TILEMAP_LOD + 1};
  render_map_tx.ensure_mip_views();

  View shadow_multi_view = {"ShadowMultiView", 64, true};
  StorageBuffer<DispatchCommand> clear_dispatch_buf;
  StorageArrayBuffer<uint, SHADOW_MAX_PAGE> clear_page_buf = {"clear_page_buf"};

  GPUShader *sh = GPU_shader_create_from_info_name("eevee_shadow_tilemap_finalize");

  PassSimple pass("Test");
  pass.shader_set(sh);
  pass.bind_ssbo("tilemaps_buf", tilemaps_data);
  pass.bind_ssbo("tiles_buf", tiles_data);
  pass.bind_ssbo("pages_infos_buf", pages_infos_data);
  pass.bind_image("tilemaps_img", tilemap_tx);
  pass.bind_ssbo("view_infos_buf", shadow_multi_view.matrices_ubo_get());
  pass.bind_ssbo("clear_dispatch_buf", clear_dispatch_buf);
  pass.bind_ssbo("clear_page_buf", clear_page_buf);
  pass.bind_image("render_map_lod0_img", render_map_tx.mip_view(0));
  pass.bind_image("render_map_lod1_img", render_map_tx.mip_view(1));
  pass.bind_image("render_map_lod2_img", render_map_tx.mip_view(2));
  pass.bind_image("render_map_lod3_img", render_map_tx.mip_view(3));
  pass.bind_image("render_map_lod4_img", render_map_tx.mip_view(4));
  pass.dispatch(int3(1, 1, tilemaps_data.size()));

  Manager manager;
  manager.submit(pass);
  GPU_finish();

  GPU_memory_barrier(GPU_BARRIER_TEXTURE_UPDATE);

  {
    uint *pixels = tilemap_tx.read<uint32_t>(GPU_DATA_UINT);

    std::string result = "";
    for (auto y : IndexRange(SHADOW_TILEMAP_RES)) {
      for (auto x : IndexRange(SHADOW_TILEMAP_RES)) {
        result += std::to_string(shadow_tile_unpack(pixels[y * SHADOW_TILEMAP_RES + x]).page.x);
      }
    }

    MEM_SAFE_FREE(pixels);

    /** The layout of these expected strings is Y down. */
    StringRefNull expected_pages =
        "1233444465555555"
        "2233444455555555"
        "3333444455555555"
        "3333444455555555"
        "4444444455555555"
        "4444444455555555"
        "4444444455555555"
        "4444444455555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555"
        "5555555555555555";

    EXPECT_EQ(expected_pages, result);
  }

  {
    auto stringify_lod = [](Span<uint> data) -> std::string {
      std::string result = "";
      for (auto x : data) {
        /* Take only X page location since we allocated manually. */
        result += std::to_string(x & 0x7u);
      }
      return result;
    };

    /** The layout of these expected strings is Y down. */
    StringRefNull expected_lod0 =
        "7777777767777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777"
        "7777777777777777";

    StringRefNull expected_lod1 =
        "77777777"
        "77777777"
        "77777777"
        "77777777"
        "77777777"
        "77777777"
        "77777777"
        "77777777";

    StringRefNull expected_lod2 =
        "3777"
        "7777"
        "7777"
        "7777";

    StringRefNull expected_lod3 =
        "77"
        "77";

    StringRefNull expected_lod4 = "5";

    uint *pixels_lod0 = render_map_tx.read<uint32_t>(GPU_DATA_UINT, 0);
    uint *pixels_lod1 = render_map_tx.read<uint32_t>(GPU_DATA_UINT, 1);
    uint *pixels_lod2 = render_map_tx.read<uint32_t>(GPU_DATA_UINT, 2);
    uint *pixels_lod3 = render_map_tx.read<uint32_t>(GPU_DATA_UINT, 3);
    uint *pixels_lod4 = render_map_tx.read<uint32_t>(GPU_DATA_UINT, 4);

    EXPECT_EQ(stringify_lod(Span<uint>(pixels_lod0, lod0_len)), expected_lod0);
    EXPECT_EQ(stringify_lod(Span<uint>(pixels_lod1, lod1_len)), expected_lod1);
    EXPECT_EQ(stringify_lod(Span<uint>(pixels_lod2, lod2_len)), expected_lod2);
    EXPECT_EQ(stringify_lod(Span<uint>(pixels_lod3, lod3_len)), expected_lod3);
    EXPECT_EQ(stringify_lod(Span<uint>(pixels_lod4, 1)), expected_lod4);

    MEM_SAFE_FREE(pixels_lod0);
    MEM_SAFE_FREE(pixels_lod1);
    MEM_SAFE_FREE(pixels_lod2);
    MEM_SAFE_FREE(pixels_lod3);
    MEM_SAFE_FREE(pixels_lod4);
  }

  pages_infos_data.read();
  EXPECT_EQ(pages_infos_data.page_free_count, 0);
  EXPECT_EQ(pages_infos_data.view_count, 1);

  GPU_shader_free(sh);
  DRW_shaders_free();
}
DRAW_TEST(eevee_shadow_finalize)

}  // namespace blender::draw
