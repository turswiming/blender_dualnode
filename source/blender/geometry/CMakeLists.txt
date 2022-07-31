# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright 2006 Blender Foundation. All rights reserved.

set(INC
  .
  ../blenkernel
  ../blenlib
  ../blentranslation
  ../functions
  ../makesdna
  ../makesrna
  ../../../intern/eigen
  ../../../intern/guardedalloc
  ${CMAKE_BINARY_DIR}/source/blender/makesdna/intern
)

set(SRC
  intern/add_curves_on_mesh.cc
  intern/constraint_solver.cc
  intern/fillet_curves.cc
  intern/mesh_merge_by_distance.cc
  intern/mesh_primitive_cuboid.cc
  intern/mesh_to_curve_convert.cc
  intern/mesh_to_volume.cc
  intern/point_merge_by_distance.cc
  intern/realize_instances.cc
  intern/resample_curves.cc
  intern/reverse_uv_sampler.cc
  intern/set_curve_type.cc
  intern/subdivide_curves.cc
  intern/uv_parametrizer.c

  GEO_add_curves_on_mesh.hh
  GEO_constraint_solver.hh
  GEO_fillet_curves.hh
  GEO_mesh_merge_by_distance.hh
  GEO_mesh_primitive_cuboid.hh
  GEO_mesh_to_curve.hh
  GEO_mesh_to_volume.hh
  GEO_point_merge_by_distance.hh
  GEO_realize_instances.hh
  GEO_resample_curves.hh
  GEO_reverse_uv_sampler.hh
  GEO_set_curve_type.hh
  GEO_subdivide_curves.hh
  GEO_uv_parametrizer.h
)

set(LIB
  bf_blenkernel
  bf_blenlib
)

if(WITH_OPENVDB)
  list(APPEND INC
    ../../../intern/openvdb
  )
  list(APPEND INC_SYS
    ${OPENVDB_INCLUDE_DIRS}
  )
  list(APPEND LIB
    bf_intern_openvdb
    ${OPENVDB_LIBRARIES}
  )
  add_definitions(-DWITH_OPENVDB ${OPENVDB_DEFINITIONS})
endif()

if(WITH_TBB)
  add_definitions(-DWITH_TBB)

  list(APPEND INC_SYS
    ${TBB_INCLUDE_DIRS}
  )

  list(APPEND LIB
    ${TBB_LIBRARIES}
  )
endif()

blender_add_lib(bf_geometry "${SRC}" "${INC}" "${INC_SYS}" "${LIB}")

if(WITH_GTESTS)
  set(TEST_SRC
    tests/GEO_constraint_solver_test.cc
  )
  set(TEST_INC
  )
  set(TEST_LIB
    bf_geometry
  )
  include(GTestTesting)
  blender_add_test_executable(geometry "${TEST_SRC}" "${INC};${TEST_INC}" "${INC_SYS}" "${LIB};${TEST_LIB}")

  add_subdirectory(tests/performance)
endif()
