# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright 2014 Blender Foundation. All rights reserved.

set(INC
  .
  ..
)

include_directories(${INC})

BLENDER_TEST_PERFORMANCE(GEO_constraint_solver_performance "bf_geometry")
