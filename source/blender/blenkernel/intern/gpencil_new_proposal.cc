/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include <algorithm>

#include "BKE_curves.hh"
#include "BKE_gpencil.h"

#include "BLI_index_mask_ops.hh"
#include "BLI_math_vec_types.hh"

#include "DNA_gpencil_types.h"

#include "gpencil_new_proposal.hh"

namespace blender::bke {

Span<float3> GPStroke::points_positions() const
{
  return {geometry_->positions().begin() + offset_, points_num_};
}

MutableSpan<float3> GPStroke::points_positions_for_write() const
{
  return {geometry_->positions_for_write().begin() + offset_, points_num_};
}

void GPStroke::transform(float4x4 matrix)
{
  threading::parallel_for(
      points_positions_for_write().index_range(), 512, [&](const IndexRange range) {
        for (float3 &position : points_positions_for_write().slice(range)) {
          position = matrix * position;
        }
      });
}

}  // namespace blender::bke