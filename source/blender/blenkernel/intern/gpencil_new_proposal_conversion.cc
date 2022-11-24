/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

#include "DNA_gpencil_types.h"
#include "gpencil_new_proposal.hh"

namespace blender::bke {

GreasePencil convert_old_to_new_gpencil_data(const bGPdata *old_gpd)
{
  GreasePencil new_gpd;

  return new_gpd;
}

bGPdata *convert_new_to_old_gpencil_data(const GreasePencil &new_gpd)
{
  bGPdata *old_gpd = reinterpret_cast<bGPdata *>(MEM_mallocN(sizeof(bGPdata), __func__));

  return old_gpd;
}

}  // namespace blender::bke