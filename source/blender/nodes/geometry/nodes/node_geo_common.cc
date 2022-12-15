/* SPDX-License-Identifier: GPL-2.0-or-later */

#include "BKE_node.h"

#include "NOD_geometry.h"
#include "NOD_node_declaration.hh"

#include "NOD_common.h"
#include "node_common.h"
#include "node_geometry_util.hh"

namespace blender::nodes {

static bool node_declare(const bNodeTree &node_tree,
                         const bNode &node,
                         NodeDeclaration &r_declaration)
{
  if (!node_group_declare_dynamic(node_tree, node, r_declaration)) {
    return false;
  }

  const bNodeTree &group = reinterpret_cast<const bNodeTree &>(*node.id);
  const FieldInferencingInterface field_interface = field_inferencing::calculate_field_inferencing(
      group);
  for (const int i : r_declaration.inputs_.index_range()) {
    r_declaration.inputs_[i]->input_field_type_ = field_interface.inputs[i];
  }
  for (const int i : r_declaration.outputs_.index_range()) {
    r_declaration.outputs_[i]->output_field_dependency_ = field_interface.outputs[i];
  }

  return true;
}

}  // namespace blender::nodes

void register_node_type_geo_group()
{
  static bNodeType ntype;

  node_type_base_custom(&ntype, "GeometryNodeGroup", "Group", NODE_CLASS_GROUP);
  ntype.type = NODE_GROUP;
  ntype.poll = geo_node_poll_default;
  ntype.poll_instance = node_group_poll_instance;
  ntype.insert_link = node_insert_link_default;
  ntype.rna_ext.srna = RNA_struct_find("GeometryNodeGroup");
  BLI_assert(ntype.rna_ext.srna != nullptr);
  RNA_struct_blender_type_set(ntype.rna_ext.srna, &ntype);

  node_type_size(&ntype, 140, 60, 400);
  ntype.labelfunc = node_group_label;
  ntype.declare_dynamic = blender::nodes::node_declare;

  nodeRegisterType(&ntype);
}

void register_node_type_geo_custom_group(bNodeType *ntype)
{
  /* These methods can be overridden but need a default implementation otherwise. */
  if (ntype->poll == nullptr) {
    ntype->poll = geo_node_poll_default;
  }
  if (ntype->insert_link == nullptr) {
    ntype->insert_link = node_insert_link_default;
  }
}
