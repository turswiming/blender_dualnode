/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2007 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup nodes
 */

#include <cstddef>
#include <cstring>

#include "DNA_node_types.h"

#include "BLI_listbase.h"
#include "BLI_map.hh"
#include "BLI_multi_value_map.hh"
#include "BLI_set.hh"
#include "BLI_stack.hh"
#include "BLI_string.h"
#include "BLI_string_ref.hh"
#include "BLI_utildefines.h"

#include "BLT_translation.h"

#include "BKE_node.h"
#include "BKE_node_runtime.hh"
#include "BKE_node_tree_update.h"

#include "RNA_types.h"

#include "MEM_guardedalloc.h"

#include "NOD_common.h"
#include "NOD_node_declaration.hh"
#include "NOD_register.hh"
#include "NOD_socket_declarations.hh"
#include "NOD_socket_declarations_geometry.hh"
#include "node_common.h"
#include "node_util.h"

using blender::Map;
using blender::MultiValueMap;
using blender::Set;
using blender::Stack;
using blender::StringRef;
using blender::Vector;

/* -------------------------------------------------------------------- */
/** \name Node Group
 * \{ */

static bNodeSocket *find_matching_socket(ListBase &sockets, StringRef identifier)
{
  LISTBASE_FOREACH (bNodeSocket *, socket, &sockets) {
    if (socket->identifier == identifier) {
      return socket;
    }
  }
  return nullptr;
}

bNodeSocket *node_group_find_input_socket(bNode *groupnode, const char *identifier)
{
  return find_matching_socket(groupnode->inputs, identifier);
}

bNodeSocket *node_group_find_output_socket(bNode *groupnode, const char *identifier)
{
  return find_matching_socket(groupnode->outputs, identifier);
}

void node_group_label(const bNodeTree * /*ntree*/, const bNode *node, char *label, int maxlen)
{
  BLI_strncpy(label, (node->id) ? node->id->name + 2 : IFACE_("Missing Data-Block"), maxlen);
}

bool node_group_poll_instance(bNode *node, bNodeTree *nodetree, const char **disabled_hint)
{
  if (node->typeinfo->poll(node->typeinfo, nodetree, disabled_hint)) {
    bNodeTree *grouptree = (bNodeTree *)node->id;
    if (grouptree) {
      return nodeGroupPoll(nodetree, grouptree, disabled_hint);
    }

    return true; /* without a linked node tree, group node is always ok */
  }

  return false;
}

bool nodeGroupPoll(const bNodeTree *nodetree,
                   const bNodeTree *grouptree,
                   const char **r_disabled_hint)
{
  /* unspecified node group, generally allowed
   * (if anything, should be avoided on operator level)
   */
  if (grouptree == nullptr) {
    return true;
  }

  if (nodetree == grouptree) {
    if (r_disabled_hint) {
      *r_disabled_hint = TIP_("Nesting a node group inside of itself is not allowed");
    }
    return false;
  }
  if (nodetree->type != grouptree->type) {
    if (r_disabled_hint) {
      *r_disabled_hint = TIP_("Node group has different type");
    }
    return false;
  }

  LISTBASE_FOREACH (const bNode *, node, &grouptree->nodes) {
    if (node->typeinfo->poll_instance &&
        !node->typeinfo->poll_instance(
            const_cast<bNode *>(node), const_cast<bNodeTree *>(nodetree), r_disabled_hint)) {
      return false;
    }
  }
  return true;
}

namespace blender::nodes {

static SocketDeclarationPtr declataion_for_interface_socket(const bNodeSocket &io_socket)
{
  SocketDeclarationPtr dst;
  switch (io_socket.type) {
    case SOCK_FLOAT: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueFloat>();
      std::unique_ptr<decl::Float> decl = std::make_unique<decl::Float>();
      decl->default_value_ = value.value;
      decl->soft_min_value_ = value.min;
      decl->soft_max_value_ = value.max;
      dst = std::move(decl);
      break;
    }
    case SOCK_VECTOR: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueVector>();
      std::unique_ptr<decl::Vector> decl = std::make_unique<decl::Vector>();
      decl->default_value_ = value.value;
      decl->soft_min_value_ = value.min;
      decl->soft_max_value_ = value.max;
      dst = std::move(decl);
      break;
    }
    case SOCK_RGBA: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueRGBA>();
      std::unique_ptr<decl::Color> decl = std::make_unique<decl::Color>();
      decl->default_value_ = value.value;
      dst = std::move(decl);
      break;
    }
    case SOCK_SHADER: {
      std::unique_ptr<decl::Color> decl = std::make_unique<decl::Color>();
      dst = std::move(decl);
      break;
    }
    case SOCK_BOOLEAN: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueBoolean>();
      std::unique_ptr<decl::Bool> decl = std::make_unique<decl::Bool>();
      decl->default_value_ = value.value;
      dst = std::move(decl);
      break;
    }
    case SOCK_INT: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueInt>();
      std::unique_ptr<decl::Int> decl = std::make_unique<decl::Int>();
      decl->default_value_ = value.value;
      decl->soft_min_value_ = value.min;
      decl->soft_max_value_ = value.max;
      dst = std::move(decl);
      break;
    }
    case SOCK_STRING: {
      const auto &value = *io_socket.default_value_typed<bNodeSocketValueString>();
      std::unique_ptr<decl::String> decl = std::make_unique<decl::String>();
      decl->default_value_ = value.value;
      dst = std::move(decl);
      break;
    }
    case SOCK_OBJECT:
      /* TODO: What happens to default values of data-block sockets? */
      dst = std::make_unique<decl::Object>();
      break;
    case SOCK_IMAGE:
      dst = std::make_unique<decl::Image>();
      break;
    case SOCK_GEOMETRY:
      dst = std::make_unique<decl::Geometry>();
      break;
    case SOCK_COLLECTION:
      dst = std::make_unique<decl::Collection>();
      break;
    case SOCK_TEXTURE:
      dst = std::make_unique<decl::Texture>();
      break;
    case SOCK_MATERIAL:
      dst = std::make_unique<decl::Material>();
      break;
  }
  dst->name_ = io_socket.name;
  dst->identifier_ = io_socket.identifier;
  dst->in_out_ = eNodeSocketInOut(io_socket.in_out);
  dst->description_ = io_socket.description;
}

bool node_group_declare_dynamic_fn(const bNodeTree & /*node_tree*/,
                                   const bNode &node,
                                   NodeDeclaration &r_declaration)
{
  if (!node.id) {
    return false;
  }
  else if (ID_IS_LINKED(node.id) && (node.id->tag & LIB_TAG_MISSING)) {
    /* TODO: Restore the behavior that keeps the sockets until the ID is found. */
    return false;
  }
  const bNodeTree &group = *reinterpret_cast<const bNodeTree *>(node.id);

  Vector<SocketDeclarationPtr> inputs;
  Vector<SocketDeclarationPtr> outputs;

  /* TODO: Specialize for geometry nodes and fields. */
  /* TODO: Figure out how this should work for custom node trees / #SOCK_CUSTOM. */
  LISTBASE_FOREACH (const bNodeSocket *, input, &group.inputs) {
    inputs.append(declataion_for_interface_socket(*input));
  }
  LISTBASE_FOREACH (const bNodeSocket *, output, &group.outputs) {
    outputs.append(declataion_for_interface_socket(*output));
  }

  if (!ID_IS_LINKED(node.id)) {
    /* TODO: The if statement is a fun possibility, but maybe not worth it right now? */
    std::make_unique<decl::Extend>();
    b.add_input<decl::Extend>("__extend__");
    b.add_output<decl::Extend>("__extend__");
  }

  return true;
}

}  // namespace blender::nodes

/** \} */

/* -------------------------------------------------------------------- */
/** \name Node Frame
 * \{ */

static void node_frame_init(bNodeTree * /*ntree*/, bNode *node)
{
  NodeFrame *data = MEM_cnew<NodeFrame>("frame node storage");
  node->storage = data;

  data->flag |= NODE_FRAME_SHRINK;

  data->label_size = 20;
}

void register_node_type_frame()
{
  /* frame type is used for all tree types, needs dynamic allocation */
  bNodeType *ntype = MEM_cnew<bNodeType>("frame node type");
  ntype->free_self = (void (*)(bNodeType *))MEM_freeN;

  node_type_base(ntype, NODE_FRAME, "Frame", NODE_CLASS_LAYOUT);
  ntype->initfunc = node_frame_init;
  node_type_storage(ntype, "NodeFrame", node_free_standard_storage, node_copy_standard_storage);
  node_type_size(ntype, 150, 100, 0);
  ntype->flag |= NODE_BACKGROUND;

  nodeRegisterType(ntype);
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Node Re-Route
 * \{ */

static void node_reroute_init(bNodeTree *ntree, bNode *node)
{
  /* NOTE: Cannot use socket templates for this, since it would reset the socket type
   * on each file read via the template verification procedure.
   */
  nodeAddStaticSocket(ntree, node, SOCK_IN, SOCK_RGBA, PROP_NONE, "Input", "Input");
  nodeAddStaticSocket(ntree, node, SOCK_OUT, SOCK_RGBA, PROP_NONE, "Output", "Output");
}

void register_node_type_reroute()
{
  /* frame type is used for all tree types, needs dynamic allocation */
  bNodeType *ntype = MEM_cnew<bNodeType>("frame node type");
  ntype->free_self = (void (*)(bNodeType *))MEM_freeN;

  node_type_base(ntype, NODE_REROUTE, "Reroute", NODE_CLASS_LAYOUT);
  ntype->initfunc = node_reroute_init;

  nodeRegisterType(ntype);
}

static void propagate_reroute_type_from_start_socket(
    bNodeSocket *start_socket,
    const MultiValueMap<bNodeSocket *, bNodeLink *> &links_map,
    Map<bNode *, const bNodeSocketType *> &r_reroute_types)
{
  Stack<bNode *> nodes_to_check;
  for (bNodeLink *link : links_map.lookup(start_socket)) {
    if (link->tonode->type == NODE_REROUTE) {
      nodes_to_check.push(link->tonode);
    }
    if (link->fromnode->type == NODE_REROUTE) {
      nodes_to_check.push(link->fromnode);
    }
  }
  const bNodeSocketType *current_type = start_socket->typeinfo;
  while (!nodes_to_check.is_empty()) {
    bNode *reroute_node = nodes_to_check.pop();
    BLI_assert(reroute_node->type == NODE_REROUTE);
    if (r_reroute_types.add(reroute_node, current_type)) {
      for (bNodeLink *link : links_map.lookup((bNodeSocket *)reroute_node->inputs.first)) {
        if (link->fromnode->type == NODE_REROUTE) {
          nodes_to_check.push(link->fromnode);
        }
      }
      for (bNodeLink *link : links_map.lookup((bNodeSocket *)reroute_node->outputs.first)) {
        if (link->tonode->type == NODE_REROUTE) {
          nodes_to_check.push(link->tonode);
        }
      }
    }
  }
}

void ntree_update_reroute_nodes(bNodeTree *ntree)
{
  /* Contains nodes that are linked to at least one reroute node. */
  Set<bNode *> nodes_linked_with_reroutes;
  /* Contains all links that are linked to at least one reroute node. */
  MultiValueMap<bNodeSocket *, bNodeLink *> links_map;
  /* Build acceleration data structures for the algorithm below. */
  LISTBASE_FOREACH (bNodeLink *, link, &ntree->links) {
    if (link->fromsock == nullptr || link->tosock == nullptr) {
      continue;
    }
    if (link->fromnode->type != NODE_REROUTE && link->tonode->type != NODE_REROUTE) {
      continue;
    }
    if (link->fromnode->type != NODE_REROUTE) {
      nodes_linked_with_reroutes.add(link->fromnode);
    }
    if (link->tonode->type != NODE_REROUTE) {
      nodes_linked_with_reroutes.add(link->tonode);
    }
    links_map.add(link->fromsock, link);
    links_map.add(link->tosock, link);
  }

  /* Will contain the socket type for every linked reroute node. */
  Map<bNode *, const bNodeSocketType *> reroute_types;

  /* Propagate socket types from left to right. */
  for (bNode *start_node : nodes_linked_with_reroutes) {
    LISTBASE_FOREACH (bNodeSocket *, output_socket, &start_node->outputs) {
      propagate_reroute_type_from_start_socket(output_socket, links_map, reroute_types);
    }
  }

  /* Propagate socket types from right to left. This affects reroute nodes that haven't been
   * changed in the loop above. */
  for (bNode *start_node : nodes_linked_with_reroutes) {
    LISTBASE_FOREACH (bNodeSocket *, input_socket, &start_node->inputs) {
      propagate_reroute_type_from_start_socket(input_socket, links_map, reroute_types);
    }
  }

  /* Actually update reroute nodes with changed types. */
  for (const auto item : reroute_types.items()) {
    bNode *reroute_node = item.key;
    const bNodeSocketType *socket_type = item.value;
    bNodeSocket *input_socket = (bNodeSocket *)reroute_node->inputs.first;
    bNodeSocket *output_socket = (bNodeSocket *)reroute_node->outputs.first;

    if (input_socket->typeinfo != socket_type) {
      nodeModifySocketType(ntree, reroute_node, input_socket, socket_type->idname);
    }
    if (output_socket->typeinfo != socket_type) {
      nodeModifySocketType(ntree, reroute_node, output_socket, socket_type->idname);
    }
  }
}

bool BKE_node_is_connected_to_output(const bNodeTree *ntree, const bNode *node)
{
  ntree->ensure_topology_cache();
  Stack<const bNode *> nodes_to_check;
  for (const bNodeSocket *socket : node->output_sockets()) {
    for (const bNodeLink *link : socket->directly_linked_links()) {
      nodes_to_check.push(link->tonode);
    }
  }
  while (!nodes_to_check.is_empty()) {
    const bNode *next_node = nodes_to_check.pop();
    for (const bNodeSocket *socket : next_node->output_sockets()) {
      for (const bNodeLink *link : socket->directly_linked_links()) {
        if (link->tonode->typeinfo->nclass == NODE_CLASS_OUTPUT &&
            link->tonode->flag & NODE_DO_OUTPUT) {
          return true;
        }
        nodes_to_check.push(link->tonode);
      }
    }
  }

  return false;
}

/** \} */

/* -------------------------------------------------------------------- */
/** \name Node #GROUP_INPUT / #GROUP_OUTPUT
 * \{ */

bNodeSocket *node_group_input_find_socket(bNode *node, const char *identifier)
{
  bNodeSocket *sock;
  for (sock = (bNodeSocket *)node->outputs.first; sock; sock = sock->next) {
    if (STREQ(sock->identifier, identifier)) {
      return sock;
    }
  }
  return nullptr;
}

void register_node_type_group_input()
{
  /* used for all tree types, needs dynamic allocation */
  bNodeType *ntype = MEM_cnew<bNodeType>("node type");
  ntype->free_self = (void (*)(bNodeType *))MEM_freeN;

  node_type_base(ntype, NODE_GROUP_INPUT, "Group Input", NODE_CLASS_INTERFACE);
  node_type_size(ntype, 140, 80, 400);
  /* TODO: Update declaration when linking to the extension sockets. */
  // ntype->declare_dynamic =

  nodeRegisterType(ntype);
}

bNodeSocket *node_group_output_find_socket(bNode *node, const char *identifier)
{
  bNodeSocket *sock;
  for (sock = (bNodeSocket *)node->inputs.first; sock; sock = sock->next) {
    if (STREQ(sock->identifier, identifier)) {
      return sock;
    }
  }
  return nullptr;
}

void register_node_type_group_output()
{
  /* used for all tree types, needs dynamic allocation */
  bNodeType *ntype = MEM_cnew<bNodeType>("node type");
  ntype->free_self = (void (*)(bNodeType *))MEM_freeN;

  node_type_base(ntype, NODE_GROUP_OUTPUT, "Group Output", NODE_CLASS_INTERFACE);
  node_type_size(ntype, 140, 80, 400);
  /* TODO: Update declaration when linking to the extension sockets. */
  // ntype->declare_dynamic = //;

  ntype->no_muting = true;

  nodeRegisterType(ntype);
}

/** \} */
