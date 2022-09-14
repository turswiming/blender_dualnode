/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright Blender Foundation. All rights reserved. */

/** \file
 * \ingroup bke
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CLG_log.h"

#include "MEM_guardedalloc.h"

/* types */
#include "DNA_collection_types.h"
#include "DNA_curve_types.h"
#include "DNA_lattice_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_object_force_types.h"
#include "DNA_object_types.h"
#include "DNA_scene_types.h"

#include "BLI_ghash.h"
#include "BLI_listbase.h"
#include "BLI_math.h"
#include "BLI_threads.h"
#include "BLI_utildefines.h"

#include "BLI_vector.hh"
#include "BLI_map.hh"
#include "BLI_set.hh"
#include <bits/stdc++.h>

#include "BKE_collection.h"
#include "BKE_collision.h"
#include "BKE_curve.h"
#include "BKE_deform.h"
#include "BKE_effect.h"
#include "BKE_global.h"
#include "BKE_layer.h"
#include "BKE_mesh.h"
#include "BKE_modifier.h"
#include "BKE_pointcache.h"
#include "BKE_scene.h"
#include "BKE_softbody.h"

#include "DEG_depsgraph.h"
#include "DEG_depsgraph_query.h"

#include "PIL_time.h"

#include "../../simulation/intern/xpbd.h"

using blender::Set;
// using blender::Map;
using blender::Vector;
using namespace std;

struct VectorHasher {
    int operator()(const Vector<int> &V) const {
        int hash = V.size();
        for(auto &i : V) {
            hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

static CLG_LogRef LOG = {"bke.softbody"};

/* callbacks for errors and interrupts and some goo */
static int (*SB_localInterruptCallBack)(void) = NULL;

int tet_face[4][4] = {{2,1,0,3}, {0,1,3,2}, {1,2,3,0}, {2,0,3,1}}; 

SoftBody *init_softbody()
{
  SoftBody *sb;

  sb = (SoftBody *)MEM_callocN(sizeof(SoftBody), "softbody");

  sb->totpoint = 0;
  sb->tottet = 0;
  sb->tot_surface_point = 0;
  sb->tot_surface_tet = 0;
  sb->bpoint = NULL;
  sb->btet = NULL;
  sb->surface_points = NULL;
  sb->surface_tets = NULL;

  sb->dt = 1.0/30.0;
  sb->alpha_edge = 0.1;
  sb->alpha_vol = 0;
  sb->substep_count = 50;
  sb->grav = 9.8f;

  sb->shared = (SoftBody_Shared *)MEM_callocN(sizeof(*sb->shared), "SoftBody_Shared");
  sb->shared->pointcache = BKE_ptcache_add(&sb->shared->ptcaches);

  sb->last_frame = MINFRAME - 1;

  return sb;
}

void free_softbody_intern(SoftBody *sb){
  if(sb->bpoint){
    MEM_freeN(sb->bpoint);
  }
  if(sb->bedge){
    MEM_freeN(sb->bedge);
  }
  if(sb->btet){
    MEM_freeN(sb->btet);
  }

  sb->bpoint = NULL;
  sb->bedge = NULL;
  sb->btet = NULL;

  sb->totpoint = 0;
  sb->tottet = 0;

  if(sb->surface_points){
    MEM_freeN(sb->surface_points);
  }
  if(sb->surface_tets){
    MEM_freeN(sb->surface_tets);
  }
  sb->surface_points = NULL;
  sb->surface_tets = NULL;
  sb->tot_surface_point = 0;
  sb->tot_surface_tet = 0;
}

void sbFree(Object *ob)
{
  SoftBody *sb = ob->soft;
  if (sb == NULL) {
    return;
  }

  const bool is_orig = (ob->id.tag & LIB_TAG_COPIED_ON_WRITE) == 0;

  free_softbody_intern(sb);

  if (is_orig) {
    /* Only free shared data on non-CoW copies */
    BKE_ptcache_free_list(&sb->shared->ptcaches);
    sb->shared->pointcache = NULL;
    MEM_freeN(sb->shared);
  }
  // if (sb->effector_weights) {
  //   MEM_freeN(sb->effector_weights);
  // }
  MEM_freeN(sb);

  ob->soft = NULL;
}

void mesh_to_softbody(Object *ob, float (*vertexCos)[3], int numVerts){
  SoftBody *sb;
  Mesh *me = (Mesh *)ob->data;
  MLoop *mloop = me->mloop;
  BodyPoint *bpoint;
  BodyTet *btet;
  int numTets = me->totpoly;

  if(ob->soft == NULL){
    ob->soft = init_softbody();
  }
  else{
    free_softbody_intern(ob->soft);
  }
  sb = ob->soft;

  // Set bodypoints
  sb->bpoint = (BodyPoint *)MEM_mallocN(numVerts * sizeof(BodyPoint), "bodypoint");
  for(int i = 0; i<numVerts; i++){
    
    float temp[3];
    copy_v3_v3(temp, vertexCos[i]);

    copy_v3_v3(sb->bpoint[i].x, vertexCos[i]);
    mul_m4_v3(ob->obmat, sb->bpoint[i].x);
    copy_v3_v3(sb->bpoint[i].x_ini, sb->bpoint[i].x);
    sb->bpoint[i].mass_inv = 0.0f;

    copy_v3_fl(sb->bpoint[i].v, 0.0);
    copy_v3_fl(sb->bpoint[i].v_ini, 0.0);
    copy_v3_fl(sb->bpoint[i].a, 0.0);
    copy_v3_fl(sb->bpoint[i].a_ini, 0.0);
  }
  sb->totpoint = numVerts;
  bpoint = sb->bpoint;

  // Set tets
  sb->btet = (BodyTet *)MEM_mallocN(numTets * sizeof(BodyTet), "bodytet");
  btet = sb->btet;
  for(int i = 0; i<numTets; i++){
    MPoly *face = &me->mpoly[i];

    for(int vert = 0; vert<4; vert++){
      btet[i].verts[vert] = mloop[face->loopstart + vert].v;
    }

    BodyTet *curr_tet = &sb->btet[i];
    curr_tet->initial_volume = get_tet_volume(bpoint, curr_tet);

    for(int vert = 0; vert<4; vert++){
      bpoint[btet[i].verts[vert]].mass_inv += 4.0/curr_tet->initial_volume; 
    }
  }
  sb->tottet = numTets;

  // Set edges
  Set<pair<int,int>> edge_set;
  for(int tetnr = 0; tetnr < sb->tottet; tetnr++){
    // Brute force, iterating over all vertex pairs for a tet
    for(int verti = 0; verti<4; verti++){
        for(int vertj = verti+1; vertj<4; vertj++){
            int vert1 = btet[tetnr].verts[verti];
            int vert2 = btet[tetnr].verts[vertj];

            int temp = min(vert1, vert2);
            vert2 = max(vert1, vert2);
            vert1 = temp;

            if(edge_set.contains({vert1, vert2})){
                continue;
            }

            edge_set.add({vert1, vert2});
        }
    }
  }
  sb->bedge = (BodyEdge *)MEM_mallocN(edge_set.size() * sizeof(BodyEdge), "bodyedge");
  int i = 0; 
  for(auto it : edge_set){
    sb->bedge[i].u = it.first;
    sb->bedge[i].v = it.second;
    i++;
  }
  sb->totedge = edge_set.size();

  // Set boundary points and tets
  unordered_map<Vector<int>, int, VectorHasher> face_tet;
  for(int tetnr = 0; tetnr < sb->tottet; tetnr++){
    for(int facenr = 0; facenr<4; facenr++){
      Vector<int> temp;
      temp.append(btet[tetnr].verts[tet_face[facenr][0]]);
      temp.append(btet[tetnr].verts[tet_face[facenr][1]]);
      temp.append(btet[tetnr].verts[tet_face[facenr][2]]);

      sort(temp.begin(), temp.end());

      if(face_tet.find(temp) != face_tet.end()){
        face_tet.erase(temp);
        continue;
      }

      face_tet[temp] = tetnr;
    }
  }

  Set<int> surface_point_set;
  Set<int> surface_tet_set;
  // sb->surface_tets = (int *)MEM_mallocN(face_tet.size() * sizeof(int), "surface_tets");
  // sb->tot_surface_tet = face_tet.size();

  int surface_tet_nr = 0;
  for(auto it : face_tet){
    for(int i = 0; i<3; i++){
      surface_point_set.add(it.first[i]);
    }

    // sb->surface_tets[surface_tet_nr++] = it.second;
    surface_tet_set.add(it.second);
  }

  sb->surface_tets = (int *)MEM_mallocN(face_tet.size() * sizeof(int), "surface_tets");
  sb->tot_surface_tet = surface_tet_set.size();

  for(auto it : surface_tet_set){
    sb->surface_tets[surface_tet_nr++] = it;
  }

  sb->surface_points = (int *)MEM_mallocN(surface_point_set.size() * sizeof(int), "surface points");
  sb->tot_surface_point = surface_point_set.size();
  int surface_point_nr = 0;
  for(auto it : surface_point_set){
    sb->surface_points[surface_point_nr++] = it;
  }
}

void softbody_to_object(Object *ob, float (*vertexCos)[3], int numVerts){
  invert_m4_m4(ob->imat, ob->obmat);
  for(int i = 0; i<numVerts; i++){
    float temp[3];
    copy_v3_v3(temp, ob->soft->bpoint[i].x);

    copy_v3_v3(vertexCos[i], ob->soft->bpoint[i].x);
    mul_m4_v3(ob->imat, vertexCos[i]);
  }
}

void sb_store_last_frame(struct Depsgraph *depsgraph, Object *object, float framenr)
{
  if (!DEG_is_active(depsgraph)) {
    return;
  }
  Object *object_orig = DEG_get_original_object(object);
  object->soft->last_frame = framenr;
  object_orig->soft->last_frame = framenr;
}

void sbObjectStep(struct Depsgraph *depsgraph,
                  Scene *scene,
                  Object *ob,
                  float cfra,
                  float (*vertexCos)[3],
                  int numVerts){
  SoftBody *sb = ob->soft;
  PointCache *cache;
  PTCacheID pid;
  float dtime, timescale;
  int framedelta, framenr, startframe, endframe;
  int cache_result;
  cache = sb->shared->pointcache;

  framenr = (int)cfra;
  framedelta = framenr - cache->simframe;

  BKE_ptcache_id_from_softbody(&pid, ob, sb);
  BKE_ptcache_id_time(&pid, scene, framenr, &startframe, &endframe, &timescale);

  /* check for changes in mesh, should only happen in case the mesh
   * structure changes during an animation */
  if (sb->bpoint && numVerts != sb->totpoint) {
    BKE_ptcache_invalidate(cache);
    return;
  }

  /* clamp frame ranges */
  if (framenr < startframe) {
    BKE_ptcache_invalidate(cache);
    return;
  }
  if (framenr > endframe) {
    framenr = endframe;
  }

  /* verify if we need to create the softbody data */
  if(sb == NULL || sb->bpoint == NULL){
    mesh_to_softbody(ob, vertexCos, numVerts);
  }

  /* still no points? go away */
  if (sb->totpoint == 0) {
    return;
  }
  if(framenr == startframe){
    BKE_ptcache_id_reset(scene, &pid, PTCACHE_RESET_OUTDATED);

    // shift code that updates sb data and reinitialize it here

    BKE_ptcache_validate(cache, framenr);
    cache->flag &= ~PTCACHE_REDO_NEEDED;

    sb_store_last_frame(depsgraph, ob, framenr);

    return;
  }

  /* try to read from cache */
  bool can_write_cache = DEG_is_active(depsgraph);
  bool can_simulate = (framenr == sb->last_frame + 1) && !(cache->flag & PTCACHE_BAKED) &&
                      can_write_cache;

  cache_result = BKE_ptcache_read(&pid, (float)framenr + scene->r.subframe, can_simulate);

  if (cache_result == PTCACHE_READ_EXACT || cache_result == PTCACHE_READ_INTERPOLATED ||
      (!can_simulate && cache_result == PTCACHE_READ_OLD)) {
    softbody_to_object(ob, vertexCos, numVerts);

    BKE_ptcache_validate(cache, framenr);

    if (cache_result == PTCACHE_READ_INTERPOLATED && cache->flag & PTCACHE_REDO_NEEDED &&
        can_write_cache) {
      BKE_ptcache_write(&pid, framenr);
    }

    sb_store_last_frame(depsgraph, ob, framenr);

    return;
  }

  /* if on second frame, write cache for first frame */
  if (cache->simframe == startframe &&
      (cache->flag & PTCACHE_OUTDATED || cache->last_exact == 0)) {
    BKE_ptcache_write(&pid, startframe);
  }

  for(int i = 0; i<sb->substep_count; i++){
    xpbd_position_update(sb);

    xpbd_enforce_constraints(sb);

    xpbd_velocity_update(sb);
  }

  softbody_to_object(ob, vertexCos, numVerts);

  BKE_ptcache_validate(cache, framenr);
  BKE_ptcache_write(&pid, framenr);

  sb_store_last_frame(depsgraph, ob, framenr);
}