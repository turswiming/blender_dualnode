/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2005 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup modifiers
 *
 * Weld modifier: Remove doubles.
 */

/* TODOs:
 * - Review weight and vertex color interpolation.;
 */

#include "MEM_guardedalloc.h"

#include "BLI_utildefines.h"

#include "BLI_array.hh"
#include "BLI_index_range.hh"
#include "BLI_span.hh"

#include "BLI_math_vector.h"
#include "BLI_utildefines.h"

#include "BLI_vector.hh"
#include "BLI_math_vec_types.hh"
#include "BLI_math_vector.hh"

// #include "bmesh_construct.h"
#include "bmesh.h"
#include "bmesh_tools.h"

#include "BLT_translation.h"

#include "DNA_defaults.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"
#include "DNA_screen_types.h"

#include "BKE_bvhutils.h"
#include "BKE_context.h"
#include "BKE_deform.h"
#include "BKE_modifier.h"
#include "BKE_screen.h"
#include "BKE_mesh.h"

#include "UI_interface.h"
#include "UI_resources.h"

#include "RNA_access.h"
#include "RNA_prototypes.h"

#include "DEG_depsgraph.h"

#include "MOD_modifiertypes.h"
#include "MOD_ui_common.h"

#include "GEO_mesh_merge_by_distance.hh"

using blender::Array;
using blender::IndexMask;
using blender::Span;
using blender::Vector;
using namespace blender;
using namespace blender::math;
// using namespace std;


Vector<Vector<int>> tetFaces({{2,1,0}, {0,1,3}, {1,2,3}, {2,0,3}});
float directions[6][3] = {{1,0,0}, {-1,0,0}, {0,1,0}, {0,-1,0}, {0,0,1}, {0,0,-1}};

float inf = FLT_MAX;// C has no FLOAT_MAX :/
float eps = 1e-4;

static float randomEps(){
  float eps = 1e-6 ;
  return -eps + 2.0 * (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * eps;
}

static bool isInside(float3 vert, BVHTreeFromMesh *treedata){
  int count =  0;
  float min_dist = 0.0f;
  for(auto dir : directions){
    float radius = 0.0f;
    
    float max_length = FLT_MAX;

    BVHTreeRayHit rayhit = {0};
    rayhit.index = -1;
    rayhit.dist = max_length;
    BLI_bvhtree_ray_cast(treedata->tree, vert, dir, radius, &rayhit, treedata->raycast_callback, treedata);

    if (rayhit.index != -1 && rayhit.dist <= max_length) {
      if(dot_v3v3(rayhit.no, dir) > eps){
        count++;
      }

      if((min_dist > eps) && (min_dist - rayhit.dist) > eps)
        return false;
    }
  }

  return count > 3;
}

// Used to populate a tets face normals and planesD
static void setTetProperties(Vector<float3> &verts, 
                      Vector<int> &tetVertId,
                      Vector<float3> &faceNormals, 
                      Vector<float> &planesD,
                      int tetNr){
  for(int i = 0; i<4; i++){ 
    float3 p0 = verts[tetVertId[4*tetNr + tetFaces[i][0]]];
    float3 p1 = verts[tetVertId[4*tetNr + tetFaces[i][1]]];
    float3 p2 = verts[tetVertId[4*tetNr + tetFaces[i][2]]];

    float3 normal = cross(p1 - p0, p2 - p0);
    normal = normalize(normal);

    faceNormals[4*tetNr + i] = normal;
    planesD[4*tetNr + i] = dot(p0, normal);
  }
}

static float3 getCircumCenter(float3 p0, float3 p1, float3 p2, float3 p3){
  float3 b = p1 - p0;
  float3 c = p2 - p0;
  float3 d = p3 - p0;

  float det = 2.0 * (b.x*(c.y*d.z - c.z*d.y) - b.y*(c.x*d.z - c.z*d.x) + b.z*(c.x*d.y - c.y*d.x));
  if (det <= eps){
    return p0;
  }
  else{
    float3 v = cross(c, d)*dot(b, b) + cross(d, b)*dot(c, c) + cross(b, c)*dot(d, d);
    v /= det;
    return p0 + v;
  }
      
}

// bool isSameSide(float3 p0, float3 p1, float3 p2, float3 p4, float3 vert){

// }

// bool isInsideTet(float3 p0, float3 p1, float3 p2, float3 p3, float3 vert){

// }

static int findContainingTet(Vector<float3> &verts, 
                    Vector<int> &tetVertId, 
                    Vector<int> &tetFaceNeighbors, 
                    Vector<float3> &faceNormals, 
                    Vector<float> &planesD, 
                    int tetMarkId, 
                    Vector<int> &tetMarks,
                    float3 currVert){
  
  /* --------------------
  Matthias Method
  ----------------------*/

  // bool found = false;
  // int tetNr = 0;

  // while(tetNr < tetVertId.size()/4 && tetVertId[4*tetNr]<0)
  //   tetNr++;

  // float3 center(0.0, 0.0, 0.0);
  // while(!found){
  //   if(tetNr < 0 || tetMarks[tetNr] == tetMarkId){
  //     break;
  //   }
    
  //   tetMarks[tetNr] = tetMarkId;

  //   center = {0.0,0.0,0.0};
  //   for(int i = 0; i<4; i++){
  //     center += verts[tetVertId[4*tetNr + i]];
  //   }
  //   center*=0.25;

  //   float minT = inf; //
  //   int minFaceNr = -1;

  //   for(int i = 0; i<4; i++){
  //     float3 normal = faceNormals[4*tetNr + i];
  //     float d = planesD[4*tetNr + i];

  //     float hp = dot(normal, currVert) - d;
  //     float hc = dot(normal, center) - d;

  //     float t = hp - hc;
  //     if(t <= eps)
  //       continue;

  //     t = -hc/t;

  //     if((t >= eps) && ((t - minT) > eps)){
  //       minT = t;
  //       minFaceNr = i;
  //     }
  //   }

  //   if(minT >= 1.0){
  //     found = true;
  //   }
  //   else{
  //     tetNr = tetFaceNeighbors[4*tetNr + minFaceNr];
  //   }
  // }

  // if(found)
  //   return tetNr;
  // else
  //   return -1;

  /* --------------------
  Brute force finding first violating tet
  ----------------------*/

  for(int currTet = 0; currTet<tetVertId.size()/4; currTet++){
    if(tetVertId[4*currTet] == -1)
      continue;
    
    float3 p0 = verts[tetVertId[4*currTet + 0]];
    float3 p1 = verts[tetVertId[4*currTet + 1]];
    float3 p2 = verts[tetVertId[4*currTet + 2]];
    float3 p3 = verts[tetVertId[4*currTet + 3]];

    float3 circumCenter = getCircumCenter(p0, p1, p2, p3);
    float circumRadius = length(p0 - circumCenter);
    
    if((circumRadius - length(currVert - circumCenter)) >= eps){
      return currTet;
    }
  }

  /* -----------------------
  My method, checks if point is inside tet, if not finds the face that is connecting to the point and center
  ------------------------*/

  // bool found = false;
  // int tetNr = 0;
  // while(tetNr < tetVertId.size()/4 && tetVertId[4*tetNr] < 0)
  //   tetNr++;

  // while(!found){
  //   float3 p0 = verts[tetVertId[4*currTet + 0]];
  //   float3 p1 = verts[tetVertId[4*currTet + 1]];
  //   float3 p2 = verts[tetVertId[4*currTet + 2]];
  //   float3 p3 = verts[tetVertId[4*currTet + 3]];
    
  //   if(isInsideTet(p0, p1, p2, p3, currVert))
  //     return tetNr;
    
  //   for(int i = 0; i<4; i++)
  // }

  return -1;
}

/*
The basic assumption employed is that 2 violating tets cannot not be neighbors. 
Simple BFS approach that checks neighbors of all violating tets and adds them to stack

If a violating tet shares a face with a non-violating tet, that's a boundary face and for each violating tet, list of its boundary faces is returned

Note - BFS can be written better
*/
static Vector<std::pair<int, Vector<int>>> getViolatingTets(Vector<float3> &verts, 
                            Vector<int> &tetVertId, 
                            Vector<int> &tetFaceNeighbors, 
                            int tetMarkId, 
                            Vector<int> &tetMarks,
                            float3 currVert,
                            int containingTetNr
                            ){
  
  Vector< std::pair<int,Vector<int>> > violatingTets;
  Vector<int> stack;

  stack.append(containingTetNr);  

  while(stack.size()){
    int currTet = stack.last();
    stack.remove_last();

    if(tetMarks[currTet] == tetMarkId){
      continue;
    }
    tetMarks[currTet] = tetMarkId; 
    Vector<int> currTetBorderFaces;

    for(int i = 0; i<4; i++){
      int neighborTet = tetFaceNeighbors[4*currTet + i];
      
      if(neighborTet<0){
        currTetBorderFaces.append(i);
        continue;
      }
      if(tetMarks[neighborTet]==tetMarkId){
        continue;
      }
      
      float3 p0 = verts[tetVertId[4*neighborTet + 0]];
      float3 p1 = verts[tetVertId[4*neighborTet + 1]];
      float3 p2 = verts[tetVertId[4*neighborTet + 2]];
      float3 p3 = verts[tetVertId[4*neighborTet + 3]];

      float3 circumCenter = getCircumCenter(p0, p1, p2, p3);
      float circumRadius = length(p0 - circumCenter);

      if((circumRadius - length(currVert - circumCenter)) >= eps){
        stack.append(neighborTet);
        // tetMarks[neighborTet] = tetMarkId;
      }
      else{
        currTetBorderFaces.append(i);
      }
    }
    violatingTets.append({currTet, currTetBorderFaces});
  }

  return violatingTets;
}

static Vector<int> createTets(Vector<float3> verts, BVHTreeFromMesh *treedata, float minTetQuality){
  Vector<int> tetVertId; // Stores indices of vertex that form a tet. Every tet is stored as 4 indices 
  Vector<int> tetFaceNeighbors; // Stores index of tet that shares the face as per tetFaces order. 4 neighbors per tet, 1 for each face

  int firstFreeTet = -1;

  Vector<float3> faceNormals; // Stores the normal of each of face of the tet
  Vector<float> planesD;

  int tetMarkId = 0;
  Vector<int> tetMarks; 
  /*Used to keep track of visited tets for various processess. In each step where visited tets need to avoided, 
  on visiting them, a unique number is stored representing that particular step
  */
 
  int bigTet = verts.size() - 4;  

  for(int i = 0; i<4; i++){
    tetVertId.append(bigTet + i);
    tetFaceNeighbors.append(-1);
    faceNormals.append({0.0, 0.0, 0.0});
    planesD.append(0.0);
  }
  tetMarks.append(0);
  setTetProperties(verts, tetVertId, faceNormals, planesD, 0);

  for(int vertNr = 0; vertNr<bigTet; vertNr++){
    float3 currVert = verts[vertNr];

    // std::cout << "Adding vert " << vertNr << std::endl;

    // Find containing tet (need to understand) center can lie outside the tet?
    tetMarkId += 1;
    int containingTetNr = -1;
    containingTetNr = findContainingTet(verts, tetVertId, tetFaceNeighbors, faceNormals, planesD, tetMarkId, tetMarks, currVert);
    
    if(containingTetNr == -1){
      std::cout << "Couldn't add vert " << vertNr << std::endl;
      continue;
    }

    // Find Violating Tets - Returns all tets that are violating the Delaunay condition along with a list of face indices that are boundary faces
    // Boundary faces may be outermost face of the structure or the face shared by a violating tet and a non violating tet
    tetMarkId += 1;
    Vector<std::pair<int, Vector<int>>> violatingTets = getViolatingTets(verts, tetVertId, tetFaceNeighbors, tetMarkId, tetMarks, currVert, containingTetNr);

    //Create new tets centered at the new point and including the boundary faces
    Vector<int> newTets; // Stores tet indices of the new tets formed. Used for making neighbors of the new tets formed
    for(int violatingTetNr = 0; violatingTetNr<violatingTets.size(); violatingTetNr++){
      int violatingTet = violatingTets[violatingTetNr].first;
      Vector<int> boundaryFaces = violatingTets[violatingTetNr].second;

      // Copying old tet information
      Vector<int> currTetVerts;
      Vector<int> currTetNeighbors;
      for(int i = 0; i<4; i++){
        currTetVerts.append(tetVertId[4*violatingTet + i]);
        currTetNeighbors.append(tetFaceNeighbors[4*violatingTet + i]);
      }

      // Deleting old tet. For each deleted tet, index 1 stores the index of next free tet
      tetVertId[4*violatingTet] = -1;
      tetVertId[4*violatingTet + 1] = firstFreeTet;
      firstFreeTet = violatingTet;

      for(int i = 0; i<boundaryFaces.size(); i++){
        Vector<int> faceVerts;
        for(int j = 2; j>=0; j--){
          faceVerts.append(currTetVerts[tetFaces[boundaryFaces[i]][j]]); 
        }

        // Make new tet
        int newTetNr = -1;
        if(firstFreeTet == -1){
          newTetNr = tetVertId.size()/4;
          for(int j = 0; j<4; j++){
            tetVertId.append(j<3 ? faceVerts[j] : vertNr);
            tetFaceNeighbors.append(-1);
            faceNormals.append({0.0,0.0,0.0});
            planesD.append(0.0);
          }

          tetMarks.append(0);
        }
        else{
          newTetNr = firstFreeTet;
          firstFreeTet = tetVertId[4*firstFreeTet + 1];
          for(int j = 0; j<3; j++){
            tetVertId[4*newTetNr + j] = faceVerts[j];
            tetFaceNeighbors[4*newTetNr + j] = -1;
          }
          tetVertId[4*newTetNr + 3] = vertNr;
          tetFaceNeighbors[4*newTetNr + 3] = -1;
          tetMarks[newTetNr] = 0;
        }
        newTets.append(newTetNr);

        tetFaceNeighbors[4*newTetNr] = currTetNeighbors[boundaryFaces[i]];
        // If the boundary face has no neighboring tet
        if(currTetNeighbors[boundaryFaces[i]] != -1){
          // Else correcting the neighbors for the shared face
          // tetFaceNeighbors[4*newTetNr] = currTetNeighbors[boundaryFaces[i]];
          for(int j = 0; j<4; j++){
            if(tetFaceNeighbors[4*currTetNeighbors[boundaryFaces[i]] + j] == violatingTet){
              tetFaceNeighbors[4*currTetNeighbors[boundaryFaces[i]] + j] = newTetNr;
              break;
            }
          }
        }

        setTetProperties(verts, tetVertId, faceNormals, planesD, newTetNr);
      }
    }

    // Setting the neighbors of internal faces of new tets
    for(int i = 0; i<newTets.size(); i++){
      for(int j = 0; j<newTets.size(); j++){
        
        if(j == i){
          continue;
        }

        for(int facei = 0; facei<4; facei++){
          int vertsI[3];
          for(int k = 0; k<3; k++){
            vertsI[k] = tetVertId[4*newTets[i] + tetFaces[facei][k]];
          }

          int count = 0;
          for(int k = 0; k<3; k++){
            for(int l = 0; l<4; l++){
              if(vertsI[k] == tetVertId[4*newTets[j] + l]){
                count++;
                break;
              }
            }
          }

          if(count == 3){
            tetFaceNeighbors[4*newTets[i] + facei] = newTets[j];
            // tetFaceNeighbors[4*newTets[i] + facei] = newTets[j];
          }
        }
      }
    }
  }

  // Remove empty tets and tets that lie outside the structure
  int emptyTet = 0;
  int tetLen = tetVertId.size()/4;
  for(int tetNr = 0; tetNr<tetLen; tetNr++){
    
    int flag = 1;
    float3 center(0.0f, 0.0f, 0.0f);
    for(int i = 0; i<4; i++){
      center += verts[tetVertId[4*tetNr + i]];
      if(tetVertId[4*tetNr + i]<0 || tetVertId[4*tetNr + i]>=bigTet){
        flag = 0;
        break;
      }
    }
    center*=0.25;

    if(!flag || !isInside(center, treedata)){
      continue;
    }

    for(int i = 0; i<4; i++){
      tetVertId[4*emptyTet + i] = tetVertId[4*tetNr + i];
    }
    emptyTet++;
  }

    tetVertId.remove(4*emptyTet, 4*(tetLen - emptyTet));

    return tetVertId;
}

static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *UNUSED(ctx), Mesh *mesh)
{
  // Parameters of tetrahedralization, need to be taken as user input, being defined here as placeholders
  float interiorResolution = 0;
  float minTetQuality = 0.001; // Exp goes from -4 to 0
  bool oneFacePerTet = false;
  float tetScale = 0.8;

  BVHTreeFromMesh treedata = {NULL};
  BKE_bvhtree_from_mesh_get(&treedata, mesh, BVHTREE_FROM_LOOPTRI, 2);

  Mesh *result;
  BMesh *bm;
  bool use_operators = true;
  
  BMeshCreateParams bmcParam;
  bmcParam.use_toolflags = true;
  bm = BM_mesh_create(&bm_mesh_allocsize_default,
                      &bmcParam);
  
  Vector<float3> tetVerts;

  float3 center(0,0,0);
  float3 bmin(inf, inf, inf);
  float3 bmax(-inf, -inf, -inf);

  // Copying surface nodes to tetVerts
  for(int i = 0; i<mesh->totvert; i++){
    tetVerts.append({0.0, 0.0, 0.0});
    tetVerts[i][0] = mesh->mvert[i].co[0] + randomEps();
    tetVerts[i][1] = mesh->mvert[i].co[1] + randomEps();
    tetVerts[i][2] = mesh->mvert[i].co[2] + randomEps();

    center += tetVerts[i]; // Can cause overflow?
    for(int axis = 0; axis<3; axis++){
      bmin[axis] = min(bmin[axis], tetVerts[i][axis]);
      bmax[axis] = max(bmax[axis], tetVerts[i][axis]);
    }
    
  } 
  center /= (float)tetVerts.size();

  // Computing radius of bounding sphere (max dist of node from center)
  float radius = 0.0;
  for(int i = 0; i<tetVerts.size(); i++){
    float dist = length(tetVerts[i] - center);
    radius = max(dist, radius);
  }
  
  // Interior sampling
  if(interiorResolution > 0.0){
    float boundLen[3];
    sub_v3_v3v3(boundLen, bmax, bmin);
    float maxBoundLen = max_axis_v3(boundLen);
    float sampleLen = maxBoundLen/interiorResolution;

    for(int xi = 0; xi<(int)(boundLen[0]/sampleLen); xi++){
      float x = bmin[0] + xi*sampleLen + randomEps();
      for(int yi = 0; yi<(int)(boundLen[1]/sampleLen); yi++){
        float y = bmin[1] + yi*sampleLen + randomEps();
        for(int zi = 0; zi<(int)(boundLen[2]/sampleLen); zi++){
          float z = bmin[2] + zi*sampleLen + randomEps(); 

          if(isInside({x,y,z}, &treedata)){
            tetVerts.append({x,y,z});  
          }
        }
      }
    }
  }

  float bigTetSize = radius*5.0;

  tetVerts.append({-bigTetSize, 0.0, -bigTetSize});
  tetVerts.append({bigTetSize, 0.0, -bigTetSize});
  tetVerts.append({0.0, bigTetSize, bigTetSize});
  tetVerts.append({0.0, -bigTetSize, bigTetSize});

  Vector<int> tetVertId = createTets(tetVerts, &treedata, minTetQuality);

  Vector<BMVert *> bmverts;
  if(oneFacePerTet){
    for(int i = 0; i<mesh->totvert; i++){
      bmverts.append(BM_vert_create(bm, mesh->mvert[i].co, NULL, BM_CREATE_NOP));
    }
    for(int i = mesh->totvert; i<tetVerts.size(); i++){
      bmverts.append(BM_vert_create(bm, tetVerts[i], NULL, BM_CREATE_NOP));
    }
  }
  else{
    // Add vertices multiple times to make the tets distinctly visible
    for(int tetNr = 0; tetNr<tetVertId.size()/4; tetNr++){
      float3 center(0,0,0);
      for(int j = 0; j<4; j++){
        center += tetVerts[tetVertId[4*tetNr + j]];
      }
      center *= 0.25;

      for(int faceNr = 0; faceNr < 4; faceNr++){
        for(int faceVertNr = 0; faceVertNr < 3; faceVertNr++){
          float3 vert = tetVerts[tetVertId[4*tetNr + tetFaces[faceNr][faceVertNr]]];
          vert = center + (vert - center)*tetScale;
          bmverts.append(BM_vert_create(bm, vert, NULL, BM_CREATE_NOP));
        }
      }
    }
  }

  int numTets = tetVertId.size()/4;
  int nr = 0;
  for(int i = 0; i<numTets; i++){
    if(oneFacePerTet){
      if(tetVertId[4*i] < 0)
        continue;
      BM_face_create_quad_tri(bm, bmverts[tetVertId[4*i + 0]], bmverts[tetVertId[4*i + 1]], bmverts[tetVertId[4*i + 2]], bmverts[tetVertId[4*i + 3]], NULL, BM_CREATE_NO_DOUBLE);
    }
    else{
      for(int faceNr = 0; faceNr<4; faceNr++){
        BM_face_create_quad_tri(bm, bmverts[nr], bmverts[nr+1], bmverts[nr+2], NULL, NULL, BM_CREATE_NO_DOUBLE);
        nr+=3;
      }
    }
  }

  CustomData_MeshMasks cd_mask_extra = {
    .vmask = CD_MASK_ORIGINDEX, .emask = CD_MASK_ORIGINDEX, .pmask = CD_MASK_ORIGINDEX};
  
  result = BKE_mesh_from_bmesh_for_eval_nomain(bm, &cd_mask_extra, mesh);
  BM_mesh_free(bm);

  return result;    
}

// --------------------------------------
// --------------------------------------

// static Span<MDeformVert> get_vertex_group(const Mesh &mesh, const int defgrp_index)
// {
//   if (defgrp_index == -1) {
//     return {};
//   }
//   const MDeformVert *vertex_group = static_cast<const MDeformVert *>(
//       CustomData_get_layer(&mesh.vdata, CD_MDEFORMVERT));
//   if (!vertex_group) {
//     return {};
//   }
//   return {vertex_group, mesh.totvert};
// }

// static Vector<int64_t> selected_indices_from_vertex_group(Span<MDeformVert> vertex_group,
//                                                           const int index,
//                                                           const bool invert)
// {
//   Vector<int64_t> selected_indices;
//   for (const int i : vertex_group.index_range()) {
//     const bool found = BKE_defvert_find_weight(&vertex_group[i], index) > 0.0f;
//     if (found != invert) {
//       selected_indices.append(i);
//     }
//   }
//   return selected_indices;
// }

// static Array<bool> selection_array_from_vertex_group(Span<MDeformVert> vertex_group,
//                                                      const int index,
//                                                      const bool invert)
// {
//   Array<bool> selection(vertex_group.size());
//   for (const int i : vertex_group.index_range()) {
//     const bool found = BKE_defvert_find_weight(&vertex_group[i], index) > 0.0f;
//     selection[i] = (found != invert);
//   }
//   return selection;
// }

// static std::optional<Mesh *> calculate_weld(const Mesh &mesh, const WeldModifierData &wmd)
// {
//   const int defgrp_index = BKE_id_defgroup_name_index(&mesh.id, wmd.defgrp_name);
//   Span<MDeformVert> vertex_group = get_vertex_group(mesh, defgrp_index);
//   const bool invert = (wmd.flag & MOD_WELD_INVERT_VGROUP) != 0;

//   if (wmd.mode == MOD_WELD_MODE_ALL) {
//     if (!vertex_group.is_empty()) {
//       Vector<int64_t> selected_indices = selected_indices_from_vertex_group(
//           vertex_group, defgrp_index, invert);
//       return blender::geometry::mesh_merge_by_distance_all(
//           mesh, IndexMask(selected_indices), wmd.merge_dist);
//     }
//     return blender::geometry::mesh_merge_by_distance_all(
//         mesh, IndexMask(mesh.totvert), wmd.merge_dist);
//   }
//   if (wmd.mode == MOD_WELD_MODE_CONNECTED) {
//     const bool only_loose_edges = (wmd.flag & MOD_WELD_LOOSE_EDGES) != 0;
//     if (!vertex_group.is_empty()) {
//       Array<bool> selection = selection_array_from_vertex_group(
//           vertex_group, defgrp_index, invert);
//       return blender::geometry::mesh_merge_by_distance_connected(
//           mesh, selection, wmd.merge_dist, only_loose_edges);
//     }
//     Array<bool> selection(mesh.totvert, true);
//     return blender::geometry::mesh_merge_by_distance_connected(
//         mesh, selection, wmd.merge_dist, only_loose_edges);
//   }

//   BLI_assert_unreachable();
//   return nullptr;
// }

// static Mesh *modifyMesh(ModifierData *md, const ModifierEvalContext *UNUSED(ctx), Mesh *mesh)
// {
// const WeldModifierData &wmd = reinterpret_cast<WeldModifierData &>(*md);

//   std::optional<Mesh *> result = calculate_weld(*mesh, wmd);
//   if (!result) {
//     return mesh;
//   }
//   return *result;
// }

static void initData(ModifierData *md)
{
  WeldModifierData *wmd = (WeldModifierData *)md;

  BLI_assert(MEMCMP_STRUCT_AFTER_IS_ZERO(wmd, modifier));

  MEMCPY_STRUCT_AFTER(wmd, DNA_struct_default_get(WeldModifierData), modifier);
}

static void requiredDataMask(Object *UNUSED(ob),
                             ModifierData *md,
                             CustomData_MeshMasks *r_cddata_masks)
{
  WeldModifierData *wmd = (WeldModifierData *)md;

  /* Ask for vertexgroups if we need them. */
  if (wmd->defgrp_name[0] != '\0') {
    r_cddata_masks->vmask |= CD_MASK_MDEFORMVERT;
  }
}

static void panel_draw(const bContext *UNUSED(C), Panel *panel)
{
  uiLayout *layout = panel->layout;

  PointerRNA ob_ptr;
  PointerRNA *ptr = modifier_panel_get_property_pointers(panel, &ob_ptr);
  int weld_mode = RNA_enum_get(ptr, "mode");

  uiLayoutSetPropSep(layout, true);

  uiItemR(layout, ptr, "mode", 0, nullptr, ICON_NONE);
  uiItemR(layout, ptr, "merge_threshold", 0, IFACE_("Distance"), ICON_NONE);
  if (weld_mode == MOD_WELD_MODE_CONNECTED) {
    uiItemR(layout, ptr, "loose_edges", 0, nullptr, ICON_NONE);
  }
  modifier_vgroup_ui(layout, ptr, &ob_ptr, "vertex_group", "invert_vertex_group", nullptr);

  modifier_panel_end(layout, ptr);
}

static void panelRegister(ARegionType *region_type)
{
  modifier_panel_register(region_type, eModifierType_Weld, panel_draw);
}

ModifierTypeInfo modifierType_Weld = {
    /* name */ "Weld",
    /* structName */ "WeldModifierData",
    /* structSize */ sizeof(WeldModifierData),
    /* srna */ &RNA_WeldModifier,
    /* type */ eModifierTypeType_Constructive,
    /* flags */
    (ModifierTypeFlag)(eModifierTypeFlag_AcceptsMesh | eModifierTypeFlag_SupportsMapping |
                       eModifierTypeFlag_SupportsEditmode | eModifierTypeFlag_EnableInEditmode |
                       eModifierTypeFlag_AcceptsCVs),
    /* icon */ ICON_AUTOMERGE_OFF, /* TODO: Use correct icon. */

    /* copyData */ BKE_modifier_copydata_generic,

    /* deformVerts */ nullptr,
    /* deformMatrices */ nullptr,
    /* deformVertsEM */ nullptr,
    /* deformMatricesEM */ nullptr,
    /* modifyMesh */ modifyMesh,
    /* modifyGeometrySet */ nullptr,

    /* initData */ initData,
    /* requiredDataMask */ requiredDataMask,
    /* freeData */ nullptr,
    /* isDisabled */ nullptr,
    /* updateDepsgraph */ nullptr,
    /* dependsOnTime */ nullptr,
    /* dependsOnNormals */ nullptr,
    /* foreachIDLink */ nullptr,
    /* foreachTexLink */ nullptr,
    /* freeRuntimeData */ nullptr,
    /* panelRegister */ panelRegister,
    /* blendWrite */ nullptr,
    /* blendRead */ nullptr,
};
