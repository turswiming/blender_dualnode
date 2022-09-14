#include <bits/stdc++.h>

#include "BKE_softbody.h"
#include "DNA_object_force_types.h"

#include "BLI_math.h"   

#include "BLI_set.hh"

#include "xpbd.h"

using blender::Set;
using namespace std;

int tet_faces[4][4] = {{2,1,0,3}, {0,1,3,2}, {1,2,3,0}, {2,0,3,1}}; 
// Why does tet_faces variable name not work?
// Stores indices permutations that make up a face in a cyclic order. The 4th indice of each face is the point not in the face

float get_tet_volume(BodyPoint *bpoint, BodyTet *curr_tet){
  float diff1[3], diff2[3], diff3[3];

  float vert0[3];
  float vert1[3];
  float vert2[3];
  float vert3[3];

  copy_v3_v3(vert0, bpoint[curr_tet->verts[0]].x);
  copy_v3_v3(vert1, bpoint[curr_tet->verts[1]].x);
  copy_v3_v3(vert2, bpoint[curr_tet->verts[2]].x);
  copy_v3_v3(vert3, bpoint[curr_tet->verts[3]].x);

  sub_v3_v3v3(diff1, bpoint[curr_tet->verts[1]].x, bpoint[curr_tet->verts[0]].x);
  sub_v3_v3v3(diff2, bpoint[curr_tet->verts[2]].x, bpoint[curr_tet->verts[0]].x);
  sub_v3_v3v3(diff3, bpoint[curr_tet->verts[3]].x, bpoint[curr_tet->verts[0]].x);
  
  float cross[3];
  cross_v3_v3v3(cross, diff1, diff2);

  return dot_v3v3(cross, diff3);
}

void xpbd_position_update(SoftBody *sb){
    BodyPoint *bpoint_arr = sb->bpoint;
    BodyTet *btet_arr = sb->btet;

    int totpoint = sb->totpoint;
    int tottet = sb->tottet;

    float sdt = sb->dt/sb->substep_count;

    for(int i = 0; i<totpoint; i++){
        if(bpoint_arr[i].mass_inv == 0.0f){ 
            continue;
        }

        copy_v3_v3(bpoint_arr[i].x_prev, bpoint_arr[i].x);
        float temp[3];
        copy_v3_fl3(temp, 0.0, 0.0, -sb->grav);
        // add_v3_v3(temp, bpoint_arr[i].a);
        mul_v3_fl(temp, sdt);
        add_v3_v3(bpoint_arr[i].v, temp);

        copy_v3_v3(temp, bpoint_arr[i].v);
        mul_v3_fl(temp, sdt);
        add_v3_v3(bpoint_arr[i].x, temp);

        if(bpoint_arr[i].x[2] < 0.0){
            // copy_v3_v3(bpoint_arr[i].x, bpoint_arr[i].x_prev);
            bpoint_arr[i].x[2] = 0.0;
        }
    }
}

void xpbd_velocity_update(SoftBody *sb){
    int totpoint = sb->totpoint;
    BodyPoint *bpoint_arr = sb->bpoint;

    float sdt = sb->dt/sb->substep_count;

    for(int i = 0; i<totpoint; i++){
        if(bpoint_arr[i].mass_inv == 0.0f){
            continue;
        }

        sub_v3_v3v3(bpoint_arr[i].v, bpoint_arr[i].x, bpoint_arr[i].x_prev);
        mul_v3_fl(bpoint_arr[i].v, 1.0/sdt);
    }
}

void xpbd_distance_constraint_for_edge(BodyPoint *p1, BodyPoint *p2, float compliance){
    float ini_length = len_v3v3(p1->x_ini, p2->x_ini);
    float curr_length = len_v3v3(p1->x, p2->x);

    float lambda = (curr_length - ini_length)/(curr_length * (p1->mass_inv + p2->mass_inv + compliance));
    float dx_dir[3];
    sub_v3_v3v3(dx_dir, p1->x, p2->x);
    mul_v3_fl(dx_dir, lambda);  

    float dx_p1[3] = {1.0f, 1.0f, 1.0f};
    mul_v3_v3fl(dx_p1, dx_dir, p1->mass_inv);
    sub_v3_v3(p1->x, dx_p1);
                    
    float dx_p2[3] = {1.0f, 1.0f, 1.0f};
    mul_v3_v3fl(dx_p2, dx_dir, -p2->mass_inv);
    sub_v3_v3(p2->x, dx_p2); 

    float new_length = len_v3v3(p1->x, p2->x);
}

void xpbd_enforce_distance_constraint(SoftBody *sb){
    BodyPoint *bpoint_arr = sb->bpoint;
    BodyEdge *bedge_arr = sb->bedge;
    int totedge = sb->totedge;

    float sdt = sb->dt/sb->substep_count;

    for(int i = 0; i<totedge; i++){
        xpbd_distance_constraint_for_edge(&bpoint_arr[bedge_arr[i].u], &bpoint_arr[bedge_arr[i].v], sb->alpha_edge/(sdt*sdt));
    }
}

void xpbd_enforce_volume_constraint(SoftBody *sb){
    BodyPoint *bpoint_arr = sb->bpoint;
    BodyTet *btet_arr = sb->btet;

    float sdt = sb->dt/sb->substep_count;

    int totpoint = sb->totpoint;
    int tottet = sb->tottet;

    float diff10[3], diff20[3], diff30[3], diff21[3], diff31[3], diff32[3];
    float delC0[3], delC1[3], delC2[3], delC3[3];

    for(int tetnr = 0; tetnr<tottet; tetnr++){
        float lambda;

        float curr_volume = get_tet_volume(bpoint_arr, &btet_arr[tetnr]);
        float ini_volume = btet_arr[tetnr].initial_volume;

        int vert0 = btet_arr[tetnr].verts[0];
        int vert1 = btet_arr[tetnr].verts[1];
        int vert2 = btet_arr[tetnr].verts[2];
        int vert3 = btet_arr[tetnr].verts[3];

        sub_v3_v3v3(diff10, bpoint_arr[vert1].x, bpoint_arr[vert0].x);
        sub_v3_v3v3(diff20, bpoint_arr[vert2].x, bpoint_arr[vert0].x);
        sub_v3_v3v3(diff30, bpoint_arr[vert3].x, bpoint_arr[vert0].x);
        sub_v3_v3v3(diff21, bpoint_arr[vert2].x, bpoint_arr[vert1].x);
        sub_v3_v3v3(diff31, bpoint_arr[vert3].x, bpoint_arr[vert1].x);
        sub_v3_v3v3(diff32, bpoint_arr[vert3].x, bpoint_arr[vert2].x);

        cross_v3_v3v3(delC0, diff31, diff21);
        cross_v3_v3v3(delC1, diff20, diff30);
        cross_v3_v3v3(delC2, diff30, diff10);
        cross_v3_v3v3(delC3, diff10, diff20);

        lambda = -1*(curr_volume - ini_volume)/(bpoint_arr[vert0].mass_inv * len_squared_v3(delC0) + 
                                                bpoint_arr[vert1].mass_inv * len_squared_v3(delC1) + 
                                                bpoint_arr[vert2].mass_inv * len_squared_v3(delC2) +
                                                bpoint_arr[vert3].mass_inv * len_squared_v3(delC3) +
                                                sb->alpha_vol/(sdt*sdt));

        mul_v3_fl(delC0, lambda*bpoint_arr[vert0].mass_inv);
        mul_v3_fl(delC1, lambda*bpoint_arr[vert1].mass_inv);
        mul_v3_fl(delC2, lambda*bpoint_arr[vert2].mass_inv);
        mul_v3_fl(delC3, lambda*bpoint_arr[vert3].mass_inv);
        
        add_v3_v3(bpoint_arr[vert0].x, delC0);
        add_v3_v3(bpoint_arr[vert1].x, delC1);
        add_v3_v3(bpoint_arr[vert2].x, delC2);
        add_v3_v3(bpoint_arr[vert3].x, delC3);
    }
}

void z_axis_collision_constraint(SoftBody *sb){
    for(int i = 0; i<sb->totpoint; i++){
        if(sb->bpoint[i].x[2] < 0.0)
            sb->bpoint[i].x[2] = 0.0; 
    }
}

static bool point_part_of_tet(int pointnr, BodyTet *tet){
    for(int i = 0; i<4; i++){
        if(tet->verts[i] == pointnr){
            return true;
        }
    }

    return false;
}

static bool same_side_of_tri(float *p1, float *p2, float *p3, float *p4, float *q){
    float normal[3];
    
    float temp1[3];
    sub_v3_v3v3(temp1, p2, p1);
    float temp2[3];
    sub_v3_v3v3(temp2, p3, p1);
    cross_v3_v3v3(normal, temp1, temp2);

    float diff_p4p1[3], diff_qp1[3];
    sub_v3_v3v3(diff_p4p1, p4, p1);
    sub_v3_v3v3(diff_qp1, q, p1);
    float dotp4 = dot_v3v3(normal, diff_p4p1);
    float dotq = dot_v3v3(normal, diff_qp1);

    return (dotq/abs(dotq)) == (dotp4/abs(dotp4));
    // Could floating point error cause problems here?
}

static bool point_inside_tet(BodyPoint *bpoint_arr, float *curr_point, BodyTet *btet){
    for(int facenr = 0; facenr < 4; facenr++){
        BodyPoint *face_point0 = &bpoint_arr[btet->verts[tet_faces[facenr][0]]];
        BodyPoint *face_point1 = &bpoint_arr[btet->verts[tet_faces[facenr][1]]];
        BodyPoint *face_point2 = &bpoint_arr[btet->verts[tet_faces[facenr][2]]];

        BodyPoint *opposite_point = &bpoint_arr[btet->verts[tet_faces[facenr][3]]];

        if(!same_side_of_tri(face_point0->x, face_point1->x, face_point2->x, opposite_point->x, curr_point)){
            return false;
        }
    }

    BodyPoint *tet_point0 = &bpoint_arr[btet->verts[0]];
    BodyPoint *tet_point1 = &bpoint_arr[btet->verts[1]];
    BodyPoint *tet_point2 = &bpoint_arr[btet->verts[2]];
    BodyPoint *tet_point3 = &bpoint_arr[btet->verts[3]];

    return true;
}

static int find_intersecting_face(BodyPoint *bpoint_arr, BodyPoint *curr_point, BodyTet *curr_tet){
    float min_dist = FLT_MAX;
    float min_face = -1;

    float translation_vec[3];
    sub_v3_v3v3(translation_vec, curr_point->x, curr_point->x_prev);
    normalize_v3(translation_vec);

    bool ahead_flag = !point_inside_tet(bpoint_arr, curr_point->x_prev, curr_tet);

    for(int facenr = 0; facenr < 4; facenr++){
        BodyPoint *face_point0 = &bpoint_arr[curr_tet->verts[tet_faces[facenr][0]]];
        BodyPoint *face_point1 = &bpoint_arr[curr_tet->verts[tet_faces[facenr][1]]];
        BodyPoint *face_point2 = &bpoint_arr[curr_tet->verts[tet_faces[facenr][2]]];
        
        float coeff[4];
        float diff10[3], diff20[3];
        sub_v3_v3v3(diff10, face_point1->x, face_point0->x);
        sub_v3_v3v3(diff20, face_point2->x, face_point0->x);
        cross_v3_v3v3(coeff, diff10, diff20);
        coeff[3] = -dot_v3v3(coeff, face_point0->x);
        
        float lambda = -(dot_v3v3(coeff, curr_point->x_prev) + coeff[3])/dot_v3v3(coeff, translation_vec);
        if(!ahead_flag){
            lambda = -lambda;
        }

        if(lambda > 0.0f && lambda<min_dist){
            min_dist = lambda;
            min_face = facenr;
        }
    }

    return min_face;

    // Checking which is the closest face
    // float min_dist = FLT_MAX;
    // float min_face = -1;
    
}

void xpbd_solve_self_collision(SoftBody *sb){
    BodyPoint *bpoint_arr = sb->bpoint;
    BodyTet *btet_arr = sb->btet;

    int totpoint = sb->totpoint;
    int tottet = sb->tottet;

    // for(int pointnr = 0; pointnr<sb->tot_surface_point; pointnr++){
    for(int pointid = 0; pointid < sb->tot_surface_point; pointid++){
        int pointnr = sb->surface_points[pointid];
        BodyPoint *curr_point = &bpoint_arr[pointnr];

        // for(int tetnr = 0; tetnr < tottet; tetnr++){
        for(int tetid = 0; tetid < sb->tot_surface_tet; tetid++){
            int tetnr = sb->surface_tets[tetid];
            BodyTet *curr_tet = &btet_arr[tetnr];
            if(point_part_of_tet(pointnr, curr_tet)){
                continue;
            }

            if(!point_inside_tet(bpoint_arr, curr_point->x, curr_tet)){
                continue;
            }

            // Check if x_prev of point is also inside tet?

            int intersecting_face = find_intersecting_face(bpoint_arr, curr_point, curr_tet);

            BodyPoint *p0 = &bpoint_arr[curr_tet->verts[tet_faces[intersecting_face][0]]];
            BodyPoint *p1 = &bpoint_arr[curr_tet->verts[tet_faces[intersecting_face][1]]];
            BodyPoint *p2 = &bpoint_arr[curr_tet->verts[tet_faces[intersecting_face][2]]];
            BodyPoint *p3 = &bpoint_arr[curr_tet->verts[tet_faces[intersecting_face][3]]];
            BodyPoint *q = curr_point;

            float diff10[3], diff20[3], diffq0[3], diff30[3];
            sub_v3_v3v3(diff10, p1->x, p0->x);
            sub_v3_v3v3(diff20, p2->x, p0->x);
            sub_v3_v3v3(diffq0, q->x, p0->x);
            

            float delC[3]; // Points outside the tet
            cross_v3_v3v3(delC, diff10, diff20);

            sub_v3_v3v3(diff30, p0->x, p3->x);
            if(dot_v3v3(diff30, delC) < 0){
                mul_v3_fl(delC, -1);
            }
            normalize_v3(delC);

            float lambda = dot_v3v3(diffq0, delC)/(p0->mass_inv + p1->mass_inv + p2->mass_inv + q->mass_inv);
            lambda /= len_squared_v3(delC);

            mul_v3_fl(delC, lambda);

            float dx_q[3];
            copy_v3_v3(dx_q, delC);
            mul_v3_fl(dx_q, -q->mass_inv);  
            add_v3_v3(q->x, dx_q);

            for(int i = 0; i<3; i++){
                BodyPoint *point = &bpoint_arr[curr_tet->verts[tet_faces[intersecting_face][i]]];

                float dx_point[3];
                copy_v3_v3(dx_point, delC);
                mul_v3_fl(dx_point, point->mass_inv);
                add_v3_v3(point->x, dx_point);
            }
        }
    }
}

void xpbd_enforce_constraints(SoftBody *sb){
    // z_axis_collision_constraint(sb);
    xpbd_enforce_distance_constraint(sb);
    xpbd_enforce_volume_constraint(sb);
    xpbd_solve_self_collision(sb);
}