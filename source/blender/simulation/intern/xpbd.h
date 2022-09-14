// #include "BKE_softbody.h"
// #include "DNA_object_force_types.h"

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

struct BodyPoint;
struct BodyTet;
struct SoftBody;

void xpbd_position_update(struct SoftBody *sb);
void xpbd_distance_constraint_for_edge(struct BodyPoint *p1, struct BodyPoint *p2, float compliance);
void xpbd_enforce_distance_constraint(struct SoftBody *sb);
void xpbd_enforce_volume_constraint(struct SoftBody *sb);
void z_axis_collision_constraint(struct SoftBody *sb);
void xpbd_solve_self_collision(struct SoftBody *sb);
void xpbd_enforce_constraints(struct SoftBody *sb);
void xpbd_velocity_update(struct SoftBody *sb);

#ifdef __cplusplus
}
#endif