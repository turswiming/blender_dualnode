/* SPDX-License-Identifier: GPL-2.0-or-later */

#pragma once

#include "BKE_bvhutils.h"
#include "BKE_curves.hh"

/** \file
 * \ingroup bke
 * \brief Constraint solver for curve deformations.
 */

namespace blender::bke::curves {

class ConstraintSolver {
 public:
  enum ErrorType {
    Ok,
    NotConverged,
  };

  struct Params {
    /* Keep the distance between points constant. */
    bool use_length_constraints = true;
    /* Root point is fixed to the surface. */
    bool use_root_constraints = true;
    /* Points do not penetrate the surface. */
    bool use_collision_constraints = true;

    /* Number of substeps to perform.
     * More substeps can be faster overall because of reduced search radius for collisions. */
    int substep_count = 20;

    /* Maximum overall distance a particle can move during a step.
     * Divide by substep count to get max substep travel.
     * This determines a the search radius for collisions.
     * A larger travel distance means the point can move faster,
     * but it can take longer to find collisions. */
    float max_travel_distance = 0.1f;

    /* Maximum number of simultaneous contacts to record per point. */
    int max_contacts_per_point = 4;

    /* Number of iterations to satisfy constraints. */
    int max_solver_iterations = 5;

    /* Acceptable error threshold for convergence, in length units. */
    float error_threshold = 1.0e-4f;
  };

  struct Result {
    ErrorType error = ErrorType::Ok;

    /* True if any point's original travel was larger than the allowed maximum.
     * If clamping is enabled the point's travel will be shorter than the input. */
    bool max_travel_exceeded = false;

    /* Root-mean-square of the constraint residuals, indicating numerical quality
     * of the solution. For a positional solver this value is in length units
     * and describes how far constraints are violated on average. */
    double rms_residual = 0.0;

    /* Sum of squared errors (in length units). */
    double error_squared_sum = 0.0;
    /* Largest squared error. */
    double max_error_squared = 0.0;
    /* Total number of constraints solved. */
    int constraint_count = 0;

    struct {
      /* Time in seconds for the entire step. */
      double step_total = 0.0;
      /* Time in seconds to build the BVH tree of the surface.
       * This is zero if the surface BVH was cached. */
      double build_bvh = 0.0;
      /* Time in seconds to find contact points, cumulative over substeps. */
      double find_contacts = 0.0;
      /* Time in seconds to solve constraints, cumulative over substeps. */
      double solve_constraints = 0.0;
    } timing;
  };

 private:
  Params params_;

  /** Length of each segment indexed by the index of the first point in the segment. */
  Array<float> segment_lengths_cu_;

  struct Contact {
    float dist_sq_;
    float3 normal_;
    float3 point_;
  };

  Array<int> contacts_num_;
  Array<Contact> contacts_;

  /** Information about the most recent step solution. */
  mutable Result result_;

 public:
  const Params &params() const;

  Span<float> segment_lengths() const;

  const Result &result() const;
  void clear_result();

  /* Initialize the solver for a given set of curves.
   * The solver must be reinitialized if the curve set changes. */
  void initialize(const Params &params, const CurvesGeometry &curves);

  /* Solve constraints for an independent subset of curves. */
  void step_curves(CurvesGeometry &curves,
                   const Mesh *surface,
                   const CurvesSurfaceTransforms &transforms,
                   Span<float3> start_positions,
                   VArray<int> changed_curves);

 private:
  void find_contact_points(const CurvesGeometry &curves,
                           const Mesh *surface,
                           const CurvesSurfaceTransforms &transforms,
                           float max_dist,
                           VArray<int> changed_curves);

  /**
   * Satisfy constraints on curve points based on initial deformation.
   */
  void solve_constraints(CurvesGeometry &curves, VArray<int> changed_curves) const;

  /**
   * Satisfy constraints on curve points based on initial deformation.
   */
  inline void solve_constraints(CurvesGeometry &curves, IndexRange changed_curves) const
  {
    solve_constraints(curves,
                      VArray<int>::ForFunc(changed_curves.size(),
                                           [changed_curves](int64_t i) { return (int)i; }));
  }
};

}  // namespace blender::bke::curves
