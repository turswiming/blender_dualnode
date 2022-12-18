#pragma once

#include "BLI_compiler_attrs.h"
#include "BLI_compiler_compat.h"

#include "BLI_math.h"
#include "BLI_math_vec_types.hh"
#include "BLI_vector.hh"

#include <algorithm>
#include <cstdio>
#include <type_traits>
#include <utility>

//#define FINITE_DIFF

/*
 * Arc length parameterized spline library.
 */
namespace blender {
/*
 Abstract curve interface. Curves must be
 arc length parameterized (evenly spaced,
 length of first derivative vector is 1).

 We get a nice set of properties by using
 arc length parameterizion.

 * Even spacing.
 * Second derivative is the curve normal multiplied by the curvature
 * Easy to calculate torsion for space curves (though this hasn't been
   implemented yet).

 Look up the Frenet-Serret formulas for more info.

 Abstract class interface:

    template<typename Float> class Curve {
      using Vector = vec_base<Float, 2>;

     public:
      Float length;

      Vector evaluate(Float s);
      Vector derivative(Float s);
      Vector derivative2(Float s);

      Float curvature(float s);
      Float dcurvature(float s); // Derivative of curvature.

      void update();
    };

*/

/** Cubic curves */

template<typename Float, int axes = 2, int table_size = 512> class CubicBezier {
  using Vector = vec_base<Float, axes>;

 public:
  Vector ps[4];

  static const int TableSize = table_size;

  CubicBezier(Vector a, Vector b, Vector c, Vector d)
  {
    ps[0] = a;
    ps[1] = b;
    ps[2] = c;
    ps[3] = d;

    deleted = false;
    _arc_to_t = new Float[table_size];
  }

  ~CubicBezier()
  {
    deleted = true;

    if (_arc_to_t) {
      delete[] _arc_to_t;
      _arc_to_t = nullptr;
    }
  }

  CubicBezier()
  {
    deleted = false;
    _arc_to_t = new Float[table_size];
  }

  CubicBezier(const CubicBezier &b)
  {
    _arc_to_t = new Float[table_size];
    *this = b;
    deleted = false;
  }

  CubicBezier &operator=(const CubicBezier &b)
  {
    ps[0] = b.ps[0];
    ps[1] = b.ps[1];
    ps[2] = b.ps[2];
    ps[3] = b.ps[3];

    length = b.length;

    if (!_arc_to_t) {
      _arc_to_t = new Float[table_size];
    }

    if (b._arc_to_t) {
      for (int i = 0; i < table_size; i++) {
        _arc_to_t[i] = b._arc_to_t[i];
      }
    }

    return *this;
  }

  CubicBezier(CubicBezier &&b)
  {
    *this = b;
  }

  CubicBezier &operator=(CubicBezier &&b)
  {
    ps[0] = b.ps[0];
    ps[1] = b.ps[1];
    ps[2] = b.ps[2];
    ps[3] = b.ps[3];

    length = b.length;

    if (b._arc_to_t) {
      _arc_to_t = std::move(b._arc_to_t);
      b._arc_to_t = nullptr;
    }
    else {
      _arc_to_t = new Float[table_size];
    }

    return *this;
  }

  Float length;

  /* Update arc length -> parameterization table. */
  void update()
  {
    Float t = 0.0, dt = 1.0 / (Float)table_size;
    Float s = 0.0;

    if (!_arc_to_t) {
      _arc_to_t = new Float[table_size];
    }

    for (int i = 0; i < table_size; i++) {
      _arc_to_t[i] = -1.0;
    }

    length = 0.0;

    for (int i = 0; i < table_size; i++, t += dt) {
      Float dlen = 0.0;
      for (int j = 0; j < axes; j++) {
        float dv = dcubic(ps[0][j], ps[1][j], ps[2][j], ps[3][j], t);

        dlen += dv * dv;
      }

      dlen = sqrt(dlen) * dt;

      length += dlen;
    }

    const int samples = table_size;
    dt = 1.0 / (Float)samples;

    t = 0.0;
    s = 0.0;

    for (int i = 0; i < samples; i++, t += dt) {
      Float dlen = 0.0;
      for (int j = 0; j < axes; j++) {
        float dv = dcubic(ps[0][j], ps[1][j], ps[2][j], ps[3][j], t);

        dlen += dv * dv;
      }

      dlen = sqrt(dlen) * dt;

      int j = (int)((s / length) * (Float)table_size * 0.999999);
      j = min_ii(j, table_size - 1);

      _arc_to_t[j] = t;

      s += dlen;
    }

    _arc_to_t[0] = 0.0;
    _arc_to_t[table_size - 1] = 1.0;

    /* Interpolate gaps in table. */
    for (int i = 0; i < table_size - 1; i++) {
      if (_arc_to_t[i] == -1.0 || _arc_to_t[i + 1] != -1.0) {
        continue;
      }

      int i1 = i;
      int i2 = i + 1;

      while (_arc_to_t[i2] == -1.0) {
        i2++;
      }

      Float start = _arc_to_t[i1];
      Float end = _arc_to_t[i2];
      Float dt2 = 1.0 / (i2 - i1);

      for (int j = i1 + 1; j < i2; j++) {
        Float factor = (Float)(j - i1) * dt2;
        _arc_to_t[j] = start + (end - start) * factor;
      }

      i = i2 - 1;
    }
  }

  inline Vector evaluate(Float s)
  {
    Float t = arc_to_t(s);
    Vector r;

    for (int i = 0; i < axes; i++) {
      r[i] = cubic(ps[0][i], ps[1][i], ps[2][i], ps[3][i], t);
    }

    return r;
  }

  Vector derivative(Float s, bool exact = true)
  {
    Float t = arc_to_t(s);
    Vector r;

    for (int i = 0; i < axes; i++) {
      r[i] = dcubic(ps[0][i], ps[1][i], ps[2][i], ps[3][i], t) * length;
    }

    /* Real arc length parameterized tangent has unit length. */
    if (exact) {
      Float len = sqrt(_dot(r, r));

      /* Use FLT_EPSILON here? */
      if (len > 0.00001) {
        r = r / len;
      }
    }

    return r;
  }

  Vector derivative2(Float s)
  {
#ifdef FINITE_DIFF
    const Float df = 0.0005;
    Float s1, s2;

    if (s >= 1.0 - df) {
      s1 = s - df;
      s2 = s;
    }
    else {
      s1 = s;
      s2 = s + df;
    }

    Vector a = derivative(s1);
    Vector b = derivative(s2);

    return (b - a) / df;
#else
    Float t = arc_to_t(s);
    Vector r;

    Float dx = dcubic(ps[0][0], ps[1][0], ps[2][0], ps[3][0], t);
    Float d2x = d2cubic(ps[0][0], ps[1][0], ps[2][0], ps[3][0], t);
    Float dy = dcubic(ps[0][1], ps[1][1], ps[2][1], ps[3][1], t);
    Float d2y = d2cubic(ps[0][1], ps[1][1], ps[2][1], ps[3][1], t);

    /* reduce algebra script
    comment: arc length second derivative;

    comment: build arc length from abstract derivative operators;

    operator x, y, z, dx, dy, dz, d2x, d2y, d2z;

    forall t let df(x(t), t) = dx(t);
    forall t let df(y(t), t) = dy(t);
    forall t let df(z(t), t) = dz(t);
    forall t let df(dx(t), t) = d2x(t);
    forall t let df(dy(t), t) = d2y(t);
    forall t let df(dz(t), t) = d2z(t);

    comment: arc length first derivative is just the normalized tangent;

    comment: 2d case;
    dlen := sqrt(df(x(t), t)**2 + df(y(t), t)**2);

    comment: final derivatives;
    df(df(x(t), t) / dlen, t);
    df(df(y(t), t) / dlen, t);

    comment: 3d case;
    dlen := sqrt(df(x(t), t)**2 + df(y(t), t)**2 + df(z(t), t)**2);

    comment: final derivatives;
    df(df(x(t), t) / dlen, t);
    df(df(y(t), t) / dlen, t);
    df(df(z(t), t) / dlen, t);
    */

    if constexpr (axes == 2) {
      /* Basically the 2d perpidicular normalized tangent multiplied by the curvature. */
      Float div = sqrt(dx * dx + dy * dy) * (dx * dx + dy * dy);

      r[0] = ((d2x * dy - d2y * dx) * dy) / div;
      r[1] = (-(d2x * dy - d2y * dx) * dx) / div;
    }
    else if constexpr (axes == 3) {
      Float dz = dcubic(ps[0][2], ps[1][2], ps[2][2], ps[3][2], t);
      Float d2z = d2cubic(ps[0][2], ps[1][2], ps[2][2], ps[3][2], t);

      Float div = sqrt(dx * dx + dy * dy + dz * dz) * (dy * dy + dz * dz + dx * dx);

      r[0] = (d2x * dy * dy + d2x * dz * dz - d2y * dx * dy - d2z * dx * dz) / div;
      r[1] = (-(d2x * dx * dy - d2y * dx * dx - d2y * dz * dz + d2z * dy * dz)) / div;
      r[2] = (-(d2x * dx * dz + d2y * dy * dz - d2z * dx * dx - d2z * dy * dy)) / div;
    }
    else {
      for (int i = 0; i < axes; i++) {
        r[i] = d2cubic(ps[0][i], ps[1][i], ps[2][i], ps[3][i], t) * length;
      }
    }

    return r;
#endif
  }

  Float curvature(Float s)
  {
    Vector dv2 = derivative2(s);

    if constexpr (axes == 2) {
      Vector dv = derivative(s, true);

      /* Calculate signed curvature. Remember that dv is normalized. */
      return dv[0] * dv2[1] - dv[1] * dv2[0];
    }

    return sqrt(_dot(dv2, dv2));
  }

  /* First derivative of curvature. */
  Float dcurvature(Float s)
  {
    const Float ds = 0.0001;
    Float s1, s2;

    if (s > 1.0 - ds) {
      s1 = s - ds;
      s2 = s;
    }
    else {
      s1 = s;
      s2 = s + ds;
    }

    Float a = curvature(s1);
    Float b = curvature(s2);

    return (b - a) / ds;
  }

 private:
  Float *_arc_to_t;
  bool deleted = false;

  /* Bernstein/bezier polynomial.*/
  Float cubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -(((3.0 * (t - 1.0) * k3 - k4 * t) * t - 3.0 * (t - 1.0) * (t - 1.0) * k2) * t +
             (t - 1) * (t - 1) * (t - 1) * k1);
  }

  /* First derivative. */
  Float dcubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -3.0 * ((t - 1.0) * (t - 1.0) * k1 - k4 * t * t + (3.0 * t - 2.0) * k3 * t -
                   (3.0 * t - 1.0) * (t - 1.0) * k2);
  }

  /* Second derivative. */
  Float d2cubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -6.0 * (k1 * t - k1 - 3.0 * k2 * t + 2.0 * k2 + 3.0 * k3 * t - k3 - k4 * t);
  }

  /* Inlinable dot product. */
  Float _dot(Vector a, Vector b)
  {
    Float sum = 0.0;

    for (int i = 0; i < axes; i++) {
      sum += a[i] * b[i];
    }

    return sum;
  }

  Float clamp_s(Float s)
  {
    s = s < 0.0 ? 0.0 : s;
    s = s >= length ? length * 0.999999 : s;

    return s;
  }

  /* Convert a unit distance along the curve to parameterization t
   * using a linearly-interpolated lookup table.
   */
  Float arc_to_t(Float s)
  {
    if (length == 0.0) {
      return 0.0;
    }

    s = clamp_s(s);

    Float t = s * (Float)(table_size - 1) / length;

    int i1 = floorf(t);
    int i2 = min_ii(i1 + 1, table_size - 1);

    t -= (Float)i1;

    Float s1 = _arc_to_t[i1];
    Float s2 = _arc_to_t[i2];

    return s1 + (s2 - s1) * t;
  }
};

template<typename Float, int axes = 2, typename BezierType = CubicBezier<Float, axes>>
class EvenSpline {
  using Vector = vec_base<Float, axes>;
  struct Segment {
    BezierType bezier;
    Float start = 0.0;

    Segment(const BezierType &bez)
    {
      bezier = bez;
    }

    Segment(const Segment &b)
    {
      *this = b;
    }

    Segment &operator=(const Segment &b)
    {
      bezier = b.bezier;
      start = b.start;

      return *this;
    }

    Segment()
    {
    }
  };

 public:
  Float length = 0.0;
  bool deleted = false;
  blender::Vector<Segment> segments;
  blender::Vector<Float> inflection_points;

  void clear()
  {
    segments.clear();
  }

  EvenSpline()
  {
  }

  ~EvenSpline()
  {
    deleted = true;
  }

  void add(BezierType &bez)
  {
    need_update = true;

    Segment seg;

    seg.bezier = bez;
    segments.append(seg);

    update();
  }

  void update()
  {
    need_update = false;

    length = 0.0;
    for (Segment &seg : segments) {
      seg.start = length;
      length += seg.bezier.length;
    }

    update_inflection_points();
  }

  /* Find inflection points, these are used to speed
   * up closest point test.
   */
  void update_inflection_points()
  {
    inflection_points.clear();

    inflection_points.append(0.0);

    const int steps = segments.size() * 5;
    Float s = 0.0, ds = length / (steps - 1);
    Float dk, lastdk;

    for (int i = 0; i < steps; i++, s += ds, lastdk = dk) {
      dk = dcurvature(s);

      if (i == 0) {
        continue;
      }

      if (dk < 0.0 != lastdk < 0.0) {
        inflection_points.append(s - ds * 0.5);
      }
    }

    inflection_points.append(1.0);
  }

  /* Number of control points inside a curve segment. */
  int order() noexcept
  {
    return sizeof(segments[0].bezier.ps) / sizeof(*segments[0].bezier.ps);
  }

  inline Vector evaluate(Float s)
  {
    if (s == 0.0) {
      return segments[0].bezier.ps[0];
    }

    if (s >= length) {
      return segments[segments.size() - 1].bezier.ps[order() - 1];
    }

    Segment *seg = get_segment(s);

    return seg->bezier.evaluate(s - seg->start);
  }

  Vector derivative(Float s, bool exact = true)
  {
    if (segments.size() == 0) {
      return Vector();
    }

    s = clamp_s(s);
    Segment *seg = get_segment(s);

    return seg->bezier.derivative(s - seg->start, exact);
  }

  /* Second derivative. */
  Vector derivative2(Float s)
  {
    if (segments.size() == 0) {
      return Vector();
    }

    s = clamp_s(s);
    Segment *seg = get_segment(s);

    return seg->bezier.derivative2(s - seg->start);
  }

  Float curvature(Float s)
  {
    if (segments.size() == 0) {
      return 0.0;
    }

    s = clamp_s(s);
    Segment *seg = get_segment(s);

    return seg->bezier.curvature(s - seg->start);
  }

  /* First derivative of curvature. */
  Float dcurvature(Float s)
  {
    if (segments.size() == 0) {
      return 0.0;
    }

    s = clamp_s(s);
    Segment *seg = get_segment(s);

    return seg->bezier.dcurvature(s - seg->start);
  }

  /* Find the closest point on the spline.  Uses a bisecting root
   * finding approach.
   */
  Vector closest_point(const Vector p, Float &r_s, Vector &r_tan, Float &r_dis)
  {
    if (segments.size() == 0) {
      return Vector();
    }

    Float mindis = FLT_MAX;
    Vector minp;
    Float mins = 0.0;
    bool found = false;

    Vector lastdv, lastp;
    Vector b, dvb;

    for (int i = 0; i < inflection_points.size(); i++, lastp = b, lastdv = dvb) {
      Float s = inflection_points[i];
      Float ds = s - inflection_points[i - 1];

      b = evaluate(s);

      /* The extra false parameter signals we can
       * accept a non-normalized first derivative.
       */
      dvb = derivative(s, false);

      if (i == 0) {
        continue;
      }

      Vector dva = lastdv;
      Vector a = lastp;

      Vector vec1 = a - p;
      Vector vec2 = b - p;

      Float sign1 = _dot(vec1, dva);
      Float sign2 = _dot(vec2, dvb);

      if ((sign1 < 0.0) == (sign2 < 0.0)) {
        found = true;

        Float len = _dot(vec1, vec1);

        if (len < mindis) {
          mindis = len;
          mins = s;
          minp = evaluate(s);
        }
        continue;
      }

      found = true;

      Float start = s - ds;
      Float end = s;
      Float mid = (start + end) * 0.5;
      const int binary_steps = 10;

      /* Main binary search loop. */
      for (int j = 0; j < binary_steps; j++) {
        Vector dvmid = derivative(mid, false);
        Vector vecmid = evaluate(mid) - p;
        Float sign_mid = _dot(vecmid, dvmid);

        if ((sign_mid < 0.0) == (sign1 < 0.0)) {
          start = mid;
        }
        else {
          end = mid;
        }
        mid = (start + end) * 0.5;
      }

      Vector p2 = evaluate(mid);
      Vector vec_mid = p2 - p;
      Float len = _dot(vec_mid, vec_mid);

      if (len < mindis) {
        mindis = len;
        minp = p2;
        mins = mid;
      }
    }

    if (!found) {
      mins = 0.0;
      minp = evaluate(mins);
      Vector vec = minp - p;
      mindis = _dot(vec, vec);
    }

    r_tan = derivative(mins, true);
    r_s = mins;
    r_dis = sqrtf(mindis);

    return minp;
  }

  void pop_front(int n = 1)
  {
    for (int i = 0; i < segments.size() - n; i++) {
      segments[i] = segments[i + n];
    }

    segments.resize(segments.size() - n);
    update();
  }

 private:
  bool need_update;

  Float _dot(Vector a, Vector b)
  {
    Float sum = 0.0;

    for (int i = 0; i < axes; i++) {
      sum += a[i] * b[i];
    }

    return sum;
  }

  Float clamp_s(Float s)
  {
    s = s < 0.0 ? 0.0 : s;
    s = s >= length ? length * (1.0 - FLT_EPSILON) : s;

    return s;
  }

  Segment *get_segment(Float s)
  {
    for (Segment &seg : segments) {
      if (s >= seg.start && s < seg.start + seg.bezier.length) {
        return &seg;
      }
    }

    return nullptr;
  }
};

using BezierSpline2f = EvenSpline<float, 2>;
using BezierSpline3f = EvenSpline<float, 3>;
}  // namespace blender
