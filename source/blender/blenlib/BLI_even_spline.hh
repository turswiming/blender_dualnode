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
 Abstract curve interface.

template<typename Float> class Curve {
  using Vector = vec_base<Float, 2>;

 public:
  Float length;

  Vector evaluate(Float s);
  Vector derivative(Float s);
  Vector derivative2(Float s);
  Float curvature(float s);

  void update();
};
*/

/** Quadratic curves */

/*
comment: Reduce algebra script;

on factor;
off period;

procedure bez(a, b);
  a + (b - a) * t;

lin := bez(k1, k2);
quad := bez(lin, sub(k2=k3, k1=k2, lin));

dquad := df(quad, t);
iquad := int(quad, t);

x1 := 0;
y1 := 0;

dx := sub(k1=x1, k2=x2, k3=x3, dquad);
dy := sub(k1=y1, k2=y2, k3=y3, dquad);
darc := sqrt(dx**2 + dy**2);

arcstep := darc*dt + 0.5*df(darc, t)*dt*dt;

d2x := df(dx / darc, t);
d2y := df(dy / darc, t);

gentran
begin
declare <<
x1,x2,x3,x4 : float;
y1,y2,y3,y4 : float;
dt,t : float;
>>;
return eval(arcstep)
end;

on fort;
quad;
dquad;
iquad;
arcstep;
d2x;
d2y;
off fort;

*/
template<typename Float, int axes = 2, int table_size = 512> class QuadBezier {
  using Vector = vec_base<Float, axes>;

 public:
  Vector ps[3];

  static const int TableSize = table_size;

  QuadBezier(Vector a, Vector b, Vector c)
  {
    ps[0] = a;
    ps[1] = b;
    ps[2] = c;

    deleted = false;
    _arc_to_t = new Float[table_size];
    _t_to_arc = new Float[table_size];
  }

  ~QuadBezier()
  {
    deleted = true;

    if (_arc_to_t) {
      delete[] _arc_to_t;
      _arc_to_t = nullptr;
    }

    if (_t_to_arc) {
      delete[] _t_to_arc;
      _t_to_arc = nullptr;
    }
  }

  QuadBezier()
  {
    deleted = false;
    _arc_to_t = new Float[table_size];
    _t_to_arc = new Float[table_size];
  }

  QuadBezier(const QuadBezier &b)
  {
    _arc_to_t = new Float[table_size];
    _t_to_arc = new Float[table_size];

    *this = b;
    deleted = false;
  }

  QuadBezier &operator=(const QuadBezier &b)
  {
    ps[0] = b.ps[0];
    ps[1] = b.ps[1];
    ps[2] = b.ps[2];

    length = b.length;

    if (!_arc_to_t) {
      _arc_to_t = new Float[table_size];
    }
    if (!_arc_to_t) {
      _t_to_arc = new Float[table_size];
    }

    if (b._arc_to_t) {
      for (int i = 0; i < table_size; i++) {
        _arc_to_t[i] = b._arc_to_t[i];
      }
    }

    if (b._t_to_arc) {
      for (int i = 0; i < table_size; i++) {
        _t_to_arc[i] = b._t_to_arc[i];
      }
    }

    return *this;
  }

#if 1
  QuadBezier(QuadBezier &&b)
  {
    *this = b;
  }

  Vector closest_point(const Vector p, Float &r_s, Vector &r_tan, Float &r_dis)
  {
#  if 0
    r_dis = 0.01;
    r_s = 0.0;
    r_tan = Vector(0.0, 1.0, 0.0);

    return p;
#  endif

#  if 0
    const int steps = 16;
    Float s = 0.0, ds = length / (steps - 1);
    Vector retco;

    r_dis = FLT_MAX;

    for (int i = 0; i < steps; i++, s += ds) {
      Vector co = evaluate(s);
      Vector vec = co - p;

      Float dist = _dot(vec, vec);

      if (dist < r_dis) {
        retco = co;
        r_dis = dist;
        r_s = s;
        r_tan = derivative(s);
      }
    }

    r_dis = std::sqrt(r_dis);

    return retco;
#elif 0
    float weights[3];
    float n[3];

    //normal_tri_v3(ps[0], ps[1], ps[2]);
    float p2[3];

    closest_on_tri_to_point_v3(p2, p, ps[0], ps[1], ps[2]);

    //barycentric_weights(ps[0], ps[1], ps[2], p, n, weights);
    resolve_tri_uv_v3(weights, p2, ps[0], ps[1], ps[2]);

    Float u = weights[0];
    Float v = weights[1];

    Float t1 = v / (2 * u + v);
    Float t2 = (2 * (v - 1 + u)) / (2 * u + v - 2);

    //v = min_ff(max_ff(v, 0.0f), 1.0f);

    //v = std::sqrt(v);


    t1 = min_ff(max_ff(t1, 0.0), 1.0);
    t2 = min_ff(max_ff(t2, 0.0), 1.0);

    // t1 = std::min(std::max(t1, 0.0), 1.0);
    // t2 = std::min(std::max(t2, 0.0), 1.0);

    //t1 = t2 = 0.5;

    Float s1 = t_to_arc(t1);
    Float s2 = t_to_arc(t2);

    Vector a = evaluate(s1);
    Vector b = evaluate(s2);

    a -= p;
    b -= p;

    Float dis1 = _dot(a, a);
    Float dis2 = _dot(b, b);

    if (dis1 < dis2) {
      r_s = s1;
      r_dis = std::sqrt(dis1);
      r_tan = derivative(s1);

      return a + p;
    }
    else {
      r_s = s2;
      r_dis = std::sqrt(dis2);
      r_tan = derivative(s2);

      return b + p;
    }
    //r_dis = v;

    //closest_on_tri_to_point_v3
    //resolve_tri_uv_v3
#  else
    /*
    on factor;
    off period;

    procedure bez(a, b);
      a + (b - a) * t;

    lin := bez(k1, k2);
    quad := bez(lin, sub(k2=k3, k1=k2, lin));

    x := sub(k1=x1, k2=x2, k3=x3, quad);
    y := sub(k1=y1, k2=y2, k3=y3, quad);
    z := sub(k1=z1, k2=z2, k3=z3, quad);

    comment: origin is at x1, y1, z1;

    x1 := 0;
    y1 := 0;
    z1 := 0;

    comment: rotate to align with control triangle normal;
    comment: z2 := 0;
    comment: z3 := 0;

    dx := df(x, t);
    dy := df(y, t);
    dz := df(z, t);

    tx := x - px;
    ty := y - py;
    tz := z - pz;

    x1 := 0.0;
    y1 := 0.0;
    z1 := 0.0;

    x2 := 0.5;
    y2 := 0.0;
    z2 := 0.0;

    x3 := 0.5;
    y3 := 0.5;
    z3 := 0.0;

    f := df(tx**2 + ty**2 + tz**2, t, 1);
    f := tx*dx + ty*dy + tz*dz;

    ff := solve(f, t);

    on fort;
    part(ff, 1, 2);
    part(ff, 2, 2);
    off fort;

    Computer Aided Geometric Design
    Thomas W. Sederberg
    page 203
    https://scholarsarchive.byu.edu/cgi/viewcontent.cgi?article=1000&context=facpub

    xs := {0.0, 1.0, 0.0};
    ys := {0.0, 0.0, 1.0};

    degree := 2;
    procedure l(i, j, x, y);
      binomial(degree, i)*binomial(degree, j)*(
        det(mat(
          (x, y, 1),
          (part(xs, i+1), part(ys, i+1), 1),
          (part(xs, j+1), part(ys, j+1), 1)
        ))
       );


    det mat(
      (l(2, 1, x, y), l(2, 0, x, y)),
      (l(2, 0, x, y), l(1, 0, x, y))
     );

     t1 := l(2, 0, x, y) / (l(2, 0, x, y) - l(2, 1, x, y));
     t2 := l(1, 0, x, y) / (l(1, 0, x, y) - l(2, 0, x, y));

     t1 := sub(x=part(xs, 1)*u + part(xs, 2)*v + part(xs, 3)*(1.0 - u - v), t1);
     t1 := sub(y=part(ys, 1)*u + part(ys, 2)*v + part(ys, 3)*(1.0 - u - v), t1);

     t2 := sub(x=part(xs, 1)*u + part(xs, 2)*v + part(xs, 3)*(1.0 - u - v), t2);
     t2 := sub(y=part(ys, 1)*u + part(ys, 2)*v + part(ys, 3)*(1.0 - u - v), t2);
    */

    
    Float x1 = ps[0][0], y1 = ps[0][1], z1 = ps[0][2];

    Float px = p[0] - x1;
    Float py = p[1] - y1;
    Float pz = p[2] - z1;

    Float x2 = ps[1][0] - x1, y2 = ps[1][1] - y1, z2 = ps[1][2] - z1;
    Float x3 = ps[2][0] - x1, y3 = ps[2][1] - y1, z3 = ps[2][2] - z1;

    Float x22 = x2 * x2, y22 = y2 * y2, z22 = z2 * z2;
    Float x32 = x3 * x3, y32 = y3 * y3, z32 = z3 * z3;

    Float B =
        (((2 * (4 * z22 - 2 * z2 * z3 - z32 - y32) + x32 + 4 * (2 * y2 - y3) * y2) * x2 -
          2 * (((2 * y2 - 3 * y3) * y2 + (2 * z2 - 3 * z3) * z2) * x3 - 2 * (x2 - x3) * x22)) *
             x2 -
         (x32 + y32 + (2 * z2 - z3) * (2 * z2 - z3) + 4 * (y2 - y3) * y2 + 4 * (x2 - x3) * x2) *
             ((2 * x2 - x3) * px + (2 * y2 - y3) * py) -
         (x32 + y32 + (2 * z2 - z3) * (2 * z2 - z3) + 4 * (y2 - y3) * y2 + 4 * (x2 - x3) * x2) *
             (2 * z2 - z3) * pz +
         2 * (2 * (y2 - y3) * y22 - (2 * z2 - 3 * z3) * y3 * z2) * y2 +
         (2 * (4 * z22 - 2 * z2 * z3 - z32) + y32) * y22 +
         ((2 * z2 - z3) * (2 * z2 - z3) - 2 * y32) * z22 - 2 * (y22 + z22) * x32);

    Float C = (3 * (x32 + y32 + (2 * z2 - z3) * (2 * z2 - z3) + 4 * (y2 - y3) * y2) +
               12 * (x2 - x3) * x2);
    const Float eps = 0.000001;

    if (B < 0.0 && B > -eps) {
      B = 0.0;
    }

    if (C == 0.0 || B < 0.0) {
      Vector p2 = evaluate(0.0);
      Vector vec = p2 - p;

      r_dis = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      r_tan = derivative(0.0);
      r_s = 0.0;

      return p2;
    }

    B = sqrt(B);

    Float t1 = (3 * (2 * z2 - z3) * z2 + B * sqrt(3) + 3 * (2 * y2 - y3) * y2 +
                3 * (2 * x2 - x3) * x2) /
               C;

    Float t2 = (3 * (2 * z2 - z3) * z2 - B * sqrt(3) + 3 * (2 * y2 - y3) * y2 +
                3 * (2 * x2 - x3) * x2) /
               C;

    t1 = min_ff(max_ff(t1, 0.0), 1.0);
    t2 = min_ff(max_ff(t2, 0.0), 1.0);

    // t1 = std::min(std::max(t1, 0.0), 1.0);
    // t2 = std::min(std::max(t2, 0.0), 1.0);

    Float s1 = t_to_arc(t1);
    Float s2 = t_to_arc(t2);

    Vector a = evaluate(s1);
    Vector b = evaluate(s2);

    a -= p;
    b -= p;

    Float dis1 = _dot(a, a);
    Float dis2 = _dot(b, b);

    if (dis1 < dis2) {
      r_s = s1;
      r_dis = std::sqrt(dis1);
      r_tan = derivative(s1);

      return a + p;
    }
    else {
      r_s = s2;
      r_dis = std::sqrt(dis2);
      r_tan = derivative(s2);

      return b + p;
    }
#  endif
  }

  QuadBezier &operator=(QuadBezier &&b)
  {
    ps[0] = b.ps[0];
    ps[1] = b.ps[1];
    ps[2] = b.ps[2];

    length = b.length;

    if (b._arc_to_t) {
      _arc_to_t = std::move(b._arc_to_t);
      _t_to_arc = std::move(b._t_to_arc);

      b._arc_to_t = nullptr;
      b._t_to_arc = nullptr;
    }
    else {
      _arc_to_t = new Float[table_size];
      _t_to_arc = new Float[table_size];
    }

    return *this;
  }
#endif

  Float length;

  void update()
  {
    Float t = 0.0, dt = 1.0 / (Float)table_size;
    Float s = 0.0;

    if (!_arc_to_t) {
      _arc_to_t = new Float[table_size];
      _t_to_arc = new Float[table_size];
    }

    for (int i : IndexRange(table_size)) {
      _arc_to_t[i] = -1.0;
    }

    length = 0.0;

    for (int i = 0; i < table_size; i++, t += dt) {
      Float dlen = 0.0;
      for (int j = 0; j < axes; j++) {
        float dv = dquad(ps[0][j], ps[1][j], ps[2][j], t);

        dlen += dv * dv;
      }

      dlen = sqrt(dlen) * dt;

      length += dlen;
    }

    t = 0.0;
    s = 0.0;

    for (int i = 0; i < table_size; i++, t += dt) {
      Float dlen = 0.0;
      for (int j = 0; j < axes; j++) {
        float dv = dquad(ps[0][j], ps[1][j], ps[2][j], t);

        dlen += dv * dv;
      }

      dlen = sqrt(dlen) * dt;

      int j = (int)((s / length) * (Float)table_size * 0.999999);
      j = min_ii(j, table_size - 1);

      _t_to_arc[i] = s;
      _arc_to_t[j] = t;

      s += dlen;
    }

    _arc_to_t[0] = 0.0;
    _arc_to_t[table_size - 1] = 1.0;

#if 1
    /* Interpolate gaps in table. */
    for (int i = 0; i < table_size - 1; i++) {
      if (_arc_to_t[i] == -1.0 || _arc_to_t[i + 1] != -1.0) {
        continue;
      }

      int i1 = i;
      int i2 = i + 1;

      while (i2 < table_size - 1 && _arc_to_t[i2] == -1.0) {
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

#  if 0
    for (int i = 0; i < table_size; i++) {
      printf("%.3f ", _arc_to_t[i]);
      if (_arc_to_t[i] == -1.0) {
        printf("BLI_even_spline.hh: error!\n");
      }
    }

    printf("\n\n");
#  endif
#endif
  }

  inline Vector evaluate(Float s)
  {
    Float t = arc_to_t(s);
    Vector r;

    for (int i = 0; i < axes; i++) {
      r[i] = quad(ps[0][i], ps[1][i], ps[2][i], t);
    }

    return r;
  }

  Vector derivative(Float s, bool exact = true)
  {
    Float t = arc_to_t(s);
    Vector r;

    for (int i = 0; i < axes; i++) {
      r[i] = dquad(ps[0][i], ps[1][i], ps[2][i], t) * length;
    }

    /* Real arc length parameterized tangent has unit length. */
    if (exact) {
      Float len = sqrt(_dot(r, r));

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

    Float dx = dquad(ps[0][0], ps[1][0], ps[2][0], t);
    Float d2x = d2quad(ps[0][0], ps[1][0], ps[2][0], t);
    Float dy = dquad(ps[0][1], ps[1][1], ps[2][1], t);
    Float d2y = d2quad(ps[0][1], ps[1][1], ps[2][1], t);

    /*
    comment: arc length second derivative;

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
      Float dz = dquad(ps[0][2], ps[1][2], ps[2][2], t);
      Float d2z = d2quad(ps[0][2], ps[1][2], ps[2][2], t);

      Float div = sqrt(dx * dx + dy * dy + dz * dz) * (dy * dy + dz * dz + dx * dx);

      r[0] = (d2x * dy * dy + d2x * dz * dz - d2y * dx * dy - d2z * dx * dz) / div;
      r[1] = (-(d2x * dx * dy - d2y * dx * dx - d2y * dz * dz + d2z * dy * dz)) / div;
      r[2] = (-(d2x * dx * dz + d2y * dy * dz - d2z * dx * dx - d2z * dy * dy)) / div;
    }
    else {
      for (int i = 0; i < axes; i++) {
        r[i] = d2quad(ps[0][i], ps[1][i], ps[2][i], t) * length;
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

 private:
  Float *_arc_to_t;
  Float *_t_to_arc;
  bool deleted = false;

  Float quad(Float k1, Float k2, Float k3, Float t)
  {
    return -((k1 - k2 + (k2 - k3) * t - (k1 - k2) * t) * t + (k1 - k2) * t - k1);
  }

  Float dquad(Float k1, Float k2, Float k3, Float t)
  {
    return -((k1 - k2 + (k2 - k3) * t - (k1 - k2) * t) * t + (k1 - k2) * t - k1);
  }

  Float d2quad(Float k1, Float k2, Float k3, Float t)
  {
    return -2 * (2 * k2 - k3 - k1);
  }

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

  Float arc_to_t(Float s)
  {
    if (length == 0.0) {
      return 0.0;
    }

    s = clamp_s(s);

    Float t = s * (Float)(table_size - 1) / length;

    int i1 = floorf(t);

    i1 = max_ii(i1, 0);

    int i2 = min_ii(i1 + 1, table_size - 1);

    t -= (Float)i1;

    Float s1 = _arc_to_t[i1];
    Float s2 = _arc_to_t[i2];

    return s1 + (s2 - s1) * t;
  }

  Float t_to_arc(Float s)
  {
    // return s * length; //XXX

    if (length == 0.0) {
      return 0.0;
    }

    s = clamp_s(s);

    Float t = s * (Float)(table_size - 1) / length;

    int i1 = floorf(t);
    int i2 = min_ii(i1 + 1, table_size - 1);

    t -= (Float)i1;

    Float s1 = _t_to_arc[i1];
    Float s2 = _t_to_arc[i2];

    return s1 + (s2 - s1) * t;
  }
};

/** Cubic curves */

/*
comment: Reduce algebra script;

on factor;
off period;

procedure bez(a, b);
  a + (b - a) * t;

lin := bez(k1, k2);
quad := bez(lin, sub(k2=k3, k1=k2, lin));

cubic := bez(quad, sub(k3=k4, k2=k3, k1=k2, quad));
dcubic := df(cubic, t);
icubic := int(cubic, t);

x1 := 0;
y1 := 0;

dx := sub(k1=x1, k2=x2, k3=x3, k4=x4, dcubic);
dy := sub(k1=y1, k2=y2, k3=y3, k4=y4, dcubic);
darc := sqrt(dx**2 + dy**2);

arcstep := darc*dt + 0.5*df(darc, t)*dt*dt;

d2x := df(dx / darc, t);
d2y := df(dy / darc, t);

gentran
begin
declare <<
x1,x2,x3,x4 : float;
y1,y2,y3,y4 : float;
dt,t : float;
>>;
return eval(arcstep)
end;

on fort;
cubic;
dcubic;
icubic;
arcstep;
d2x;
d2y;
off fort;

*/
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

#if 1
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
#endif

  Float length;

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

#if 1
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

#  if 0
    for (int i = 0; i < table_size; i++) {
      printf("%.3f ", _arc_to_t[i]);
    }
    printf("\n\n");
#  endif
#endif
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

    /*
    comment: arc length second derivative;

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

 private:
  Float *_arc_to_t;
  bool deleted = false;

  Float cubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -(((3.0 * (t - 1.0) * k3 - k4 * t) * t - 3.0 * (t - 1.0) * (t - 1.0) * k2) * t +
             (t - 1) * (t - 1) * (t - 1) * k1);
  }

  Float dcubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -3.0 * ((t - 1.0) * (t - 1.0) * k1 - k4 * t * t + (3.0 * t - 2.0) * k3 * t -
                   (3.0 * t - 1.0) * (t - 1.0) * k2);
  }

  Float d2cubic(Float k1, Float k2, Float k3, Float k4, Float t)
  {
    return -6.0 * (k1 * t - k1 - 3.0 * k2 * t + 2.0 * k2 + 3.0 * k3 * t - k3 - k4 * t);
  }

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

  EvenSpline<Float, axes, QuadBezier<Float, axes>> split_to_quadratics()
  {
    EvenSpline<Float, axes, QuadBezier<Float, axes>> spline;

    for (auto &seg : segments) {
      if (seg.bezier.length == 0.0) {
        continue;
      }

      const int steps = 32;
      const Float df = 0.00001;
      Float dt = seg.bezier.length / steps, t = dt;

      blender::Vector<Float> points;

      points.append(0.0);

#if 1
      /* Find local minima/maxima of second derivatives. */
      for (int i = 1; i < steps - 1; i++, t += dt) {
        Float k1 = seg.bezier.curvature(t);
        Float k2 = seg.bezier.curvature(t + df);
        Float dk1 = (k2 - k1) / df;

        k1 = seg.bezier.curvature(t + dt);
        k2 = seg.bezier.curvature(t + dt + df);
        Float dk2 = (k2 - k1) / df;

        if (dk1 <= 0.0 == dk2 <= 0.0) {
          continue;
        }

        points.append(t);
      }
#endif

      points.append(seg.bezier.length);

      blender::Vector<Float> points2;

      for (int i : IndexRange(points.size())) {
        points2.append(points[i]);

        if (i < points.size() - 1) {
          //points2.append(points[i] * 0.5 + points[i + 1] * 0.5);
        }
      }
      points = points2;

      /*
      on factor;
      off period;

      procedure bez(a, b);
        a + (b - a) * t;

      lin := bez(k1, k2);
      quad := bez(lin, sub(k2=k3, k1=k2, lin));

      dquad := df(quad, t);

      f1 := sub(t=0.0, dquad) - dv1;
      f2 := sub(t=1.0, dquad) - dv2;

      k2new := part(solve(f1, k2), 1, 2)*0.5 + part(solve(f2, k2), 1, 2)*0.5;

      */

      Float tscale = 1.0 / (points.size() - 1);

      for (int i : IndexRange(points.size() - 1)) {
        Float t1 = points[i], t2 = points[i + 1] * 0.99999;

        Vector dv1 = seg.bezier.derivative(t1) * seg.bezier.length * tscale;
        Vector dv2 = seg.bezier.derivative(t2) * seg.bezier.length * tscale;
        Vector p1 = seg.bezier.evaluate(t1);
        Vector p2 = seg.bezier.evaluate(t2);

        Vector b = (2.0 * (p1 + p2) - dv2 + dv1) / 4.0;

        // b = (p1 + p2) * 0.5;

        QuadBezier<Float, axes> quad(p1, b, p2);
        quad.update();

        spline.add(quad);
      }
    }

    spline.update();
    return spline;
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
  }

  int components() noexcept
  {
    return sizeof(segments[0].bezier.ps) / sizeof(*segments[0].bezier.ps);
  }

  inline Vector evaluate(Float s)
  {
    if (s == 0.0) {
      return segments[0].bezier.ps[0];
    }

    if (s >= length) {
      return segments[segments.size() - 1].bezier.ps[components() - 1];
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

  Vector closest_point(const Vector p, Float &r_s, Vector &r_tan, Float &r_dis)
  {
    if constexpr (std::is_same_v<BezierType, QuadBezier<Float, axes, BezierType::TableSize>>) {
      r_dis = FLT_MAX;
      Vector retco;

      for (auto &seg : segments) {
        Float s, dis;
        Vector tan;
        Vector co = seg.bezier.closest_point(p, s, tan, dis);

        if (dis < r_dis) {
          r_dis = dis;
          r_s = s + seg.start;
          r_tan = tan;
          retco = co;
        }
      }

      return retco;
    }
    else {
      return closest_point_generic(p, r_s, r_tan, r_dis);
    }
  }

  /* Find the closest point on the spline.  Uses a bisecting root finding approach.
   * Note: in thoery we could split the spline into quadratic segments and solve
   * for the closest point directy.
   */
  Vector closest_point_generic(const Vector p, Float &r_s, Vector &r_tan, Float &r_dis)
  {
    if (segments.size() == 0) {
      return Vector();
    }

    const int steps = 12;
    Float s = 0.0, ds = length / steps;
    Float mindis = FLT_MAX;
    Vector minp;
    Float mins = 0.0;
    bool found = false;

    Vector lastdv, lastp;
    Vector b, dvb;

    for (int i = 0; i < steps + 1; i++, s += ds, lastp = b, lastdv = dvb) {
      b = evaluate(s);
      dvb = derivative(s, false); /* We don't need real normalized derivative here. */

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
    s = s >= length ? length * 0.999999 : s;

    return s;
  }

  Segment *get_segment(Float s)
  {
    // printf("\n");

    for (Segment &seg : segments) {
      // printf("s: %f %f\n", seg.start, seg.start + seg.bezier.length);

      if (s >= seg.start && s < seg.start + seg.bezier.length) {
        return &seg;
      }
    }

    // printf("\n");

    return nullptr;
  }
};

using BezierSpline2f = EvenSpline<float, 2>;
using BezierSpline3f = EvenSpline<float, 3>;
}  // namespace blender
