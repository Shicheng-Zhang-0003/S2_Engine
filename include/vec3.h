#ifndef VEC3_H
#define VEC3_H

#include <math.h>
#include <stdio.h>

/*
 * vec3.h
 * Inline 3D vector arithmetic in double precision.
 * All operations are branchless and compile to efficient SIMD on modern GCC/Clang.
 */

typedef struct { double x, y, z; } Vec3;

/* ── Construction ────────────────────────────────────────────────────────── */
static inline Vec3 vec3(double x, double y, double z) {
    return (Vec3){x, y, z};
}
static inline Vec3 vec3_zero(void) { return (Vec3){0.0, 0.0, 0.0}; }
static inline Vec3 vec3_ones(void) { return (Vec3){1.0, 1.0, 1.0}; }

/* ── Arithmetic ──────────────────────────────────────────────────────────── */
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z};
}
static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z};
}
static inline Vec3 vec3_scale(Vec3 a, double s) {
    return (Vec3){a.x*s, a.y*s, a.z*s};
}
static inline Vec3 vec3_negate(Vec3 a) {
    return (Vec3){-a.x, -a.y, -a.z};
}
static inline Vec3 vec3_mul(Vec3 a, Vec3 b) {   /* element-wise product */
    return (Vec3){a.x*b.x, a.y*b.y, a.z*b.z};
}

/* ── Dot / cross / norm ──────────────────────────────────────────────────── */
static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline double vec3_norm2(Vec3 a) {        /* squared length           */
    return a.x*a.x + a.y*a.y + a.z*a.z;
}
static inline double vec3_norm(Vec3 a) {
    return sqrt(vec3_norm2(a));
}
static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

/* ── Normalise (returns zero vector if near-zero input) ──────────────────── */
static inline Vec3 vec3_normalize(Vec3 a) {
    double n = vec3_norm(a);
    if (n < 1.0e-300) return vec3_zero();
    return vec3_scale(a, 1.0/n);
}

/* ── Distance between two points ─────────────────────────────────────────── */
static inline double vec3_dist(Vec3 a, Vec3 b) {
    return vec3_norm(vec3_sub(b, a));
}
static inline double vec3_dist2(Vec3 a, Vec3 b) {
    Vec3 d = vec3_sub(b, a);
    return vec3_norm2(d);
}

/* ── Minimum-image convention for periodic boundaries ───────────────────── */
/* Legacy all-axes wrapper: only correct when all three axes share the
 * same periodicity. Prefer vec3_pbc_box() below for mixed PBC.
 * A degenerate (zero or negative) box length is passed through
 * UNWRAPPED rather than producing NaN: 0 * round(dr/0) = 0 * inf = NaN,
 * which used to silently poison every coordinate downstream. */
static inline Vec3 vec3_pbc(Vec3 dr, Vec3 box) {
    Vec3 r;
    r.x = (box.x > 1e-12) ? dr.x - box.x * round(dr.x / box.x) : dr.x;
    r.y = (box.y > 1e-12) ? dr.y - box.y * round(dr.y / box.y) : dr.y;
    r.z = (box.z > 1e-12) ? dr.z - box.z * round(dr.z / box.z) : dr.z;
    return r;
}

/* Per-axis minimum image: wraps only axes flagged periodic[3].
 * Non-periodic axes pass through unwrapped. Guards degenerate box. */
static inline Vec3 vec3_pbc_box(Vec3 dr, Vec3 box, const int periodic[3]) {
    Vec3 r = dr;
    if (periodic) {
        if (periodic[0] && box.x > 1e-12) r.x = dr.x - box.x * round(dr.x / box.x);
        if (periodic[1] && box.y > 1e-12) r.y = dr.y - box.y * round(dr.y / box.y);
        if (periodic[2] && box.z > 1e-12) r.z = dr.z - box.z * round(dr.z / box.z);
    }
    return r;
}

/* ── Angle between two vectors (radians) ─────────────────────────────────── */
/* Full-audit P8: zero-vector guard. na=|a|, nb=|b|; if either <1e-300 the
 * angle is undefined (0/0=NaN survives the clamp and acos(NaN)=NaN, which
 * sim_rebuild_angles_geometric stored as theta0). Return 0.0 for degenerate
 * input rather than NaN. */
static inline double vec3_angle(Vec3 a, Vec3 b) {
    double na = vec3_norm(a), nb = vec3_norm(b);
    if (!(na > 1e-300) || !(nb > 1e-300) || !isfinite(na) || !isfinite(nb)) return 0.0;
    double c = vec3_dot(a, b) / (na * nb);
    /* clamp to [-1,1] for numerical safety */
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    if (!isfinite(c)) return 0.0;
    return acos(c);
}

/* ── Dihedral angle (radians), properly SIGNED ───────────────────────────── */
/*
 * Computes the signed dihedral angle for 4 sequential atoms p1-p2-p3-p4,
 * given as 3 bond vectors b1=p2-p1, b2=p3-p2, b3=p4-p3.
 *
 * Standard atan2-based formula (used throughout computational chemistry
 * precisely because it is numerically robust and, critically, SIGNED -
 * unlike an acos-based unsigned angle, which cannot distinguish +60 deg
 * from -60 deg. That distinction is not a refinement: the sign of a
 * dihedral IS the chirality/handedness of the twist it describes, and
 * any real torsion potential (which must, for example, favor a right-
 * handed helix over a left-handed one) is meaningless without it.
 *
 * Verified against 4 independent, hand-constructed test cases with
 * known angles (+90, -90, 0, 180 degrees) before use in any force
 * calculation - see conversation record for the numerical verification
 * performed in Python before this C implementation was written.
 *
 * n1 = b1 x b2  (normal to the plane containing p1,p2,p3)
 * n2 = b2 x b3  (normal to the plane containing p2,p3,p4)
 * m1 = n1 x b2_hat  (a vector in the n1 plane, perpendicular to b2)
 * phi = atan2(m1 . n2, n1 . n2)
 */
/* FULL-AUDIT C32: degenerate quads (zero/collinear legs) make phi
 * undefined; direct callers outside forces.c got atan2(NaN). Guard here. */
static inline double vec3_dihedral(Vec3 b1, Vec3 b2, Vec3 b3) {
    if (!isfinite(b1.x + b1.y + b1.z + b2.x + b2.y + b2.z + b3.x + b3.y + b3.z)) return 0.0;
    if (vec3_norm(b1) < 1.0e-10 || vec3_norm(b2) < 1.0e-10 || vec3_norm(b3) < 1.0e-10) return 0.0;
    if (vec3_norm(vec3_cross(b1, b2)) < 1.0e-12 || vec3_norm(vec3_cross(b2, b3)) < 1.0e-12) return 0.0;
    Vec3 n1 = vec3_cross(b1, b2);
    Vec3 n2 = vec3_cross(b2, b3);
    Vec3 b2_hat = vec3_normalize(b2);
    Vec3 m1 = vec3_cross(n1, b2_hat);

    double x = vec3_dot(n1, n2);
    double y = vec3_dot(m1, n2);
    return atan2(y, x);
}

/* ── Accumulate (v += a) ─────────────────────────────────────────────────── */
static inline void vec3_iadd(Vec3 *v, Vec3 a) {
    v->x += a.x; v->y += a.y; v->z += a.z;
}
static inline void vec3_isub(Vec3 *v, Vec3 a) {
    v->x -= a.x; v->y -= a.y; v->z -= a.z;
}
static inline void vec3_iscale(Vec3 *v, double s) {
    v->x *= s; v->y *= s; v->z *= s;
}

/* ── Printing ─────────────────────────────────────────────────────────────── */
static inline void vec3_print(const char *label, Vec3 v) {
    printf("%s: (%.6f, %.6f, %.6f)\n", label, v.x, v.y, v.z);
}

/* ── Rodrigues' rotation formula ─────────────────────────────────────────── */
/*
 * Rotates vector v around unit axis `axis` by `angle` radians (right-hand
 * rule). Standard formula: v_rot = v*cos(t) + (axis x v)*sin(t)
 *                                   + axis*(axis.v)*(1-cos(t))
 * `axis` need not be pre-normalised; this function normalises it.
 */
static inline Vec3 vec3_rotate_axis_angle(Vec3 v, Vec3 axis, double angle) {
    /* FULL-AUDIT C32: degenerate axis previously returned v*cos(t) (silent
     * wrong rotation). Return v unchanged: the honest identity. */
    if (!isfinite(axis.x + axis.y + axis.z) || vec3_norm(axis) < 1.0e-12) return v;
    axis = vec3_normalize(axis);
    double c = cos(angle), s = sin(angle);
    Vec3 term1 = vec3_scale(v, c);
    Vec3 term2 = vec3_scale(vec3_cross(axis, v), s);
    Vec3 term3 = vec3_scale(axis, vec3_dot(axis, v) * (1.0 - c));
    return vec3_add(vec3_add(term1, term2), term3);
}

#endif /* VEC3_H */
