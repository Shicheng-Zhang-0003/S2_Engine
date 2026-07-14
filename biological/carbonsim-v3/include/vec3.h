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
static inline Vec3 vec3_pbc(Vec3 dr, Vec3 box) {
    Vec3 r;
    r.x = dr.x - box.x * round(dr.x / box.x);
    r.y = dr.y - box.y * round(dr.y / box.y);
    r.z = dr.z - box.z * round(dr.z / box.z);
    return r;
}

/* ── Angle between two vectors (radians) ─────────────────────────────────── */
static inline double vec3_angle(Vec3 a, Vec3 b) {
    double c = vec3_dot(a, b) / (vec3_norm(a) * vec3_norm(b));
    /* clamp to [-1,1] for numerical safety */
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    return acos(c);
}

/* ── Dihedral angle (radians) ────────────────────────────────────────────── */
static inline double vec3_dihedral(Vec3 b1, Vec3 b2, Vec3 b3) {
    Vec3 n1 = vec3_cross(b1, b2);
    Vec3 n2 = vec3_cross(b2, b3);
    double cos_phi = vec3_dot(n1, n2) / (vec3_norm(n1) * vec3_norm(n2));
    if (cos_phi >  1.0) cos_phi =  1.0;
    if (cos_phi < -1.0) cos_phi = -1.0;
    return acos(cos_phi);
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

#endif /* VEC3_H */
