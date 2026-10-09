/* test_eht.c — EHT oracle suite (7th gate).
 *
 * Every closed form here is checked against an INDEPENDENT implementation:
 * brute-force 3D Cartesian grid integration of separately-coded STOs
 * (not the engine's radial code, not the prolate machinery). A wrong
 * formula misses by O(0.1); the grid resolves ~1e-4. Jacobi is checked
 * against 2x2 closed forms. SCF properties (signs, neutrality, BO
 * dissolution, SCC convergence, rotation invariance) are law-level.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/qm_eht.h"
#include "../include/sim.h"
#include "../include/vec3.h"

static int npass = 0, nfail = 0;
static void ok (const char *name, int cond, const char *detail) {
    if (cond) {
        npass++;
        printf ("  PASS %s\n", name);
    } else {
        nfail++;
        printf ("  FAIL %s %s\n", name, detail ? detail : "");
    }
}
static void okrel (const char *name, double got, double ref, double tol) {
    char d[160];
    snprintf (d, sizeof d, "(got %.9g ref %.9g tol %.1e)", got, ref, tol);
    ok (name, fabs (got - ref) <= tol * (1.0 + fabs (ref)), d);
}

/* Correct Cartesian STOs, derived from chi = N_n r^{n-1} e^{-zr} Y_lm:
 * N_1 = (2z)^{3/2}/sqrt(2), N_2 = (2z)^{5/2}/sqrt(24). */
static double n1 (double z) {
    double t = 2.0 * z;
    return t * sqrt (t) / sqrt (2.0);
}
static double n2 (double z) {
    double t = 2.0 * z;
    return t * t * sqrt (t) / sqrt (24.0);
}
static double o1s (double z, double x, double y, double w) {
    double r = sqrt (x * x + y * y + w * w);
    return n1 (z) * exp (-z * r) * 0.28209479177387814;
}
static double o2s (double z, double x, double y, double w) {
    double r = sqrt (x * x + y * y + w * w);
    return n2 (z) * r * exp (-z * r) * 0.28209479177387814;
}
static double o2pz (double z, double x, double y, double w) {
    double r = sqrt (x * x + y * y + w * w);
    if (r < 1e-300) return 0.0;
    return n2 (z) * r * exp (-z * r) * 0.4886025119029199 * (w / r);
}
static double o2px (double z, double x, double y, double w) {
    double r = sqrt (x * x + y * y + w * w);
    if (r < 1e-300) return 0.0;
    return n2 (z) * r * exp (-z * r) * 0.4886025119029199 * (x / r);
}

/* midpoint grid overlap of fa (centered -R/2) and fb (+R/2) */
typedef double (*orb_t)(double, double, double, double);
static double grid_overlap (orb_t fa, double za, orb_t fb, double zb,
                            double R, double L, int N) {
    double h = 2.0 * L / N;
    double s = 0.0;
    for (int i = 0; i < N; i++) {
        double x = -L + (i + 0.5) * h;
        for (int j = 0; j < N; j++) {
            double y = -L + (j + 0.5) * h;
            for (int k = 0; k < N; k++) {
                double z = -L + (k + 0.5) * h;
                s += fa (za, x, y, z + R / 2) * fb (zb, x, y, z - R / 2);
            }
        }
    }
    return s * h * h * h;
}

static void test_overlap_grid (void) {
    printf ("== EHT overlaps vs brute-force 3D grid ==\n");
    const double zH = 1.0 / 0.529177210903;   /* H 1s, independent 1/a0 */
    const double z2 = 1.0 / (2.0 * 0.529177210903); /* H-like 2p test exp */
    char d[160];
    /* 1s-1s at R=1.0: closed S(rho) = e^-rho(1+rho+rho^2/3), rho=zR */
    {
        double R = 1.0;
        double eng = eht_overlap_ss (1, zH, 1, zH, R);
        double grd = grid_overlap (o1s, zH, o1s, zH, R, 7.0, 110);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", eng, grd);
        ok ("ss/H2R1 vs grid", fabs (eng - grd) < 2e-3, d);
        double rho = zH * R;
        double cf = exp (-rho) * (1 + rho + rho * rho / 3.0);
        okrel ("ss/H2R1 closed form", eng, cf, 1e-9);
    }
    /* R=0 normalization */
    okrel ("ss R=0 is 1", eht_overlap_ss (1, zH, 1, zH, 0.0), 1.0, 1e-12);
    /* large R vanishes */
    ok ("ss R=20 ~ 0", fabs (eht_overlap_ss (1, zH, 1, zH, 20.0)) < 1e-9,
        "");
    /* 1s-2s */
    {
        double R = 1.5;
        double eng = eht_overlap_ss (1, zH, 2, z2, R);
        double grd = grid_overlap (o1s, zH, o2s, z2, R, 9.0, 120);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", eng, grd);
        ok ("1s2s vs grid", fabs (eng - grd) < 3e-3, d);
    }
    /* s-psigma both orientations, antisymmetric pair check */
    {
        double R = 1.5;
        double e1 = eht_overlap_sps (1, zH, 2, z2, R, 1);
        double g1 = grid_overlap (o1s, zH, o2pz, z2, R, 9.0, 120);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", e1, g1);
        ok ("s-ps(p on B) vs grid", fabs (e1 - g1) < 3e-3, d);
        double e0 = eht_overlap_sps (1, zH, 2, z2, R, 0);
        double g0 = grid_overlap (o2pz, z2, o1s, zH, R, 9.0, 120);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", e0, g0);
        ok ("s-ps(p on A) vs grid", fabs (e0 - g0) < 3e-3, d);
    }
    /* p-p sigma and pi */
    {
        double R = 1.5;
        double es = eht_overlap_pps (z2, z2, R);
        double gs = grid_overlap (o2pz, z2, o2pz, z2, R, 9.0, 120);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", es, gs);
        ok ("pp-sigma vs grid", fabs (es - gs) < 3e-3, d);
        double ep = eht_overlap_ppi (z2, z2, R);
        double gp = grid_overlap (o2px, z2, o2px, z2, R, 9.0, 120);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", ep, gp);
        ok ("pp-pi vs grid", fabs (ep - gp) < 3e-3, d);
        /* pi closed form spot: e^-r(1+r+2r^2/5+r^3/15) with r=zR */
        double r = z2 * R;
        double cf = exp (-r) * (1 + r + 0.4 * r * r + r * r * r / 15.0);
        snprintf (d, sizeof d, "(eng %.6f closed %.6f)", ep, cf);
        ok ("pp-pi closed form", fabs (ep - cf) < 1e-6, d);
    }
    /* heteronuclear exactness (za != zb): grid decides */
    {
        double R = 1.2, zb = 1.3 * zH;
        double eng = eht_overlap_ss (1, zH, 1, zb, R);
        double grd = grid_overlap (o1s, zH, o1s, zb, R, 7.0, 110);
        snprintf (d, sizeof d, "(eng %.6f grid %.6f)", eng, grd);
        ok ("hetero ss vs grid", fabs (eng - grd) < 2e-3, d);
    }
}

static Simulation *mk_h2 (double R) {
    Simulation *s = sim_create (8, 8);
    if (!s) return NULL;
    sim_add_atom (s, 1, vec3 (-R / 2, 0, 0), 0.0);
    sim_add_atom (s, 1, vec3 (R / 2, 0, 0), 0.0);
    s->dt = 0.5;
    s->cutoff = 12.0;
    return s;
}

static void test_h2_laws (void) {
    printf ("== H2 EHT laws ==\n");
    char d[160];
    /* 2x2 closed form: e+- = (H11 +- H12)/(1 +- S) */
    {
        Simulation *s = mk_h2 (0.74);
        static qm_eht_t E;
        int rc = qm_eht_solve (s, &E);
        snprintf (d, sizeof d, "(rc=%d)", rc);
        ok ("H2 solves", rc == 0, d);
        if (rc == 0) {
            double zH = 1.0 / 0.529177210903;   /* H 1s exponent */
            double S = eht_overlap_ss (1, zH, 1, zH, 0.74);
            double H11 = -13.6, H12 = 1.75 * S * H11;
            double ep = (H11 + H12) / (1 + S), em = (H11 - H12) / (1 - S);
            snprintf (d, sizeof d, "(e0 %.6f e1 %.6f)", E.eps[0], E.eps[1]);
            ok ("H2 MO zn closed-2x2",
                fabs (E.eps[0] - (ep < em ? ep : em)) < 1e-6 &&
                    fabs (E.eps[1] - (ep < em ? em : ep)) < 1e-6,
                d);
            double bo = qm_eht_bond_order (&E, 0, 1);
            snprintf (d, sizeof d, "(BO %.4f)", bo);
            ok ("H2 BO ~ 1", bo > 0.8 && bo < 1.2, d);
            double q = E.charges[0] + E.charges[1];
            snprintf (d, sizeof d, "(sum %.3e)", q);
            ok ("H2 neutrality", fabs (q) < 1e-9, d);
            snprintf (d, sizeof d, "(e0 %.6f e1 %.6f gap %.6f)",
                      E.eps[0], E.eps[1], E.gap);
            /* E7: minimal-basis EHT gaps run large (big-S x K=1.75);
             * band pins finiteness/determinism + ballpark, not spectroscopy. */
            ok ("H2 gap physical", E.gap > 20.0 && E.gap < 50.0, d);
            snprintf (d, sizeof d, "(scc %d resid %.1e)", E.scc_iters,
                      E.scc_resid);
            ok ("H2 SCC converged", E.scc_resid < 1e-6, d);
        }
        sim_destroy (s);
    }
    /* dissociation: BO dissolves (the anti-harmonic gate) */
    {
        Simulation *s = mk_h2 (4.0);
        static qm_eht_t E;
        int rc = qm_eht_solve (s, &E);
        double bo = rc == 0 ? qm_eht_bond_order (&E, 0, 1) : -1.0;
        snprintf (d, sizeof d, "(rc=%d BO %.4f)", rc, bo);
        ok ("H2(4A) solves", rc == 0, d);
        /* E8/E12: at floor Fermi kT the near-degenerate stretched pair
         * fractionalizes (Mermin finite-T ensemble), so Mayer falls with
         * the switch together here; at zero T the restricted determinant
         * would keep B~1 (RHF static-correlation failure, inherited and
         * documented in qm_eht.h). Either way the printed BO is pure
         * Mayer with its notes, and the switch carries dissociation. */
        ok ("H2(4A) Mayer dissolved", bo >= 0.0 && bo < 0.6, d);
        if (rc == 0) {
            snprintf (d, sizeof d, "(Eband %.4f vs 2*H %.4f)", E.e_band,
                      2.0 * -13.6);
            ok ("H2(4A) band dissociates", fabs (E.e_band + 27.2) < 0.5, d);
        }
        sim_destroy (s);
    }
    /* overlap-gated switch: unity at equilibrium, dissolved at 4 A.
     * THIS is the reactive milestone (harmonic wall removed). */
    {
        Simulation *se = mk_h2 (0.74);
        static qm_eht_t Ee;
        double fe = -1.0;
        if (qm_eht_solve (se, &Ee) == 0)
            fe = qm_eht_bond_factor (&Ee, se, 0, 1);
        snprintf (d, sizeof d, "(switch %.4f)", fe);
        ok ("H2(eq) switch ~ 1", fe > 0.9 && fe <= 1.0, d);
        sim_destroy (se);
        Simulation *sd = mk_h2 (4.0);
        static qm_eht_t Ed;
        double fd = -1.0;
        if (qm_eht_solve (sd, &Ed) == 0)
            fd = qm_eht_bond_factor (&Ed, sd, 0, 1);
        snprintf (d, sizeof d, "(switch %.4f)", fd);
        ok ("H2(4A) switch dissolved", fd >= 0.0 && fd < 0.1, d);
        sim_destroy (sd);
    }
    /* E11: stretched heteronuclear level-crossing flip has no reachable
     * fixed point (annealing tried, reverted — roots drift with kT).
     * Gate pins FAIL-LOUD (nonzero rc), not silence: a future level-shift
     * that truly converges will trip this pin and must update it. */
    {
        Simulation *sx = sim_create (8, 8);
        sim_add_atom (sx, 8, vec3 (0, 0, 0), 0.0);
        sim_add_atom (sx, 1, vec3 (0, 0, 1.9144), 0.0);
        sim_add_atom (sx, 1, vec3 (0, 0, -1.9144), 0.0);
        sx->dt = 0.5;
        sx->cutoff = 12.0;
        static qm_eht_t Ex;
        int rcx = qm_eht_solve (sx, &Ex);
        snprintf (d, sizeof d, "(rc=%d)", rcx);
        ok ("stretched heteronuclear fails loud (E11)", rcx != 0, d);
        sim_destroy (sx);
    }
    /* unsupported / over-cap fail loud */
    {
        Simulation *s = sim_create (8, 8);
        sim_add_atom (s, 2, vec3 (0, 0, 0), 0.0);   /* He: no EHT basis */
        static qm_eht_t E;
        ok ("He refused", qm_eht_solve (s, &E) != 0, "");
        sim_destroy (s);
    }
}

static void test_h2o_laws (void) {
    printf ("== H2O EHT laws ==\n");
    char d[160];
    Simulation *s = sim_create (8, 8);
    if (!s) {
        ok ("h2o alloc", 0, "");
        return;
    }
    sim_add_atom (s, 8, vec3 (0, 0, 0), 0.0);
    double ang = 104.52 * 3.14159265358979323846 / 180.0;
    sim_add_atom (s, 1, vec3 (0.9572 * sin (ang / 2), 0.9572 * cos (ang / 2), 0), 0.0);
    sim_add_atom (s, 1, vec3 (-0.9572 * sin (ang / 2), 0.9572 * cos (ang / 2), 0), 0.0);
    s->dt = 0.5;
    s->cutoff = 12.0;
    static qm_eht_t E;
    int rc = qm_eht_solve (s, &E);
    snprintf (d, sizeof d, "(rc=%d)", rc);
    ok ("H2O solves", rc == 0, d);
    if (rc == 0) {
        snprintf (d, sizeof d, "(qO %.4f qH %.4f %.4f)", E.charges[0],
                  E.charges[1], E.charges[2]);
        ok ("H2O O negative, H positive",
            E.charges[0] < -0.15 && E.charges[1] > 0.05 &&
                E.charges[2] > 0.05,
            d);
        double sum = E.charges[0] + E.charges[1] + E.charges[2];
        snprintf (d, sizeof d, "(sum %.3e)", sum);
        ok ("H2O neutrality", fabs (sum) < 1e-9, d);
        snprintf (d, sizeof d, "(gap %.3f)", E.gap);
        ok ("H2O gap physical", E.gap > 10.0 && E.gap < 45.0, d);
        /* rotation invariance (Slater-Koster assembly check) */
        for (int i = 0; i < s->num_atoms; i++) {
            Vec3 p = s->atoms[i].position;
            s->atoms[i].position = vec3 (-p.y, p.x, p.z + 0.5);
        }
        static qm_eht_t E2;
        int rc2 = qm_eht_solve (s, &E2);
        snprintf (d, sizeof d, "(rc=%d)", rc2);
        ok ("H2O rotated solves", rc2 == 0, d);
        if (rc2 == 0) {
            double dq = fabs (E2.charges[0] - E.charges[0]) +
                        fabs (E2.charges[1] - E.charges[1]) +
                        fabs (E2.charges[2] - E.charges[2]);
            snprintf (d, sizeof d, "(dq %.3e dgap %.3e)", dq,
                      fabs (E2.gap - E.gap));
            ok ("rotation-invariant charges+gap",
                dq < 1e-6 && fabs (E2.gap - E.gap) < 1e-6, d);
        }
    }
    sim_destroy (s);
    /* VSIP table spots (Hoffmann 1963) */
    ok ("VSIP table", eht_vsip (1, 1, 0) == 13.6 &&
            eht_vsip (6, 2, 0) == 21.4 && eht_vsip (6, 2, 1) == 11.4 &&
            eht_vsip (7, 2, 0) == 26.0 && eht_vsip (7, 2, 1) == 13.4 &&
            eht_vsip (8, 2, 0) == 32.3 && eht_vsip (8, 2, 1) == 14.8 &&
            eht_vsip (2, 1, 0) == 0.0,
        "");
}

int main (void) {
    test_overlap_grid ();
    test_h2_laws ();
    test_h2o_laws ();
    printf ("EHT: %d passed, %d failed\n", npass, nfail);
    return nfail ? 1 : 0;
}
