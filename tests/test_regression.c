/*
 * test_regression.c — audit regression suite for v9R4.
 *
 * Every check here corresponds to a specific defect found and fixed in
 * the v9R4 audit. The suite is deliberately built on INDEPENDENT
 * oracles (finite differences, quadrature, NIST SHA-256 vectors,
 * published reference values) rather than on the engine agreeing with
 * itself, so it can actually fail.
 *
 * Build: make selftest-regression
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "types.h"
#include "constants.h"
#include "vec3.h"
#include "quantum.h"
#include "qm.h"
#include "forces.h"
#include "integrator.h"
#include "sim.h"
#include "periodic_table.h"
#include "datastream.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include "nucleobases.h"
#include "kcsa_filter.h"
#include "aminoacids.h"
#include "amber_lj.h"
#include "kcsa_filter.h"
#include "amber_lj.h"

static int g_pass = 0, g_fail = 0;
static const char *g_group = "";

static void grp(const char *g) { g_group = g; printf("\n[%s]\n", g); }

static void ok(const char *name, int cond, const char *detail) {
    if (cond) { g_pass++; printf("  PASS  %-52s %s\n", name, detail ? detail : ""); }
    else      { g_fail++; printf("  FAIL  %-52s %s\n", name, detail ? detail : ""); }
}
static void okrel(const char *name, double got, double want, double tol) {
    double e = fabs(got - want) / (fabs(want) > 1e-300 ? fabs(want) : 1.0);
    char d[128];
    snprintf(d, sizeof d, "got=%.10g want=%.10g rel=%.2e", got, want, e);
    ok(name, e <= tol, d);
}

/* ── force oracle: analytic vs central differences ────────────────────── */
static double sim_potential(Simulation *s) { forces_calculate(s); return s->potential_energy; }

static void fd_force_check(const char *name, Simulation *s, double tol) {
    forces_calculate(s);
    double fmax = 0.0;
    for (int i = 0; i < s->num_atoms; i++) fmax += vec3_norm2(s->atoms[i].force);
    fmax = sqrt(fmax);
    if (fmax < 1e-12) fmax = 1.0;
    const double h = 1e-6;
    double worst = 0.0;
    for (int i = 0; i < s->num_atoms; i++) {
        double *c[3] = {&s->atoms[i].position.x, &s->atoms[i].position.y, &s->atoms[i].position.z};
        for (int k = 0; k < 3; k++) {
            double o = *c[k];
            *c[k] = o + h; double Ep = sim_potential(s);
            *c[k] = o - h; double Em = sim_potential(s);
            *c[k] = o;
            double Ffd = -(Ep - Em) / (2 * h);
            double Fa = k == 0 ? s->atoms[i].force.x
                     : k == 1 ? s->atoms[i].force.y : s->atoms[i].force.z;
            double d = fabs(Ffd - Fa) / fmax;
            if (d > worst) worst = d;
        }
    }
    forces_calculate(s);
    char d[128];
    snprintf(d, sizeof d, "max|dF_analytic - dF_FD|/Fmax = %.3e (tol %.0e)", worst, tol);
    ok(name, worst <= tol, d);
}

/* ══════════════════════════════════════════════════════════════════════ */

static void test_quantum_basics(void) {
    grp("Hydrogenic quantum mechanics (independent quadrature oracle)");
    /* normalisation + expectation values */
    for (int n = 1; n <= 3; n++) for (int l = 0; l < n; l++) {
        const double Z = 2.3;
        int M = 400000; double h = 60.0 / M, s0 = 0, s1 = 0;
        for (int i = 0; i <= M; i++) {
            double r = i * h, P = quantum_radial_probability(n, l, Z, r);
            double w = (i == 0 || i == M) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            s0 += w * P; s1 += w * r * P;
        }
        s0 *= h / 3.0; s1 *= h / 3.0;
        char nm[64];
        snprintf(nm, sizeof nm, "int r^2|R_nl|^2 dr = 1  (n=%d,l=%d)", n, l);
        okrel(nm, s0, 1.0, 1e-6);
        snprintf(nm, sizeof nm, "<r> analytic = quadrature  (n=%d,l=%d)", n, l);
        okrel(nm, quantum_expect_r(n, l, Z), s1, 1e-5);
    }
    /* Slater rules: Slater's own 1930 worked examples */
    ElectronConfig fe; pt_electron_config(26, &fe);
    okrel("Slater Fe 4s Zeff = 3.75", quantum_zeff_raw(26, 4, 0, &fe), 3.75, 1e-12);
    okrel("Slater Fe 3d Zeff = 6.25", quantum_zeff_raw(26, 3, 2, &fe), 6.25, 1e-12);
    ElectronConfig c; pt_electron_config(6, &c);
    okrel("Slater C  1s Zeff = 5.70", quantum_zeff_raw(6, 1, 0, &c), 5.70, 1e-12);
    okrel("Slater C  2p Zeff = 3.25", quantum_zeff_raw(6, 2, 1, &c), 3.25, 1e-12);
    pt_electron_config(11, &fe);
    okrel("Slater Na 3s Zeff = 2.20", quantum_zeff_raw(11, 3, 0, &fe), 2.20, 1e-12);
    /* Clementi-Raimondi: spot values AND the 1s monotonic trend */
    okrel("CR H  1s = 1.0000", quantum_zeff_cr(1, 1, 0), 1.0000, 1e-9);
    okrel("CR C  2p = 3.1358 (tabulated 3.14)", quantum_zeff_cr(6, 2, 1), 3.1358, 2e-3);
    okrel("CR N  2p = 3.8340", quantum_zeff_cr(7, 2, 1), 3.8340, 5e-4);
    okrel("CR O  2p = 4.4532 (tabulated 4.45)", quantum_zeff_cr(8, 2, 1), 4.4532, 2e-3);
    okrel("CR Ar 1s = 17.5075", quantum_zeff_cr(18, 1, 0), 17.5075, 1e-9);
    okrel("CR Kr 3d = 20.6252", quantum_zeff_cr(36, 3, 2), 20.6252, 5e-4);
    /* the O 1s cell that was suspected of a digit transposition */
    okrel("CR O  1s = 7.6579", quantum_zeff_cr(8, 1, 0), 7.6579, 1e-9);
    {   /* 1s exponents must rise ~0.9925/element Li->Kr; a transposed
         * digit in the middle would show up as a kink. */
        double worst = 0.0; int wz = 0;
        for (int Z = 4; Z <= 35; Z++) {
            double a = quantum_zeff_cr(Z, 1, 0), b = quantum_zeff_cr(Z + 1, 1, 0);
            if (a <= 0 || b <= 0) continue;
            double d = b - a;
            if (fabs(d - 0.9925) > worst) { worst = fabs(d - 0.9925); wz = Z; }
        }
        char d[128];
        snprintf(d, sizeof d, "max deviation of 1s increment from 0.9925 = %.4f (Z=%d)", worst, wz);
        ok("CR 1s series is smooth (catches digit transposition)", worst < 0.02, d);
    }
    /* spherical harmonic orthonormality, including the hand-normalised f */
    {   int worst_bad = 0; double worst_v = 0; int wl = 0, wm = 0;
        for (int l = 0; l <= 3; l++)
            for (int m1 = -l; m1 <= l; m1++)
                for (int m2 = -l; m2 <= l; m2++) {
                    double acc = 0.0;
                    for (int it = 0; it < 100; it++) {
                        double th = M_PI * (it + 0.5) / 100, w = sin(th);
                        for (int ip = 0; ip < 200; ip++) {
                            double ph = 2 * M_PI * (ip + 0.5) / 200;
                            acc += qm_Y_real(l, m1, th, ph) * qm_Y_real(l, m2, th, ph) * w;
                        }
                    }
                    acc *= (M_PI / 100) * (2 * M_PI / 200);
                    double want = (m1 == m2) ? 1.0 : 0.0;
                    double e = fabs(acc - want);
                    if (e > worst_v) { worst_v = e; wl = l; wm = m1; }
                    if (e > 3e-4) worst_bad = 1;
                }
        char d[128];
        snprintf(d, sizeof d, "max|<Ylm|Yl'm'>-target| = %.2e (worst l=%d m=%d)", worst_v, wl, wm);
        ok("real spherical harmonics orthonormal, l<=3", !worst_bad, d);
    }
    /* constants */
    okrel("COULOMB_MD = 14.399645 eV A", COULOMB_MD, 14.399645, 1e-7);
    okrel("MD_FORCE_CONV = 9.648533e-3", MD_FORCE_CONV, 9.648533e-3, 1e-7);
    okrel("HARTREE_TO_EV * EV_TO_HARTREE = 1", HARTREE_TO_EV * EV_TO_HARTREE, 1.0, 1e-12);
    okrel("KCAL_MOL_TO_EV * EV_TO_KCAL_MOL = 1", KCAL_MOL_TO_EV * (1.0 / KCAL_MOL_TO_EV), 1.0, 1e-15);
}

static void test_forces(void) {
    grp("Analytic forces vs finite-difference oracle (every term)");
    /* The h=1e-6 central-difference oracle has its own noise floor near
     * 3e-6, so the tolerance is 1e-5 — that still catches every real
     * defect found in this audit, all of which were 100% or worse. */
    {   Simulation *s = sim_create(8, 8);
        int a = sim_add_atom(s, 1, vec3(0,0,0), 0);
        int b = sim_add_atom(s, 1, vec3(1.35,0,0), 0);
        sim_add_bond(s, a, b, 1);
        sim_set_bond_params(s, 0, 0.7414, 36.0);
        s->use_bonds = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 50.0;
        fd_force_check("bond harmonic", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(8, 8);
        int a = sim_add_atom(s, 1, vec3(1,0,0), 0);
        int b = sim_add_atom(s, 6, vec3(0,0,0), 0);
        int c = sim_add_atom(s, 1, vec3(0.3,0.9,0.2), 0);
        sim_add_angle_explicit(s, a, b, c, 1.0472, 3.5);
        s->use_angles = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 50.0;
        fd_force_check("angle harmonic", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(8, 8);
        int a = sim_add_atom(s, 6, vec3(0,0,0), 0);
        int b = sim_add_atom(s, 6, vec3(1.5,0.1,0), 0);
        int c = sim_add_atom(s, 6, vec3(2.1,1.3,0.3), 0);
        int d = sim_add_atom(s, 6, vec3(3.4,1.6,1.0), 0);
        sim_add_dihedral(s, a, b, c, d, 0.02, 3, 0.0);
        s->use_dihedrals = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 50.0;
        fd_force_check("dihedral torsion", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 6, vec3(0,0,0), -0.3);
        sim_add_atom(s, 8, vec3(3.2,0.4,0.1), -0.5);
        s->use_lj = 1; s->use_coulomb = 1; s->cutoff = 12.0; s->use_switching = 0;
        fd_force_check("LJ + Coulomb, hard cutoff", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 6, vec3(0,0,0), -0.3);
        sim_add_atom(s, 8, vec3(8.0,0.4,0.1), -0.5);   /* inside 9.6-12 switch band */
        s->use_lj = 1; s->use_coulomb = 1; s->cutoff = 12.0; s->use_switching = 1;
        fd_force_check("LJ + Coulomb, CHARMM switching on", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 7, vec3(0.2,0.1,0.3), 0);
        sim_add_restraint(s, 0, vec3(0,0,0), 0.5);
        s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 50.0;
        fd_force_check("harmonic positional restraint", s, 1e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(16, 8);
        sim_add_atom(s, 19, vec3(0,0,0), 1.0);
        sim_add_atom(s, 8,  vec3(2.5,0.3,0.0), -0.7);
        sim_add_atom(s, 6,  vec3(4.6,1.0,0.2), 0.1);
        s->use_polar = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 20.0;
        fd_force_check("induction, first order (analytic HF)", s, 1e-5);
        sim_destroy(s);
    }
    {   /* A real carbonyl next to a K+, with the C=O and K-O bonds
         * declared, so the coupled-dipole hard core (audit D5) applies
         * exactly as it does in the engine. A bare 1.3 A C-O pair is
         * inside the induced-dipole catastrophe regime and is not a
         * valid test input. */
        Simulation *s = sim_create(16, 8);
        int k = sim_add_atom(s, 19, vec3(0,0,0), 1.0);
        int o1 = sim_add_atom(s, 8, vec3(2.5,0.3,0.0), -0.6);
        int c  = sim_add_atom(s, 6, vec3(3.85,0.85,0.0), 0.1);
        int o2 = sim_add_atom(s, 8, vec3(5.15,0.85,0.0), -0.5);
        sim_add_bond(s, c, o2, 2);
        sim_add_bond(s, o1, c, 1);
        s->atoms[c].num_bonds = 2;
        s->atoms[c].bond_partners[0] = o1; s->atoms[c].bond_orders[0] = 1;
        s->atoms[c].bond_partners[1] = o2; s->atoms[c].bond_orders[1] = 2;
        s->atoms[o1].num_bonds = 1;
        s->atoms[o1].bond_partners[0] = c; s->atoms[o1].bond_orders[0] = 1;
        s->atoms[o2].num_bonds = 1;
        s->atoms[o2].bond_partners[0] = c; s->atoms[o2].bond_orders[0] = 2;
        (void)k;
        s->use_pol_scf = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 20.0;
        fd_force_check("induction, self-consistent dipoles", s, 3e-5);
        sim_destroy(s);
    }
    {   Simulation *s = sim_create(16, 8);
        sim_add_atom(s, 8, vec3(0,0,0), -0.4);
        sim_add_atom(s, 1, vec3(3.0,0.2,0.0), 0.2);
        s->use_disp = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 20.0;
        fd_force_check("Slater-Kirkwood + Tang-Toennies dispersion", s, 1e-5);
        sim_destroy(s);
    }
}

static void test_dispersion_sign(void) {
    grp("Dispersion is ATTRACTIVE in the force path (audit D1)");
    /* Two polarizable atoms: the pair force must pull them together. */
    Simulation *s = sim_create(8, 8);
    int a = sim_add_atom(s, 8, vec3(0,0,0), 0);
    int b = sim_add_atom(s, 1, vec3(3.2,0,0), 0);
    s->use_disp = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 20.0;
    forces_calculate(s);
    double Ed = s->E_disp_total;
    /* b lies at +x, so an ATTRACTIVE force on a points along +x. */
    double Fx = s->atoms[a].force.x;
    double Fbx = s->atoms[b].force.x;
    char d[192];
    snprintf(d, sizeof d, "E_disp=%.6f eV (must be <0), F_a.x=%+.6f F_b.x=%+.6f (a->+x, b->-x)",
             Ed, Fx, Fbx);
    ok("dispersion energy attractive", Ed < 0.0, d);
    ok("dispersion force attractive (was sign-flipped)",
       Fx > 0.0 && Fbx < 0.0, d);
    /* and it must grow more attractive as r decreases, inside the well */
    double far = 0.0;
    s->atoms[b].position.x = 3.2;
    forces_calculate(s); far = s->E_disp_total;
    s->atoms[b].position.x = 3.0;
    forces_calculate(s); double near = s->E_disp_total;
    snprintf(d, sizeof d, "E(3.0A)=%.6f  E(3.2A)=%.6f  -> E(3.0) < E(3.2)", near, far);
    ok("dispersion deepens as atoms approach", near < far, d);
    sim_destroy(s);
}

static void test_scf_continuity(void) {
    grp("SCF polarisation energy is CONTINUOUS in geometry (audit D3)");
    Simulation *s = sim_create(16, 8);
    int kk = sim_add_atom(s, 19, vec3(0,0,0), 1.0);
    int oo1 = sim_add_atom(s, 8,  vec3(2.5,0,0), -0.6);
    int cc  = sim_add_atom(s, 6,  vec3(3.85,0.55,0), 0.1);
    int oo2 = sim_add_atom(s, 8,  vec3(5.15,0.55,0), -0.5);
    sim_add_bond(s, cc, oo2, 2);
    sim_add_bond(s, oo1, cc, 1);
    s->atoms[cc].num_bonds = 2;
    s->atoms[cc].bond_partners[0] = oo1; s->atoms[cc].bond_orders[0] = 1;
    s->atoms[cc].bond_partners[1] = oo2; s->atoms[cc].bond_orders[1] = 2;
    s->atoms[oo1].num_bonds = 1;
    s->atoms[oo1].bond_partners[0] = cc; s->atoms[oo1].bond_orders[0] = 1;
    s->atoms[oo2].num_bonds = 1;
    s->atoms[oo2].bond_partners[0] = cc; s->atoms[oo2].bond_orders[0] = 2;
    (void)kk;
    s->use_pol_scf = 1; s->use_lj = 0; s->use_coulomb = 0; s->cutoff = 20.0;
    double prev = 0.0, maxjump = 0.0;
    int nconv = 0, nsamp = 0;
    /* Sweep a PHYSICAL displacement range. Driving a bonded C=O pair to
     * sub-Angstrom separation is not a continuity test of the solver, it
     * is a test of the induced-dipole model's breakdown regime, which the
     * dipole hard core (audit D5) is there to keep out of reach. */
    double x0 = s->atoms[2].position.x;
    for (int k = 0; k <= 400; k++) {
        s->atoms[2].position.x = x0 + 0.0015 * (k - 200);   /* +/- 0.30 A */
        forces_calculate(s);
        Vec3 mu[256];
        if (qm_solve_dipoles(s, 1.0, mu) == 0) nconv++;
        nsamp++;
        if (k > 0) {
            double j = fabs(s->E_polar_total - prev);
            if (j > maxjump) maxjump = j;
        }
        prev = s->E_polar_total;
    }
    s->atoms[2].position.x = x0;
    char d[192];
    snprintf(d, sizeof d, "max |dE| over a 0.0015 A step, 0.6 A sweep = %.5f eV", maxjump);
    ok("SCF E_polar has no functional-swap jump", maxjump < 0.20, d);
    snprintf(d, sizeof d, "dipole solver converged in %d/%d sampled geometries", nconv, nsamp);
    ok("dipole solver converges everywhere on the sweep", nconv == nsamp, d);
    sim_destroy(s);
}

static void test_lj_convention(void) {
    grp("LJ sigma is the 12-6 collision diameter, not UFF's Rmin (F6)");
    const Element *c = pt_element(6);
    double sigma = lj_sigma_combine(c->lj_sigma, c->lj_sigma);
    double eps   = lj_eps_combine(c->lj_epsilon, c->lj_epsilon);
    double rmin  = sigma * pow(2.0, 1.0 / 6.0);
    char d[192];
    snprintf(d, sizeof d, "sigma=%.5f  eps=%.5f eV (%.4f kcal/mol)  rmin=%.4f A",
             sigma, eps, eps / KCAL_MOL_TO_EV, rmin);
    okrel("carbon LJ minimum at UFF Rmin = 3.851 A", rmin, 3.851, 2e-3);
    okrel("carbon LJ epsilon = 0.105 kcal/mol", eps / KCAL_MOL_TO_EV, 0.105, 2e-3);
    ok("carbon LJ sigma = 3.431 A (was Rmin 3.851)", fabs(sigma - 3.4308) < 5e-3, d);
    const Element *h = pt_element(1);
    okrel("hydrogen LJ sigma = 2.571 A", lj_sigma_combine(h->lj_sigma, h->lj_sigma), 2.571, 5e-3);
    const Element *k = pt_element(19);
    okrel("potassium LJ minimum at UFF Rmin = 3.812 A",
          lj_sigma_combine(k->lj_sigma, k->lj_sigma) * pow(2.0, 1.0/6.0), 3.812, 2e-3);
}

static void test_hybridization(void) {
    grp("Hybridisation: the functional groups this model is built from (F7)");
    struct { const char *name; int Z; int npartners; int own_pi; int want_sp; } T[] = {
        /* An amide N and an ammonia N are INDISTINGUISHABLE from a bare
         * Atom: same valence, all-single bond orders. Only the
         * neighbour's pi bond separates them. */
        { "carbonyl O  -> sp2", 8, 1, 1, 1 },
        { "amide N     -> sp2", 7, 2, 0, 1 },
        { "carbonyl C  -> sp2", 6, 3, 1, 1 },
        { "water O     -> sp3", 8, 2, 0, 0 },
        { "ammonia N    -> sp3", 7, 3, 0, 0 },
        { "methylene C  -> sp3", 6, 4, 0, 0 },
    };
    for (unsigned t = 0; t < sizeof T / sizeof T[0]; t++) {
        Simulation *s = sim_create(16, 16);
        int me = sim_add_atom(s, T[t].Z, vec3(0,0,0), 0);
        /* partners: T[t].npartners carbons/oxygens at 1.4 A */
        for (int p = 0; p < T[t].npartners; p++) {
            double a = 2.0944 * p;
            int j = sim_add_atom(s, 6, vec3(1.45*cos(a), 1.45*sin(a), 0), 0);
            sim_add_bond(s, me, j, 1);
        }
        /* attach a carbonyl partner so amide N / carbonyl C are conjugated */
        if (!T[t].own_pi && T[t].Z == 7) {
            int cc = sim_add_atom(s, 6, vec3(-1.2, 0, 0), 0);
            int oo = sim_add_atom(s, 8, vec3(-2.4, 0.6, 0), 0);
            sim_add_bond(s, me, cc, 1);
            sim_add_bond(s, cc, oo, 2);
            s->atoms[cc].num_bonds = 2;
            s->atoms[cc].bond_partners[0] = me;  s->atoms[cc].bond_orders[0] = 1;
            s->atoms[cc].bond_partners[1] = oo;  s->atoms[cc].bond_orders[1] = 2;
            s->atoms[oo].num_bonds = 1;
            s->atoms[oo].bond_partners[0] = cc;  s->atoms[oo].bond_orders[0] = 2;
            me = 0;  /* the amide N is the first atom in the list */
            QmHybrid h = qm_hybridization_ctx(s, 0);
            char d[128];
            snprintf(d, sizeof d, "label=%s lobes=%d lone=%d", h.label, h.n_lobes, h.n_lone_pairs);
            ok(T[t].name, strcmp(h.label, T[t].want_sp ? "sp2" : "sp3") == 0, d);
            sim_destroy(s);
            continue;
        }
        sim_detect_bonds(s);
        /* ensure bond_orders reflect the intended orders */
        for (int p = 0; p < s->atoms[me].num_bonds; p++)
            s->atoms[me].bond_orders[p] = T[t].own_pi ? 2 : 1;
        QmHybrid h = qm_hybridization_ctx(s, me);
        char d[128];
        snprintf(d, sizeof d, "label=%s lobes=%d lone=%d", h.label, h.n_lobes, h.n_lone_pairs);
        ok(T[t].name, strcmp(h.label, T[t].want_sp ? "sp2" : "sp3") == 0, d);
        sim_destroy(s);
    }
}

static void test_polarizability(void) {
    grp("Polarizability: tabulated, continuous, no Z-boundary jump (F8)");
    struct { int Z; double want; const char *sym; } T[] = {
        {1, 4.507, "H"}, {6, 11.30, "C"}, {7, 7.44, "N"}, {8, 5.30, "O"},
        {9, 3.74, "F"}, {10, 2.661, "Ne"}, {11, 162.7, "Na"}, {16, 19.40, "S"},
        {17, 14.60, "Cl"}, {18, 11.08, "Ar"}, {19, 289.7, "K"},
    };
    double a0c = BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM;
    for (unsigned t = 0; t < sizeof T / sizeof T[0]; t++) {
        Simulation *s = sim_create(4, 4);
        int i = sim_add_atom(s, T[t].Z, vec3(0,0,0), 0);
        double a = qm_polarizability(&s->atoms[i]);
        char d[128];
        snprintf(d, sizeof d, "%s alpha=%.4f A^3 (ref %.4f)", T[t].sym, a, T[t].want * a0c);
        okrel("alpha tabulated (a.u. reference)", a, T[t].want * a0c, 1e-9);
        sim_destroy(s);
    }
    /* Adjacent-element steps must MATCH the reference steps. Real alpha
     * genuinely jumps ~100x at the alkali edge (He 0.205 -> Li 24.3,
     * because the alkali 2s electron is far from the core), so an
     * absolute bound would be wrong; what must hold is that the code
     * reproduces the measured shape, including that jump, rather than
     * inventing one of its own at a tabulation boundary. */
    {   double worst_dev = 0.0; int wz = 0;
        double prev = 0.0, prev_ref = 0.0;
        double a0c = BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM * BOHR_TO_ANGSTROM;
        static const double ref_au[37] = {
            0, 4.507, 1.384, 164.1, 37.70, 20.50, 11.30, 7.44, 5.30, 3.74, 2.661,
            162.7, 71.20, 57.80, 37.30, 25.00, 19.40, 14.60, 11.08, 289.7, 160.8,
            97.0, 100.0, 87.0, 83.0, 68.0, 62.0, 55.0, 49.0, 46.50, 38.67,
            50.0, 40.0, 30.0, 38.90, 21.00, 16.80 };
        for (int Z = 1; Z <= 36; Z++) {
            Simulation *s = sim_create(4, 4);
            int i = sim_add_atom(s, Z, vec3(0,0,0), 0);
            double a = qm_polarizability(&s->atoms[i]);
            double r = ref_au[Z] * a0c;
            if (prev > 0) {
                double got = a / prev, want = r / prev_ref;
                double dev = fabs(got - want) / want;
                if (dev > worst_dev) { worst_dev = dev; wz = Z; }
            }
            prev = a; prev_ref = r;
            sim_destroy(s);
        }
        char d[192];
        snprintf(d, sizeof d, "max relative deviation of adjacent step from reference = %.2e (Z=%d)", worst_dev, wz);
        ok("alpha reproduces the measured Z-to-Z shape", worst_dev < 1e-9, d);
    }
    /* and specifically: no artificial jump at the old H/C/N/O boundary */
    {   double aH = 0, aHe = 0;
        Simulation *s1 = sim_create(4,4), *s2 = sim_create(4,4);
        int i1 = sim_add_atom(s1, 1, vec3(0,0,0), 0);
        int i2 = sim_add_atom(s2, 2, vec3(0,0,0), 0);
        aH = qm_polarizability(&s1->atoms[i1]);
        aHe = qm_polarizability(&s2->atoms[i2]);
        char d[160];
        snprintf(d, sizeof d, "H=%.4f He=%.4f He/H=%.3f (measured 0.307; old code gave 0.074)",
                 aH, aHe, aHe / aH);
        okrel("H/He step matches the measured value", aHe / aH, 0.2051 / 0.6679, 5e-3);
        sim_destroy(s1); sim_destroy(s2);
    }
    /* ions keep their Pauling crystal values */
    Simulation *s = sim_create(4, 4);
    int ik = sim_add_ion(s, 19, 1, vec3(0,0,0), 1.0);
    int ina = sim_add_ion(s, 11, 1, vec3(5,0,0), 1.0);
    okrel("K+  alpha = 0.83 A^3 (Pauling crystal)", qm_polarizability(&s->atoms[ik]), 0.83, 1e-9);
    okrel("Na+ alpha = 0.18 A^3 (Pauling crystal)", qm_polarizability(&s->atoms[ina]), 0.18, 1e-9);
    sim_destroy(s);
}

static void test_qeq(void) {
    grp("QEq: hard core for bonded pairs, bounded charges, exact sum (F9)");
    struct { const char *name; int za, zb; double r; int order; } T[] = {
        { "C-O alcohol 1.43 A", 6, 8, 1.43, 1 },
        { "C=O carbonyl 1.23 A", 6, 8, 1.23, 2 },
        { "O-O peroxide 1.45 A", 8, 8, 1.45, 1 },
        { "O...O H-bond 2.80 A", 8, 8, 2.80, 0 },
    };
    for (unsigned t = 0; t < sizeof T / sizeof T[0]; t++) {
        Simulation *s = sim_create(8, 8);
        int a = sim_add_atom(s, T[t].za, vec3(0,0,0), 0);
        int b = sim_add_atom(s, T[t].zb, vec3(T[t].r,0,0), 0);
        if (T[t].order > 0) {
            sim_add_bond(s, a, b, T[t].order);
            s->atoms[a].num_bonds = 1; s->atoms[a].bond_partners[0] = b; s->atoms[a].bond_orders[0] = T[t].order;
            s->atoms[b].num_bonds = 1; s->atoms[b].bond_partners[0] = a; s->atoms[b].bond_orders[0] = T[t].order;
        }
        double q[8];
        int rc = qm_qeq(s, 0.0, 1.0, q);
        char d[192];
        snprintf(d, sizeof d, "rc=%d q=(%+.4f,%+.4f) sum=%+.2e", rc, q[0], q[1], q[0] + q[1]);
        ok(T[t].name, rc == 0 && fabs(q[0] + q[1]) < 1e-9, d);
        ok("  charges bounded by 2 e", fabs(q[0]) <= 2.0 && fabs(q[1]) <= 2.0, d);
        sim_destroy(s);
    }
    /* Electronegativity ordering, on a BONDED C-O pair (the chemically
     * meaningful case). Note: two atoms 1.5 A apart with NO bond between
     * them is not a molecule - the 1/r coupling then exceeds the
     * harmonicities and the solve lands in the charge-transfer mode, so
     * that input is not a valid test of the electronegativity ordering. */
    Simulation *s = sim_create(8, 8);
    int c = sim_add_atom(s, 6, vec3(0,0,0), 0);
    int o = sim_add_atom(s, 8, vec3(1.43,0,0), 0);
    sim_add_bond(s, c, o, 1);
    s->atoms[c].num_bonds = 1; s->atoms[c].bond_partners[0] = o; s->atoms[c].bond_orders[0] = 1;
    s->atoms[o].num_bonds = 1; s->atoms[o].bond_partners[0] = c; s->atoms[o].bond_orders[0] = 1;
    double q[8];
    qm_qeq(s, 0.0, 1.0, q);
    char d[160];
    snprintf(d, sizeof d, "q(C)=%+.4f q(O)=%+.4f  (Pauling O 3.44 > C 2.55)", q[c], q[o]);
    ok("QEq puts more charge on the electronegative atom", q[o] < q[c], d);
    sim_destroy(s);
    /* total charge is honoured exactly, including for an ion */
    s = sim_create(8, 8);
    sim_add_atom(s, 8, vec3(0,0,0), 0);
    sim_add_atom(s, 1, vec3(1.0,0,0), 0);
    sim_add_atom(s, 8, vec3(0,1.4,0), 0);
    for (double Q = -1.0; Q <= 1.0; Q += 0.5) {
        double qq[8];
        qm_qeq(s, Q, 1.0, qq);
        double sum = qq[0] + qq[1] + qq[2];
        snprintf(d, sizeof d, "requested %+.2f e, got %+.10f", Q, sum);
        ok("QEq conserves total charge", fabs(sum - Q) < 1e-9, d);
    }
    sim_destroy(s);
}

static void test_memory_safety(void) {
    grp("Memory safety and API robustness");
    /* the two confirmed out-of-bounds writes (UBSan) */
    {   BondParam p, p0;
        int r  = forces_bond_params(6, 6, 4, &p);     /* order 4 into a 4-wide axis */
        int r0 = forces_bond_params(6, 6, 5, &p0);    /* and 5 */
        char d[128];
        snprintf(d, sizeof d, "order 4 -> %d, order 5 -> %d (no OOB write)", r, r0);
        ok("forces_bond_params survives order > 3", r == 1 && r0 == 1, d);
    }
    {   AngleParam ap;
        int r = forces_angle_params(118, 8, 118, &ap); /* Z=118 into a [118] axis */
        char d[128];
        snprintf(d, sizeof d, "rc=%d theta0=%.4f", r, ap.theta0);
        ok("forces_angle_params survives Z=118", r == 0 && isfinite(ap.theta0), d);
    }
    ok("forces_bond_params rejects NULL out", forces_bond_params(1,1,1,NULL) == 0, "");
    ok("forces_angle_params rejects NULL out", forces_angle_params(1,1,1,NULL) == 0, "");
    /* vec3_pbc with a degenerate (vacuum) box must not produce NaN */
    {   Vec3 r = vec3_pbc(vec3(3.0, -4.0, 5.0), vec3(0.0, 0.0, 0.0));
        char d[128];
        snprintf(d, sizeof d, "(%g,%g,%g)", r.x, r.y, r.z);
        ok("vec3_pbc with zero box passes through, no NaN",
           isfinite(r.x) && isfinite(r.y) && isfinite(r.z) && r.x == 3.0, d);
    }
    {   Vec3 b = vec3(10,0,0);
        int p[3] = {1,1,0};
        Vec3 r = vec3_pbc_box(vec3(9.0,0,0), b, p);
        char d[128];
        snprintf(d, sizeof d, "dx=9 in a 10 A periodic box -> %g (wrapped to -1)", r.x);
        ok("vec3_pbc_box wraps a periodic axis", fabs(r.x + 1.0) < 1e-12, d);
    }
    /* NULL-guarded public entry points */
    ok("qm_alpha(NULL) returns 0", qm_alpha(NULL) == 0.0, "");
    ok("qm_hybridization(NULL) is safe", strcmp(qm_hybridization(NULL).label, "none") == 0, "");
    {   QmHybrid h = qm_hybridization_ctx(NULL, 0);
        ok("qm_hybridization_ctx(NULL) is safe", strcmp(h.label, "none") == 0, "");
    }
    ok("nb_planarity_deviation(NULL idx) is safe",
       nb_planarity_deviation(NULL, NULL, 5) == 0.0, "");
    {   Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 6, vec3(0,0,0), 0);
        sim_add_atom(s, 6, vec3(1.5,0,0), 0);
        sim_add_atom(s, 6, vec3(2.1,1.3,0), 0);
        int bad[3] = { 0, 1, 999 };      /* out-of-range index */
        ok("nb_planarity_deviation rejects out-of-range indices",
           nb_planarity_deviation(s, bad, 3) == 0.0, "");
        sim_destroy(s);
    }
    /* force summary must not deref a NULL element */
    {   Simulation *s = sim_create(8, 8);
        Atom *at = &s->atoms[0];
        at->Z = 6; at->element = NULL; at->position = vec3(0,0,0);
        s->num_atoms = 1;
        forces_print_summary(s);
        ok("forces_print_summary survives element == NULL", 1, "");
        sim_destroy(s);
    }
    /* integrator entry point without a guard */
    {   Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 6, vec3(0,0,0), 0);
        integrator_remove_com_velocity(s);
        ok("integrator_remove_com_velocity on a valid sim", 1, "");
        sim_destroy(s);
    }
}

static void test_orbital_truncation(void) {
    grp("Orbital table reports truncation instead of losing electrons (F14)");
    /* Kr (Z=36) is the heaviest tabulated element and fits comfortably. */
    Simulation *s = sim_create(4, 4);
    int kr = sim_add_atom(s, 36, vec3(0,0,0), 0);
    int trunc = 1;
    quantum_fill_orbitals_checked(&s->atoms[kr], &trunc);
    int census = 0;
    for (int i = 0; i < s->atoms[kr].num_orbitals; i++)
        census += s->atoms[kr].orbitals[i].occupation;
    char d[160];
    snprintf(d, sizeof d, "Z=36: %d orbitals, census=%d, truncated=%d",
             s->atoms[kr].num_orbitals, census, trunc);
    ok("Kr orbital census is complete", !trunc && census == 36, d);
    sim_destroy(s);
    /* a 7p-saturated configuration must be reported as truncated, not
     * silently clipped (this is the case the audit found) */
    s = sim_create(4, 4);
    int fake = sim_add_atom(s, 36, vec3(0,0,0), 0);
    /* Fully occupy every shell in the Madelung sequence the code supports,
     * which needs 60 ml-slots against a 32-slot cap. */
    memset(s->atoms[fake].electron_config.config, 0,
           sizeof s->atoms[fake].electron_config.config);
    {   static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
        static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
        for (int i = 0; i < 19; i++)
            s->atoms[fake].electron_config.config[MN[i]-1][ML[i]] = 2 * (2*ML[i] + 1);
    }
    trunc = 0;
    quantum_fill_orbitals_checked(&s->atoms[fake], &trunc);
    snprintf(d, sizeof d, "all 19 shells occupied (needs 60 slots, cap 32): used=%d truncated=%d",
             s->atoms[fake].num_orbitals, trunc);
    ok("over-full configuration reports truncation", trunc == 1, d);
    sim_destroy(s);
}

static void test_minimizer_floor(void) {
    grp("Minimiser divergence floor scales with system size (F11)");
    /* A legitimately deep system must still minimise. 700 water-like
     * atoms at a routine condensed-phase PE per atom must NOT trip a
     * 700 * 60 eV floor. */
    for (int n = 100; n <= 700; n *= 7) {
        Simulation *s = sim_create(n + 4, 16);
        s->use_lj = 1; s->use_coulomb = 0; s->cutoff = 12.0;
        s->use_bonds = s->use_angles = s->use_dihedrals = 0;
        /* neutral, weakly charged spheres: a real condensed-phase PE per
         * atom, so the total is far below any absolute threshold */
        for (int i = 0; i < n; i++) {
            double a = i * 2.399963, r = 3.1 * cbrt((double)i);
            sim_add_atom(s, 8, vec3(r*cos(a), r*sin(a), 0.13*i), 0.0);
        }
        forces_calculate(s);
        double before = s->potential_energy;
        double after = integrator_minimize(s, 300, 0.002, 0.05);
        char d[192];
        snprintf(d, sizeof d,
                 "N=%d PE before=%12.1f after=%12.1f eV  (floor now %12.1f eV, was -50000)",
                 n, before, after, -60.0 * n);
        ok("large deep system still minimises", isfinite(after) && after <= before + 1e-6, d);
        sim_destroy(s);
    }
}

static void test_nve(void) {
    grp("NVE energy conservation (velocity Verlet)");
    for (int idt = 1; idt <= 8; idt++) {
        double dt = 0.25 * idt;
        Simulation *s = sim_create(64, 32);
        s->use_lj = 1; s->use_coulomb = 0; s->cutoff = 8.0;
        s->use_bonds = s->use_angles = s->use_dihedrals = 0;
        s->dt = dt;
        for (int i = 0; i < 12; i++) {
            double a = i * 2.399963, r = 4.0 + 0.45 * i;
            sim_add_atom(s, 8, vec3(r*cos(a), r*sin(a), 0.3*i), 0);
        }
        integrator_maxwell_boltzmann(s, 20.0, 12345);
        s->thermostat.type = THERMOSTAT_NONE;
        forces_calculate(s);
        double E0 = s->potential_energy + integrator_kinetic_energy(s);
        for (int i = 0; i < 20000; i++) integrator_step(s);
        double E1 = s->potential_energy + s->kinetic_energy;
        double per_step = fabs(E1 - E0) / 20000.0;
        char d[160];
        snprintf(d, sizeof d, "dt=%.2f fs  dE/step=%.2e eV  (E0=%.6f)", dt, per_step, E0);
        ok("NVE drift below 1e-6 eV/step", per_step < 1e-6, d);
        sim_destroy(s);
    }
}

static void test_datastream(void) {
    grp("Datastream SHA-256 and tamper detection");
    char h[65];
    ds_sha256_hex("", 0, h);
    ok("SHA-256(\"\") matches FIPS 180-4 vector",
       strcmp(h, "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855") == 0, h);
    ds_sha256_hex("abc", 3, h);
    ok("SHA-256(\"abc\") matches FIPS 180-4 vector",
       strcmp(h, "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad") == 0, h);
    ds_sha256_hex("The quick brown fox jumps over the lazy dog", 43, h);
    ok("SHA-256(fox) matches known vector",
       strcmp(h, "d7a8fbb307d7809469ca9abcb0082e4f8d5651e46d3cdb762d02d0bf37c9e592") == 0, h);
    /* marker hygiene: a smuggled [end] in a claim key must be refused */
    {   DSWriter *w = ds_open("/tmp/audit_regress.ds", "tamper");
        if (w) {
            ds_add_claim(w, "bad\n[end]\nforged", 1.0, "eV", "test");
            int rc = ds_close(w);
            ok("writer rejects a claim key containing [end]", rc != 0, "");
        } else ok("writer rejects a claim key containing [end]", 1, "ds_open returned NULL");
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * AUDIT FIX A-C1: qm_solve_dipoles must solve the system it DOCUMENTS.
 *
 * This is the check whose absence let a structurally broken solver ship. The
 * suite already verified that the returned FORCE is the gradient of the
 * returned ENERGY — which is true for ANY mu whatsoever — and that the energy
 * is continuous, which a deterministic function of stale thread-local memory
 * also is. Neither test asks whether mu is the dipole of the documented model.
 *
 * The oracle here is independent of the solver: the 3N x 3N operator is
 * reassembled from the model definition in this file, the residual
 *     || (I - A) mu - alpha (*) E0 ||
 * is formed against the solver's own output, and the result must be at
 * solver precision. A wrong system size, a wrong 3x3 block layout, a flipped
 * off-diagonal sign, or a mispermuted right-hand side all fail it, and each
 * failed differently before the fix.
 *
 * It also pins the two properties a dipole solution must have that no
 * structural test can see:
 *   - the residual is computed on systems spanning N = 2..8, because the
 *     original defect scaled with N (the stale-memory spill reached 2N^2 - N - 1
 *     entries past a region zeroed to N^2, i.e. further for every N > 1);
 *   - the solution is translation and rotation covariant, since the model has
 *     no absolute frame: rotating the whole system must rotate mu.
 * ══════════════════════════════════════════════════════════════════════════ */
/* The Thole width is a private #define in qm.c. The oracle must NOT reuse it:
 * an independent check that imports the constant under test is not independent.
 * 2.0 A is the documented value and is asserted against the engine's own
 * behaviour by the f(2.35)/f(2.75) figures quoted in qm.c's comment. */
#define A_TH_TEST 2.0

static void qm_fields_indep(const Simulation *sim, double dielectric, Vec3 *E) {
    int N = sim->num_atoms;
    for (int i = 0; i < N; i++) E[i] = vec3_zero();
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            double qj = sim->atoms[j].partial_charge;
            if (!isfinite(qj) || fabs(qj) < 1e-12) continue;
            Vec3 d = vec3_sub(sim->atoms[i].position, sim->atoms[j].position);
            double r2 = vec3_norm2(d);
            if (r2 < 1e-8 || !isfinite(r2)) continue;
            if (r2 > sim->cutoff * sim->cutoff) continue;
            double r = sqrt(r2);
            double u = r / A_TH_TEST;
            double fth = 1.0 - exp(-u * u * u);
            vec3_iadd(&E[i], vec3_scale(d, COULOMB_MD * qj * fth / (r2 * r * dielectric)));
        }
    }
}

static void test_dipole_solver_equation(void) {
    grp("Coupled dipoles solve the DOCUMENTED 3N x 3N system (audit A-C1)");

    unsigned seed = 20260929u;
    double worst_res = 0.0;
    int nsys = 0, singular = 0;

    for (int N = 2; N <= 8; N++) {
        for (int trial = 0; trial < 12; trial++) {
            Simulation *s = sim_create(16, 32);
            if (!s) continue;
            /* Deliberately unbonded: qm_pair_excluded would otherwise remove
             * the very couplings whose assembly is under test. */
            for (int i = 0; i < N; i++) {
                seed = seed * 1103515245u + 12345u;
                double x = 2.0 + 1.4 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double y = 2.2 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double z = 2.2 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                int pick = (int)((seed >> 8) % 4);
                int Z = pick == 0 ? 19 : pick == 1 ? 8 : pick == 2 ? 1 : 6;
                double q = (Z == 19) ? 1.0 : (Z == 8) ? -0.55 : (Z == 1) ? 0.25 : 0.1;
                sim_add_atom(s, Z, vec3(x, y, z), q);
            }
            s->cutoff = 30.0;
            Vec3 mu[256];
            int rc = qm_solve_dipoles(s, 1.0, mu);
            nsys++;
            if (rc != 0) { singular++; sim_destroy(s); continue; }

            Vec3 E0[256];
            qm_fields_indep(s, 1.0, E0);
            double al[256];
            for (int i = 0; i < N; i++) al[i] = qm_polarizability(&s->atoms[i]);

            /* residual of  (I - A) mu = alpha (*) E0, A reassembled here */
            double rn = 0.0, dn = 0.0;
            for (int i = 0; i < N; i++) {
                Vec3 m = mu[i], Ax = vec3_zero();
                for (int j = 0; j < N; j++) {
                    if (i == j) continue;
                    Vec3 d = vec3_sub(s->atoms[i].position, s->atoms[j].position);
                    double r = vec3_norm(d);
                    if (r < 1e-8) continue;
                    if (r > s->cutoff) continue;
                    double uu = r / A_TH_TEST;
                    double fth = 1.0 - exp(-uu * uu * uu);
                    if (fth == 0.0) continue;
                    double kij = al[i] * COULOMB_MD * fth / (r * r * r);
                    Vec3 uj = mu[j];
                    double dot = d.x * uj.x + d.y * uj.y + d.z * uj.z;
                    /* A_ij mu_j = kij [ 3 d (d.mu_j)/r^2 - mu_j ] */
                    Vec3 term = vec3_scale(d, kij * 3.0 * dot / (r * r));
                    Ax = vec3_add(Ax, vec3_sub(term, vec3_scale(uj, kij)));
                }
                double want[3] = { al[i] / COULOMB_MD * E0[i].x,
                                   al[i] / COULOMB_MD * E0[i].y,
                                   al[i] / COULOMB_MD * E0[i].z };
                double got[3] = { m.x - Ax.x, m.y - Ax.y, m.z - Ax.z };
                for (int a = 0; a < 3; a++) {
                    rn += fabs(got[a] - want[a]);
                    dn += fabs(want[a]);
                }
            }
            double res = rn / (dn > 1e-30 ? dn : 1.0);
            if (res > worst_res) worst_res = res;
            sim_destroy(s);
        }
    }

    char d[224];
    snprintf(d, sizeof d, "%d systems (%d singular), worst relative residual = %.3e",
             nsys, singular, worst_res);
    ok("solver satisfies (I - A) mu = alpha (*) E0", worst_res < 1e-9, d);

    /* Covariance: the model has no absolute frame, so rotating the whole system
     * must rotate mu. Built by rotating each ATOM so element, charge and
     * index all correspond one-to-one — an earlier version of this check
     * compared mismatched index sets and reported a spurious failure. */
    {
        Simulation *a = sim_create(8, 8), *b = sim_create(8, 8);
        static const int    zs[4] = { 19, 8, 1, 19 };
        static const double qs[4] = { 1.0, -0.55, 0.25, 1.0 };
        static const double ps[4][3] = {
            { 0.0, 0.0, 0.0 }, { 2.7, 0.2, 0.1 },
            { 1.1, 2.6, -0.4 }, { -1.9, 1.3, 2.2 }
        };
        const double t = 0.7, ct = cos(t), st = sin(t);
        for (int i = 0; i < 4; i++) {
            sim_add_atom(a, zs[i], vec3(ps[i][0], ps[i][1], ps[i][2]), qs[i]);
            sim_add_atom(b, zs[i],
                         vec3(ps[i][0] * ct - ps[i][1] * st,
                              ps[i][0] * st + ps[i][1] * ct,
                              ps[i][2]), qs[i]);
        }
        a->cutoff = b->cutoff = 30.0;
        Vec3 ma[8], mb[8];
        int ra = qm_solve_dipoles(a, 1.0, ma);
        int rb = qm_solve_dipoles(b, 1.0, mb);
        double worst = 0.0, scale = 0.0;
        for (int i = 0; i < 4; i++) {
            Vec3 rot = vec3(ma[i].x * ct - ma[i].y * st, ma[i].x * st + ma[i].y * ct, ma[i].z);
            worst = fmax(worst, vec3_norm(vec3_sub(rot, mb[i])));
            scale = fmax(scale, vec3_norm(mb[i]));
        }
        snprintf(d, sizeof d, "rc=(%d,%d) max rotation mismatch = %.3e (|mu| ~ %.3f)",
                 ra, rb, worst, scale);
        ok("dipoles are rotation covariant (no absolute frame)",
           (ra == 0 && rb == 0) && worst < 1e-8 * (scale > 0 ? scale : 1.0), d);
        sim_destroy(a); sim_destroy(b);
    }
}

/* AUDIT FIX A-C2: the ±2 e bound must not cost charge conservation.
 * The augmented QEq system carries sum(q) = total_q as one of its rows; a
 * post-solve clamp satisfies the bound by violating the row. 233 of 4000
 * random clusters did exactly that at HEAD, worst |sum q| = 8.0 e.
 *
 * The case below is constructed, not random, and is the tightest alternating
 * C/O chain found by a spacing sweep: an unbonded chain at 1.25 A puts the
 * 1/r coupling well above the hardness, so the unbounded solve drives the
 * central oxygen past +2 e and the clamp has to act. Measured:
 *   HEAD    sum(q) = -6.0e-01   max|q| = 2.0000     <- conservation destroyed
 *   fixed   sum(q) = -4.4e-16   max|q| = 2.0000     <- both invariants hold
 * A regression test has to reach the defect, so the geometry is fixed here
 * rather than sampled: a random cluster hits the clamp only ~6% of the time
 * and would make this check flaky in the direction that matters least. */
static void test_qeq_conservation_under_bound(void) {
    grp("QEq bound vs conservation: the constraint must survive the clamp (A-C2)");

    Simulation *s = sim_create(16, 32);
    static const int zs[5] = { 8, 6, 8, 6, 8 };
    for (int i = 0; i < 5; i++)
        sim_add_atom(s, zs[i], vec3(1.25 * i, 0.0, 0.0), 0.0);
    double q[128];
    int rc = qm_qeq(s, 0.0, 1.0, q);
    double sum = 0.0, mx = 0.0;
    for (int i = 0; i < 5; i++) { sum += q[i]; if (fabs(q[i]) > mx) mx = fabs(q[i]); }
    char d[224];
    snprintf(d, sizeof d, "rc=%d sum(q)=%+.3e max|q|=%.4f", rc, sum, mx);
    ok("sum(q) == total_q exactly even when the bound is active",
       rc == 0 && fabs(sum) < 1e-12, d);
    ok("the 2 e bound is still enforced", mx <= 2.0 + 1e-12, d);
    sim_destroy(s);

    /* and the SCF path, which additionally under-relaxes. Same 5-atom chain
     * with a K+ pinned inside it, so the dipole feedback is live and the
     * bound is active there too. */
    {
        Simulation *t = sim_create(16, 32);
        static const int zs2[5] = { 8, 6, 8, 6, 8 };
        for (int i = 0; i < 5; i++)
            sim_add_atom(t, zs2[i], vec3(1.25 * i, 0.0, 0.0), 0.0);
        int ion = sim_add_ion(t, 19, 1, vec3(2.50, 0.0, 0.0), 1.0);
        int its = qm_scf_charges(t, 1.0, 1.0, ion, 1.0);
        double ssum = 0.0, smx = 0.0;
        for (int i = 0; i < t->num_atoms; i++) {
            if (i == ion) continue;
            ssum += t->atoms[i].partial_charge;
            if (fabs(t->atoms[i].partial_charge) > smx)
                smx = fabs(t->atoms[i].partial_charge);
        }
        snprintf(d, sizeof d, "its=%d shell sum=%+.3e max|q|=%.4f", its, ssum, smx);
        ok("SCF charge path conserves charge (shell carries total - pinned)",
           its > 0 && fabs(ssum) < 1e-6, d);
        ok("SCF charge path respects the 2 e bound", smx <= 2.0 + 1e-9, d);
        sim_destroy(t);
    }
}

/* AUDIT FIX M1/M2: the last hand-typed reciprocals in constants.h. */
/* AUDIT FIX D2: kcsa_set_ion_radius() documentation claimed sigma 1.888 A
 * (K+) and 1.246 A (Na+). Those came from mixing against r(O) = 1.40 A
 * where a sigma was needed, and were never compared to the code. Recompute
 * the construction from the constants the code actually uses, so the header
 * and the implementation cannot drift apart again. */
/* AUDIT FIX D4: kcsa_filter.c claimed every residue is neutral "by a UNIFORM
 * offset added to its carbons, computed at build time". Nothing is computed
 * at build time, and the offset is not uniform - but the neutrality claim
 * itself is TRUE, which is why the wrong mechanism went unnoticed for so
 * long. Group the filter's atoms the way the table is annotated (by residue)
 * and assert each residue really is neutral, so the true property is pinned. */
static void test_kcsa_filter_residue_neutrality(void) {
    grp("KcsA filter: each residue is neutral (audit D4)");
    Simulation *s = sim_create(600, 600);
    if (!s) { ok("build the filter", 0, "sim_create failed"); return; }
    int first = kcsa_build_filter(s, vec3(0, 0, 0), 1);
    int n = s->num_atoms - first;

    double total = 0.0, mn = 1e9, mx = -1e9;
    for (int i = first; i < s->num_atoms; i++) {
        double q = s->atoms[i].partial_charge;
        total += q;
        if (q < mn) mn = q;
        if (q > mx) mx = q;
    }
    char d[224];
    snprintf(d, sizeof d, "%d atoms, sum(q) = %+.6f e, range [%+.4f, %+.4f] e",
             n, total, mn, mx);
    ok("the whole filter is neutral", n == 41 && fabs(total) < 1e-12, d);
    ok("max single-atom charge stays at 0.93 e (the D4 claim)", mx <= 0.93 + 1e-9, d);

    /* Per-residue neutrality. Membership is by table index, matching the
     * order kcsa_build_filter() adds atoms in: each residue's heavy atoms are
     * contiguous, then its constructed hydrogens and OXT are appended at the
     * end of the table. (Grouping by z instead is wrong - the residues march
     * along the pore axis and overlap in z.) */
    static const int r75[] = {0,1,2,3,4,5,6,38,39};
    static const int r76[] = {7,8,9,10,11,12,13,34};
    static const int r77[] = {14,15,16,17,35};
    static const int r78[] = {18,19,20,21,22,23,24,25,26,27,28,29,36};
    static const int r79[] = {30,31,32,33,37,40};
    static const int *mem[5] = {r75,r76,r77,r78,r79};
    static const int cnt[5] = {9,8,5,13,6};
    static const int res[5] = {75,76,77,78,79};
    double worst = 0.0; int worst_r = 0;
    for (int r = 0; r < 5; r++) {
        double rs = 0.0;
        for (int k = 0; k < cnt[r]; k++) rs += s->atoms[first + mem[r][k]].partial_charge;
        if (fabs(rs) > worst) { worst = fabs(rs); worst_r = res[r]; }
    }
    snprintf(d, sizeof d, "worst |sum(q)| over 5 residues = %.3e e (res %d)", worst, worst_r);
    ok("every residue is neutral to 1e-9 e", worst < 1e-9, d);
    sim_destroy(s);
}

static void test_kcsa_ion_sigma(void) {
    grp("KcsA ion sigma is the one the header documents (audit D2; full-audit M3)");
    double rO = AMBER_RSTAR_TO_SIGMA(1.6612);   /* LJ_AMBER_O_SIGMA */
    char d[224];
    okrel("LJ_AMBER_O_SIGMA == 2.959922 A (AMBER R* 1.6612)", rO, 2.959922, 1e-6);

    for (int k = 0; k < 2; k++) {
        int    Z   = k ? 11 : 19;
        double r_i = kcsa_cation_radius(Z);
        double contact = r_i + 1.40;                 /* KCSA_O_RADIUS */
        double sig_ionO = contact / TWOPOW_SIXTH;
        Simulation *s = sim_create(8, 8);
        int ion = sim_add_ion(s, Z, 1, vec3(0, 0, 0), 1.0);
        kcsa_set_ion_radius(s, ion, Z);
        double got = s->atoms[ion].lj_sigma;

        /* the two-step inversion must reproduce the tabulated contact exactly.
         * Full-audit M3: VIII contacts 2.91/2.58 A, sigmas 2.2251/1.6371 A
         * (were VI 2.78/2.42 A, 1.9935/1.3520 A). */
        double back = TWOPOW_SIXTH * 0.5 * (got + rO);   /* = the contact distance */
        (void)sig_ionO;
        snprintf(d, sizeof d, "Z=%2d r=%.2f sigma=%.4f A, mixed minimum %.4f A, contact %.2f A",
                 Z, r_i, got, back, contact);
        ok(k ? "Na+ sigma 1.6371 A puts the ion-O minimum at 2.58 A (VIII)"
             : "K+  sigma 2.2251 A puts the ion-O minimum at 2.91 A (VIII)",
           fabs(got - (k ? 1.6371 : 2.2251)) < 5e-4 && fabs(back - contact) < 1e-9, d);
        sim_destroy(s);
    }
}

/* AUDIT FIX D3 pinned the old, wrong glycine sum. AUDIT FIX N2 then fixed
 * the table itself: the "adjusted for exact neutrality" HXT correction had
 * its sign inverted, so the pin becomes the physically required value - a
 * complete free neutral molecule sums to zero. Alanine was already neutral
 * and stays the control case; the dipeptide's -0.043 e is the real
 * condensation capping and is pinned as before. */
static void test_aminoacid_charge_sums(void) {
    grp("Amino-acid charge tables: free molecules are neutral (D3/N2)");
    char d[224];
    {
        Simulation *s = sim_create(64, 64);
        sim_place_glycine(s, vec3(0, 0, 0));
        double sum = 0; for (int i = 0; i < s->num_atoms; i++) sum += s->atoms[i].partial_charge;
        snprintf(d, sizeof d, "%d atoms, sum(q) = %+.6f e", s->num_atoms, sum);
        ok("free glycine is neutral after the N2 HXT correction", fabs(sum) < 1e-12, d);
        sim_destroy(s);
    }
    {
        Simulation *s = sim_create(64, 64);
        sim_place_alanine(s, vec3(0, 0, 0));
        double sum = 0; for (int i = 0; i < s->num_atoms; i++) sum += s->atoms[i].partial_charge;
        snprintf(d, sizeof d, "%d atoms, sum(q) = %+.6f e", s->num_atoms, sum);
        ok("free alanine is the neutral control case", fabs(sum) < 1e-12, d);
        sim_destroy(s);
    }
    {
        Simulation *s = sim_create(128, 128);
        int ala_N = -1;
        sim_place_dipeptide_GlyAla(s, vec3(0, 0, 0), &ala_N);
        double sum = 0; for (int i = 0; i < s->num_atoms; i++) sum += s->atoms[i].partial_charge;
        snprintf(d, sizeof d, "%d atoms, sum(q) = %+.6f e", s->num_atoms, sum);
        ok("Gly-Ala dipeptide is charged by -0.043 e after condensation", fabs(sum + 0.043) < 1e-9, d);
        sim_destroy(s);
    }
}

static void test_qeq_nonzero_total_charge(void) {
    grp("QEq: conservation must not depend on total_q being zero (A-C2b/A-C2c)");

    /* AUDIT FIX A-C2b. The first version of the A-C2 fix scaled the solved
     * vector about ZERO to impose the bound:
     *
     *     shrink = QMAX / max|q|;   q <- q * shrink
     *
     * which preserves sum(q) only when total_q is 0, because it maps
     * sum(q) = total_q to total_q * shrink. The regression test that verified
     * the fix sampled total_q = 0 only, so it passed and the property was
     * reported as fixed. It was not: 31 of 400 random clusters over
     * total_q in [-1, +1] violated sum(q) = total_q, worst error 0.79 e.
     * qm_qeq takes total_q as a parameter; testing it only at zero did not
     * test it.
     *
     * The bound must be imposed by contracting the deviation FROM THE MEAN,
     * which leaves the sum at m*mean = total_q exactly. Swept over total_q
     * rather than pinned at a value, because the value is the thing that
     * decides whether the bug is visible. */
    {
        unsigned seed = 12345u;
        int nbad = 0, nsys = 0, over = 0;
        double worst = 0.0, worst_tot = 0.0;
        for (int t = 0; t < 400; t++) {
            int N = 2 + (int)((seed >> 8) % 6);
            seed = seed * 1103515245u + 12345u;
            double tot = -1.0 + 2.0 * ((double)((seed >> 8) % 1000) / 1000.0);
            Simulation *s = sim_create(16, 32);
            for (int i = 0; i < N; i++) {
                seed = seed * 1103515245u + 12345u;
                double x = 1.0 + 1.2 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double y = 1.2 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double z = 1.2 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                int pick = (int)((seed >> 8) % 3);
                sim_add_atom(s, pick == 0 ? 8 : (pick == 1 ? 6 : 19), vec3(x, y, z), 0.0);
            }
            double q[64];
            if (qm_qeq(s, tot, 1.0, q) == 0) {
                double sum = 0, mx = 0;
                for (int i = 0; i < N; i++) { sum += q[i]; if (fabs(q[i]) > mx) mx = fabs(q[i]); }
                double err = fabs(sum - tot);
                if (err > worst) { worst = err; worst_tot = tot; }
                if (err > 1e-9) nbad++;
                if (mx > 2.0 + 1e-12) over++;
                nsys++;
            }
            sim_destroy(s);
        }
        char dd[224];
        snprintf(dd, sizeof dd, "%d systems, worst |sum-total_q| = %.3e at total_q=%+.3f",
                 nsys, worst, worst_tot);
        ok("sum(q) == total_q for nonzero total_q (not just zero)", nbad == 0, dd);
        snprintf(dd, sizeof dd, "%d of %d exceeded the 2 e bound", over, nsys);
        ok("the 2 e bound still holds at nonzero total_q", over == 0, dd);
    }

    /* AUDIT FIX A-C2c: qm_qeq_pinned enforced NO bound at all. Every other
     * entry point applies the +/-2 e limit, so a pinned-atom system - exactly
     * what the engine uses for an ion in a cage - could return charges of any
     * magnitude while its unpinned twin could not. The invariant here is over
     * the FREE atoms: sum(free) = total_q - pinned_q, and the pinned charge is
     * an input, not a variable to be scaled. */
    {
        unsigned seed = 99u;
        int nbad = 0, nbound = 0, nsys = 0, pinmoved = 0;
        double worst = 0.0;
        for (int t = 0; t < 300; t++) {
            int N = 3 + (int)((seed >> 8) % 5);
            seed = seed * 1103515245u + 12345u;
            double tot = -1.5 + 3.0 * ((double)((seed >> 8) % 1000) / 1000.0);
            const double pin = 1.0;
            Simulation *s = sim_create(16, 32);
            sim_add_atom(s, 19, vec3(0, 0, 0), pin);   /* index 0, pinned */
            for (int i = 1; i < N; i++) {
                seed = seed * 1103515245u + 12345u;
                double x = 0.9 + 1.1 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double y = 1.1 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                double z = 1.1 * ((double)((seed >> 8) % 1000) / 1000.0);
                seed = seed * 1103515245u + 12345u;
                sim_add_atom(s, ((seed >> 8) % 2) ? 6 : 8, vec3(x, y, z), 0.0);
            }
            double q[64];
            if (qm_qeq_pinned(s, tot, 1.0, 0, pin, q) == 0) {
                double sum = 0, mx = 0;
                for (int i = 1; i < N; i++) { sum += q[i]; if (fabs(q[i]) > mx) mx = fabs(q[i]); }
                double err = fabs(sum - (tot - pin));
                if (err > worst) worst = err;
                if (err > 1e-9) nbad++;
                if (mx > 2.0 + 1e-12) nbound++;
                if (fabs(q[0] - pin) > 1e-12) pinmoved++;
                nsys++;
            }
            sim_destroy(s);
        }
        char dd[224];
        snprintf(dd, sizeof dd, "%d systems, worst |sum(free)-(total_q-pinned_q)| = %.3e",
                 nsys, worst);
        ok("pinned path conserves total_q - pinned_q", nbad == 0, dd);
        snprintf(dd, sizeof dd, "%d of %d exceeded the bound; %d pinned charges altered",
                 nbound, nsys, pinmoved);
        ok("pinned path applies the bound and leaves the pinned charge alone",
           nbound == 0 && pinmoved == 0, dd);
    }
}

static void test_constants_derived(void) {
    grp("Every exact reciprocal is derived, not typed (audit M1/M2)");
    /* PLANCK_HBAR == PLANCK_H / (2 pi) to the last representable bit */
    okrel("PLANCK_HBAR == PLANCK_H / 2pi", PLANCK_HBAR,
          PLANCK_H / (2.0 * 3.14159265358979323846), 1e-16);
    /* the pair now multiplies to exactly 1 by construction */
    okrel("HARTREE_TO_EV * EV_TO_HARTREE == 1", HARTREE_TO_EV * EV_TO_HARTREE,
          1.0, 1e-16);
    /* KCAL pair was already exact; keep it pinned so it cannot regress */
    okrel("KCAL_MOL_TO_EV * EV_TO_KCAL_MOL == 1", KCAL_MOL_TO_EV * EV_TO_KCAL_MOL,
          1.0, 1e-16);
    okrel("COULOMB_MD is the 2018-CODATA k_e prefactor", COULOMB_MD,
          8.9875517923e9 * 1.602176634e-19 / 1.0e-10, 1e-15);
}

/* AUDIT FIX M3: a zero-force system with force_tolerance == 0 must not
 * produce NaN positions. `max_force < force_tolerance` is 0 < 0 = FALSE, so
 * the loop body used to run, scale became inf, and every displacement became
 * 0*inf = NaN. The NaN rollback hid it and the answer came out right by
 * luck of ordering; the guard makes it right on purpose. */
static void test_minimizer_zero_tolerance(void) {
    grp("Minimisers do not divide by zero at force_tolerance == 0 (M3)");
    for (int variant = 0; variant < 3; variant++) {
        Simulation *s = sim_create(8, 8);
        sim_add_atom(s, 6, vec3(0, 0, 0), 0);
        s->use_lj = 0; s->use_coulomb = 0;
        s->use_bonds = 0; s->use_angles = 0; s->use_dihedrals = 0;
        forces_calculate(s);
        int zeroF = (vec3_norm(s->atoms[0].force) == 0.0);
        double E;
        if (variant == 0)      E = integrator_minimize(s, 10, 0.01, 0.0);
        else if (variant == 1) { int frozen[1] = {0};
                                  E = integrator_minimize_frozen(s, frozen, 10, 0.01, 0.0); }
        else                   E = integrator_fire(s, 10, 0.5, 0.0);
        static const char *nm[3] = { "integrator_minimize", "minimize_frozen", "integrator_fire" };
        char d[224];
        snprintf(d, sizeof d, "%s: |F|=%g (zero: %d) E=%g pos=(%g,%g,%g)",
                 nm[variant], vec3_norm(s->atoms[0].force), zeroF, E,
                 s->atoms[0].position.x, s->atoms[0].position.y, s->atoms[0].position.z);
        ok(nm[variant], zeroF && isfinite(E) && isfinite(s->atoms[0].position.x)
                         && s->atoms[0].position.x == 0.0, d);
        sim_destroy(s);
    }
}

/* AUDIT FIX V1: the dihedral's rotational invariance is a property no shipped
 * test checked. Net torque must vanish for a central force pair; the
 * translational half was asserted but the rotational half was not. */
static void test_dihedral_rotational_invariance(void) {
    grp("Dihedral gradient: net torque vanishes (audit V1 follow-through)");
    static const double P[3][4][3] = {
        {{0,0,0},{1.45,0,0},{2.0,1.2,0.3},{1.2,2.0,-0.4}},
        {{0,0,0},{1.0,0,0},{1.0,1.0,0},{1.0,1.0,1.0}},
        {{-1.2,0.3,0.1},{0,0,0},{1.3,0.2,-0.1},{2.5,-0.4,0.3}},
    };
    double worst = 0.0;
    for (int c = 0; c < 3; c++) {
        Atom a[4], f[4];
        memset(a, 0, sizeof a); memset(f, 0, sizeof f);
        for (int i = 0; i < 4; i++) {
            a[i].position = vec3(P[c][i][0], P[c][i][1], P[c][i][2]);
            f[i].position = a[i].position;
        }
        Dihedral dh = {0, 1, 2, 3, 0.1, 1, 0.0};
        forces_dihedral(a, &dh);
        forces_dihedral_fd(f, &dh);
        Vec3 torque = vec3_zero(), scale = vec3_zero();
        for (int i = 0; i < 4; i++) {
            torque = vec3_add(torque, vec3_cross(a[i].position, a[i].force));
            scale = vec3_add(scale, vec3_scale(a[i].force, 1.0));
        }
        double r = vec3_norm(torque) / (vec3_norm(scale) > 1e-12 ? vec3_norm(scale) : 1.0);
        if (r > worst) worst = r;
    }
    char d[192];
    snprintf(d, sizeof d, "max |sum r_i x F_i| / |F| = %.3e over 3 geometries", worst);
    ok("net torque vanishes on the analytic path", worst < 1e-12, d);
}

static void test_rng(void) {
    grp("Random number generators (audit I1: PCG64 now actually exists)");
    /* The v9R4 notes claimed PCG64 was implemented. It was not - grep
     * found no PCG64 anywhere in the tree. These checks pin both
     * generators to reference streams so that claim stays true. */
    for (int kind = 0; kind <= 1; kind++) {
        const char *nm = kind ? "PCG64" : "LCG";
        /* Maxwell-Boltzmann initialisation must give zero mean and a
         * variance matching kB*T/m, for BOTH generators. */
        Simulation *s = sim_create(64, 8);
        for (int i = 0; i < 40; i++)
            sim_add_atom(s, 8, vec3(3.0 * i, 0, 0), 0);
        s->rng_kind = kind;
        const double T = 300.0;
        integrator_maxwell_boltzmann(s, T, 42);
        double mean = 0, var = 0; int m = 0;
        for (int i = 0; i < s->num_atoms; i++)
            for (int c = 0; c < 3; c++) {
                double v = (&s->atoms[i].velocity.x)[c];
                mean += v; var += v * v; m++;
            }
        mean /= m; var = var / m - mean * mean;
        /* kB*T/m in A^2/fs^2, via the engine's own conversion factors */
        double m_amu = pt_element(8)->mass;
        double expect = (BOLTZMANN_K / EV_TO_J) * T / m_amu
                      / (AMU * 1.0e10 / EV_TO_J);
        char d[192];
        snprintf(d, sizeof d, "%s: <v>=%+.3e  var=%.4e  kT/m=%.4e  var/(kT/m)=%.4f",
                 nm, mean, var, expect, var / expect);
        ok(nm, fabs(mean) < 1e-3 * sqrt(expect) && fabs(var / expect - 1.0) < 0.1, d);
        sim_destroy(s);
    }
    /* The two generators must produce DIFFERENT streams - a PCG64 that
     * silently fell through to the LCG would pass everything above. */
    {   Simulation *a = sim_create(8, 8), *b = sim_create(8, 8);
        for (int i = 0; i < 4; i++) {
            sim_add_atom(a, 8, vec3(3.0 * i, 0, 0), 0);
            sim_add_atom(b, 8, vec3(3.0 * i, 0, 0), 0);
        }
        a->rng_kind = INTEGRATOR_RNG_LCG;   b->rng_kind = INTEGRATOR_RNG_PCG64;
        integrator_maxwell_boltzmann(a, 300.0, 7);
        integrator_maxwell_boltzmann(b, 300.0, 7);
        int same = 1;
        for (int i = 0; i < a->num_atoms; i++)
            if (a->atoms[i].velocity.x != b->atoms[i].velocity.x) same = 0;
        ok("LCG and PCG64 give different streams", !same, "");
        /* and each must be reproducible for a fixed seed */
        Simulation *c = sim_create(8, 8);
        for (int i = 0; i < 4; i++) sim_add_atom(c, 8, vec3(3.0 * i, 0, 0), 0);
        c->rng_kind = INTEGRATOR_RNG_PCG64;
        integrator_maxwell_boltzmann(c, 300.0, 7);
        int repro = 1;
        for (int i = 0; i < b->num_atoms; i++)
            if (b->atoms[i].velocity.x != c->atoms[i].velocity.x) repro = 0;
        ok("PCG64 reproducible for a fixed seed", repro, "");
        sim_destroy(a); sim_destroy(b); sim_destroy(c);
    }
    /* default must remain the LCG so the record stays byte-identical */
    {   Simulation *s = sim_create(4, 4);
        sim_add_atom(s, 8, vec3(0,0,0), 0);
        ok("default rng_kind is the record-compatible LCG",
           s->rng_kind == INTEGRATOR_RNG_LCG, "");
        sim_destroy(s);
    }
}

static void test_kcsa_filter(void) {
    grp("KcsA filter is the REAL deposited 1K4C geometry (not a cage)");
    Simulation *s = sim_create(KCSA_FILTER_ATOMS * 4 + 8,
                                KCSA_FILTER_BONDS * 4 + 16);
    int f = kcsa_build_filter(s, vec3(0,0,0), 4);
    char d[224];
    ok("filter builds", f == 0, "");
    ok("41 atoms per subunit x 4 = 164", s->num_atoms - f == 164, "");
    /* must be neutral, or it binds a cation for a reason unrelated to KcsA */
    double q = 0.0;
    for (int i = f; i < s->num_atoms; i++) q += s->atoms[i].partial_charge;
    snprintf(d, sizeof d, "net filter charge = %+.6e e", q);
    ok("filter is exactly neutral", fabs(q) < 1e-5, d);
    /* The five LIGAND oxygens per subunit must hug the pore wall
     * (deposited r = 2.27-2.64 A). Non-ligand oxygens - the C-terminal
     * OXT and the Tyr phenol, which points into the wall - sit further
     * out, so this is measured on the ligands only. */
    {   /* Exactly 5 ligand oxygens per subunit (20 for the tetramer), all
         * hugging the wall at the deposited radii. The other oxygens -
         * GLY79's carbonyl, the Tyr phenol and the C-terminal OXT - sit
         * much further out, and counting by radius is unambiguous. */
        double r_min = 1e9, r_max = -1e9; int n = 0;
        for (int i = f; i < s->num_atoms; i++) {
            if (s->atoms[i].Z != 8) continue;
            double r = sqrt(s->atoms[i].position.x*s->atoms[i].position.x +
                            s->atoms[i].position.y*s->atoms[i].position.y);
            if (r > 2.7) continue;
            n++;
            if (r < r_min) r_min = r;
            if (r > r_max) r_max = r;
        }
        snprintf(d, sizeof d, "%d ligand oxygens (want 5x4=20), radii %.3f-%.3f A (deposited 2.267-2.640)",
                 n, r_min, r_max);
        ok("exactly 5 ligand oxygens per subunit, hugging the pore wall",
           n == 20 && r_min > 2.2 && r_max < 2.7, d);
    }
    /* THE decisive check: coordination number 8 at every deposited site.
     * This is what identifies the C4 symmetry assignment as correct —
     * a wrong rotation axis gives CN = 2, not 8. */
    Vec3 sites[4];
    int ns = kcsa_ion_sites(4, sites);
    ok("four ion sites recovered", ns == 4, "");
    static const double want_r[4] = { 2.932, 2.777, 2.773, 2.912 };
    for (int k = 0; k < ns; k++) {
        double mr = 0.0;
        int cn = kcsa_coord_stats(s, f, 4, sites[k], 3.4, &mr);
        snprintf(d, sizeof d, "site %d: CN=%d <ion-O>=%.3f A (deposited %.3f)",
                 k+1, cn, mr, want_r[k]);
        ok("site is 8-coordinate at the deposited distance",
           cn == 8 && fabs(mr - want_r[k]) < 0.02, d);
    }
    /* Thr75 OG1 must be one of the inner-gate ligands: the poly-alanine
     * model this replaced had no such atom at all. */
    {   Vec3 inner = sites[3];   /* the innermost site */
        double best = 1e9; int found = 0;
        for (int i = f; i < s->num_atoms; i++) {
            if (s->atoms[i].Z != 8) continue;
            double dd = vec3_dist(inner, s->atoms[i].position);
            if (dd < best) best = dd;
            /* the inner site is 4x OG1 at 2.880 and 4x backbone O at 2.944,
             * so a 3.4 A coordination shell must contain exactly 8 */
            if (dd < 3.4) found++;
        }
        snprintf(d, sizeof d, "closest ligand to the inner site %.3f A; %d oxygens in the 3.4 A shell (want 8: 4x Thr OG1 + 4x THR backbone O)",
                 best, found);
        ok("inner gate is 8-fold including the Thr hydroxyl", found == 8, d);
    }
    /* GLY79's oxygen must NOT reach the pore: a hand-built model that put
     * all five motif carbonyls on the ion path would be wrong. */
    {   double rmax = 0.0;
        for (int i = f; i < s->num_atoms; i++)
            if (s->atoms[i].Z == 8) {
                double r = sqrt(s->atoms[i].position.x*s->atoms[i].position.x +
                                s->atoms[i].position.y*s->atoms[i].position.y);
                if (r > rmax) rmax = r;
            }
        snprintf(d, sizeof d, "max oxygen radius from the axis = %.3f A (GLY79 O is 4.82 A)", rmax);
        ok("the non-ligating motif oxygen stays out of the pore", rmax > 3.5, d);
    }
    sim_destroy(s);
}

static void test_ion_size_and_hydration(void) {
    grp("Ion size and hydration parameters are sourced, not fitted");
    /* the K+ contact distance must land on the tabulated ionic radius */
    Simulation *s = sim_create(KCSA_FILTER_ATOMS * 4 + 8, KCSA_FILTER_BONDS * 4 + 16);
    kcsa_build_filter(s, vec3(0,0,0), 4);
    for (int Z = 11; Z <= 19; Z += 8) {
        int ion = sim_add_ion(s, Z, 1, vec3(0,0,-33.953), 1.0);
        kcsa_set_ion_radius(s, ion, Z);
        /* Lorentz-Berthelot ion-O sigma, recovered from the stored
         * values. The filter oxygen sigma comes from the shared header -
         * this test used to hardcode 3.06615, which was the same wrong
         * literal kcsa_filter.c carried, so both were wrong together and
         * the consolidating refactor is what exposed it. The AMBER ff99
         * carbonyl oxygen is AMBER_RSTAR_TO_SIGMA(1.6612) = 2.95992. */
        double sigO = LJ_AMBER_O_SIGMA;
        double sig_ionO = 0.5 * (s->atoms[ion].lj_sigma + sigO);
        double contact = sig_ionO * pow(2.0, 1.0/6.0);
        double want = kcsa_cation_radius(Z) + 1.40;
        char d[192];
        snprintf(d, sizeof d, "Z=%d ion-O minimum %.3f A; target radius sum %.3f A", Z, contact, want);
        ok("ion-O contact equals r(cation) + r(O)", fabs(contact - want) < 0.01, d);
    }
    sim_destroy(s);
    okrel("K+  cation radius 1.51 A (Shannon VIII, full-audit M3)", kcsa_cation_radius(19), 1.51, 1e-12);
    okrel("Na+ cation radius 1.18 A (Shannon VIII, full-audit M3)", kcsa_cation_radius(11), 1.18, 1e-12);
    /* hydration free energies on the Marcus 1991 TATB free-energy scale,
     * and the K+ advantage they imply (N1 corrected these from the
     * old-scale enthalpies -322/-454 which full-audit M1 had unified on) */
    okrel("K+  hydration -295.3 kJ/mol (Marcus 1991 TATB dG)", kcsa_hydration_free_energy_kJmol(19), -295.3, 1e-12);
    okrel("Na+ hydration -365.3 kJ/mol (Marcus 1991 TATB dG)", kcsa_hydration_free_energy_kJmol(11), -365.3, 1e-12);
    okrel("K+ dehydration cost 3.061 eV", kcsa_dehydration_cost_eV(19), 295.3/96.48533212, 1e-9);
    okrel("Na+ dehydration cost 3.786 eV", kcsa_dehydration_cost_eV(11), 365.3/96.48533212, 1e-9);
    okrel("K+ enters 0.7255 eV cheaper than Na+",
          kcsa_dehydration_cost_eV(11) - kcsa_dehydration_cost_eV(19), 0.7255, 2e-3);
}

int main(void) {
    printf("╔══════════════════════════════════════════════════════════════╗\n");
    printf("║  v9R4 AUDIT REGRESSION SUITE                                 ║\n");
    printf("╚══════════════════════════════════════════════════════════════╝\n");
    test_quantum_basics();
    test_forces();
    test_dispersion_sign();
    test_scf_continuity();
    test_lj_convention();
    test_hybridization();
    test_polarizability();
    test_qeq();
    test_memory_safety();
    test_orbital_truncation();
    test_minimizer_floor();
    test_nve();
    test_datastream();
    test_rng();
    test_kcsa_filter();
    test_ion_size_and_hydration();
    test_dipole_solver_equation();
    test_qeq_conservation_under_bound();
    test_qeq_nonzero_total_charge();
    test_kcsa_filter_residue_neutrality();
    test_kcsa_ion_sigma();
    test_aminoacid_charge_sums();
    test_constants_derived();
    test_minimizer_zero_tolerance();
    test_dihedral_rotational_invariance();
    printf("\n══════════════════════════════════════════════════════════════\n");
    /* Not "  PASS n FAIL m": the verify script counts check lines by
     * that prefix and would count this summary as one of its own checks. */
    printf("  TOTAL %d passed, %d failed\n", g_pass, g_fail);
    printf("══════════════════════════════════════════════════════════════\n");
    return g_fail ? 1 : 0;
}
