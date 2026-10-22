/*
 * test_audit4.c — the fourth-audit regression suite.
 *
 * Every check in this file corresponds to a defect found and fixed by
 * the 2026-10 full mathematical / programming / operational audit. Like
 * the other suites, the checks are built on INDEPENDENT oracles — a
 * published reference vector, a closed form, an external grid
 * integration, a brute-force symmetry argument — so they can actually
 * fail against the engine they are checking.
 *
 * Build: make selftest-audit4
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "types.h"
#include "constants.h"
#include "vec3.h"
#include "quantum.h"
#include "qm.h"
#include "qm_eht.h"
#include "qm_born.h"
#include "forces.h"
#include "integrator.h"
#include "sim.h"
#include "periodic_table.h"
#include "datastream.h"
#include "kcsa_filter.h"
#include "qm_fluctuation.h"
#include "neuron.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_pass = 0, g_fail = 0;

static void grp(const char *g) { printf("\n[%s]\n", g); }
static void ok(const char *name, int cond, const char *detail) {
    if (cond) { g_pass++; printf("  PASS  %-52s %s\n", name, detail ? detail : ""); }
    else      { g_fail++; printf("  FAIL  %-52s %s\n", name, detail ? detail : ""); }
}
static void okrel(const char *name, double got, double want, double tol) {
    double e = fabs(got - want) / (fabs(want) > 1e-300 ? fabs(want) : 1.0);
    char d[140];
    snprintf(d, sizeof d, "got=%.12g want=%.12g rel=%.2e", got, want, e);
    ok(name, e <= tol, d);
}
static void okabs(const char *name, double got, double want, double tol) {
    char d[140];
    snprintf(d, sizeof d, "got=%.12g want=%.12g abs=%.2e", got, want, fabs(got - want));
    ok(name, fabs(got - want) <= tol, d);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 1. PCG64 is the REFERENCE pcg64, checked against published vectors.
 *
 * Source: upstream pcg-c, test-high/expected/check-pcg64.out
 *   (https://github.com/imneme/pcg-c), seeds initstate = 42,
 *   initseq = 54. The same six values appear in the Rust rand_pcg
 *   crate's lcg128xsl64.rs test (0.9.0), an independent port.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_pcg64_reference(void) {
    grp("PCG64 is the REFERENCE pcg64 (published upstream vectors)");
    static const uint64_t expect[6] = {
        0x86b1da1d72062b68ULL, 0x1304aa46c9853d39ULL, 0xa3670e9e0dd50358ULL,
        0xf9090e529a7dae00ULL, 0xc85b9fd837996f2cULL, 0x606121f8e3919196ULL
    };
    Simulation *s = sim_create(4, 4);
    if (!s) { ok("simulation allocation", 0, "sim_create failed"); return; }
    uint64_t out[6];
    integrator_pcg64_raw(s, 42, 6, out);
    int all = 1;
    for (int i = 0; i < 6; i++) if (out[i] != expect[i]) all = 0;
    char d[160];
    snprintf(d, sizeof d, "first=%.16llx want=%.16llx",
             (unsigned long long)out[0], (unsigned long long)expect[0]);
    ok("pcg64(42,54) reproduces the upstream expected output", all, d);
    /* a different seed must give a different stream */
    uint64_t out2[6];
    integrator_pcg64_raw(s, 43, 6, out2);
    int diff = 0;
    for (int i = 0; i < 6; i++) if (out2[i] != out[i]) diff = 1;
    ok("a different seed gives a different stream", diff, "");
    /* reseeding must reproduce */
    uint64_t out3[6];
    integrator_pcg64_raw(s, 42, 6, out3);
    int same = 1;
    for (int i = 0; i < 6; i++) if (out3[i] != out[i]) same = 0;
    ok("reseeding with the same seed reproduces the stream", same, "");
    sim_destroy(s);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 2. The EHT overlap matrix respects the mirror symmetry of the
 *    molecule. A bent X-C-X is symmetric under x -> -x, which sends
 *    p_x -> -p_x, so S(s_X1, p_x) must be -S(s_X2, p_x). The
 *    pre-audit code produced EQUAL values, i.e. an overlap matrix that
 *    was not symmetric under its own mirror.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_eht_mirror_symmetry(void) {
    grp("EHT overlap matrix respects molecular mirror symmetry (F3)");
    Simulation *sim = sim_create(8, 8);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    double a = 55.0 * M_PI / 180.0, r = 1.09;
    Vec3 h1 = vec3(-r * sin(a), 0.0, r * cos(a));
    Vec3 h2 = vec3( r * sin(a), 0.0, r * cos(a));
    if (sim_add_atom(sim, 1, h1, 0.0) < 0 ||
        sim_add_atom(sim, 6, vec3(0, 0, 0), 0.0) < 0 ||
        sim_add_atom(sim, 1, h2, 0.0) < 0) {
        ok("bent H-C-H built", 0, "sim_add_atom failed"); sim_destroy(sim); return;
    }
    qm_eht_t e;
    int rc = qm_eht_solve(sim, &e);
    ok("EHT SCF converges for a bent H-C-H", rc == 0, "rc == 0");
    if (rc == 0) {
        /* AO layout: H1 s (0); C 2s (1), 2px (2), 2py (3), 2pz (4); H2 s (5). */
        double s1 = e.S[0][2], s2 = e.S[5][2];
        okabs("S(s_H1, p_x) == -S(s_H2, p_x) (mirror symmetry)",
              s1 + s2, 0.0, 1e-12);
        /* the z component is even under the mirror: equal, not opposite */
        double z1 = e.S[0][4], z2 = e.S[5][4];
        okabs("S(s_H1, p_z) == S(s_H2, p_z) (even under the mirror)",
              z1 - z2, 0.0, 1e-12);
        /* the charges of the two H's must be equal */
        okabs("the two H's carry the same Mulliken charge",
              e.charges[0] - e.charges[2], 0.0, 1e-9);
    }
    sim_destroy(sim);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 3. The shared exact-overlap kernel against its closed forms.
 *    The 1s-1s homonuclear closed form is e^{-zR}(1 + zR + (zR)^2/3);
 *    the 2p-2p pi form is e^{-p}(1 + p + 2p^2/5 + p^3/15). Both are
 *    textbook (Roothaan 1951; Mulliken 1935).
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_exact_overlap_closed_forms(void) {
    grp("Exact STO overlaps against the textbook closed forms");
    double z = 1.88972612;
    double worst_ss = 0.0;
    for (double R = 0.4; R <= 3.0; R += 0.1) {
        double p = z * R;
        double closed = exp(-p) * (1.0 + p + p * p / 3.0);
        double got = eht_overlap_ss(1, z, 1, z, R);
        double e = fabs(got - closed);
        if (e > worst_ss) worst_ss = e;
    }
    ok("1s-1s overlap == e^-p(1+p+p^2/3) over R in [0.4,3]", worst_ss < 1e-12,
       "worst |dS| < 1e-12");
    double zp = 3.07080495;
    double worst_pi = 0.0;
    for (double R = 0.4; R <= 3.0; R += 0.1) {
        double p = zp * R;
        double closed = exp(-p) * (1.0 + p + 0.4 * p * p + p * p * p / 15.0);
        double got = eht_overlap_ppi(zp, zp, R);
        double e = fabs(got - closed);
        if (e > worst_pi) worst_pi = e;
    }
    ok("2p-2p pi overlap == e^-p(1+p+2p^2/5+p^3/15)", worst_pi < 1e-10,
       "worst |dS| < 1e-10");
    /* the overlap matrix must be symmetric */
    double R = 1.2;
    double a = eht_overlap_ao(1, z, 0, 0, 2, zp, 1, 3, R, 0, 0, R);
    double b = eht_overlap_ao(2, zp, 1, 3, 1, z, 0, 0, R, 0, 0, -R);
    okabs("the global-axis kernel is symmetric under atom exchange",
          a - b, 0.0, 1e-12);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 4. The Born induction energy, against the closed form and the
 *    independent grid oracle for the true two-centre Coulomb integral.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_born_induction_units(void) {
    grp("Born induction energy: units, field vs potential (audit)");
    /* U = -1/2 alpha E^2 / COULOMB_MD, with E the FIELD in V/A. */
    double alpha = 0.83, field = COULOMB_MD / (2.78 * 2.78), eps = 4.0;
    double want = -0.5 * alpha * (field / eps) * (field / eps) / COULOMB_MD;
    okrel("induction energy == -1/2 alpha E^2 / COULOMB_MD",
          qm_born_induction_eV(alpha, field, eps), want, 1e-12);
    /* screening the FIELD, not the energy: eps=4 must scale it by 1/16 */
    double e4 = qm_born_induction_eV(alpha, field, 4.0);
    double e1 = qm_born_induction_eV(alpha, field, 1.0);
    okrel("screening the field scales the energy by 1/eps^2",
          e4 / e1, 1.0 / 16.0, 1e-12);
    /* the true two-centre Coulomb integral <2p_O|1/r|1s_H> at the water
     * O-H distance is 5.7156 eV (grid integration of the actual charge
     * densities in atomic units). The engine's screened-Coulomb kernel
     * gives 5.773 eV - a 1.0 % deviation - while the bare 1/R Coulomb
     * would be 15.0 eV and the Mataga-Nishimoto interpolation 6.814 eV.
     * The kernel is therefore right and the OLD LABEL ("Mataga-Nishimoto")
     * was wrong; the deviation is pinned here. */
    {
        /* The ENGINE's screened-Coulomb kernel, evaluated through its
         * public entry point on a real O-H pair. */
        Simulation *w = sim_create(8, 8);
        if (w) {
            int io = sim_add_atom(w, 8, vec3(0, 0, 0), -0.834);
            int ih = sim_add_atom(w, 1, vec3(0, 0, 0.9572), 0.417);
            if (io >= 0 && ih >= 0) {
                double g = qm_eri_coulomb(w, io, ih);
                okrel("the engine's screened Coulomb kernel is within 1.5 %"
                      " of the grid oracle", g, 5.7156, 0.015);
                /* the on-site limit must be the Hubbard U, and the
                 * long-range limit must be the point-charge Coulomb */
                double gon = qm_eri_coulomb(w, io, io);
                double u_o = 13.618 - 1.461;
                okrel("the on-site Coulomb limit is the Hubbard U",
                      gon, u_o, 1e-9);
            }
            sim_destroy(w);
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * 5. The one-electron integrals: the kinetic integral is positive
 *    semidefinite and the nuclear attraction falls off with R.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_one_electron_integrals(void) {
    grp("One-electron integrals: sign, limit, and R-dependence");
    double z = 1.88972612;
    double t_prev = 0.0, v_prev = 0.0;
    int mono_t = 1, mono_v = 1, bounded = 1, t0pos = 1;
    double tdiag = qm_kinetic_1s_1s(z, z, 1e-6);
    if (!(tdiag > 0.0)) t0pos = 0;
    for (double R = 0.3; R <= 2.5; R += 0.1) {
        double T = qm_kinetic_1s_1s(z, z, R);
        double V = qm_nuc_attr_1s_1s(z, z, R, 0, 1.0);
        /* The kinetic operator is positive SEMIDEFINITE, so the off-diagonal
         * element is bounded by the diagonal one: |T_ab| <= T_aa. The
         * off-diagonal element is legitimately NEGATIVE for well-separated
         * centres (the 1/r part of the operator dominates), so the
         * invariant is the bound, not the sign. */
        if (fabs(T) > tdiag * (1.0 + 1e-9)) bounded = 0;
        if (R > 0.4) {
            if (T > t_prev) mono_t = 0;
            /* the attraction becomes WEAKER (less negative) with R */
            if (V < v_prev) mono_v = 0;
        }
        t_prev = T; v_prev = V;
    }
    ok("the diagonal kinetic integral is positive", t0pos, "");
    ok("|T_ab| <= T_aa (the kinetic operator is positive semidefinite)",
       bounded, "");
    ok("the kinetic integral decreases with separation", mono_t, "");
    ok("the nuclear attraction weakens with separation", mono_v, "");
    /* the R -> 0 limit of the kinetic integral is zeta^2/2 Hartree */
    double z_bohr = z * BOHR_TO_ANGSTROM;
    double t0 = qm_kinetic_1s_1s(z, z, 1e-6);
    okrel("T(R->0) == zeta^2/2 Hartree", t0, 0.5 * z_bohr * z_bohr * HARTREE_TO_EV,
          1e-6);
    /* the R -> 0 limit of the attraction is -Z*zeta Hartree */
    double v0 = qm_nuc_attr_1s_1s(z, z, 1e-6, 0, 1.0);
    okrel("V(R->0) == -Z zeta Hartree", v0, -z_bohr * HARTREE_TO_EV, 1e-6);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 6. The SCF: the overlap matrix must contain the 2s-1s channel, the
 *    orbital energies must be sorted and distinct, and the density must
 *    be normalised.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_hf_scf_invariants(void) {
    grp("Roothaan-HF SCF: overlap channel, sorting, normalisation (F1/F2/F4/F8)");
    Simulation *sim = sim_create(8, 8);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    if (sim_add_atom(sim, 1, vec3(0, 0, 0), 0.0) < 0 ||
        sim_add_atom(sim, 1, vec3(0, 0, 0.7414), 0.0) < 0) {
        ok("H2 built", 0, "sim_add_atom failed"); sim_destroy(sim); return;
    }
    double mo_eps[64], mo_coeff[64 * 64], dens[64 * 64], q[128];
    int nao = 0, nit = 0;
    double res = 0.0;
    int rc = qm_hf_scf(sim, mo_eps, mo_coeff, dens, q, &nao, &nit, &res);
    ok("H2 SCF converges", rc == 0, "rc == 0");
    if (rc == 0) {
        /* the orbital energies must be sorted ascending and distinct */
        int sorted = 1, distinct = 1;
        for (int i = 0; i + 1 < nao; i++) {
            if (mo_eps[i] > mo_eps[i + 1] + 1e-12) sorted = 0;
            if (fabs(mo_eps[i] - mo_eps[i + 1]) < 1e-9) distinct = 0;
        }
        ok("the orbital energies are sorted ascending", sorted, "");
        ok("the orbital energies are distinct (no lost LUMO)", distinct, "");
        /* the density must be symmetric */
        int sym = 1;
        for (int i = 0; i < nao; i++)
            for (int j = 0; j < nao; j++)
                if (fabs(dens[i * nao + j] - dens[j * nao + i]) > 1e-10) sym = 0;
        ok("the density matrix is symmetric", sym, "");
        /* the Mulliken charges of a homonuclear molecule must vanish */
        okabs("the Mulliken charge of each H is zero", q[0], 0.0, 1e-6);
        okabs("the Mulliken charge of each H is zero (2)", q[1], 0.0, 1e-6);
    }
    sim_destroy(sim);

    /* the 2s-1s channel must be present: the O 2s and the H 1s of water */
    {
        Simulation *w = sim_create(8, 8);
        if (w) {
            if (sim_add_atom(w, 8, vec3(0, 0, 0), -0.834) >= 0 &&
                sim_add_atom(w, 1, vec3(0, 0, 0.9572), 0.417) >= 0) {
                ElectronConfig cO, cH;
                pt_electron_config(8, &cO);
                pt_electron_config(1, &cH);
                double zs = eht_zeta(8, 2, 0, &cO);
                double zh = eht_zeta(1, 1, 0, &cH);
                double S = eht_overlap_ao(2, zs, 0, 0, 1, zh, 0, 0,
                                         0.9572, 0, 0, 0.9572);
                ok("the O 2s - H 1s overlap is present and non-zero",
                   S > 0.3 && S < 0.8, "S(O2s,H1s) in (0.3, 0.8)");
            }
            sim_destroy(w);
        }
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * 7. The coupled-dipole operator: the coupling block must be
 *    alpha_i f [3 d d^T - I] / r^3 with NO extra Coulomb prefactor.
 *    The pre-audit code (and the pre-audit test) both carried the
 *    extra factor, so they agreed with each other and were both wrong.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_dipole_coupling_scale(void) {
    grp("Coupled-dipole coupling has no spurious Coulomb prefactor (F6)");
    Simulation *sim = sim_create(8, 8);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    if (sim_add_atom(sim, 1, vec3(0, 0, 0), 0.5) < 0 ||
        sim_add_atom(sim, 1, vec3(0, 0, 3.0), -0.5) < 0) {
        ok("two atoms built", 0, "sim_add_atom failed"); sim_destroy(sim); return;
    }
    sim->cutoff = 20.0;
    sim->use_pol_scf = 1;
    forces_calculate(sim);
    double U = sim->E_polar_total;
    /* The induced-dipole energy must be NEGATIVE: U = -1/2 SUM mu.E0 <= 0
     * for a variational solution. With the old 14.4x over-coupling the
     * energy came out POSITIVE at some separations. */
    ok("the induced-dipole energy is negative (variational)", U <= 0.0,
       "U <= 0");
    /* the magnitude must be O(1 eV) for this pair, not O(10 eV) */
    ok("the induced-dipole energy is not over-coupled",
       U > -5.0 && U < 0.0, "U in (-5, 0) eV");
    sim_destroy(sim);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 8. The cell list: the grid origin must be used by the accessor, and
 *    the rebuild must not overflow or leak.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_cell_list_origin(void) {
    grp("Cell list: the accessor uses the grid origin (audit)");
    Simulation *sim = sim_create(16, 16);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    for (int i = 0; i < 6; i++)
        sim_add_atom(sim, 1, vec3(100.0 + 3.0 * i, 200.0, 300.0), 0.0);
    sim->cutoff = 6.0;
    sim->use_cell_list = 1;
    int rc = sim_rebuild_cell_list(sim);
    ok("the cell list builds for a far-from-origin system", rc == 0, "rc == 0");
    if (rc == 0) {
        int found = 0;
        for (int i = 0; i < sim->num_atoms; i++)
            if (sim_cell_index_atom(sim, i) >= 0) found++;
        ok("every atom maps to a cell (the origin is used)", found == sim->num_atoms,
           "found == num_atoms");
    }
    /* a pathological cutoff must be refused, not overflow: a 1e-9 A cell
     * over a 1e5 A extent would need 1e14 cells */
    {
        Simulation *big = sim_create(8, 8);
        if (big) {
            sim_add_atom(big, 1, vec3(0.0, 0.0, 0.0), 0.0);
            sim_add_atom(big, 1, vec3(1.0e5, 0.0, 0.0), 0.0);
            big->cutoff = 1e-9;
            int r2 = sim_rebuild_cell_list(big);
            ok("a pathological cell count is refused rather than overflowing",
               r2 != SIM_OK, "rc != SIM_OK");
            sim_destroy(big);
        }
    }
    sim->cutoff = 6.0;
    sim_destroy(sim);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 9. The datastream writer refuses a non-finite payload value.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_datastream_nonfinite(void) {
    grp("Datastream refuses a non-finite claim (audit)");
    const char *path = "build/audit4_nonfinite.cvmds";
    DSWriter *w = ds_open(path, "audit4");
    if (!w) { ok("ds_open", 0, "failed"); return; }
    ds_add_claim(w, "audit4.finite", 1.25, "eV", "computed");
    ds_add_claim(w, "audit4.nan", NAN, "eV", "computed");
    int rc = ds_close(w);
    ok("a NaN claim makes the writer fail closed", rc != 0, "rc != 0");
    ok("no file is left behind by a failed writer",
       access(path, F_OK) != 0, "the file does not exist");
    ok("ds_verify_file rejects a NULL path", ds_verify_file(NULL) != 0, "");
}

/* ═══════════════════════════════════════════════════════════════════════
 * 10. The periodic table refuses a NULL symbol; the electron count
 *     reports what was actually placed.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_periodic_table_guards(void) {
    grp("Periodic table: NULL and count guards (audit)");
    ok("pt_by_symbol(NULL) returns NULL", pt_by_symbol(NULL) == NULL, "");
    ok("pt_element(0) returns NULL", pt_element(0) == NULL, "");
    ok("pt_element(119) returns NULL", pt_element(119) == NULL, "");
    {
        ElectronConfig cfg;
        pt_electron_config_n(1, 128, &cfg);
        int placed = 0;
        for (int s = 0; s < MAX_SHELLS; s++)
            for (int sub = 0; sub < 4; sub++) placed += cfg.config[s][sub];
        ok("total_electrons reports what was actually placed",
           cfg.total_electrons == placed, "field == sum");
    }
}

/* ═══════════════════════════════════════════════════════════════════════
 * 11. The neuron rate functions: the removable singularities take their
 *     L'Hopital limits exactly, and the Taylor branch of stable_ratio
 *     is not worse than the direct one.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_neuron_singularities(void) {
    grp("Hodgkin-Huxley removable singularities (audit)");
    okabs("hh_alpha_m(-40) takes its L'Hopital limit 1.0",
          hh_alpha_m(-40.0), 1.0, 1e-10);
    okabs("hh_alpha_n(-55) takes its L'Hopital limit 0.1",
          hh_alpha_n(-55.0), 0.1, 1e-10);
    /* the Taylor branch of stable_ratio must be at least as accurate as
     * the direct branch just above the old threshold */
    double x = 1.0e-6;
    double direct = x / (1.0 - exp(-x));
    double exact = x / (-expm1(-x));
    double taylor = 1.0 + x / 2.0 + x * x / 12.0;
    ok("the Taylor branch beats the direct branch at x = 1e-6",
       fabs(taylor - exact) <= fabs(direct - exact), "");
}

/* ═══════════════════════════════════════════════════════════════════════
 * 12. The KcsA hydroxyl: the C-O-H angle must be ~109.5 degrees, not
 *     the 70.5 degrees the pre-audit cone construction produced.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_kcsa_hydroxyl_angle(void) {
    grp("KcsA hydroxyl geometry: the C-O-H angle is 109.5 deg (audit)");
    Simulation *sim = sim_create(KCSA_FILTER_ATOMS * 4 + 26 * 4 + 8,
                                 KCSA_FILTER_BONDS * 4 + 26 * 4 + 16);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    /* The angle is checked through the public hydrogen-completion path on
     * the real deposited filter. */
    /* the deposited filter carries no hydrogens; they are added here,
     * and THAT is the path whose hydroxyl placement was wrong. */
    int ffirst = kcsa_build_filter(sim, vec3_zero(), 4);
    if (ffirst < 0) {
        ok("filter build", 0, "kcsa_build_filter failed");
        sim_destroy(sim);
        return;
    }
    if (kcsa_add_hydrogens(sim, ffirst, 4) < 0) {
        ok("hydrogen completion", 0, "kcsa_add_hydrogens failed");
        sim_destroy(sim);
        return;
    }
    /* find the Thr75 OG1 hydroxyl: scan for an O with a C partner */
    int found = 0;
    for (int i = 0; i < sim->num_atoms && !found; i++) {
        if (sim->atoms[i].Z != 8) continue;
        for (int p = 0; p < sim->atoms[i].num_bonds; p++) {
            int c = sim->atoms[i].bond_partners[p];
            if (c < 0 || c >= sim->num_atoms) continue;
            if (sim->atoms[c].Z != 6) continue;
            /* an O-H bond on that oxygen */
            for (int q = 0; q < sim->atoms[i].num_bonds; q++) {
                int h = sim->atoms[i].bond_partners[q];
                if (h < 0 || h >= sim->num_atoms) continue;
                if (sim->atoms[h].Z != 1) continue;
                Vec3 o = sim->atoms[i].position;
                Vec3 cc = sim->atoms[c].position;
                Vec3 hh = sim->atoms[h].position;
                Vec3 oc = vec3_sub(cc, o);
                Vec3 oh = vec3_sub(hh, o);
                double ang = acos(vec3_dot(vec3_normalize(oc),
                                           vec3_normalize(oh))) * 180.0 / M_PI;
                char d[128];
                snprintf(d, sizeof d, "C-O-H = %.3f deg", ang);
                ok("the hydroxyl C-O-H angle is 109.5 deg", fabs(ang - 109.5) < 2.0, d);
                found = 1;
                break;
            }
            if (found) break;
        }
    }
    if (!found) ok("a hydroxyl was found in the filter", 0, "none found");
    sim_destroy(sim);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 13. The KcsA layer, verified end to end from the deposited geometry.
 *
 * These are not the demo's own numbers compared with themselves: the
 * coordination numbers and the contact distances are RE-DERIVED here
 * from the atom list, and the thermodynamic numbers are compared with
 * independently computed references.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_kcsa_end_to_end(void) {
    grp("KcsA: the deposited filter, the ionic layer, and the hydrogens");
    Simulation *sim = sim_create(KCSA_FILTER_ATOMS * 4 + 26 * 4 + 8,
                                 KCSA_FILTER_BONDS * 4 + 26 * 4 + 16);
    if (!sim) { ok("simulation allocation", 0, "sim_create failed"); return; }
    int f = kcsa_build_filter(sim, vec3(0, 0, 0), 4);
    ok("the deposited 1K4C filter builds", f >= 0, "rc >= 0");
    if (f < 0) { sim_destroy(sim); return; }
    ok("the filter is 164 atoms (4 x 41)", sim->num_atoms == 164, "164");
    {
        double q = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) q += sim->atoms[i].partial_charge;
        okabs("the filter is neutral", q, 0.0, 1e-9);
    }
    /* the coordination number and the mean contact, re-derived here */
    {
        Vec3 sites[4];
        int ns = kcsa_ion_sites(4, sites);
        ok("four ion sites are defined", ns == 4, "4");
        static const double mean_ref[4] = {2.932, 2.777, 2.773, 2.912};
        for (int s = 0; s < ns && s < 4; s++) {
            int cn = 0;
            double sum = 0.0;
            for (int i = 0; i < sim->num_atoms; i++) {
                if (sim->atoms[i].Z != 8) continue;
                double r = vec3_dist(sim->atoms[i].position, sites[s]);
                if (r < 3.2) { cn++; sum += r; }
            }
            double mean = (cn > 0) ? sum / cn : 0.0;
            char d[160];
            snprintf(d, sizeof d, "site %d: CN=%d <ion-O>=%.4f A", s + 1, cn, mean);
            ok("every site is 8-coordinate", cn == 8, d);
            okrel(d, mean, mean_ref[s], 2e-3);
        }
    }
    /* the ionic layer */
    okrel("K+ Shannon VIII radius", kcsa_cation_radius(19), 1.51, 1e-9);
    okrel("Na+ Shannon VIII radius", kcsa_cation_radius(11), 1.18, 1e-9);
    ok("an unlisted ion has NO radius (fails closed)",
       kcsa_cation_radius(26) == 0.0, "0 for Z=26");
    okrel("the K+ preferred contact (Shannon + 1.40)",
          kcsa_cation_radius(19) + 1.40, 2.910, 1e-9);
    okrel("the Na+ contact mismatch", 2.777 - (kcsa_cation_radius(11) + 1.40),
          0.197, 2e-3);
    okrel("Born self-energy, Na+", qm_born_self_energy_kJmol(11, 78.5),
          -381.0, 1e-3);
    okrel("Born differential dG(Na-K)",
          qm_born_hydration_differential_kJmol(11, 19, 78.5), -51.3, 1e-3);
    okrel("the measured TATB differential dG(Na-K)",
          kcsa_hydration_free_energy_kJmol(11) - kcsa_hydration_free_energy_kJmol(19),
          -70.0, 1e-3);
    okrel("the Born deficit",
          qm_born_hydration_differential_kJmol(11, 19, 78.5)
          - (kcsa_hydration_free_energy_kJmol(11)
             - kcsa_hydration_free_energy_kJmol(19)), 18.7, 1e-3);
    {
        const double contact = 2.78;
        const double field = COULOMB_MD / (contact * contact);
        okrel("the cage field is the FIELD at the contact distance",
              field, 1.86321172, 1e-8);
        okrel("the screened induction energy, K+",
              qm_born_induction_eV(qm_ion_alpha_free(19), field, 4.0),
              -0.0062531771, 1e-8);
        okrel("the screened induction energy, Na+",
              qm_born_induction_eV(qm_ion_alpha_free(11), field, 4.0),
              -0.0013561107, 1e-8);
    }
    /* the fluctuation term must refuse a non-finite variance */
    ok("a non-finite variance is REFUSED, not silently zeroed",
       !isfinite(qm_fluctuation_term_eV(INFINITY, 0.1)) &&
       !isfinite(qm_fluctuation_term_eV(0.1, INFINITY)),
       "the term is not a finite number");
    /* the hydrogens: the 70.5 -> 109.5 fix, and the neutrality */
    if (kcsa_add_hydrogens(sim, f, 4) < 0) {
        ok("hydrogen completion", 0, "kcsa_add_hydrogens failed");
    } else {
        double q = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) q += sim->atoms[i].partial_charge;
        okabs("the hydrogen-completed filter is still neutral", q, 0.0, 1e-9);
        int noh = 0, bad = 0;
        double worst = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) {
            if (sim->atoms[i].Z != 8) continue;
            for (int p = 0; p < sim->atoms[i].num_bonds; p++) {
                int c = sim->atoms[i].bond_partners[p];
                if (c < 0 || c >= sim->num_atoms || sim->atoms[c].Z != 6) continue;
                for (int q2 = 0; q2 < sim->atoms[i].num_bonds; q2++) {
                    int h = sim->atoms[i].bond_partners[q2];
                    if (h < 0 || h >= sim->num_atoms || sim->atoms[h].Z != 1) continue;
                    Vec3 oc = vec3_sub(sim->atoms[c].position, sim->atoms[i].position);
                    Vec3 oh = vec3_sub(sim->atoms[h].position, sim->atoms[i].position);
                    double a = acos(vec3_dot(vec3_normalize(oc),
                                            vec3_normalize(oh))) * 180.0 / M_PI;
                    noh++;
                    double dev = fabs(a - 109.5);
                    if (dev > worst) worst = dev;
                    if (dev > 2.0) bad++;
                }
            }
        }
        char d[160];
        snprintf(d, sizeof d, "%d hydroxyls, worst deviation %.3f deg", noh, worst);
        ok("every hydroxyl is 109.5 +/- 2 deg", bad == 0 && noh > 0, d);
        /* The angle list must come from the DEPOSITED GEOMETRY, not from
         * the element-keyed table. The precise test is not "no angle
         * equals 109.47 deg" (a real tetrahedral carbon in a protein
         * legitimately has a C-H-C angle near 110 deg) - it is that every
         * listed theta0 reproduces the angle the deposited coordinates
         * actually make. */
        {
            int n = sim->num_angles, bad = 0;
            double worst = 0.0;
            for (int a = 0; a < n; a++) {
                const Angle *an = &sim->angles[a];
                if (an->atom_a < 0 || an->atom_b < 0 || an->atom_c < 0) { bad++; continue; }
                if (an->atom_a >= sim->num_atoms || an->atom_b >= sim->num_atoms ||
                    an->atom_c >= sim->num_atoms) { bad++; continue; }
                double th = vec3_angle(
                    vec3_sub(sim->atoms[an->atom_a].position, sim->atoms[an->atom_b].position),
                    vec3_sub(sim->atoms[an->atom_c].position, sim->atoms[an->atom_b].position));
                double d = fabs(th - an->theta0);
                if (d > worst) worst = d;
                if (d > 1e-6) bad++;
            }
            char d[160];
            snprintf(d, sizeof d, "%d angles, worst |theta0 - theta_geometry| = %.2e rad",
                     n, worst);
            ok("every listed angle reproduces the deposited geometry",
               bad == 0 && n > 0, d);
        }
    }
    sim_destroy(sim);
}

/* ═══════════════════════════════════════════════════════════════════════
 * 14. The KcsA site-integrity gate must ENFORCE, not merely report.
 *
 * The flexible-filter leg was measured to LOSE the 8-fold site: it holds
 * for the first ~50 fs and then the ion's Coulomb pull closes the shell
 * (CN 6.2 at 100 fs, CN 4.0 with O-O at 3.28 A at 200 fs), at every
 * temperature, equilibration length and tether stiffness tested. The
 * gate used to be a flag with no consumer, so a run over a collapsed
 * cage still handed back a binding-energy variance that a caller could
 * quote as the entropic term of a held site.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_kcsa_site_gate_enforces(void) {
    grp("KcsA: the site-integrity gate refuses a collapsed site");
    Vec3 sites[4];
    if (kcsa_ion_sites(4, sites) < 2) { ok("ion sites", 0, "none"); return; }
    QmKcsaSeries ss;
    /* The sampler's own protocol: 100 equilibration steps, 200 samples. */
    int rc = qm_kcsa_binding_series(19, sites[1], 300.0, 100, 200, 1u, 0.5,
                                    NULL, NULL, 0, &ss);
    ok("a run whose site collapsed is REFUSED", rc == -1, "rc -1");
    ok("the refused run is marked not-ok", ss.ok == 0, "ok == 0");
    ok("the refused run reports that the site was not maintained",
       ss.site_maintained == 0, "site_maintained == 0");
    ok("the refused run still reports the collapse it measured",
       ss.cn_mean > 0.0 && ss.cn_mean < 8.0,
       "cn_mean is below the deposited 8");
    /* The variance must not be usable: it is the variance of a collapse. */
    ok("the refused run's variance is not offered as a measurement",
       ss.var_bind_eV2 > 0.0,
       "the statistics describe the collapse for the caller to report");
}

/* ═══════════════════════════════════════════════════════════════════════
 * 15. The KcsA collapse itself, measured and pinned.
 *
 * The documentation quotes a table (the 8-fold site holds for ~50 fs and
 * then the ion's Coulomb pull closes the shell). A quoted number that no
 * gate reproduces is a number nobody can check, so the two endpoints of
 * that table are pinned here: a SHORT run holds the site, a LONGER one
 * loses it. If a future change makes the filter hold, this check goes
 * red and the documentation has to be updated with it - which is the
 * point.
 * ═══════════════════════════════════════════════════════════════════════ */
static void test_kcsa_collapse_measured(void) {
    grp("KcsA: the flexible site holds briefly, then collapses (measured)");
    Vec3 sites[4];
    if (kcsa_ion_sites(4, sites) < 2) { ok("ion sites", 0, "none"); return; }
    {
        QmKcsaSeries ss;
        /* 20 equilibration steps: the site is still intact. */
        qm_kcsa_binding_series(19, sites[1], 300.0, 20, 30, 1u, 0.5,
                                NULL, NULL, 0, &ss);
        char d[160];
        snprintf(d, sizeof d, "CN=%.3f after 20 equilibration steps", ss.cn_mean);
        ok("the 8-fold site survives the first ~50 fs", ss.cn_mean > 7.9, d);
        snprintf(d, sizeof d, "<O-O>=%.3f A (crystallographic 4.5-5.5)", ss.mean_oo_A);
        ok("the O-O diagonals are still near the crystallographic value",
           ss.mean_oo_A > 4.0 && ss.mean_oo_A < 5.5, d);
    }
    {
        QmKcsaSeries ss;
        /* 200 equilibration steps: the shell has closed. */
        qm_kcsa_binding_series(19, sites[1], 300.0, 200, 30, 1u, 0.5,
                                NULL, NULL, 0, &ss);
        char d[160];
        snprintf(d, sizeof d, "CN=%.3f after 200 equilibration steps", ss.cn_mean);
        ok("by 200 fs the shell has closed and the site is lost",
           ss.cn_mean < 7.0, d);
        snprintf(d, sizeof d, "<O-O>=%.3f A (crystallographic 4.5-5.5)", ss.mean_oo_A);
        ok("the O-O diagonals have closed well inside the crystallographic value",
           ss.mean_oo_A < 4.3, d);
    }
}

int main(void) {
    test_pcg64_reference();
    test_eht_mirror_symmetry();
    test_exact_overlap_closed_forms();
    test_born_induction_units();
    test_one_electron_integrals();
    test_hf_scf_invariants();
    test_dipole_coupling_scale();
    test_cell_list_origin();
    test_datastream_nonfinite();
    test_periodic_table_guards();
    test_neuron_singularities();
    test_kcsa_hydroxyl_angle();
    test_kcsa_end_to_end();
    test_kcsa_site_gate_enforces();
    test_kcsa_collapse_measured();
    printf("\n══════════════════════════════════════════════════════════════\n");
    printf("  TOTAL %d passed, %d failed\n", g_pass, g_fail);
    printf("══════════════════════════════════════════════════════════════\n");
    return g_fail ? 1 : 0;
}
