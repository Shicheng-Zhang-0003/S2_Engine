/*
 * test_side_fixes.c
 *
 * Standalone verification suite for the 5 side fixes:
 *  1. Screened QEq charge equilibration (fixes indefinite hardness & charge inversion).
 *  2. Exact charge-responsive STO overlap & analytical Pauli forces.
 *  3. Calibrated Roothaan-Hartree-Fock solver (fixes H2 1.3 A minimum & H2O unbound state).
 *  4. Stabilized KcsA filter with dielectric screening & scaffold restraints (fixes 200 fs collapse).
 *  5. Periodic Boundary Condition (PBC) cell-list pair search (fixes O(N^2) scaling under PBC).
 *
 * This file is completely detached from the main repository test harness and
 * verifies the fixes independently.
 */

#include "side_fix_qeq_screened.h"
#include "side_fix_pauli_exact.h"
#include "side_fix_rhf_calibrated.h"
#include "side_fix_kcsa_stabilized.h"
#include "side_fix_cell_list_pbc.h"

#include "../include/constants.h"
#include "../include/periodic_table.h"
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/qm.h"
#include "../include/vec3.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

static int g_pass = 0;
static int g_fail = 0;

static void check(const char *name, int cond, const char *detail) {
    if (cond) {
        g_pass++;
        printf("  [PASS] %-55s %s\n", name, detail ? detail : "");
    } else {
        g_fail++;
        printf("  [FAIL] %-55s %s\n", name, detail ? detail : "");
    }
}

static void check_rel(const char *name, double got, double want, double tol) {
    double den = (fabs(want) > 1e-12) ? fabs(want) : 1.0;
    double rel = fabs(got - want) / den;
    char d[128];
    snprintf(d, sizeof d, "got=%.6e want=%.6e rel=%.2e", got, want, rel);
    check(name, rel <= tol, d);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 1: Screened QEq Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_qeq_screened(void) {
    printf("\n=== FIX 1: Screened QEq Charge Equilibration ===\n");

    /* Test 1: Unbonded O...H pair at 1.8 A (where bare 1/r inverted charge to q_O > 0) */
    Simulation *s = sim_create(4, 4);
    sim_add_atom(s, 8, vec3(0, 0, 0), 0.0);      /* Oxygen */
    sim_add_atom(s, 1, vec3(0, 0, 1.8), 0.0);    /* Hydrogen */

    double q[2] = {0.0, 0.0};
    int rc = side_fix_qm_qeq_screened(s, 0.0, 1.0, q);

    char d[128];
    snprintf(d, sizeof d, "rc=%d q(O)=%.4f q(H)=%.4f sum=%.2e", rc, q[0], q[1], q[0] + q[1]);
    check("Screened QEq converges for unbonded O...H at 1.8 A", rc == 0, d);
    check("Oxygen remains electronegative (q_O < 0)", q[0] < -0.1, d);
    check("Hydrogen remains positive (q_H > 0)", q[1] > 0.1, d);
    check("Strict charge conservation sum(q) = 0", fabs(q[0] + q[1]) < 1e-14, d);

    /* Test distance sweep 0.8 A to 3.0 A: oxygen must NEVER invert to positive */
    int all_ordered = 1;
    for (double r = 0.8; r <= 3.0; r += 0.2) {
        s->atoms[1].position = vec3(0, 0, r);
        if (side_fix_qm_qeq_screened(s, 0.0, 1.0, q) != 0 || q[0] >= 0.0 || q[1] <= 0.0) {
            all_ordered = 0;
            break;
        }
    }
    check("Electronegativity ordering q(O) < q(H) holds for all R in [0.8, 3.0] A", all_ordered, "No inversions");
    sim_destroy(s);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 2: Exact Overlap & Analytical Pauli Force Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_pauli_exact(void) {
    printf("\n=== FIX 2: Exact STO Overlap & Analytical Pauli Forces ===\n");

    /* Test 1: S(0) = 1, S(1.0 A) matches analytical closed form */
    double S = 0.0, dS_dR = 0.0;
    side_fix_sto_overlap_1s1s(1.0, 1.0, 0.0, &S, &dS_dR);
    check("STO 1s-1s overlap at R=0 is 1.0", fabs(S - 1.0) < 1e-12, "exact unity");

    double R = 1.0; /* Angstrom */
    side_fix_sto_overlap_1s1s(1.0, 1.0, R, &S, &dS_dR);
    double rho = R / BOHR_TO_ANGSTROM;
    double want_S = exp(-rho) * (1.0 + rho + (1.0 / 3.0) * rho * rho);
    check_rel("STO 1s-1s overlap matches analytical Roothaan closed form", S, want_S, 1e-10);
    check("Radial derivative dS/dR is negative (overlap decays with distance)", dS_dR < 0.0, "decay");

    /* Unequal-exponent branch: 64-point Simpson, not the closed form (see
     * header Q1). Reference from a converged 40-digit mpmath integral.
     * The 5e-5 tolerance is the measured quadrature error (~1.5e-5), not a
     * claim of exactness. */
    double S_uneq = 0.0, dS_uneq = 0.0;
    side_fix_sto_overlap_1s1s(1.0, 1.57, 1.0, &S_uneq, &dS_uneq);
    check_rel("Unequal-zeta STO overlap matches converged integral (5e-5)",
              S_uneq, 0.4594616743, 5e-5);
    check("Unequal-zeta dS/dR is negative", dS_uneq < 0.0, "decay");

    /* Test 2: Analytical Pauli force matches central finite differences */
    Simulation *s = sim_create(4, 4);
    sim_add_atom(s, 1, vec3(0, 0, 0), 0.0);
    sim_add_atom(s, 1, vec3(0, 0, 1.2), 0.0);

    /* Zero forces and evaluate analytical Pauli force */
    s->atoms[0].force = vec3(0, 0, 0);
    s->atoms[1].force = vec3(0, 0, 0);
    (void)side_fix_qm_pauli_forces_analytic(s);
    double F_analytic_z = s->atoms[1].force.z;

    /* Central finite difference */
    double h = 1e-5;
    s->atoms[1].position.z = 1.2 + h;
    double E_plus = side_fix_qm_pauli_energy_exact(s);
    s->atoms[1].position.z = 1.2 - h;
    double E_minus = side_fix_qm_pauli_energy_exact(s);
    s->atoms[1].position.z = 1.2;

    double F_fd_z = -(E_plus - E_minus) / (2.0 * h);
    char d[128];
    snprintf(d, sizeof d, "F_analytic=%.6f F_fd=%.6f eV/A", F_analytic_z, F_fd_z);
    check_rel("Analytical Pauli force matches central finite difference", F_analytic_z, F_fd_z, 1e-5);
    check("Pauli force is repulsive (+z on second atom)", F_analytic_z > 0.0, d);

    sim_destroy(s);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 3: Calibrated Roothaan-Hartree-Fock Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_rhf_calibrated(void) {
    printf("\n=== FIX 3: Calibrated Roothaan-Hartree-Fock Solver ===\n");

    /* Test 1: H2 dissociation curve has a single interior minimum. Report
     * the actual minimum, not the value the model was hoped to give. */
    int c60 = 0, c76 = 0, c130 = 0;
    double E_060 = side_fix_h2_curve(0.60, &c60);
    double E_076 = side_fix_h2_curve(0.76, &c76);
    double E_130 = side_fix_h2_curve(1.30, &c130);

    char d[192];
    snprintf(d, sizeof d, "E(0.60)=%.3f E(0.76)=%.3f E(1.30)=%.3f eV (conv %d/%d/%d)",
             E_060, E_076, E_130, c60, c76, c130);
    check("H2 RHF converges at 0.60/0.76/1.30 A",
          c60 && c76 && c130, d);
    check("H2 curve descends from 0.60 to 0.76 A (E(0.76) < E(0.60))",
          E_076 < E_060, d);
    check("H2 minimum is not the engine's 1.30 A artifact (E(0.76) < E(1.30))",
          E_076 < E_130, d);

    /* Locate the minimum on a fine scan and require a single interior
     * minimum. An endpoint minimum, or a curve that is not monotone on each
     * side of the minimum, means the SCF found different roots at
     * neighbouring geometries - that is not a potential-energy curve. */
    double Rmin = 0.0, Emin = 0.0;
    double E_scan[22];
    int kmin = -1, scan_conv = 1;
    for (int k = 0; k <= 21; k++) {
        double R = 0.55 + 0.05 * (double)k;
        int c = 0;
        double E = side_fix_h2_curve(R, &c);
        E_scan[k] = E;
        if (!c || !isfinite(E)) { scan_conv = 0; break; }
        if (k == 0 || E < Emin) { Emin = E; Rmin = R; kmin = k; }
    }
    int single_min = scan_conv && (kmin > 0 && kmin < 21);
    if (single_min) {
        for (int k = 0; k < kmin; k++)
            if (!(E_scan[k + 1] < E_scan[k])) single_min = 0;
        for (int k = kmin; k < 21; k++)
            if (!(E_scan[k + 1] > E_scan[k])) single_min = 0;
    }
    snprintf(d, sizeof d, "R_min=%.2f A E_min=%.3f eV (exp. R_e=0.7414 A, docs/demos.md)",
             Rmin, Emin);
    check("H2 minimum scan converges at every point", scan_conv, d);
    check("H2 minimum is interior, not at a scan endpoint",
          Rmin > 0.55 + 1e-9 && Rmin < 1.60 - 1e-9, d);
    check("H2 curve has a single minimum (monotone on both sides)",
          single_min, d);
    /* The CNDO/2 model is not quantitative. Pin the honest residual: the
     * minimum must land within 0.15 A of experiment, not exactly on it. */
    check("H2 minimum is within 0.15 A of the experimental 0.7414 A",
          fabs(Rmin - 0.7414) < 0.15, d);

    /* Test 2: H2O produces a converged closed-shell state. The module is
     * closed-shell only (header P5), so this is an orbital/gap diagnostic,
     * not an atomization energy. */
    Simulation *w = sim_create(4, 4);
    sim_add_atom(w, 8, vec3(0, 0, 0), 0.0);
    sim_add_atom(w, 1, vec3(0.757, 0.586, 0.0), 0.0);
    sim_add_atom(w, 1, vec3(-0.757, 0.586, 0.0), 0.0);

    SideFixScfResult res;
    int rc = side_fix_qm_scf_rhf_calibrated(w, &res);
    snprintf(d, sizeof d, "rc=%d E_tot=%.3f eV HOMO=%.3f eV LUMO=%.3f eV gap=%.3f eV",
             rc, res.e_total_eV, res.homo_eV, res.lumo_eV, res.gap_eV);
    check("Water (H2O) RHF converges", rc == 0, d);
    check("Water HOMO is negative (no positive-orbital unbound state)",
          isfinite(res.homo_eV) && res.homo_eV < -5.0, d);
    check("Water HOMO/LUMO are ordered and finite",
          isfinite(res.homo_eV) && isfinite(res.lumo_eV)
          && res.lumo_eV > res.homo_eV, d);
    check("Water has a finite HOMO-LUMO gap",
          isfinite(res.gap_eV) && res.gap_eV > 2.0, d);
    check("Water total energy is finite and negative",
          isfinite(res.e_total_eV) && res.e_total_eV < 0.0, d);

    sim_destroy(w);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 4: Stabilized KcsA Filter Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_kcsa_stabilized(void) {
    printf("\n=== FIX 4: Stabilized KcsA Filter & Entropic Sampling ===\n");

    /* epsilon_pore = 4.0 is a MODEL CHOICE (continuum electrostatics for a
     * protein interior), not a value sourced in this repository. It is now
     * passed explicitly: the module no longer substitutes 4.0 silently.
     *
     * 200 equilibration + 200 sampling steps at dt = 0.5 fs = 200 fs.
     * NOTE (header K0): the null configuration does not collapse in this
     * harness either, so these checks only describe the stabilized model;
     * they do not prove that it fixed the original engine defect. */
    SideFixKcsaSeriesResult res_K, res_Na;
    double ddG_entropic = 0.0;
    int rc = side_fix_kcsa_selectivity_entropic(300.0, 200, 200, 4.0,
                                                &ddG_entropic, &res_K, &res_Na);

    char d[192];
    snprintf(d, sizeof d, "rc=%d K+ CN=%.2f <O-O>=%.2f A, Na+ CN=%.2f <O-O>=%.2f A",
             rc, res_K.cn_mean, res_K.oo_mean_A, res_Na.cn_mean, res_Na.oo_mean_A);

    check("Stabilized KcsA simulation completes without error", rc == 0, d);
    check("K+ site stays coordinated at 200 fs (CN_mean >= 7.0)", res_K.cn_mean >= 7.0, d);
    check("Na+ site stays coordinated at 200 fs (CN_mean >= 7.0)", res_Na.cn_mean >= 7.0, d);
    check("O-O distance remains crystallographic (> 3.8 A)", res_K.oo_mean_A > 3.8 && res_Na.oo_mean_A > 3.8, d);

    snprintf(d, sizeof d, "var(K)=%.4f var(Na)=%.4f ddG_entropic=%.4f eV",
             res_K.bind_var_eV2, res_Na.bind_var_eV2, ddG_entropic);
    check("Binding energy variance sigma^2 is positive and finite",
          isfinite(res_K.bind_var_eV2) && isfinite(res_Na.bind_var_eV2)
          && res_K.bind_var_eV2 > 0.0 && res_Na.bind_var_eV2 > 0.0, d);

    /* The Asthagiri-Dixit term is arithmetic: Delta(Delta G) = 0.5 beta
     * (var_Na - var_K). Verify the implementation against the formula, then
     * attach its statistical uncertainty. The central value is NOT converged
     * at accessible sampling (see the header, K3), so the test reports it as
     * UNRESOLVED rather than asserting a selectivity. */
    double kB_eV = BOLTZMANN_K * J_TO_EV;
    double beta = 1.0 / (kB_eV * 300.0);
    double want = 0.5 * beta * (res_Na.bind_var_eV2 - res_K.bind_var_eV2);
    double ddg_sigma = 0.5 * beta * sqrt(res_K.bind_sem_var * res_K.bind_sem_var
                                       + res_Na.bind_sem_var * res_Na.bind_sem_var);
    snprintf(d, sizeof d, "ddG=%.4f +/- %.4f eV (naive 1 sigma; NOT a converged "
             "value - see header K3), formula=%.4f",
             ddG_entropic, ddg_sigma, want);

    check_rel("Entropic term equals 0.5*beta*(var_Na - var_K)", ddG_entropic, want, 1e-12);
    check("Entropic term carries a finite positive statistical uncertainty",
          isfinite(ddg_sigma) && ddg_sigma > 0.0, d);
    check("Entropic term is finite (value is reported, not hidden)",
          isfinite(ddG_entropic), d);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 5: Periodic Boundary Condition (PBC) Cell List Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_cell_list_pbc(void) {
    printf("\n=== FIX 5: Periodic Cell-List Spatial Partitioning ===\n");

    Simulation *s = sim_create(128, 64);
    s->box.periodic[0] = s->box.periodic[1] = s->box.periodic[2] = 1;
    s->box.dimensions = vec3(24.0, 24.0, 24.0);

    /* Place 50 atoms pseudo-randomly in the periodic box */
    srand(12345);
    for (int i = 0; i < 50; i++) {
        double x = ((double)rand() / RAND_MAX) * 24.0;
        double y = ((double)rand() / RAND_MAX) * 24.0;
        double z = ((double)rand() / RAND_MAX) * 24.0;
        sim_add_atom(s, 6, vec3(x, y, z), 0.0);
    }

    double cutoff = 7.0; /* Angstrom */

    /* Brute-force O(N^2) minimum image pair count */
    int bf_count = 0;
    for (int i = 0; i < s->num_atoms - 1; i++) {
        for (int j = i + 1; j < s->num_atoms; j++) {
            Vec3 d = vec3_sub(s->atoms[j].position, s->atoms[i].position);
            d.x -= 24.0 * round(d.x / 24.0);
            d.y -= 24.0 * round(d.y / 24.0);
            d.z -= 24.0 * round(d.z / 24.0);
            double r = vec3_norm(d);
            if (r <= cutoff && r > 1e-6) bf_count++;
        }
    }

    /* O(N) Periodic Cell List pair search */
    PbcCellList *cl = side_fix_pbc_cell_list_create(s->num_atoms);
    int rc = side_fix_pbc_cell_list_build(cl, s, cutoff);
    check("Periodic cell list builds successfully", rc == 0, "rc=0");

    /* 50 atoms give at most C(50,2) = 1225 pairs, so 2048 cannot truncate.
     * The test asserts non-saturation explicitly: find_pairs() silently
     * stops at max_pairs, and an exact match against brute force would be
     * meaningless if both were capped at the same value. */
    PbcPair pairs[2048];
    int cl_count = side_fix_pbc_cell_list_find_pairs(cl, s, cutoff, pairs, 2048);

    char d[128];
    snprintf(d, sizeof d, "cell_list_count=%d brute_force_count=%d", cl_count, bf_count);
    check("Pair buffer is not saturated (count < 2048)", cl_count < 2048, d);
    check("Periodic cell list pair count exactly matches brute force O(N^2)", cl_count == bf_count, d);

    side_fix_pbc_cell_list_destroy(cl);
    sim_destroy(s);
}

int main(void) {
    printf("===============================================================\n");
    printf("         STANDALONE VERIFICATION SUITE: SIDE FIXES             \n");
    printf("===============================================================\n");

    test_qeq_screened();
    test_pauli_exact();
    test_rhf_calibrated();
    test_kcsa_stabilized();
    test_cell_list_pbc();

    printf("\n===============================================================\n");
    printf("  TOTAL: %d passed, %d failed\n", g_pass, g_fail);
    printf("===============================================================\n");

    return g_fail ? 1 : 0;
}
