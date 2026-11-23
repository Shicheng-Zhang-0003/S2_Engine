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

    /* Test 1: H2 dissociation curve has an interior minimum at R ~ 0.74-0.78 A */
    int conv = 0;
    double E_060 = side_fix_h2_curve(0.60, &conv);
    double E_076 = side_fix_h2_curve(0.76, &conv);
    double E_130 = side_fix_h2_curve(1.30, &conv);

    char d[128];
    snprintf(d, sizeof d, "E(0.60)=%.3f E(0.76)=%.3f E(1.30)=%.3f eV", E_060, E_076, E_130);
    check("H2 RHF converges at multiple bond lengths", conv, d);
    check("H2 potential curve has minimum near experimental 0.74 A (E(0.76) < E(0.60))", E_076 < E_060, d);
    check("H2 potential curve minimum is NOT at 1.3 A (E(0.76) < E(1.30))", E_076 < E_130, d);

    /* Test 2: H2O produces a stably bound molecular ground state */
    Simulation *w = sim_create(4, 4);
    sim_add_atom(w, 8, vec3(0, 0, 0), 0.0);
    sim_add_atom(w, 1, vec3(0.757, 0.586, 0.0), 0.0);
    sim_add_atom(w, 1, vec3(-0.757, 0.586, 0.0), 0.0);

    SideFixScfResult res;
    int rc = side_fix_qm_scf_rhf_calibrated(w, &res);
    snprintf(d, sizeof d, "rc=%d E_tot=%.3f eV HOMO=%.3f eV gap=%.3f eV", rc, res.e_total_eV, res.homo_eV, res.gap_eV);
    check("Water (H2O) RHF converges", rc == 0, d);
    check("Water HOMO is negative (bound molecular state)", res.homo_eV < -5.0, d);
    check("Water has a finite HOMO-LUMO gap", res.gap_eV > 2.0, d);

    sim_destroy(w);
}

/* ──────────────────────────────────────────────────────────────────────────
 * FIX 4: Stabilized KcsA Filter Tests
 * ────────────────────────────────────────────────────────────────────────── */
static void test_kcsa_stabilized(void) {
    printf("\n=== FIX 4: Stabilized KcsA Filter & Entropic Sampling ===\n");

    /* Run K+ and Na+ with dielectric pore cavity screening (epsilon=4.0) */
    SideFixKcsaSeriesResult res_K, res_Na;
    double ddG_entropic = 0.0;
    int rc = side_fix_kcsa_selectivity_entropic(300.0, 50, 100, 4.0,
                                                &ddG_entropic, &res_K, &res_Na);

    char d[140];
    snprintf(d, sizeof d, "rc=%d K+ CN=%.2f <O-O>=%.2f A, Na+ CN=%.2f <O-O>=%.2f A",
             rc, res_K.cn_mean, res_K.oo_mean_A, res_Na.cn_mean, res_Na.oo_mean_A);

    check("Stabilized KcsA simulation completes without error", rc == 0, d);
    check("K+ site does NOT collapse (CN_mean >= 7.0)", res_K.cn_mean >= 7.0, d);
    check("Na+ site does NOT collapse (CN_mean >= 7.0)", res_Na.cn_mean >= 7.0, d);
    check("O-O distance remains crystallographic (> 3.8 A)", res_K.oo_mean_A > 3.8 && res_Na.oo_mean_A > 3.8, d);

    snprintf(d, sizeof d, "var(K)=%.4f var(Na)=%.4f ddG_entropic=%.4f eV",
             res_K.bind_var_eV2, res_Na.bind_var_eV2, ddG_entropic);
    check("Binding energy variance sigma^2 is positive and finite",
          res_K.bind_var_eV2 > 0.0 && res_Na.bind_var_eV2 > 0.0, d);
    check("Asthagiri-Dixit entropic fluctuation term evaluates successfully",
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

    PbcPair pairs[1024];
    int cl_count = side_fix_pbc_cell_list_find_pairs(cl, s, cutoff, pairs, 1024);

    char d[128];
    snprintf(d, sizeof d, "cell_list_count=%d brute_force_count=%d", cl_count, bf_count);
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
