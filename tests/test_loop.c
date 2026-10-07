/*
 * tests/test_loop.c — bio/QC/QM loop-closure gate.
 *
 * Kills the "three inconsistent physics" concept by recomputing every
 * cross-track number from the SAME primaries and asserting agreement
 * with outside references:
 *
 *   Nernst (HH 1952 squid 440/50 Na, 20/400 K, 6.3 C=279.45 K)
 *   selectivity scale -kT ln(1000) at 300 K = -0.179 eV
 *   Marcus 1991 TATB free energies -295.3/-365.3 kJ/mol
 *   Shannon 1976 VIII radii K 1.51 / Na 1.18 A
 *   UFF Rmin->sigma / AMBER R*->sigma conversions
 *   CODATA-derived COULOMB_MD, MD_FORCE_CONV, kT
 *
 * Plus fail-closed NaN behaviour and the irreducible-abstraction
 * inventory (mirrors readme §7; adding a new assumption without
 * listing it here fails the gate by design).
 *
 * Build: make selftest-loop. Record impact: none (not linked into
 * carbonsim stdout path beyond shared objects).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "types.h"
#include "constants.h"
#include "amber_lj.h"
#include "kcsa_filter.h"
#include "neuron.h"
#include "loop.h"
#include "qm.h"
#include "forces.h"
#include "integrator.h"
#include "sim.h"
#include "vec3.h"
#include "periodic_table.h"

static int g_pass = 0, g_fail = 0;

static void ok(const char *name, int cond, const char *detail) {
    if (cond) { g_pass++; printf("  PASS  %-56s %s\n", name, detail ? detail : ""); }
    else      { g_fail++; printf("  FAIL  %-56s %s\n", name, detail ? detail : ""); }
}

static void okrel(const char *name, double got, double want, double tol) {
    char d[224];
    double den = fabs(want) > 1e-300 ? fabs(want) : 1.0;
    double rel = fabs(got - want) / den;
    snprintf(d, sizeof d, "got %.12g want %.12g rel %.2e", got, want, rel);
    ok(name, rel <= tol, d);
}

int main(void) {
    printf("[loop-closure: thermal scales from same k/e]\n");
    /* RT/F at squid 6.3 C (279.45 K) ~24.08 mV; at 300 K ~25.85 mV. */
    okrel("thermal-mV 279.45K", loop_thermal_mV(279.45), 24.08, 0.01);
    okrel("thermal-mV 300K", loop_thermal_mV(300.0), 25.851999, 0.005);
    okrel("thermal-eV 300K", loop_thermal_eV(300.0), 0.025851999, 0.005);
    /* kT in eV must equal BOLTZMANN_K*T/EV_TO_J exactly (same primaries). */
    okrel("kT-eV identity", loop_thermal_eV(300.0),
          (BOLTZMANN_K * 300.0) / EV_TO_J, 1e-15);

    printf("\n[loop-closure: HH reversals are Nernst of same ions]\n");
    /* Squid concentrations from HH 1952: Na_o 440/Na_i 50, K_o 20/K_i 400. */
    {
        double e_na = loop_nernst_mV(1.0, 440.0, 50.0, 279.45);
        double e_k = loop_nernst_mV(1.0, 20.0, 400.0, 279.45);
        char d1[128], d2[128];
        snprintf(d1, sizeof d1, "nernst %.3f vs HH %.1f", e_na, (double)HH_E_NA);
        snprintf(d2, sizeof d2, "nernst %.3f vs HH %.1f", e_k, (double)HH_E_K);
        /* Activity + leak mixing accounts for 2-5 mV; larger gap is incoherence. */
        ok("E_Na Nernst-consistent", fabs(e_na - (double)HH_E_NA) < 5.0, d1);
        ok("E_K Nernst-consistent", fabs(e_k - (double)HH_E_K) < 5.0, d2);
        okrel("E_Na value", e_na, 52.3, 0.05);
        okrel("E_K value", e_k, -72.1, 0.05);
    }

    printf("\n[loop-closure: selectivity scale is same kT]\n");
    /* -kT ln(1000) at 300 K = Demo 12 expt reference -0.1786 eV. */
    okrel("1000:1 scale 300K", loop_selectivity_scale_eV(1000.0, 300.0),
          -0.17857928799699513, 1e-9);
    /* Marcus TATB free energies drive the K advantage in the same units. */
    ok("Marcus K tabulated", kcsa_hydration_free_energy_kJmol(19) == -295.3, "");
    ok("Marcus Na tabulated", kcsa_hydration_free_energy_kJmol(11) == -365.3, "");
    okrel("dehyd K eV", kcsa_dehydration_cost_eV(19), 3.0605688294942199, 1e-9);
    okrel("dehyd Na eV", kcsa_dehydration_cost_eV(11), 3.7860677054325724, 1e-9);
    okrel("K enters cheaper", kcsa_dehydration_cost_eV(11) - kcsa_dehydration_cost_eV(19),
          0.72549887593835249, 1e-9);

    printf("\n[loop-closure: sizes are VIII, consistent with filter CN=8]\n");
    okrel("K VIII", kcsa_cation_radius(19), 1.51, 1e-12);
    okrel("Na VIII", kcsa_cation_radius(11), 1.18, 1e-12);
    /* Preferred contacts K-O 2.91 / Na-O 2.58 A (Shannon VIII + O 1.40). */
    okrel("K-O contact", kcsa_cation_radius(19) + 1.40, 2.91, 1e-12);
    okrel("Na-O contact", kcsa_cation_radius(11) + 1.40, 2.58, 1e-12);

    printf("\n[loop-closure: LJ + Coulomb from same primaries]\n");
    okrel("AMBER O sigma", (double)LJ_AMBER_O_SIGMA, 2.959922, 1e-5);
    /* AMBER minimum = 2*R* = 3.3224 A; UFF minimum = x1 = 3.50 A. */
    okrel("AMBER O min", loop_lj_minimum((double)LJ_AMBER_O_SIGMA),
          2.0 * 1.6612, 1e-12);
    okrel("UFF O min", loop_lj_minimum(3.11815), 3.50, 1e-4);
    okrel("COULOMB_MD derived", (double)COULOMB_MD,
          (double)(COULOMB_K * ELEM_CHARGE / 1.0e-10), 1e-15);
    ok("COULOMB_MD ~14.4", fabs((double)COULOMB_MD - 14.399645) < 0.001, "");

    printf("\n[loop-closure: fail-closed, no silent numbers]\n");
    ok("thermal bad-T NaN", isnan(loop_thermal_mV(0.0)), "");
    ok("thermal bad-T NaN eV", isnan(loop_thermal_eV(-1.0)), "");
    ok("nernst z=0 NaN", isnan(loop_nernst_mV(0.0, 440, 50, 279.45)), "");
    ok("nernst bad-c NaN", isnan(loop_nernst_mV(1.0, -1.0, 50, 279.45)), "");
    ok("selectivity bad-ratio NaN", isnan(loop_selectivity_scale_eV(0.0, 300.0)), "");
    ok("lj bad-sigma NaN", isnan(loop_lj_minimum(0.0)), "");

    printf("\n[loop-closure: abstraction inventory is explicit]\n");
    ok("11 irreducible listed", loop_abstraction_count() == 11, "");
    {
        int has_solvent = 0, has_rigid = 0, has_tracks = 0;
        for (int i = 0; i < loop_abstraction_count(); i++) {
            const char *n = loop_abstraction_name(i);
            if (!n) { ok("name non-NULL", 0, ""); break; }
            if (strstr(n, "no-bulk-solvent")) has_solvent = 1;
            if (strstr(n, "rigid-filter")) has_rigid = 1;
            if (strstr(n, "tracks-bridged")) has_tracks = 1;
        }
        ok("solvent listed", has_solvent, "");
        ok("rigid listed", has_rigid, "");
        ok("tracks listed", has_tracks, "");
        ok("OOB NULL", loop_abstraction_name(99) == 0, "");
    }

    printf("\n[best-layer: BJ + penetration damping]\n");
    okrel("BJ long-range ->1", qm_bj_f6(20.0, 3.0, 0.4289, 4.4407), 1.0, 1e-3);
    ok("BJ short-range <1", qm_bj_f6(1.0, 3.0, 0.4289, 4.4407) < 0.05, "");
    ok("BJ bad-r0 ->1", qm_bj_f6(3.0, 0.0, 0.4289, 4.4407) == 1.0, "");
    okrel("pen long-range ->1", qm_coulomb_pen(2.0, 20.0), 1.0, 1e-6);
    ok("pen short-range <1", qm_coulomb_pen(2.0, 0.5) < 0.5, "");
    ok("pen bad-b ->1", qm_coulomb_pen(0.0, 2.0) == 1.0, "");

    printf("\n[best-layer: Langevin canonical + flat-bottom + tethers]\n");
    {
        Simulation *s = sim_create(8, 8);
        int okc = (s != NULL);
        ok("alloc", okc, "");
        if (okc) {
            sim_add_atom(s, 1, vec3(0, 0, 0), 0.0);
            sim_add_atom(s, 8, vec3(2.9, 0, 0), 0.0);
            s->dt = 0.5;
            s->thermostat.type = THERMOSTAT_LANGEVIN;
            s->thermostat.target_temperature = 300.0;
            s->thermostat.tau = 100.0;
            s->thermostat.gamma = 0.01;
            forces_calculate(s);
            integrator_langevin(s);
            ok("langevin finite", isfinite(s->atoms[0].velocity.x), "");
            /* flat-bottom: free inside, harmonic outside */
            int a0 = 0;
            Vec3 an = s->atoms[a0].position;
            int r0 = sim_add_restraint_fb(s, a0, an, 0.5, 0.3);
            ok("fb add", r0 >= 0, "");
            forces_calculate(s);
            double e_in = s->E_restraint_total;
            s->atoms[a0].position = vec3(an.x + 0.1, an.y, an.z);
            forces_calculate(s);
            double e_near = s->E_restraint_total;
            s->atoms[a0].position = vec3(an.x + 1.0, an.y, an.z);
            forces_calculate(s);
            double e_out = s->E_restraint_total;
            ok("fb free inside", fabs(e_near - e_in) < 1e-12, "");
            ok("fb harmonic outside", e_out > e_in + 0.05, "");
            ok("fb bad-flat rejected", sim_add_restraint_fb(s, a0, an, 0.5, 99.0) < 0, "");
            sim_destroy(s);
        }
    }
    {
        Simulation *s = sim_create(512, 512);
        int f = -1, n = 0;
        if (s) {
            f = kcsa_build_filter(s, vec3_zero(), 4);
            if (f >= 0) n = kcsa_add_calpha_tethers(s, f, 4, 0.5);
            ok("tethers 20", n == 20, "");
            ok("tethers bad-k", kcsa_add_calpha_tethers(s, f, 4, -1.0) < 0, "");
            sim_destroy(s);
        } else ok("tether alloc", 0, "");
    }
    printf("\n[best-layer: H-complete filter holds with ions, collapses empty]\n");
    {
        Simulation *s = sim_create(2048, 4096);
        int okc = (s != NULL);
        ok("h-alloc", okc, "");
        if (okc) {
            s->dt = 0.5; s->cutoff = 12.0;
            int f = kcsa_build_filter(s, vec3_zero(), 4);
            int nh = kcsa_add_hydrogens(s, f, 4);
            ok("H count 104", nh == 104, "");
            double q = 0;
            for (int i = 0; i < s->num_atoms; i++) q += s->atoms[i].partial_charge;
            ok("H neutral", fabs(q) < 1e-9, "");
            Vec3 sites[4]; kcsa_ion_sites(4, sites);
            double mr = 0;
            ok("rigid+H CN8", kcsa_coord_stats(s, f, 4, sites[1], 3.2, &mr) == 8, "");
            /* zero-strain + geometric angles + tethers + Langevin */
            for (int b = 0; b < s->num_bonds; b++) {
                int a = s->bonds[b].atom_a, c = s->bonds[b].atom_b;
                sim_set_bond_params(s, b, vec3_dist(s->atoms[a].position, s->atoms[c].position), s->bonds[b].k);
            }
            sim_rebuild_angles_geometric(s, 4.0);
            s->use_bonds = 1; s->use_angles = 1; s->use_dihedrals = 0;
            kcsa_add_calpha_tethers(s, f, 4, 2.0);
            /* empty collapses, ions hold — 200 steps each */
            s->thermostat.type = THERMOSTAT_LANGEVIN;
            s->thermostat.target_temperature = 50.0;
            s->thermostat.tau = 100.0; s->thermostat.gamma = 0.02;
            integrator_maxwell_boltzmann(s, 50.0, 7);
            for (int i = 0; i < 200; i++) integrator_step(s);
            double mre = 0;
            int cne = kcsa_coord_stats(s, f, 4, sites[1], 3.2, &mre);
            char de[64]; snprintf(de, sizeof de, "empty CN=%d", cne);
            ok("empty collapses", cne <= 3, de);
            sim_destroy(s);
        }
    }
    {
        Simulation *s = sim_create(2048, 4096);
        if (s) {
            s->dt = 0.5; s->cutoff = 12.0;
            int f = kcsa_build_filter(s, vec3_zero(), 4);
            kcsa_add_hydrogens(s, f, 4);
            for (int b = 0; b < s->num_bonds; b++) {
                int a = s->bonds[b].atom_a, c = s->bonds[b].atom_b;
                sim_set_bond_params(s, b, vec3_dist(s->atoms[a].position, s->atoms[c].position), s->bonds[b].k);
            }
            sim_rebuild_angles_geometric(s, 4.0);
            s->use_bonds = 1; s->use_angles = 1; s->use_dihedrals = 0;
            kcsa_add_calpha_tethers(s, f, 4, 2.0);
            Vec3 sites[4]; kcsa_ion_sites(4, sites);
            for (int i = 0; i < 3; i++) { int ion = sim_add_ion(s, 19, 1, sites[i], 1.0); kcsa_set_ion_radius(s, ion, 19); }
            s->thermostat.type = THERMOSTAT_LANGEVIN;
            s->thermostat.target_temperature = 50.0;
            s->thermostat.tau = 100.0; s->thermostat.gamma = 0.02;
            integrator_maxwell_boltzmann(s, 50.0, 7);
            for (int i = 0; i < 200; i++) integrator_step(s);
            double mr = 0;
            int cn = kcsa_coord_stats(s, f, 4, sites[1], 3.2, &mr);
            char di[64]; snprintf(di, sizeof di, "ions CN=%d <r>=%.2f", cn, mr);
            ok("ions hold CN8", cn == 8, di);
            sim_destroy(s);
        } else ok("ions alloc", 0, "");
    }

    printf("\n=====================================================\n");
    printf("  LOOP: %d passed, %d failed\n", g_pass, g_fail);
    printf("=====================================================\n");
    return g_fail == 0 ? 0 : 1;
}
