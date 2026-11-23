#include "side_fix_kcsa_stabilized.h"
#include "../include/constants.h"
#include "../include/kcsa_filter.h"
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/vec3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define CN_CUTOFF 3.4
#ifndef QM_KCSA_ANGLE_K
#define QM_KCSA_ANGLE_K 0.35
#endif

static int find_shell_oxygens(const Simulation *sim, int ion_idx, int *out_idx, int max_n) {
    int count = 0;
    for (int i = 0; i < sim->num_atoms && count < max_n; i++) {
        if (i == ion_idx || sim->atoms[i].Z != 8) continue;
        double r = vec3_dist(sim->atoms[ion_idx].position, sim->atoms[i].position);
        if (r <= CN_CUTOFF) {
            out_idx[count++] = i;
        }
    }
    return count;
}

static double calc_mean_oo(const Simulation *sim, const int *idx, int n) {
    if (n < 2) return 0.0;
    double sum = 0.0;
    int pairs = 0;
    for (int a = 0; a < n; a++) {
        for (int b = a + 1; b < n; b++) {
            sum += vec3_dist(sim->atoms[idx[a]].position, sim->atoms[idx[b]].position);
            pairs++;
        }
    }
    return pairs ? sum / pairs : 0.0;
}

int side_fix_kcsa_sample_stabilized(int ion_Z, Vec3 site, double temp_K,
                                    int n_equil, int n_sample, unsigned long seed,
                                    double epsilon_pore,
                                    SideFixKcsaSeriesResult *out) {
    if (!out || n_sample < 2 || n_equil < 0) return -1;
    if (!(temp_K > 0.0) || !isfinite(temp_K)) return -1;
    if (epsilon_pore < 1.0) epsilon_pore = 4.0; /* Default physiological channel dielectric */

    memset(out, 0, sizeof *out);

    /* Allocate simulation for TVGYG tetramer + hydrogens + ion */
    const int cap_atoms = KCSA_FILTER_ATOMS * 4 + 26 * 4 + 8;
    const int cap_bonds = KCSA_FILTER_BONDS * 4 + 26 * 4 + 16;
    Simulation *sim = sim_create(cap_atoms, cap_bonds);
    if (!sim) return -1;

    int f = kcsa_build_filter(sim, vec3(0, 0, 0), 4);
    if (f < 0) { sim_destroy(sim); return -1; }

    if (kcsa_add_hydrogens(sim, f, 4) < 0) { sim_destroy(sim); return -1; }

    /* Set up physics parameters */
    sim->use_lj = 1;
    sim->use_coulomb = 1;
    sim->cutoff = 30.0;
    sim->dielectric = epsilon_pore; /* Dielectric cavity screening prevents Coulomb collapse */
    sim->use_bonds = 1;
    sim->use_angles = 1;
    sim->use_dihedrals = 1;

    if (sim_zero_strain_bonds(sim) < 1) { sim_destroy(sim); return -1; }
    if (sim_rebuild_angles_geometric(sim, QM_KCSA_ANGLE_K) < 1) { sim_destroy(sim); return -1; }

    /* Add C-alpha harmonic tethers representing the surrounding pore helix bundle */
    if (kcsa_add_calpha_tethers(sim, f, 4, 8.0) < 0) { sim_destroy(sim); return -1; }
    if (kcsa_add_dihedrals(sim, f, 4) < 0) { sim_destroy(sim); return -1; }

    /* Add the cation at the designated site (e.g. S2 site on pore axis) */
    int ion = sim_add_ion(sim, ion_Z, 1, site, 1.0);
    if (ion < 0) { sim_destroy(sim); return -1; }
    kcsa_set_ion_radius(sim, ion, ion_Z);

    /* Constrain ion to the pore axis (x=0, y=0, z=site.z) */
    sim_add_restraint(sim, ion, vec3(0, 0, site.z), 50.0);

    /* Minimize briefly to relax contact strain */
    (void)integrator_minimize(sim, 200, 0.005, 1e-4);

    /* Setup Langevin thermostat */
    sim->dt = 0.5;
    sim->thermostat.type = THERMOSTAT_LANGEVIN;
    sim->thermostat.target_temperature = temp_K;
    sim->thermostat.gamma = 0.02;

    /* Thermalize velocities to target temperature */
    integrator_maxwell_boltzmann(sim, temp_K, seed);
    integrator_remove_com_velocity(sim);

    /* Equilibration */
    for (int step = 0; step < n_equil; step++) {
        integrator_step(sim);
    }

    /* Sampling loop */
    double cn_sum = 0.0;
    double oo_sum = 0.0;
    double bind_mean = 0.0;
    double M2 = 0.0; /* Welford running variance accumulator */
    int samples_collected = 0;

    for (int step = 0; step < n_sample; step++) {
        integrator_step(sim);

        /* Evaluate binding energy: E_with_ion - E_without_ion */
        forces_calculate(sim);
        double E_live = sim->potential_energy;

        /* Temporarily turn off ion to get same-conformation baseline */
        double orig_q = sim->atoms[ion].partial_charge;
        double orig_eps = sim->atoms[ion].lj_epsilon;
        sim->atoms[ion].partial_charge = 0.0;
        sim->atoms[ion].lj_epsilon = 0.0;
        forces_calculate(sim);
        double E_neutral = sim->potential_energy;

        /* Restore ion properties */
        sim->atoms[ion].partial_charge = orig_q;
        sim->atoms[ion].lj_epsilon = orig_eps;

        double E_bind = E_live - E_neutral;

        /* Sample coordination shell */
        int shell_idx[32];
        int cn = find_shell_oxygens(sim, ion, shell_idx, 32);
        double oo = calc_mean_oo(sim, shell_idx, cn);

        cn_sum += cn;
        oo_sum += oo;
        samples_collected++;

        /* Welford update */
        double delta = E_bind - bind_mean;
        bind_mean += delta / samples_collected;
        double delta2 = E_bind - bind_mean;
        M2 += delta * delta2;
    }

    sim_destroy(sim);

    out->n_samples = samples_collected;
    out->cn_mean = cn_sum / samples_collected;
    out->oo_mean_A = oo_sum / samples_collected;
    out->bind_mean_eV = bind_mean;
    out->bind_var_eV2 = (samples_collected > 1) ? M2 / (samples_collected - 1) : 0.0;
    out->bind_sem_var = out->bind_var_eV2 * sqrt(2.0 / (double)samples_collected);

    /* Site integrity validation:
     * Does the site maintain its 8-fold coordination structure? */
    if (out->cn_mean >= 7.0 && out->oo_mean_A >= 3.8) {
        out->ok = 1;
    } else {
        out->ok = 0;
    }

    return 0;
}

int side_fix_kcsa_selectivity_entropic(double temp_K, int n_equil, int n_sample,
                                       double epsilon_pore,
                                       double *ddG_entropic_eV,
                                       SideFixKcsaSeriesResult *res_K,
                                       SideFixKcsaSeriesResult *res_Na) {
    if (!ddG_entropic_eV || !res_K || !res_Na) return -1;

    /* Get crystallographic S2 ion site */
    Vec3 sites[4];
    if (kcsa_ion_sites(4, sites) < 2) return -1;
    Vec3 site_s2 = sites[1];

    /* Run K+ (Z=19) */
    if (side_fix_kcsa_sample_stabilized(19, site_s2, temp_K, n_equil, n_sample,
                                        424242UL, epsilon_pore, res_K) != 0) {
        return -1;
    }

    /* Run Na+ (Z=11) */
    if (side_fix_kcsa_sample_stabilized(11, site_s2, temp_K, n_equil, n_sample,
                                        848484UL, epsilon_pore, res_Na) != 0) {
        return -1;
    }

    if (!res_K->ok || !res_Na->ok) {
        return -2; /* One of the sites collapsed */
    }

    /* Asthagiri-Dixit entropic fluctuation formula:
     * Delta(Delta G_fluct) = 0.5 * beta * (sigma_Na^2 - sigma_K^2)
     * beta = 1 / (k_B * T)
     */
    double kB_eV = BOLTZMANN_K * J_TO_EV;
    double beta = 1.0 / (kB_eV * temp_K);
    *ddG_entropic_eV = 0.5 * beta * (res_Na->bind_var_eV2 - res_K->bind_var_eV2);

    return 0;
}
