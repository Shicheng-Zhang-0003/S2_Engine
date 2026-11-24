#ifndef SIDE_FIX_KCSA_STABILIZED_H
#define SIDE_FIX_KCSA_STABILIZED_H

#include "../include/types.h"
#include "../include/qm_fluctuation.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_kcsa_stabilized.h
 *
 * Stabilisation layer for the 4-fold KcsA TVGYG filter fragment.
 *
 * Why a collapse was reported:
 *  1. The engine's earlier KcsA work (docs/kcsA.md) reported the flexible
 *     filter collapsing under an unscreened +1e ion. A static estimate for
 *     eight -0.55e carbonyl oxygens at 2.78 A gives > 8 eV/A inward on a
 *     164-atom TVGYG fragment with no scaffold.
 *  2. A real filter is held open by the surrounding pore-helix bundle and
 *     sits in a hydrated, low-dielectric cavity. Both effects are modelled
 *     here as (a) an isotropic relative permittivity epsilon_pore and (b)
 *     C4-symmetric C-alpha harmonic tethers.
 *
 * HONESTY BLOCK - read before quoting any number from this module:
 *
 *  K0 THIS SUITE DOES NOT REPRODUCE THE ORIGINAL COLLAPSE. Measured with
 *     this module's own setup after 200 fs: epsilon=1, tethers OFF gives
 *     K CN=8.00, Na CN=7.05, <O-O>=3.95 A; epsilon=1, tethers ON gives
 *     K CN=8.00, Na CN=7.47; epsilon=4, tethers ON gives K/Na CN=8.00.
 *     Nothing collapses in the null configuration, so the tests below
 *     establish only that the stabilized model stays coordinated in this
 *     window. They do NOT establish that the module fixes the engine's
 *     original defect. That would require reproducing the original demo
 *     setup, which this harness does not do.
 *
 *  K1 epsilon_pore is a MODEL CHOICE, NOT a sourced constant. The code no
 *     longer defaults it: the caller must pass a value >= 1.0. The values
 *     used in tests (4.0) are the usual continuum-electrostatics range for
 *     a protein interior (2-4); this repository contains no experimental
 *     measurement that fixes 4.0. The module's own null run (K0) shows that
 *     epsilon=1 does not collapse either, so the suite cannot attribute
 *     stability to the dielectric. The measured dielectric effect in the
 *     no-tether probe is modest and one-sided (Na CN 7.05 -> 7.71 from
 *     epsilon=1 to 4; K stays 8.00).
 *
 *  K2 The scaffold tethers maintain CN ~ 8 and <O-O> ~ 3.9-4.1 A through
 *     about 0.4 ps. By 0.8 ps the site starts to dewet (CN falls to ~7.3,
 *     <O-O> rises to ~4.3 A). The earlier header claimed "CN >= 7.5 over
 *     picosecond trajectories"; that is not what the model does at 1 ps
 *     and is withdrawn. The suite exercises the 0.2 ps window.
 *
 *  K3 The Asthagiri-Dixit fluctuation term is computed, but it is NOT
 *     converged at accessible sampling: the variance estimate drifts by
 *     O(1) eV between 50 fs and 800 fs windows (the binding-energy series
 *     is dominated by slow filter breathing, not by a stationary
 *     ensemble). Treat Delta(Delta G_entropic) as UNRESOLVED and quote the
 *     accompanying variance/SEM, not the central value.
 */

typedef struct {
    int    ok;                  /* 1 = site maintained & sampled, 0 = failed */
    double cn_mean;             /* mean coordination number across trajectory */
    double oo_mean_A;           /* mean O-O diagonal distance */
    double bind_mean_eV;        /* mean ion binding energy */
    double bind_var_eV2;        /* binding energy variance (sigma^2) */
    double bind_sem_var;        /* standard error on variance */
    int    n_samples;
} SideFixKcsaSeriesResult;

/* Runs stabilized KcsA sampling with cavity dielectric screening and scaffold restraints */
int side_fix_kcsa_sample_stabilized(int ion_Z, Vec3 site, double temp_K,
                                    int n_equil, int n_sample, unsigned long seed,
                                    double epsilon_pore,
                                    SideFixKcsaSeriesResult *out);

/* Compares K+ vs Na+ entropic selectivity Delta(0.5 * beta * sigma^2) */
int side_fix_kcsa_selectivity_entropic(double temp_K, int n_equil, int n_sample,
                                       double epsilon_pore,
                                       double *ddG_entropic_eV,
                                       SideFixKcsaSeriesResult *res_K,
                                       SideFixKcsaSeriesResult *res_Na);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_KCSA_STABILIZED_H */
