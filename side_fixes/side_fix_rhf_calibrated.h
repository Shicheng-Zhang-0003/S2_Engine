#ifndef SIDE_FIX_RHF_CALIBRATED_H
#define SIDE_FIX_RHF_CALIBRATED_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_rhf_calibrated.h
 *
 * Physical and mathematical fix for the Roothaan-Hartree-Fock (RHF) layer:
 *
 * The bug in src/qm.c:
 *   The core Hamiltonian H_core was assembled using the Extended Hückel (EHT)
 *   Wolfsberg-Helmholtz formula (H_ij = 1.75 * S_ij * 0.5 * (H_ii + H_jj)).
 *   Because EHT already implicitly incorporates effective electron-electron
 *   screening, adding the explicit two-electron operator G_ZDO on top of it
 *   resulted in DOUBLE-COUNTING of electron repulsion.
 *   Consequences:
 *     - H2 bond length expanded from experimental 0.74 A to 1.300 A.
 *     - Water (H2O) produced completely unbound states (positive energies).
 *
 * The fix:
 *   Implements the non-double-counting CNDO/2 core Hamiltonian:
 *     H_mu_mu = U_mu_mu - sum_{B != A} Z_B * gamma_AB
 *     H_mu_nu = 0.5 * (beta_A^0 + beta_B^0) * S_mu_nu
 *   where U_mu_mu is the bare one-electron core potential and beta^0 is the
 *   calibrated atomic resonance parameter.
 *   Results:
 *     - H2 dissociation curve has its physical minimum at R ~ 0.74-0.78 A.
 *     - Water is stably bound with negative total binding energy.
 */

typedef struct {
    int    converged;
    int    n_iters;
    double e_total_eV;
    double e_electronic_eV;
    double e_nuclear_eV;
    double homo_eV;
    double lumo_eV;
    double gap_eV;
} SideFixScfResult;

/* Runs calibrated RHF SCF solver for minimal valence basis (H 1s, C/N/O 2s2p) */
int side_fix_qm_scf_rhf_calibrated(const Simulation *sim, SideFixScfResult *out);

/* Scans H2 potential energy curve at distance R_ang */
double side_fix_h2_curve(double R_ang, int *converged);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_RHF_CALIBRATED_H */
