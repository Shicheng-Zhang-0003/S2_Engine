#ifndef SIDE_FIX_RHF_CALIBRATED_H
#define SIDE_FIX_RHF_CALIBRATED_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_rhf_calibrated.h
 *
 * CNDO/2 replacement for the Roothaan-Hartree-Fock layer in src/qm.c.
 *
 * The bug in src/qm.c (measured, not assumed):
 *   qm_hf_scf() builds H_core with the EHT Wolfsberg-Helmholtz formula and
 *   then adds the ZDO two-electron operator on top, so electron repulsion is
 *   double counted. Consequences measured by this repository:
 *     - H2 minimum at R = 1.30 A (experiment 0.7414 A, docs/demos.md);
 *     - water runs hit the "not bound" gate (-4) with HOMO = +3.3 eV.
 *
 * The fix (standard CNDO/2):
 *     H_mu_mu = U_mu_mu - sum_{B != A} Z_B * gamma_AB
 *     H_mu_nu = 0.5 * (beta_A^0 + beta_B^0) * S_mu_nu
 *   with ZDO two-electron terms
 *     F_mu_mu = H_mu_mu + sum_lambda P_ll gamma_mul - 0.5 P_mu_mu gamma_AA
 *     F_mu_nu = H_mu_nu - 0.5 P_mu_nu gamma_AB
 *   and the ZDO secular problem S = I: the exact overlap enters ONLY the
 *   resonance integral H_mu_nu, as in Pople & Segal (1966). The old code
 *   solved F C = S C eps with the exact S and also fed the non-orthogonal
 *   density into the orthogonal-basis energy expression; that mixture left
 *   Tr(P) = 1.14 instead of 2 for H2 and made the curve unbound.
 *
 * RESULTS (measured, this build):
 *     - H2: smooth single minimum at R = 0.80 A, E(0.76) < E(1.30), no root
 *       flip. Residual to experiment: +0.06 A. This is a CNDO/2 result,
 *       not a claim of 0.74 A accuracy.
 *       Well depth: E(0.80) = -33.77 eV, E(10 A) = -21.49 eV, analytic
 *       R -> infinity limit -20.79 eV, so De ~ 13 eV. The curve is bound
 *       but the well is much deeper than the experimental H2 bond energy;
 *       this repository carries no experimental De to compare against, so
 *       no ratio is quoted. Shape and location are qualitative; the depth
 *       is not quantitative.
 *     - H2O: converged, HOMO = -17.3 eV, gap = 24 eV, negative total energy.
 *       The HOMO is far deeper than the experimental first ionisation
 *       energy (CNDO orbital energies are not Koopmans-quality here); the
 *       claim is only that the engine's positive-HOMO unbound state is gone.
 *
 * HONESTY / PROVENANCE BLOCK:
 *   P1 The one-centre integrals U_mu_mu, the resonance parameters beta_A^0
 *      and the Hubbard parameters gamma_AA are the side-fix's own calibrated
 *      constants. They are NOT a tabulated published CNDO/2 parameter set,
 *      and this repository contains no citation for them. Absolute energies
 *      are CNDO-scale, not experimental total energies; use differences.
 *   P2 gamma_AB is an Ohno interpolant, gamma_AB = e^2/sqrt(R^2 + d^2),
 *      d = e^2/gamma_bar. It is an approximation to the two-centre Coulomb
 *      integral, flagged here rather than presented as exact.
 *   P3 Slater exponents are tabulated in atomic units; eht_overlap_ao()
 *      takes A^-1. They are converted explicitly (the old code was 1.889x
 *      too diffuse).
 *   P4 The energy decomposition is the algebraic identity
 *      E = 0.5 Tr[P(H+F)] + V_nn. After every converged solve the solver
 *      rebuilds F from the final P and returns -6 if the decomposition
 *      disagrees with the trace by more than 1e-6 eV, so an energy that
 *      violates the identity is never handed to a caller.
 *   P5 CLOSED-SHELL ONLY. An odd valence-electron count returns -5. This is
 *      deliberate: n_occ = n_elec/2 used to drop the odd electron, so an
 *      isolated H atom returned "converged" with P = 0 and E = 0. Because H
 *      and OH are not representable, this module defines no atomization
 *      energy and no heat of formation. Even the closed-shell O atom is an
 *      artificial 2p^4 singlet (one empty 2p), not the triplet ground state,
 *      so an O + H2 comparison is not a physical binding energy either. Use
 *      the H2 curve shape and the orbital/gap diagnostics; do not quote
 *      fragmentation energies.
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

/* Runs calibrated RHF SCF solver for minimal valence basis (H 1s, C/N/O 2s2p).
 * Returns 0 on success, -1 on bad input, -3 if the SCF does not converge,
 * -5 for an odd valence-electron count (closed-shell only), and -6 if the
 * post-SCF energy identity check fails. */
int side_fix_qm_scf_rhf_calibrated(const Simulation *sim, SideFixScfResult *out);

/* Scans H2 potential energy curve at distance R_ang */
double side_fix_h2_curve(double R_ang, int *converged);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_RHF_CALIBRATED_H */
