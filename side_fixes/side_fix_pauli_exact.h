#ifndef SIDE_FIX_PAULI_EXACT_H
#define SIDE_FIX_PAULI_EXACT_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_pauli_exact.h
 *
 * Pauli repulsion layer:
 *  1. Replaces the empirical heuristic exp(-zeta*R) in qm_overlap_q with
 *     STO overlap integrals evaluated with charge-responsive exponents.
 *  2. Implements analytical spatial force gradients for E_Pauli = A * S^2,
 *     replacing O(N^2 * 6N) finite-difference force evaluations.
 *
 * HONESTY BLOCK:
 *  Q1 "Exact" applies ONLY to the equal-exponent 1s-1s closed form
 *     S = exp(-rho)(1 + rho + rho^2/3), rho = zeta*R. The unequal-exponent
 *     branch is a 64-point Simpson quadrature over the prolate-spheroidal
 *     mu coordinate with an analytic nu integral. Measured against a
 *     converged mpmath integral it is accurate to 1e-5 relative at
 *     R <= 2 A and 1.4e-4 at R = 3 A (where S is already ~1e-2). It is an
 *     approximation, not a closed form.
 *  Q2 The base exponents (H 1.0; C 3.14, N 3.83, O 4.45 divided by n = 2)
 *     and the linear charge response zeta = zeta0 + gamma*q (gamma 0.15
 *     for H, 0.35 otherwise) are MODEL PARAMETERS. They resemble Slater /
 *     Clementi-Raimondi effective charges, but this repository carries no
 *     citation or fit for them. Do not present them as a published set.
 *  Q3 The pair prefactor A = (J_a + J_b)/2 uses qm_chi_J, the engine's
 *     (IE - EA)/2 hardness. That is a model choice for the Pauli scale,
 *     not a derived result.
 */

/* Computes exact 1s-1s STO overlap and its radial derivative dS/dR */
void side_fix_sto_overlap_1s1s(double zeta_a, double zeta_b, double R_ang,
                               double *S_out, double *dS_dR_out);

/* Computes exact charge-responsive STO valence overlap between two atoms */
double side_fix_qm_overlap_exact_q(const Atom *a, const Atom *b, double R_ang,
                                   double qa, double qb, double *dS_dR_out);

/* Evaluates total Pauli repulsion energy using exact STO overlaps */
double side_fix_qm_pauli_energy_exact(const Simulation *sim);

/* Evaluates analytical Pauli repulsion forces (accumulates into sim->atoms[i].force) */
double side_fix_qm_pauli_forces_analytic(Simulation *sim);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_PAULI_EXACT_H */
