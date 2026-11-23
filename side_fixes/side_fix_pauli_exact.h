#ifndef SIDE_FIX_PAULI_EXACT_H
#define SIDE_FIX_PAULI_EXACT_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_pauli_exact.h
 *
 * Physical and computational fix for the Pauli repulsion layer:
 *  1. Replaces the empirical heuristic exp(-zeta*R) in qm_overlap_q with
 *     exact Slater-Type Orbital (STO) overlap integrals evaluated with
 *     charge-responsive orbital exponents zeta(q).
 *  2. Implements exact analytical spatial force gradients for E_Pauli = A * S^2,
 *     replacing expensive O(N^2 * 6N) finite-difference force evaluations with
 *     direct analytical derivatives.
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
