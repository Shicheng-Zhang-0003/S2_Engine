#ifndef SIDE_FIX_QEQ_SCREENED_H
#define SIDE_FIX_QEQ_SCREENED_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_qeq_screened.h
 *
 * Physical and mathematical fix for QEq charge equilibration:
 * Replaces the bare 1/r unbonded electrostatic coupling with the
 * Ohno-Klopman / Rappé-Goddard screened Coulomb integral:
 *
 *    J_ij(R) = COULOMB_MD / sqrt(R^2 + d_ij^2)
 *    d_ij    = 2 * COULOMB_MD / (J_i + J_j)
 *
 * Mathematical properties:
 *  1. At R -> 0: J_ii(0) = J_i (exact atomic hardness self-coupling).
 *  2. For all R >= 0: J_ij(R) < (J_i + J_j)/2 <= max(J_i, J_j).
 *  3. The hardness matrix A is strictly positive-definite for all geometries.
 *  4. Eliminates the indefinite matrix bug where bare 1/r exceeded J at R < 2.2 A,
 *     which caused the unconstrained saddle-point solve to put positive charge
 *     on oxygen in O...H pairs.
 *  5. Strictly conserves total charge sum(q_i) = total_q.
 */

/* Screened Coulomb coupling between two atom sites */
double side_fix_qeq_coulomb_screened(double R_ang, double J_i_eV, double J_j_eV, double dielectric);

/* Screened QEq charge solver */
int side_fix_qm_qeq_screened(const Simulation *sim, double total_q, double dielectric, double *out_q);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_QEQ_SCREENED_H */
