#include "side_fix_pauli_exact.h"
#include "../include/constants.h"
#include "../include/quantum.h"
#include "../include/qm.h"
#include "../include/vec3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

void side_fix_sto_overlap_1s1s(double zeta_a, double zeta_b, double R_ang,
                               double *S_out, double *dS_dR_out) {
    if (R_ang < 1e-12) {
        if (S_out) *S_out = 1.0;
        if (dS_dR_out) *dS_dR_out = 0.0;
        return;
    }

    double a0 = BOHR_TO_ANGSTROM;
    /* Convert R to Bohr, zeta to a0^-1 */
    double R_bohr = R_ang / a0;
    double za = zeta_a;
    double zb = zeta_b;

    double diff = fabs(za - zb);
    if (diff < 1e-5) {
        double z = 0.5 * (za + zb);
        double rho = z * R_bohr;
        double exp_rho = exp(-rho);
        double S = exp_rho * (1.0 + rho + (1.0 / 3.0) * rho * rho);
        if (S_out) *S_out = S;
        if (dS_dR_out) {
            /* dS/dR_bohr = -(1/3) * z * rho * (1 + rho) * exp(-rho) */
            double dS_dbohr = -(1.0 / 3.0) * z * rho * (1.0 + rho) * exp_rho;
            *dS_dR_out = dS_dbohr / a0; /* Convert dS/dR to A^-1 */
        }
    } else {
        /* General unequal 1s-1s STO overlap */
        double p = 0.5 * (za + zb) * R_bohr;
        double t = (za - zb) / (za + zb);

        /* Closed form helper via A and B functions */
        double S_val = 0.0;
        /* Stable Simpson integration over prolate mu in [1, 20] and nu in [-1, 1] */
        int n_mu = 64;
        double mu_max = 20.0 / p;
        if (mu_max < 2.0) mu_max = 2.0;
        double h_mu = (mu_max - 1.0) / (double)n_mu;
        double sum = 0.0;
        for (int i = 0; i <= n_mu; i++) {
            double mu = 1.0 + i * h_mu;
            double w_mu = (i == 0 || i == n_mu) ? 1.0 : (i % 2 == 1 ? 4.0 : 2.0);
            double exp_p_mu = exp(-p * mu);
            /* Int_{-1}^1 (mu^2 - nu^2) exp(-p*t*nu) dnu */
            double pt = p * t;
            double I_nu = 0.0;
            if (fabs(pt) < 1e-4) {
                I_nu = 2.0 * mu * mu - (2.0 / 3.0);
            } else {
                double e_pt = exp(pt);
                double e_mpt = exp(-pt);
                double sinh_pt = 0.5 * (e_pt - e_mpt);
                double cosh_pt = 0.5 * (e_pt + e_mpt);
                double term0 = (2.0 / pt) * sinh_pt;
                double term2 = (2.0 / pt) * sinh_pt - (4.0 / (pt * pt)) * cosh_pt + (4.0 / (pt * pt * pt)) * sinh_pt;
                I_nu = mu * mu * term0 - term2;
            }
            sum += w_mu * exp_p_mu * I_nu;
        }
        double prefactor = pow(za * zb, 1.5) * pow(R_bohr, 3.0) / 4.0 * (h_mu / 3.0);
        S_val = prefactor * sum;
        if (S_val < 0.0) S_val = 0.0;
        if (S_val > 1.0) S_val = 1.0;
        if (S_out) *S_out = S_val;

        if (dS_dR_out) {
            /* Central finite difference for radial derivative */
            double hR = 1e-5;
            double S_plus, S_minus;
            side_fix_sto_overlap_1s1s(zeta_a, zeta_b, R_ang + hR, &S_plus, NULL);
            side_fix_sto_overlap_1s1s(zeta_a, zeta_b, R_ang - hR, &S_minus, NULL);
            *dS_dR_out = (S_plus - S_minus) / (2.0 * hR);
        }
    }
}

double side_fix_qm_overlap_exact_q(const Atom *a, const Atom *b, double R_ang,
                                   double qa, double qb, double *dS_dR_out) {
    if (!a || !b || R_ang < 1e-9) {
        if (dS_dR_out) *dS_dR_out = 0.0;
        return (R_ang < 1e-9 && a && b) ? 1.0 : 0.0;
    }

    /* Orbital exponents with charge-responsive adjustment */
    double za0 = (a->Z == 1) ? 1.0 : (a->Z == 6) ? 3.14 : (a->Z == 7) ? 3.83 : 4.45;
    double zb0 = (b->Z == 1) ? 1.0 : (b->Z == 6) ? 3.14 : (b->Z == 7) ? 3.83 : 4.45;

    /* Slope gamma: cations contract, anions expand */
    double ga = (a->Z == 1) ? 0.15 : 0.35;
    double gb = (b->Z == 1) ? 0.15 : 0.35;

    double za = za0 + ga * qa;
    double zb = zb0 + gb * qb;
    if (za < 0.5) za = 0.5;
    if (zb < 0.5) zb = 0.5;

    int na = (a->Z == 1) ? 1 : 2;
    int nb = (b->Z == 1) ? 1 : 2;
    double zeta_a = za / (double)na;
    double zeta_b = zb / (double)nb;

    double S = 0.0;
    double dS_dR = 0.0;
    side_fix_sto_overlap_1s1s(zeta_a, zeta_b, R_ang, &S, &dS_dR);

    if (dS_dR_out) *dS_dR_out = dS_dR;
    return S;
}

double side_fix_qm_pauli_energy_exact(const Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    double E_pauli = 0.0;
    for (int i = 0; i < sim->num_atoms - 1; i++) {
        for (int j = i + 1; j < sim->num_atoms; j++) {
            Vec3 r_ij = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
            double R = vec3_norm(r_ij);
            if (R < 1e-6 || (sim->cutoff > 0.0 && R > sim->cutoff)) continue;

            double qa = sim->atoms[i].partial_charge;
            double qb = sim->atoms[j].partial_charge;
            double S = side_fix_qm_overlap_exact_q(&sim->atoms[i], &sim->atoms[j], R, qa, qb, NULL);
            if (S <= 1e-12) continue;

            double chi_a, J_a, chi_b, J_b;
            qm_chi_J(sim->atoms[i].element, &chi_a, &J_a);
            qm_chi_J(sim->atoms[j].element, &chi_b, &J_b);
            double A_pair = 0.5 * (J_a + J_b);
            E_pauli += A_pair * S * S;
        }
    }
    return E_pauli;
}

double side_fix_qm_pauli_forces_analytic(Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 2) return 0.0;
    double E_pauli = 0.0;
    for (int i = 0; i < sim->num_atoms - 1; i++) {
        for (int j = i + 1; j < sim->num_atoms; j++) {
            Vec3 r_ij = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
            double R = vec3_norm(r_ij);
            if (R < 1e-6 || (sim->cutoff > 0.0 && R > sim->cutoff)) continue;

            double qa = sim->atoms[i].partial_charge;
            double qb = sim->atoms[j].partial_charge;
            double dS_dR = 0.0;
            double S = side_fix_qm_overlap_exact_q(&sim->atoms[i], &sim->atoms[j], R, qa, qb, &dS_dR);
            if (S <= 1e-12) continue;

            double chi_a, J_a, chi_b, J_b;
            qm_chi_J(sim->atoms[i].element, &chi_a, &J_a);
            qm_chi_J(sim->atoms[j].element, &chi_b, &J_b);
            double A_pair = 0.5 * (J_a + J_b);
            E_pauli += A_pair * S * S;

            /* Analytical force:
             * E = A * S^2
             * dE/dR = 2 * A * S * (dS/dR)
             * F_i = - dE/d(r_i) = + dE/dR * (r_ij / R) = 2 * A * S * (dS/dR) * (r_ij / R)
             * (Since dS/dR < 0, F_i points along -r_ij, i.e. repelling away from atom j).
             */
            Vec3 u_ij = vec3_scale(r_ij, 1.0 / R);
            double dE_dR = 2.0 * A_pair * S * dS_dR;
            Vec3 F_on_i = vec3_scale(u_ij, dE_dR);
            Vec3 F_on_j = vec3_scale(u_ij, -dE_dR);

            sim->atoms[i].force = vec3_add(sim->atoms[i].force, F_on_i);
            sim->atoms[j].force = vec3_add(sim->atoms[j].force, F_on_j);
        }
    }
    return E_pauli;
}
