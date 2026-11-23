#include "side_fix_rhf_calibrated.h"
#include "../include/constants.h"
#include "../include/sim.h"
#include "../include/qm.h"
#include "../include/qm_eht.h"
#include "../include/vec3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_AO 32
#define MAX_ITER 100
#define SCF_TOL 1e-6

static int jacobi_diag(double A[MAX_AO][MAX_AO], double V[MAX_AO][MAX_AO], double d[MAX_AO], int n) {
    for (int i = 0; i < n; i++) {
        d[i] = A[i][i];
        for (int j = 0; j < n; j++) V[i][j] = (i == j) ? 1.0 : 0.0;
    }
    for (int sweep = 0; sweep < 150; sweep++) {
        double off = 0.0;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
        if (off < 1e-24) return 0;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = A[p][q];
                if (fabs(apq) < 1e-15) continue;
                double app = A[p][p], aqq = A[q][q];
                double phi = 0.5 * (aqq - app) / apq;
                double t = ((phi >= 0.0) ? 1.0 : -1.0) / (fabs(phi) + sqrt(phi * phi + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0);
                double s = t * c;
                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = A[q][p] = 0.0;
                for (int k = 0; k < n; k++) {
                    if (k == p || k == q) continue;
                    double akp = A[k][p], akq = A[k][q];
                    A[k][p] = A[p][k] = c * akp - s * akq;
                    A[k][q] = A[q][k] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {
                    double vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
                d[p] = A[p][p];
                d[q] = A[q][q];
            }
        }
    }
    return 0;
}

static double get_beta0(int Z) {
    switch (Z) {
        case 1:  return -12.0;  /* H  */
        case 6:  return -21.0;  /* C  */
        case 7:  return -25.0;  /* N  */
        case 8:  return -31.0;  /* O  */
        default: return -15.0;
    }
}

static double get_gamma_self(int Z) {
    switch (Z) {
        case 1:  return 12.85;  /* H (1s-1s) */
        case 6:  return 11.11;  /* C (2s-2s) */
        case 7:  return 12.01;  /* N (2s-2s) */
        case 8:  return 13.00;  /* O (2s-2s) */
        default: return 10.0;
    }
}

static double get_gamma_ab(int Za, int Zb, double R_ang) {
    double g_a = get_gamma_self(Za);
    double g_b = get_gamma_self(Zb);
    double g_avg = 0.5 * (g_a + g_b);
    double d = COULOMB_MD / g_avg;
    return COULOMB_MD / sqrt(R_ang * R_ang + d * d);
}

static double get_u_core(int Z, int l) {
    if (Z == 1) return -13.606; /* H 1s */
    if (Z == 6) return (l == 0) ? -50.69 : -41.53;
    if (Z == 7) return (l == 0) ? -70.09 : -57.85;
    if (Z == 8) return (l == 0) ? -99.96 : -84.28;
    return -15.0;
}

static int get_val_e(int Z) {
    switch (Z) {
        case 1:  return 1;
        case 6:  return 4;
        case 7:  return 5;
        case 8:  return 6;
        default: return 1;
    }
}

int side_fix_qm_scf_rhf_calibrated(const Simulation *sim, SideFixScfResult *out) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return -1;
    int Nat = sim->num_atoms;

    int nao = 0;
    int ao_atom[MAX_AO], ao_n[MAX_AO], ao_l[MAX_AO], ao_m[MAX_AO];
    double ao_zeta[MAX_AO];

    for (int a = 0; a < Nat && nao < MAX_AO; a++) {
        int Z = sim->atoms[a].Z;
        if (Z == 1) {
            ao_atom[nao] = a; ao_n[nao] = 1; ao_l[nao] = 0; ao_m[nao] = 0;
            ao_zeta[nao] = 1.2; nao++;
        } else if (Z == 6 || Z == 7 || Z == 8) {
            double zeta = (Z == 6) ? 1.625 : (Z == 7) ? 1.950 : 2.275;
            ao_atom[nao] = a; ao_n[nao] = 2; ao_l[nao] = 0; ao_m[nao] = 0;
            ao_zeta[nao] = zeta; nao++;
            for (int m = 1; m <= 3 && nao < MAX_AO; m++) {
                ao_atom[nao] = a; ao_n[nao] = 2; ao_l[nao] = 1; ao_m[nao] = m;
                ao_zeta[nao] = zeta; nao++;
            }
        }
    }
    if (nao < 1) return -1;

    int n_elec = 0;
    for (int a = 0; a < Nat; a++) n_elec += get_val_e(sim->atoms[a].Z);
    int n_occ = n_elec / 2;

    double S[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++) {
        S[i][i] = 1.0;
        for (int j = i + 1; j < nao; j++) {
            int a = ao_atom[i], b = ao_atom[j];
            double s_val = 0.0;
            if (a == b) {
                s_val = (i == j) ? 1.0 : 0.0;
            } else {
                Vec3 ra = sim->atoms[a].position;
                Vec3 rb = sim->atoms[b].position;
                Vec3 d = vec3_sub(rb, ra);
                double R = vec3_norm(d);
                if (R < 1e-12) {
                    s_val = 0.0;
                } else {
                    double dx = d.x / R, dy = d.y / R, dz = d.z / R;
                    s_val = eht_overlap_ao(ao_n[i], ao_zeta[i], ao_l[i], ao_m[i],
                                           ao_n[j], ao_zeta[j], ao_l[j], ao_m[j],
                                           R, dx, dy, dz);
                }
            }
            S[i][j] = S[j][i] = s_val;
        }
    }

    /* Core Hamiltonian H */
    double H[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++) {
        int a = ao_atom[i];
        int Za = sim->atoms[a].Z;
        double h_diag = get_u_core(Za, ao_l[i]);
        for (int b = 0; b < Nat; b++) {
            if (b == a) continue;
            int Zb = sim->atoms[b].Z;
            double R = vec3_dist(sim->atoms[a].position, sim->atoms[b].position);
            h_diag -= (double)get_val_e(Zb) * get_gamma_ab(Za, Zb, R);
        }
        H[i][i] = h_diag;

        for (int j = i + 1; j < nao; j++) {
            int b = ao_atom[j];
            int Zb = sim->atoms[b].Z;
            if (a == b) {
                H[i][j] = H[j][i] = 0.0;
            } else {
                double beta_avg = 0.5 * (get_beta0(Za) + get_beta0(Zb));
                H[i][j] = H[j][i] = beta_avg * S[i][j];
            }
        }
    }

    double W[MAX_AO][MAX_AO], Vmat[MAX_AO][MAX_AO], lam[MAX_AO], Sm[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) Sm[i][j] = S[i][j];

    if (jacobi_diag(Sm, Vmat, lam, nao) != 0) return -2;
    for (int i = 0; i < nao; i++) {
        if (lam[i] <= 1e-10) return -2;
    }
    for (int i = 0; i < nao; i++) {
        for (int j = 0; j < nao; j++) {
            double sum = 0.0;
            for (int k = 0; k < nao; k++) {
                sum += Vmat[i][k] * Vmat[j][k] / sqrt(lam[k]);
            }
            W[i][j] = sum;
        }
    }

    double P[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) P[i][j] = 0.0;

    double F[MAX_AO][MAX_AO];
    double C[MAX_AO][MAX_AO];
    double eps[MAX_AO];
    int converged = 0;
    int it = 0;

    for (it = 0; it < MAX_ITER; it++) {
        for (int i = 0; i < nao; i++) {
            int a = ao_atom[i];
            int Za = sim->atoms[a].Z;
            double P_AA = 0.0;
            for (int k = 0; k < nao; k++) {
                if (ao_atom[k] == a) P_AA += P[k][k];
            }
            double f_ii = H[i][i] + (P_AA - 0.5 * P[i][i]) * get_gamma_self(Za);
            for (int b = 0; b < Nat; b++) {
                if (b == a) continue;
                int Zb = sim->atoms[b].Z;
                double R = vec3_dist(sim->atoms[a].position, sim->atoms[b].position);
                double g_ab = get_gamma_ab(Za, Zb, R);
                double P_BB = 0.0;
                for (int k = 0; k < nao; k++) {
                    if (ao_atom[k] == b) P_BB += P[k][k];
                }
                f_ii += P_BB * g_ab;
            }
            F[i][i] = f_ii;

            for (int j = i + 1; j < nao; j++) {
                int b = ao_atom[j];
                int Zb = sim->atoms[b].Z;
                if (a == b) {
                    F[i][j] = F[j][i] = -0.5 * P[i][j] * get_gamma_self(Za);
                } else {
                    double R = vec3_dist(sim->atoms[a].position, sim->atoms[b].position);
                    double g_ab = get_gamma_ab(Za, Zb, R);
                    F[i][j] = F[j][i] = H[i][j] - 0.5 * P[i][j] * g_ab;
                }
            }
        }

        double Fp[MAX_AO][MAX_AO];
        for (int i = 0; i < nao; i++) {
            for (int j = 0; j < nao; j++) {
                double s = 0.0;
                for (int k = 0; k < nao; k++)
                    for (int l = 0; l < nao; l++)
                        s += W[i][k] * F[k][l] * W[l][j];
                Fp[i][j] = s;
            }
        }

        double Cp[MAX_AO][MAX_AO];
        if (jacobi_diag(Fp, Cp, eps, nao) != 0) return -3;

        int idx[MAX_AO];
        for (int i = 0; i < nao; i++) idx[i] = i;
        for (int i = 0; i < nao - 1; i++) {
            int best = i;
            for (int j = i + 1; j < nao; j++) {
                if (eps[idx[j]] < eps[idx[best]]) best = j;
            }
            if (best != i) {
                int t = idx[i]; idx[i] = idx[best]; idx[best] = t;
            }
        }
        double eps_sorted[MAX_AO];
        for (int i = 0; i < nao; i++) eps_sorted[i] = eps[idx[i]];
        for (int i = 0; i < nao; i++) eps[i] = eps_sorted[i];

        for (int i = 0; i < nao; i++) {
            for (int j = 0; j < nao; j++) {
                double sum = 0.0;
                for (int k = 0; k < nao; k++) sum += W[i][k] * Cp[k][idx[j]];
                C[i][j] = sum;
            }
        }

        double P_new[MAX_AO][MAX_AO];
        for (int i = 0; i < nao; i++) {
            for (int j = 0; j < nao; j++) {
                double sum = 0.0;
                for (int k = 0; k < n_occ; k++) sum += 2.0 * C[i][k] * C[j][k];
                P_new[i][j] = sum;
            }
        }

        double max_dp = 0.0;
        for (int i = 0; i < nao; i++) {
            for (int j = 0; j < nao; j++) {
                double dp = fabs(P_new[i][j] - P[i][j]);
                if (dp > max_dp) max_dp = dp;
                P[i][j] = 0.5 * P[i][j] + 0.5 * P_new[i][j];
            }
        }

        if (max_dp < SCF_TOL) {
            converged = 1;
            break;
        }
    }

    /* Pople CNDO/2 total energy decomposition:
     * E_tot = sum_A E_A + sum_{A < B} E_AB
     */
    double E_tot = 0.0;
    /* One-center terms */
    for (int a = 0; a < Nat; a++) {
        int Za = sim->atoms[a].Z;
        double P_AA = 0.0;
        double sum_p_diag = 0.0;
        double e_a = 0.0;
        for (int m = 0; m < nao; m++) {
            if (ao_atom[m] != a) continue;
            P_AA += P[m][m];
            sum_p_diag += P[m][m] * P[m][m];
            e_a += P[m][m] * get_u_core(Za, ao_l[m]);
        }
        e_a += 0.5 * (P_AA * P_AA - 0.5 * sum_p_diag) * get_gamma_self(Za);
        E_tot += e_a;
    }

    /* Two-center terms */
    for (int a = 0; a < Nat; a++) {
        int Za = sim->atoms[a].Z;
        double P_AA = 0.0;
        for (int m = 0; m < nao; m++) if (ao_atom[m] == a) P_AA += P[m][m];
        double Qa = (double)get_val_e(Za) - P_AA;

        for (int b = a + 1; b < Nat; b++) {
            int Zb = sim->atoms[b].Z;
            double P_BB = 0.0;
            for (int n = 0; n < nao; n++) if (ao_atom[n] == b) P_BB += P[n][n];
            double Qb = (double)get_val_e(Zb) - P_BB;

            double R = vec3_dist(sim->atoms[a].position, sim->atoms[b].position);
            double g_ab = get_gamma_ab(Za, Zb, R);
            double beta_avg = 0.5 * (get_beta0(Za) + get_beta0(Zb));

            /* 1. Resonance (covalent bonding) energy */
            double e_res = 0.0;
            double p_ab_sq = 0.0;
            for (int m = 0; m < nao; m++) {
                if (ao_atom[m] != a) continue;
                for (int n = 0; n < nao; n++) {
                    if (ao_atom[n] != b) continue;
                    e_res += 2.0 * P[m][n] * beta_avg * S[m][n];
                    p_ab_sq += P[m][n] * P[m][n];
                }
            }

            /* 2. Exchange correlation energy */
            double e_ex = -0.25 * p_ab_sq * g_ab;

            /* 3. Electrostatic interaction between net charges */
            double e_elec = Qa * Qb * g_ab;

            /* 4. Core-core repulsion deficit: Z_A * Z_B * (14.4/R - g_ab) */
            double e_core = (double)get_val_e(Za) * (double)get_val_e(Zb) * (COULOMB_MD / R - g_ab);

            E_tot += e_res + e_ex + e_elec + e_core;
        }
    }

    if (out) {
        out->converged = converged;
        out->n_iters = it + 1;
        out->e_total_eV = E_tot;
        out->homo_eV = (n_occ > 0) ? eps[n_occ - 1] : 0.0;
        out->lumo_eV = (n_occ < nao) ? eps[n_occ] : 0.0;
        out->gap_eV = out->lumo_eV - out->homo_eV;
    }

    return converged ? 0 : -3;
}

double side_fix_h2_curve(double R_ang, int *converged) {
    Simulation *s = sim_create(4, 4);
    sim_add_atom(s, 1, vec3(0, 0, 0), 0.0);
    sim_add_atom(s, 1, vec3(0, 0, R_ang), 0.0);

    SideFixScfResult res;
    int rc = side_fix_qm_scf_rhf_calibrated(s, &res);
    sim_destroy(s);

    if (converged) *converged = (rc == 0);
    return res.e_total_eV;
}
