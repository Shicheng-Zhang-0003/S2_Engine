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
#define MAX_ITER 500
#define SCF_TOL 1e-6
/* CONVERGENCE IS NOT ONE SMALL STEP.
 *
 * A single max_dp < SCF_TOL test is unsafe. In the pre-fix hybrid model the
 * damped iteration approached an UNSTABLE antibonding fixed point, parked
 * near it for a few iterations, then departed for the physical solution. At
 * R = 2.50 A the trajectory ran
 *
 *     it=21 max_dp=7.5e-07   it=22 4.3e-07   it=23 3.6e-07
 *     it=24 max_dp=5.9e-07   it=25 1.5e-06 -> rising again
 *
 * and the single test fired at it=21 with a spurious P_01 = -1.535. The
 * orthogonal CNDO/2 secular problem no longer has that pathological root,
 * but the window is retained: it costs a few iterations and rejects exactly
 * the parked-but-about-to-depart case if a future parameter change
 * reintroduces an unstable fixed point. */
#define SCF_STABLE_WINDOW 5

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

/* Assemble the CNDO/2 Fock matrix from a density P. Factored out so the
 * post-SCF energy-identity check rebuilds F with exactly the same code path
 * that produced it during the iteration. */
static void build_fock(const Simulation *sim, int Nat, int nao,
                       const int *ao_atom,
                       const double H[MAX_AO][MAX_AO],
                       const double P[MAX_AO][MAX_AO],
                       double F[MAX_AO][MAX_AO]) {
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
            /* Slater exponents are tabulated in atomic units; eht_overlap_ao()
             * takes zeta in A^-1. The old code passed 1.2 / 1.625 / 1.950 /
             * 2.275 straight through, i.e. every orbital was 1.889x too
             * diffuse (zeta_eff = zeta_a.u. * a0). The H2 curve and the
             * water HOMO both move when this is corrected. */
            ao_zeta[nao] = 1.2 / BOHR_TO_ANGSTROM; nao++;
        } else if (Z == 6 || Z == 7 || Z == 8) {
            double zeta = ((Z == 6) ? 1.625 : (Z == 7) ? 1.950 : 2.275)
                          / BOHR_TO_ANGSTROM;
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
    /* CLOSED-SHELL ONLY. With n_occ = n_elec / 2 an odd electron count used
     * to lose its electron silently: an isolated H atom became n_occ = 0,
     * P = 0, E = 0 and returned "converged". That is a wrong answer, not a
     * converged one, and it invalidates any atomization energy built from
     * such a fragment. Refuse instead. */
    if (n_elec % 2 != 0) return -5;
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

    /* CNDO/2 SECULAR PROBLEM: S = I.
     *
     * CNDO/2 is a zero-differential-overlap model. The atomic orbitals are
     * treated as an ORTHOGONAL basis in the secular equations, and the real
     * overlap S_mn survives only in the resonance integral
     * H_mn = beta_AB S_mn (above). The previous code solved the generalized
     * problem F C = S C eps with the EXACT S and also fed the resulting
     * non-orthogonal P into an expression derived for an orthogonal basis.
     * That mixed two models: Tr(P S) = n_elec held, but the two-electron
     * terms saw Tr(P) < n_elec, the H2 symmetric solution had P_00 = 0.571
     * instead of 1, and the H2 curve came out unbound with a spurious root
     * flip at R ~ 2.5 A.
     *
     * With W = I the Roothaan step reduces to the ordinary eigenproblem,
     * the density is normalised by Tr(P) = n_elec, and the CNDO/2 energy
     * decomposition below is the matching expression. This is the standard
     * CNDO/2 convention (Pople & Segal 1966, eqs. 13-16). */
    double W[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
            W[i][j] = (i == j) ? 1.0 : 0.0;

    double P[MAX_AO][MAX_AO];
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) P[i][j] = 0.0;

    double F[MAX_AO][MAX_AO];
    double C[MAX_AO][MAX_AO];
    double eps[MAX_AO];
    int converged = 0;
    int stable = 0;   /* consecutive iterations that stayed under SCF_TOL */
    int it = 0;

    for (it = 0; it < MAX_ITER; it++) {
        build_fock(sim, Nat, nao, ao_atom, H, P, F);

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

        /* Stability window, not a single test - see SCF_STABLE_WINDOW. */
        if (max_dp < SCF_TOL) {
            if (++stable >= SCF_STABLE_WINDOW) {
                converged = 1;
                break;
            }
        } else {
            stable = 0;
        }
    }

    /* Pople CNDO/2 total energy decomposition:
     * E_tot = sum_A E_A + sum_{A < B} E_AB
     */
    double E_tot = 0.0;
    /* One-center terms.
     *
     * ENERGY IDENTITY: for a closed-shell SCF the total energy must equal
     *   E = 1/2 * sum_mn P_mn (H_mn + F_mn) + V_nn
     * for ANY density P, not just a converged one. Both terms below are
     * required for that to hold; each was verified to be missing by
     * evaluating the identity against an independently rebuilt Fock matrix
     * on the same P (see test_rhf_energy_identity).
     *
     * The second line - the SAME-ATOM off-diagonal exchange - was absent.
     * It contributes -1/2 * gamma_AA * sum_{m<n on A} P_mn^2, because
     * F_mn = -1/2 P_mn gamma_AA for m != n on one atom, so
     *   1/2 sum_{m!=n} P_mn F_mn = -1/4 gamma_AA * 2 sum_{m<n} P_mn^2.
     * For H2 it is identically zero (one AO per atom), which is why the
     * defect only surfaced on water, where it is worth 4.87 eV. */
    for (int a = 0; a < Nat; a++) {
        int Za = sim->atoms[a].Z;
        double P_AA = 0.0;
        double sum_p_diag = 0.0;
        double sum_off = 0.0;
        double e_a = 0.0;
        int ao_on_a[MAX_AO];
        int n_on_a = 0;
        for (int m = 0; m < nao; m++) {
            if (ao_atom[m] != a) continue;
            ao_on_a[n_on_a++] = m;
            P_AA += P[m][m];
            sum_p_diag += P[m][m] * P[m][m];
            e_a += P[m][m] * get_u_core(Za, ao_l[m]);
        }
        for (int x = 0; x < n_on_a; x++)
            for (int y = x + 1; y < n_on_a; y++)
                sum_off += P[ao_on_a[x]][ao_on_a[y]] * P[ao_on_a[x]][ao_on_a[y]];
        e_a += 0.5 * (P_AA * P_AA - 0.5 * sum_p_diag) * get_gamma_self(Za);
        e_a -= 0.5 * sum_off * get_gamma_self(Za);
        E_tot += e_a;
    }

    /* Two-center terms */
    double V_nn = 0.0;
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
            double z_a = (double)get_val_e(Za);
            double z_b = (double)get_val_e(Zb);
            V_nn += z_a * z_b * (COULOMB_MD / R);

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

            /* 2. Exchange correlation energy.
             *
             * COEFFICIENT IS -1/2, NOT -1/4. From the identity
             *   1/2 sum_mn P_mn (H_mn + F_mn)
             * the off-diagonal (m on A, n on B) contribution is
             *   P_mn H_mn + P_mn F_mn = 2 P_mn H_mn - 1/2 P_mn^2 gamma_AB,
             * since F_mn = H_mn - 1/2 P_mn gamma_AB. The resonance term
             * below already carries the 2 P_mn H_mn, so the exchange term
             * must carry -1/2. With -1/4 every total energy came out too
             * high by exactly 1/4 * P_mn^2 * gamma_AB (measured: 0.722582
             * eV at R = 2.50 A vs 0.72256 predicted). */
            double e_ex = -0.5 * p_ab_sq * g_ab;

            /* 3. Electrostatic interaction between net charges */
            double e_elec = Qa * Qb * g_ab;

            /* 4. Core-core repulsion deficit: Z_A * Z_B * (14.4/R - g_ab) */
            double e_core = z_a * z_b * (COULOMB_MD / R - g_ab);

            E_tot += e_res + e_ex + e_elec + e_core;
        }
    }

    /* ENERGY-IDENTITY GATE.
     *
     * For a closed-shell SCF, E = 1/2 Tr[P(H+F)] + V_nn holds for ANY
     * density P, not only at convergence. Rebuild F from the final P and
     * evaluate the trace independently of the decomposition above. If the
     * two disagree the decomposition is wrong (this is how the missing
     * same-atom exchange term and the -1/4 coefficient were found), and a
     * number that violates the identity must not be returned. */
    double F_check[MAX_AO][MAX_AO];
    build_fock(sim, Nat, nao, ao_atom, H, P, F_check);
    double E_trace = 0.0;
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++)
            E_trace += P[i][j] * (H[i][j] + F_check[i][j]);
    E_trace = 0.5 * E_trace + V_nn;
    if (!isfinite(E_trace) || fabs(E_trace - E_tot) > 1e-6) return -6;

    if (out) {
        out->converged = converged;
        out->n_iters = it + 1;
        out->e_total_eV = E_tot;
        out->e_nuclear_eV = V_nn;
        out->e_electronic_eV = E_tot - V_nn;
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
    /* An unconverged solve has no energy. Returning the last iterate's
     * number invites a caller to scan it as if it were physics, so return
     * NAN and let the arithmetic fail loudly. */
    return (rc == 0) ? res.e_total_eV : NAN;
}
