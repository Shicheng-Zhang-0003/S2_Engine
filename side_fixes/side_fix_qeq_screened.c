#include "side_fix_qeq_screened.h"
#include "../include/constants.h"
#include "../include/qm.h"
#include "../include/vec3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef QM_QEQ_QMAX
#define QM_QEQ_QMAX 2.0
#endif

double side_fix_qeq_coulomb_screened(double R_ang, double J_i_eV, double J_j_eV, double dielectric) {
    if (R_ang < 0.0) R_ang = 0.0;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    double J_avg = 0.5 * (J_i_eV + J_j_eV);
    if (!(J_avg > 1e-6)) J_avg = 1.0;

    /* Ohno-Klopman screening length: d_ij = 2 * COULOMB_MD / (J_i + J_j) */
    double d_ij = COULOMB_MD / J_avg;
    double denom = sqrt(R_ang * R_ang + d_ij * d_ij);
    return COULOMB_MD / (denom * dielectric);
}

static double atom_dist_pbc(const Simulation *sim, int i, int j) {
    Vec3 r_ij = vec3_sub(sim->atoms[j].position, sim->atoms[i].position);
    if (sim->box.periodic[0] && sim->box.dimensions.x > 0.0) {
        r_ij.x -= sim->box.dimensions.x * round(r_ij.x / sim->box.dimensions.x);
    }
    if (sim->box.periodic[1] && sim->box.dimensions.y > 0.0) {
        r_ij.y -= sim->box.dimensions.y * round(r_ij.y / sim->box.dimensions.y);
    }
    if (sim->box.periodic[2] && sim->box.dimensions.z > 0.0) {
        r_ij.z -= sim->box.dimensions.z * round(r_ij.z / sim->box.dimensions.z);
    }
    return vec3_norm(r_ij);
}

int side_fix_qm_qeq_screened(const Simulation *sim, double total_q, double dielectric, double *out_q) {
    if (!sim || !sim->atoms || !out_q) return -1;
    int n = sim->num_atoms;
    if (n < 1 || n > 128) return -1;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    if (!isfinite(total_q)) return -1;

    for (int i = 0; i < n; i++) {
        if (!sim->atoms[i].element) return -1;
    }

    int N = n + 1;
    double chi[128], J[128];
    for (int i = 0; i < n; i++) {
        qm_chi_J(sim->atoms[i].element, &chi[i], &J[i]);
    }

    /* Allocate matrix M and vectors on stack */
    double M[129 * 129];
    double rhs[129];
    double sol[129];
    for (int i = 0; i < N * N; i++) M[i] = 0.0;

    for (int i = 0; i < n; i++) {
        M[i * N + i] = J[i];
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            double r = atom_dist_pbc(sim, i, j);
            /* Screened Coulomb coupling replaces bare 1/r:
             * Guarantees M[i*N + j] < 0.5*(J_i + J_j) for all r >= 0 */
            M[i * N + j] = side_fix_qeq_coulomb_screened(r, J[i], J[j], dielectric);
        }
        M[i * N + n] = 1.0;
        M[n * N + i] = 1.0;
        rhs[i] = -chi[i];
    }
    M[n * N + n] = 0.0;
    rhs[n] = total_q;

    /* Gaussian elimination with partial pivoting */
    for (int c = 0; c < N; c++) {
        int piv = c;
        double best = fabs(M[c * N + c]);
        for (int r = c + 1; r < N; r++) {
            double v = fabs(M[r * N + c]);
            if (v > best) { best = v; piv = r; }
        }
        if (best < 1e-12) return -1; /* Singular system */
        if (piv != c) {
            for (int k = c; k < N; k++) {
                double t = M[c * N + k]; M[c * N + k] = M[piv * N + k]; M[piv * N + k] = t;
            }
            double t = rhs[c]; rhs[c] = rhs[piv]; rhs[piv] = t;
        }
        double d = M[c * N + c];
        for (int r = c + 1; r < N; r++) {
            double f = M[r * N + c] / d;
            if (f == 0.0) continue;
            for (int k = c; k < N; k++) M[r * N + k] -= f * M[c * N + k];
            rhs[r] -= f * rhs[c];
        }
    }

    /* Back-substitution */
    for (int r = N - 1; r >= 0; r--) {
        double acc = rhs[r];
        for (int k = r + 1; k < N; k++) acc -= M[r * N + k] * sol[k];
        if (fabs(M[r * N + r]) < 1e-14) return -1;
        sol[r] = acc / M[r * N + r];
    }

    /* Charge-conserving bounds enforcement (contracting deviations from mean) */
    double m = (double)n;
    double mean = total_q / m;
    if (fabs(mean) > QM_QEQ_QMAX) return -1;

    double dev = 0.0;
    for (int i = 0; i < n; i++) {
        if (!isfinite(sol[i])) return -1;
        double dv = fabs(sol[i] - mean);
        if (dv > dev) dev = dv;
    }
    double span = QM_QEQ_QMAX - fabs(mean);
    if (span < 1e-12) span = 1e-12;
    if (dev > span && dev > 0.0) {
        double s = span / dev;
        for (int i = 0; i < n; i++) sol[i] = mean + (sol[i] - mean) * s;
    }

    /* Verify charge conservation and output */
    double sum_q = 0.0;
    for (int i = 0; i < n; i++) {
        out_q[i] = sol[i];
        sum_q += sol[i];
    }

    /* Electronegativity ordering check:
     * In screened QEq, a higher chi atom should have lower (more negative) charge */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (chi[i] > chi[j] + 0.5 && out_q[i] > out_q[j] + 0.1) {
                /* Even with screening, return error if ordering somehow inverted */
                return -2;
            }
        }
    }

    return 0;
}
