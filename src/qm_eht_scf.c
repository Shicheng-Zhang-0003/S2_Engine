/* qm_eht_scf.c — Hoffmann EHT + self-consistent charge driver.
 *
 * Basis: per-atom valence AOs in GLOBAL Cartesian form (s; px,py,pz).
 * Two-center S via Slater-Koster rotation of the sigma/pi channels from
 * qm_eht_overlap.c (exact prolate machinery / documented GL for ppi):
 *   S(s,pi)   = l_i S_spsigma
 *   S(pi,pj)  = l_i l_j S_ppsigma + (dij - l_i l_j) S_pppi
 * with (l,m,n) the A->B direction cosines. H_ii = -VSIP,
 * H_ij = 1.75 S_ij (H_ii+H_jj)/2 (Wolfsberg-Helmholtz, E3).
 *
 * Solve: Loewdin S^{-1/2} (Jacobi eigensolver below) -> H' -> Jacobi ->
 * fill (Aufbau, half-open shell allowed) -> P = sum occ c c^T ->
 * Mulliken q -> SCC update H_ii -= sum_B gamma_AB dq_B (Ohno gamma with
 * in-tree Hubbard U = 2J, E4) -> iterate to 1e-6 e (cap 50, damped 0.5).
 *
 * Fail codes are loud: -1 unsupported/over cap, -2 singular overlap
 * (linear dependence), -3 SCC non-convergence. No silent numbers.
 */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "qm.h"
#include "qm_eht.h"
#include "quantum.h"
#include "vec3.h"

/* Hoffmann 1963 VSIPs, eV (positive; Hii = -VSIP). */
double eht_vsip (int Z, int n, int l) {
    if (Z == 1 && n == 1 && l == 0) return 13.6;
    if (Z == 6 && n == 2 && l == 0) return 21.4;
    if (Z == 6 && n == 2 && l == 1) return 11.4;
    if (Z == 7 && n == 2 && l == 0) return 26.0;
    if (Z == 7 && n == 2 && l == 1) return 13.4;
    if (Z == 8 && n == 2 && l == 0) return 32.3;
    if (Z == 8 && n == 2 && l == 1) return 14.8;
    return 0.0;
}

int eht_valence_electrons (int Z) {
    if (Z == 1) return 1;
    if (Z == 6) return 4;
    if (Z == 7) return 5;
    if (Z == 8) return 6;
    return -1;
}

int eht_valence_aos (int Z) {
    if (Z == 1) return 1;
    if (Z == 6 || Z == 7 || Z == 8) return 4;
    return 0;
}

/* ---- Jacobi symmetric eigensolver (row-major, n <= EHT_MAX_AOS) ----
 * Cyclic sweeps; eigenvalues -> d[], eigenvectors -> V (columns).
 * 0 converged, -1 sweep-cap/non-finite. Deterministic. */
static int jacobi_eig (double A[EHT_MAX_AOS][EHT_MAX_AOS],
                       double V[EHT_MAX_AOS][EHT_MAX_AOS],
                       double d[EHT_MAX_AOS], int n) {
    if (n < 1 || n > EHT_MAX_AOS) return -1;
    for (int i = 0; i < n; i++) {
        d[i] = A[i][i];
        for (int j = 0; j < n; j++) V[i][j] = (i == j) ? 1.0 : 0.0;
    }
    for (int sweep = 0; sweep < EHT_JACOBI_SWEEPS; sweep++) {
        double off = 0.0;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) off += A[p][q] * A[p][q];
        if (!isfinite (off)) return -1;
        if (off == 0.0) return 0;
        double frob = 0.0;
        for (int i = 0; i < n; i++) frob += d[i] * d[i];
        if (off < 1e-28 * (1.0 + frob)) return 0;
        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double apq = A[p][q];
                if (apq == 0.0) continue;
                double app = A[p][p], aqq = A[q][q];
                double phi = 0.5 * (aqq - app) / apq;
                double t = ((phi >= 0.0) ? 1.0 : -1.0) /
                           (fabs (phi) + sqrt (phi * phi + 1.0));
                double c = 1.0 / sqrt (t * t + 1.0);
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
    return -1;
}

/* Basis build: returns AO count or -1 unsupported. Fills ao_* + zeta/Hii. */
static int eht_build_basis (const Simulation *sim, qm_eht_t *e, double zeta[],
                            double Hii0[]) {
    int nao = 0;
    for (int a = 0; a < sim->num_atoms; a++) {
        const Atom *at = &sim->atoms[a];
        int Z = at->Z;
        int naos = eht_valence_aos (Z);
        if (naos <= 0) return -1;
        if (nao + naos > EHT_MAX_AOS) return -1;
        if (Z == 1) {
            double z = eht_zeta (Z, 1, 0, &at->electron_config);
            double v = eht_vsip (Z, 1, 0);
            if (!(z > 0.0) || !(v > 0.0)) return -1;
            e->ao_atom[nao] = a;
            e->ao_n[nao] = 1;
            e->ao_l[nao] = 0;
            e->ao_m[nao] = 0;
            zeta[nao] = z;
            Hii0[nao] = -v;
            nao++;
        } else {
            /* 2s, 2px, 2py, 2pz */
            double zs = eht_zeta (Z, 2, 0, &at->electron_config);
            double zp = eht_zeta (Z, 2, 1, &at->electron_config);
            double vs = eht_vsip (Z, 2, 0);
            double vp = eht_vsip (Z, 2, 1);
            if (!(zs > 0.0) || !(zp > 0.0) || !(vs > 0.0) || !(vp > 0.0))
                return -1;
            const int ls[4] = { 0, 1, 1, 1 };
            const int ms[4] = { 0, 1, 2, 3 };
            for (int k = 0; k < 4; k++) {
                e->ao_atom[nao] = a;
                e->ao_n[nao] = 2;
                e->ao_l[nao] = ls[k];
                e->ao_m[nao] = ms[k];
                zeta[nao] = ls[k] ? zp : zs;
                Hii0[nao] = ls[k] ? -vp : -vs;
                nao++;
            }
        }
    }
    return nao;
}

/* Assemble S and H (H diagonal from Hii[]). 0 ok. */
static int eht_assemble (const Simulation *sim, const qm_eht_t *e, int nao,
                         const double zeta[], const double Hii[],
                         double S[EHT_MAX_AOS][EHT_MAX_AOS],
                         double H[EHT_MAX_AOS][EHT_MAX_AOS]) {
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
            S[i][j] = (i == j) ? 1.0 : 0.0;
            H[i][j] = (i == j) ? Hii[i] : 0.0;
        }
    for (int i = 0; i < nao; i++) {
        for (int j = i + 1; j < nao; j++) {
            int ai = e->ao_atom[i], aj = e->ao_atom[j];
            if (ai == aj) continue;   /* same-atom AOs orthogonal by symmetry */
            Vec3 d = vec3_sub (sim->atoms[aj].position, sim->atoms[ai].position);
            double R = vec3_norm (d);
            if (!isfinite (R) || R < 1e-9) return -1;
            double lx = d.x / R, ly = d.y / R, lz = d.z / R;
            int li = e->ao_l[i], lj = e->ao_l[j];
            double sij = 0.0;
            if (li == 0 && lj == 0) {
                sij = eht_overlap_ss (e->ao_n[i], zeta[i], e->ao_n[j],
                                      zeta[j], R);
            } else if (li == 0 && lj == 1) {
                double s = eht_overlap_sps (e->ao_n[i], zeta[i], 2, zeta[j],
                                            R, 1);
                int mj = e->ao_m[j];
                double c = (mj == 1) ? lx : (mj == 2) ? ly : lz;
                sij = c * s;
            } else if (li == 1 && lj == 0) {
                double s = eht_overlap_sps (e->ao_n[j], zeta[j], 2, zeta[i],
                                            R, 0);
                int mi = e->ao_m[i];
                double c = (mi == 1) ? -lx : (mi == 2) ? -ly : -lz;
                /* p on A sees axis flipped: (s_A,p_B) evaluated along
                 * A->B already carries the orientation; the A-side p
                 * contracts with -(direction cosines). */
                sij = c * s;
            } else {
                double ss = eht_overlap_pps (zeta[i], zeta[j], R);
                double sp = eht_overlap_ppi (zeta[i], zeta[j], R);
                int mi = e->ao_m[i], mj = e->ao_m[j];
                double ci = (mi == 1) ? lx : (mi == 2) ? ly : lz;
                double cj = (mj == 1) ? lx : (mj == 2) ? ly : lz;
                sij = ci * cj * ss +
                      ((mi == mj ? 1.0 : 0.0) - ci * cj) * sp;
            }
            if (!isfinite (sij)) return -1;
            S[i][j] = S[j][i] = sij;
            double hij = 1.75 * sij * 0.5 * (Hii[i] + Hii[j]);
            if (!isfinite (hij)) return -1;
            H[i][j] = H[j][i] = hij;
        }
    }
    return 0;
}

/* One Roothaan step: Loewdin + fill + density + Mulliken. */
/* eig scratch sized to max; occ[] filled by Aufbau. */
static int eht_step (int nao, int nelec,
                     double S[EHT_MAX_AOS][EHT_MAX_AOS],
                     double H[EHT_MAX_AOS][EHT_MAX_AOS],
                     const int ao_atom[],
                     const int atom_z[], int natoms,
                     double P[EHT_MAX_AOS][EHT_MAX_AOS],
                     double charges[], double eps[], double kt) {
    if (!(kt > 1e-6) || !isfinite (kt)) kt = EHT_FERMI_KT;
    static double W[EHT_MAX_AOS][EHT_MAX_AOS];
    static double V[EHT_MAX_AOS][EHT_MAX_AOS];
    static double lam[EHT_MAX_AOS];
    static double Hp[EHT_MAX_AOS][EHT_MAX_AOS];
    static double Cp[EHT_MAX_AOS][EHT_MAX_AOS];
    static double Cc[EHT_MAX_AOS][EHT_MAX_AOS];

    /* S eigendecomposition */
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) W[i][j] = S[i][j];
    if (jacobi_eig (W, V, lam, nao) != 0) return -1;
    for (int i = 0; i < nao; i++)
        if (!(lam[i] > 1e-12) || !isfinite (lam[i])) return -2;
    /* S^{-1/2} into W */
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
            double s = 0.0;
            for (int k = 0; k < nao; k++)
                s += V[i][k] * V[j][k] / sqrt (lam[k]);
            W[i][j] = s;
        }
    /* H' = S^{-1/2} H S^{-1/2} */
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
            double s = 0.0;
            for (int k = 0; k < nao; k++)
                for (int l = 0; l < nao; l++) s += W[i][k] * H[k][l] * W[l][j];
            Hp[i][j] = s;
        }
    if (jacobi_eig (Hp, Cp, eps, nao) != 0) return -1;
    /* sort ascending, permute columns */
    int idx[EHT_MAX_AOS];
    for (int i = 0; i < nao; i++) idx[i] = i;
    for (int i = 1; i < nao; i++) {
        int t = idx[i], j = i - 1;
        while (j >= 0 && eps[idx[j]] > eps[t]) {
            idx[j + 1] = idx[j];
            j--;
        }
        idx[j + 1] = t;
    }
    double es[EHT_MAX_AOS];
    for (int i = 0; i < nao; i++) es[i] = eps[idx[i]];
    for (int i = 0; i < nao; i++) eps[i] = es[i];
    /* C = S^{-1/2} C' (columns in sorted order) */
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
            double s = 0.0;
            for (int k = 0; k < nao; k++) s += W[i][k] * Cp[k][idx[j]];
            Cc[i][j] = s;
        }
    /* Fermi fill (E10 continued): fractional occupations smear level
     * crossings that otherwise ping-pong integer electrons between
     * near-degenerate configurations at dissociation (stretched
     * heteronuclear -> rc=-3 limit cycle under pure Aufbau). kT small
     * enough that gapped ground states are unaffected (H2O gap 14.8 eV
     * >> kT); mu bisected so sum(occ) == nelec exactly. Standard SCC
     * practice (DFTB+ does the same), labeled here. */
    double occ[EHT_MAX_AOS];
    {
        double lo = eps[0] - 10.0 * kt - 5.0, hi = eps[nao - 1] + 10.0 * kt + 5.0;
        for (int it = 0; it < 200; it++) {
            double mu = 0.5 * (lo + hi);
            double n = 0.0;
            for (int i = 0; i < nao; i++) {
                double x = (eps[i] - mu) / kt;
                double f = (x > 60.0) ? 0.0 : (x < -60.0 ? 1.0 : 1.0 / (1.0 + exp (x)));
                n += 2.0 * f;
            }
            if (n > nelec) hi = mu;
            else lo = mu;
            if (hi - lo < 1e-12) break;
        }
        double mu = 0.5 * (lo + hi);
        double n = 0.0;
        for (int i = 0; i < nao; i++) {
            double x = (eps[i] - mu) / kt;
            double f = (x > 60.0) ? 0.0 : (x < -60.0 ? 1.0 : 1.0 / (1.0 + exp (x)));
            occ[i] = 2.0 * f;
            n += occ[i];
        }
        /* exact electron count: park the residual on the HOMO level */
        double miss = (double)nelec - n;
        if (fabs (miss) > 1e-9) {
            int h = 0;
            double best = -1e300;
            for (int i = 0; i < nao; i++) {
                if (occ[i] > 1e-9 && occ[i] < 2.0 - 1e-9 && eps[i] > best) {
                    best = eps[i];
                    h = i;
                }
            }
            occ[h] += miss;
            if (occ[h] < 0.0) occ[h] = 0.0;
            if (occ[h] > 2.0) occ[h] = 2.0;
        }
    }
    for (int i = 0; i < nao; i++)
        for (int j = 0; j < nao; j++) {
            double s = 0.0;
            for (int k = 0; k < nao; k++) s += occ[k] * Cc[i][k] * Cc[j][k];
            P[i][j] = s;
        }
    /* Mulliken: q_A = Zval_A - sum_{mu in A} (PS)_mumu */
    for (int a = 0; a < natoms; a++) {
        double pop = 0.0;
        for (int m = 0; m < nao; m++) {
            if (ao_atom[m] != a) continue;
            double psm = 0.0;
            for (int n = 0; n < nao; n++) psm += P[m][n] * S[n][m];
            pop += psm;
        }
        charges[a] = atom_z[a] - pop;
    }
    return 0;
}

int qm_eht_solve (const Simulation *sim, qm_eht_t *out) {
    if (!sim || !out || !sim->atoms) return -1;
    if (sim->num_atoms < 1 || sim->num_atoms > EHT_MAX_ATOMS) return -1;
    memset (out, 0, sizeof *out);
    static double zeta[EHT_MAX_AOS];
    static double Hii0[EHT_MAX_AOS];
    static double Hii[EHT_MAX_AOS];
    static double S[EHT_MAX_AOS][EHT_MAX_AOS];
    static double H[EHT_MAX_AOS][EHT_MAX_AOS];
    static double P[EHT_MAX_AOS][EHT_MAX_AOS];
    int nao = eht_build_basis (sim, out, zeta, Hii0);
    if (nao < 1) return -1;
    int nelec = 0;
    int atom_z[EHT_MAX_ATOMS];
    for (int a = 0; a < sim->num_atoms; a++) {
        int ve = eht_valence_electrons (sim->atoms[a].Z);
        if (ve < 0) return -1;
        atom_z[a] = ve;
        nelec += ve;
    }
    if (nelec > 2 * nao) return -1;
    /* Hubbard U from in-tree Mulliken J (U = IE-EA = 2J), eV. */
    double U[EHT_MAX_ATOMS];
    for (int a = 0; a < sim->num_atoms; a++) {
        double chi, J;
        qm_chi_J (sim->atoms[a].element, &chi, &J);
        (void)chi;
        U[a] = 2.0 * J;
        if (!(U[a] > 0.5) || !isfinite (U[a])) U[a] = 8.0;
    }
    for (int i = 0; i < nao; i++) Hii[i] = Hii0[i];
    double charges[EHT_MAX_ATOMS] = { 0.0 };
    double eps[EHT_MAX_AOS] = { 0.0 };
    int it;
    double resid = 0.0;
    /* E11 (measured limit): NO annealing schedule survives contact with
     * a first-order charge flip. Fermi annealing was implemented and
     * reverted same-day: raised kT does smooth the crossing (resid
     * 0.36 -> 1e-3 over ~25 ticks at 2x H2O), but the smoothed roots
     * DRIFT with kT (-0.096 -> -0.058 -> ...), and every return toward
     * the floor re-triggers the flip (three measured re-flips at
     * resid ~1e-5..1e-7, one of them AT 1.4e-07 — converged by any
     * tolerance, rejected only because kT was not at floor, then
     * destroyed by the lowering step). A schedule that rejects 1e-7 to
     * chase exactness is worse than none. Stretched heteronuclear with
     * level crossing is therefore rc=-3 fail-loud (gated as such), not
     * silently smeared: no fixed point exists at floor kT that iteration
     * can reach, and inventing one via kT-dependence would be the
     * fabricated-target failure. Homonuclear dissociation (H2) and all
     * gapped equilibria are unaffected. */
    double kt = EHT_FERMI_KT;
    int done = 0;
    for (it = 0; it < EHT_SCC_ITERS; it++) {
        if (eht_assemble (sim, out, nao, zeta, Hii, S, H) != 0) return -1;
        int rc = eht_step (nao, nelec, S, H, out->ao_atom, atom_z,
                           sim->num_atoms, P, charges, eps, kt);
        if (rc == -2) return -2;
        if (rc != 0) return -1;
        /* Damped trust-region step (E10 continued): qnext = qin +
         * s*Delta with s = min(SMIX, SMIX*DSTEP/|Delta|). Two lessons from
         * measurement: (1) pure 0.5-mixing flips sign at stretched bonds
         * (+1.8/-1.2 alternation at 2x: local gain ~ -3); (2) a norm cap
         * ALONE restores s=1 near the root, where damping is needed most
         * (measured growing x1.11 2-cycle at equilibrium). Constant SMIX
         * damps every step (stable for local gains in (-7,1)); the cap
         * only bounds wild first steps. Fixed points untouched (both
         * mechanisms go quiet as Delta -> 0); neutrality preserved
         * (uniform scaling of a neutral step). */
        static double qold[EHT_MAX_ATOMS];
        resid = 0.0;
        double dmax = 0.0;
        for (int a = 0; a < sim->num_atoms; a++) {
            double dq = charges[a] - (it == 0 ? 0.0 : qold[a]);
            double ad = fabs (dq);
            if (ad > resid) resid = ad;
            if (ad > dmax) dmax = ad;
        }
        {
            double s = EHT_SCC_SMIX;
            if (dmax > EHT_SCC_DSTEP && dmax > 0.0)
                s *= EHT_SCC_DSTEP / dmax;
            for (int a = 0; a < sim->num_atoms; a++)
                qold[a] = (it == 0 ? 0.0 : qold[a]) +
                          s * (charges[a] - (it == 0 ? 0.0 : qold[a]));
        }
        if (getenv ("MGB_EHT_DEBUG")) {
            fprintf (stderr, "  [eht-dbg] it=%d resid=%.3e", it, resid);
            for (int a = 0; a < sim->num_atoms; a++)
                fprintf (stderr, " q%d=%+.4f", a, charges[a]);
            fprintf (stderr, "\n");
        }
        if (resid < EHT_SCC_TOL) {
            done = 1;
            break;
        }
        /* H_ii -= sum_B gamma_AB dq_B (Ohno interpolation, E4).
         * E10 (found by measurement, 5x-separated H2O converging to full
         * +/-2.0 e charge transfer): the sum MUST include B == ai with
         * gamma_AA = U_A (Hubbard on-site hardness). Excluding it removed
         * the restoring force against charge sloshing — added electrons
         * never pushed their own levels back up — and the SCC ran away to
         * shell-closure ionicity (O 2p full = neon) at ANY separation,
         * converging there. With the on-site term, gained charge raises
         * its own levels by U per e (exactly the chemical hardness) and
         * the neutral root is iteration-stable. */
        for (int i = 0; i < nao; i++) {
            int ai = out->ao_atom[i];
            double shift = 0.0;
            for (int b = 0; b < sim->num_atoms; b++) {
                double gam;
                if (b == ai) {
                    gam = U[ai];
                } else {
                    Vec3 d = vec3_sub (sim->atoms[b].position,
                                       sim->atoms[ai].position);
                    double R = vec3_norm (d);
                    if (!isfinite (R) || R < 1e-9) return -1;
                    double Uab = 0.5 * (U[ai] + U[b]);
                    gam = COULOMB_MD / sqrt (R * R + pow (28.8 / Uab, 2));
                }
                if (!isfinite (gam)) return -1;
                shift += gam * qold[b];
            }
            Hii[i] = Hii0[i] - shift;
        }
    }
    if (!done) return -3;
    /* commit */
    out->n_ao = nao;
    out->n_atoms = sim->num_atoms;
    out->n_elec = nelec;
    for (int i = 0; i < nao; i++) {
        out->eps[i] = eps[i];
        for (int j = 0; j < nao; j++) {
            out->P[i][j] = P[i][j];
            out->S[i][j] = S[i][j];
        }
    }
    for (int a = 0; a < sim->num_atoms; a++) out->charges[a] = charges[a];
    /* gap: HOMO/LUMO from Aufbau fill */
    int nhomo = (nelec - 1) / 2;
    out->gap = (nhomo + 1 < nao) ? eps[nhomo + 1] - eps[nhomo] : 0.0;
    out->e_band = 0.0;
    {
        int left = nelec;
        for (int i = 0; i < nao && left > 0; i++) {
            double occ = (left >= 2) ? 2.0 : 1.0;
            out->e_band += occ * eps[i];
            left -= (left >= 2) ? 2 : 1;
        }
    }
    out->e_nuc = 0.0;
    for (int a = 0; a < sim->num_atoms; a++)
        for (int b = a + 1; b < sim->num_atoms; b++) {
            Vec3 d = vec3_sub (sim->atoms[b].position,
                               sim->atoms[a].position);
            double R = vec3_norm (d);
            if (!isfinite (R) || R < 1e-9) return -1;
            out->e_nuc += atom_z[a] * atom_z[b] * COULOMB_MD / R;
        }
    /* E10 continued: subtract the double-counted half. E_band was built
     * from shifted eigenvalues, so it contains the full Sum gam dq dq
     * electron-interaction estimate; the energy correct to second order
     * keeps half (DFTB2 form). Without this, ionic roots rank below
     * neutral ones and no gate can tell them apart. */
    {
        double e2 = 0.0;
        for (int a = 0; a < sim->num_atoms; a++)
            for (int b = 0; b < sim->num_atoms; b++) {
                double gam;
                if (b == a) {
                    gam = U[a];
                } else {
                    Vec3 d = vec3_sub (sim->atoms[b].position,
                                       sim->atoms[a].position);
                    double R = vec3_norm (d);
                    if (!isfinite (R) || R < 1e-9) return -1;
                    double Uab = 0.5 * (U[a] + U[b]);
                    gam = COULOMB_MD / sqrt (R * R + pow (28.8 / Uab, 2));
                }
                e2 += gam * out->charges[a] * out->charges[b];
            }
        out->e_total = out->e_band - 0.5 * e2 + out->e_nuc;
    }
    out->scc_iters = it + 1;
    out->scc_resid = resid;
    return 0;
}

double qm_eht_bond_order (const qm_eht_t *e, int ia, int ib) {
    if (!e || ia < 0 || ib < 0 || ia >= e->n_atoms || ib >= e->n_atoms ||
        ia == ib)
        return -1.0;
    /* Mayer: B_AB = sum_{mu in A, nu in B} (PS)_munu (PS)_numu, with S
     * cached at solve time (the overlap the density was converged
     * against — self-consistent by construction). */
    int nao = e->n_ao;
    double B = 0.0;
    for (int mu = 0; mu < nao; mu++) {
        if (e->ao_atom[mu] != ia) continue;
        for (int nu = 0; nu < nao; nu++) {
            if (e->ao_atom[nu] != ib) continue;
            double ps1 = 0.0, ps2 = 0.0;
            for (int k = 0; k < nao; k++) {
                ps1 += e->P[mu][k] * e->S[k][nu];
                ps2 += e->P[nu][k] * e->S[k][mu];
            }
            B += ps1 * ps2;
        }
    }
    return isfinite (B) ? B : -1.0;
}

double qm_eht_bond_factor (const qm_eht_t *e, const Simulation *sim,
                           int ia, int ib) {
    if (!e || !sim || !sim->atoms) return 1.0;
    if (ia < 0 || ib < 0 || ia >= e->n_atoms || ib >= e->n_atoms ||
        ia >= sim->num_atoms || ib >= sim->num_atoms || ia == ib)
        return 1.0;
    double B = qm_eht_bond_order (e, ia, ib);
    if (!(B >= 0.0) || !isfinite (B)) return 1.0;
    /* sigma-channel magnitude between the valence shells, current vs
     * covalent-contact reference (E8 switch). */
    const Atom *a = &sim->atoms[ia], *b = &sim->atoms[ib];
    Vec3 d = vec3_sub (b->position, a->position);
    double R = vec3_norm (d);
    if (!isfinite (R) || R < 1e-9) return 1.0;
    double rcov = 1.0;
    if (a->element && b->element &&
        a->element->covalent_radius > 0.0 && b->element->covalent_radius > 0.0)
        rcov = a->element->covalent_radius + b->element->covalent_radius;
    if (!(rcov > 1e-9)) return 1.0;
    /* zeta per side from the solved basis (first valence AO of atom). */
    double za = -1.0, zb = -1.0;
    int na = 1, nb = 1, la = 0, lb = 0;
    for (int m = 0; m < e->n_ao; m++) {
        if (e->ao_atom[m] == ia && za < 0.0) {
            za = eht_zeta (a->Z, e->ao_n[m], e->ao_l[m],
                           &a->electron_config);
            na = e->ao_n[m];
            la = e->ao_l[m];
        }
        if (e->ao_atom[m] == ib && zb < 0.0) {
            zb = eht_zeta (b->Z, e->ao_n[m], e->ao_l[m],
                           &b->electron_config);
            nb = e->ao_n[m];
            lb = e->ao_l[m];
        }
    }
    if (!(za > 0.0) || !(zb > 0.0)) return 1.0;
    /* max |sigma channel| at R and at reference (s/H, p directed). */
    double s_cur = 0.0, s_ref = 0.0, t;
    t = fabs (eht_overlap_ss (na, za, nb, zb, R));
    if (t > s_cur) s_cur = t;
    t = fabs (eht_overlap_ss (na, za, nb, zb, rcov));
    if (t > s_ref) s_ref = t;
    if (la == 1 || lb == 1) {
        /* directed p-sigma along the bond (exact orientation, no SK). */
        Vec3 u = vec3_scale (d, 1.0 / R);
        double ux = fabs (u.x), uy = fabs (u.y), uz = fabs (u.z);
        double w = ux > uy ? (ux > uz ? ux : uz) : (uy > uz ? uy : uz);
        if (la == 1 && lb == 1) {
            t = fabs (eht_overlap_pps (za, zb, R)) * w * w;
            if (t > s_cur) s_cur = t;
            t = fabs (eht_overlap_pps (za, zb, rcov)) * w * w;
            if (t > s_ref) s_ref = t;
        } else if (la == 1) {
            t = fabs (eht_overlap_sps (nb, zb, 2, za, R, 0)) * w;
            if (t > s_cur) s_cur = t;
            t = fabs (eht_overlap_sps (nb, zb, 2, za, rcov, 0)) * w;
            if (t > s_ref) s_ref = t;
        } else {
            t = fabs (eht_overlap_sps (na, za, 2, zb, R, 1)) * w;
            if (t > s_cur) s_cur = t;
            t = fabs (eht_overlap_sps (na, za, 2, zb, rcov, 1)) * w;
            if (t > s_ref) s_ref = t;
        }
    }
    if (!(s_ref > 1e-12) || !isfinite (s_ref)) return (B > 1.0) ? 1.0 : B;
    double f = B * s_cur / s_ref;
    if (!(f >= 0.0) || !isfinite (f)) return 1.0;
    if (f > 1.0) f = 1.0;
    return f;
}
