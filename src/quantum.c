#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/quantum.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"

/*
 * quantum.c
 * Atomic quantum mechanics: Slater effective nuclear charge,
 * orbital energies, and hydrogen-like wave functions.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Slater's effective principal quantum number n*
 *
 * This corrects for the fact that inner electrons penetrate the nucleus
 * less efficiently at higher n, giving an effective quantum number that
 * produces better orbital energies than the bare n.
 *
 * Table from Slater (1930), Phys. Rev. 36, 57 (n=1..6). n=7 reuses 4.2
 * as an unsourced extrapolation for completeness; no H..Kr demo reaches
 * n=7 in the ground state, so this value never affects the record.
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_nstar(int n) {
    if (n < 1) return 1.0;
    switch (n) {
        case 1: return 1.0;
        case 2: return 2.0;
        case 3: return 3.0;
        case 4: return 3.7;
        case 5: return 4.0;
        case 6: return 4.2;
        case 7: return 4.2;
        default: return (double)n;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Slater grouping helper
 *
 * Maps each (n, l) orbital to its Slater group index:
 *   Group 0 : 1s
 *   Group 1 : 2s, 2p
 *   Group 2 : 3s, 3p
 *   Group 3 : 3d
 *   Group 4 : 4s, 4p
 *   Group 5 : 4d
 *   Group 6 : 4f
 *   Group 7 : 5s, 5p
 *   Group 8 : 5d
 *   Group 9 : 5f
 *   Group 10: 6s, 6p
 *   ...
 *
 * For s/p: group = 2*(n-1) - (n>2 ? 1 : 0) ... simplest lookup:
 * ══════════════════════════════════════════════════════════════════════════ */
static int slater_group(int n, int l) {
    /* d and f orbitals get their own groups below the corresponding sp */
    if (l >= 2) {
        /* 3d=3, 4d=5, 4f=6, 5d=8, 5f=9, 6d=11, ... */
        if (n == 3 && l == 2) return 3;
        if (n == 4 && l == 2) return 5;
        if (n == 4 && l == 3) return 6;
        if (n == 5 && l == 2) return 8;
        if (n == 5 && l == 3) return 9;
        if (n == 6 && l == 2) return 11;
        if (n == 6 && l == 3) return 12;
        return 2*n;           /* fallback */
    }
    /* s and p: group by n */
    switch (n) {
        case 1: return 0;
        case 2: return 1;
        case 3: return 2;
        case 4: return 4;
        case 5: return 7;
        case 6: return 10;
        case 7: return 13;
        default: return 2*n;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Slater screening constant (raw: no clamp, no stderr — for probes)
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_zeff_raw(int Z, int n, int l, const ElectronConfig *cfg) {
    if (!cfg || Z < 1 || n < 1 || l < 0 || l >= n) return 1.0;
    int target_group = slater_group(n, l);
    double S = 0.0;

    /* Iterate over all occupied orbitals in the Madelung filling order */
    static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    static const int MLEN = 19;

    for (int i = 0; i < MLEN; i++) {
        int on = MN[i];
        int ol = ML[i];
        int occ = cfg->config[on-1][ol]; /* electrons in this orbital type */
        if (occ == 0) continue;

        int other_group = slater_group(on, ol);

        /* Count electrons in this orbital, excluding the target electron */
        int electrons = occ;
        if (on == n && ol == l) electrons -= 1; /* don't count self */
        if (electrons <= 0) continue;

        /*
         * Determine contribution to S based on Slater's rules.
         *
         * IMPORTANT: for s/p targets, the 0.85 "one shell inner" tier is
         * defined by PRINCIPAL QUANTUM NUMBER (n-1) directly - it covers
         * (n-1)s, (n-1)p, AND (n-1)d together, NOT just "the immediately
         * preceding entry in the sequential group list" (which would
         * wrongly put (n-1)d and (n-1)s,(n-1)p into DIFFERENT tiers
         * whenever (n-1)d is populated). Verified against Slater's own
         * original worked example (Phys. Rev. 36, 57, 1930) for iron's
         * 4s electron: sigma = 0.35x1 + 0.85x14 + 1.00x10, where the 14
         * is 3s^2+3p^6+3d^6 = 14 electrons combined into ONE 0.85 tier.
         * An earlier version of this code used group-list adjacency
         * here and got iron's 4s Zeff wrong by 1.2 (3.75 correct vs
         * 2.55 computed) - confirmed by reproducing Slater's own
         * example numerically before fixing this.
         */
        if (l <= 1) {
            /* Target is s or p electron */
            if (on == n) {
                /* Same shell (ns, np together): 0.35 (0.30 for 1s) */
                double contrib = (n == 1) ? 0.30 : 0.35;
                S += contrib * electrons;
            } else if (on == n - 1) {
                /* One shell inner - (n-1)s, (n-1)p, AND (n-1)d together */
                S += 0.85 * electrons;
            } else if (on <= n - 2) {
                /* Two or more shells inner */
                S += 1.00 * electrons;
            }
            /* Higher shells do not screen (impossible in ground state) */
        } else {
            /* Target is d or f electron - group-list adjacency IS correct
             * here (verified against Slater's Fe 3d example: sigma =
             * 0.35x5 + 1.00x18 = 19.75, matching exactly) */
            if (other_group == target_group) {
                S += 0.35 * electrons;
            } else if (other_group < target_group) {
                S += 1.00 * electrons;
            }
        }
    }

    double Zeff = (double)Z - S;
    return Zeff;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Slater screening constant (public: clamped, warns once on stderr)
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_zeff(int Z, int n, int l, const ElectronConfig *cfg) {
    if (!cfg || Z < 1 || n < 1 || l < 0 || l >= n) return 1.0;
    double Zeff = quantum_zeff_raw(Z, n, l, cfg);
    if (Zeff < 1.0) {
        /* Clamp honesty: S exceeding Z-1 is unphysical for the supplied
         * configuration (e.g. over-screened anion census); floor at 1.0
         * and warn once rather than silently returning a bare proton. */
        static volatile int warned = 0;
        if (!warned) {
            warned = 1;
            fprintf(stderr,
                    "quantum: WARNING: Slater Zeff=%.2f for Z=%d below 1.0; "
                    "floored at 1.0\n", Zeff, Z);
        }
        return 1.0;
    }
    return Zeff;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Orbital energy
 * E_nl = -13.6058 eV × (Z_eff / n*)²
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_orbital_energy(int Z, int n, int l, const ElectronConfig *cfg) {
    double Zeff = quantum_zeff(Z, n, l, cfg);
    double nstar = quantum_nstar(n);
    return -13.605693122994 * (Zeff / nstar) * (Zeff / nstar);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Fill atom->orbitals[] from its electron configuration
 * ══════════════════════════════════════════════════════════════════════════ */
void quantum_fill_orbitals(Atom *atom) {
    if (!atom) return;
    static const int MN[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
    static const int ML[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
    static const int MLEN = 19;

    int orb_idx = 0;

    for (int i = 0; i < MLEN; i++) {
        int n = MN[i];
        int l = ML[i];
        int total_e = atom->electron_config.config[n-1][l];
        if (total_e == 0) continue;

        double energy = quantum_orbital_energy(atom->Z, n, l,
                                               &atom->electron_config);

        /*
         * Distribute electrons across ml values using Hund's rule:
         * first fill each ml with one electron (spin up), then pair.
         * ml runs from -l to +l, so 2l+1 distinct values.
         */
        int num_ml = 2*l + 1;
        int filled[7] = {0}; /* occupation per ml slot, max 7 for f */

        /* Pass 1: one electron each slot (spin up) */
        int left = total_e;
        for (int m = 0; m < num_ml && left > 0; m++, left--)
            filled[m] = 1;
        /* Pass 2: second electron (spin down) */
        for (int m = 0; m < num_ml && left > 0; m++, left--) {
            if (filled[m] == 1) filled[m] = 2;
        }

        for (int m = 0; m < num_ml; m++) {
            if (filled[m] == 0) continue;
            if (orb_idx >= MAX_ORBITALS) break;

            Orbital *orb = &atom->orbitals[orb_idx++];
            orb->qn.n    = n;
            orb->qn.l    = l;
            orb->qn.ml   = m - l;      /* ml = -l … +l */
            /* Representative spin: one entry per ml slot. A singly-occupied
             * slot holds one electron (+0.5). A doubly-occupied slot holds
             * +0.5/-0.5; no single value equals both, so store +0.5 as the
             * representative electron spin (a valid single-electron m_s)
             * and rely on occupation==2 for the pair count. Storing 0.0
             * would violate the m_s = ±1/2 definition (types.h); net zero
             * is a property of the pair, not of either electron. No
             * energy/force path reads ms. */
            orb->qn.ms   = 0.5;
            orb->orbital_energy = energy;
            orb->occupation = filled[m];
        }
    }
    atom->num_orbitals = orb_idx;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Associated Laguerre polynomial L_p^q(x)
 * Three-term recurrence: L_0^q = 1, L_1^q = 1+q-x,
 *   L_p^q = [(2p-1+q-x)L_{p-1}^q - (p-1+q)L_{p-2}^q] / p
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_laguerre(int p, int q, double x) {
    if (p == 0) return 1.0;
    if (p == 1) return 1.0 + (double)q - x;

    double L_prev2 = 1.0;
    double L_prev1 = 1.0 + (double)q - x;
    double L_curr  = 0.0;

    for (int k = 2; k <= p; k++) {
        L_curr = ((2.0*k - 1.0 + q - x) * L_prev1
                  - (k - 1.0 + q)         * L_prev2) / (double)k;
        L_prev2 = L_prev1;
        L_prev1 = L_curr;
    }
    return L_curr;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Factorial (integer, up to 20)
 * ══════════════════════════════════════════════════════════════════════════ */
static double factorial(int k) {
    static const double F[21] = {
        1,1,2,6,24,120,720,5040,40320,362880,3628800,
        39916800,479001600,6227020800.0,87178291200.0,
        1307674368000.0,20922789888000.0,355687428096000.0,
        6402373705728000.0,121645100408832000.0,
        2432902008176640000.0
    };
    if (k < 0) return 0.0; /* invalid: caller must guard n>l; 0 flags error */
    if (k > 20) {
        double r = F[20];
        for (int i = 21; i <= k; i++) r *= i;
        return r;
    }
    return F[k];
}

/* ══════════════════════════════════════════════════════════════════════════
 * Hydrogen-like radial wave function R_nl(r)
 *
 *  ρ  = 2 Z_eff r / (n a₀)           (dimensionless)
 *  N  = -sqrt[(2Z_eff/na₀)³ × (n-l-1)! / (2n (n+l)!)]
 *  R_nl(r) = N × e^(-ρ/2) × ρ^l × L_{n-l-1}^{2l+1}(ρ)
 *
 * Return units: Å^(-3/2)
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_radial_wavefunction(int n, int l, double Z_eff, double r_ang) {
    if (r_ang < 0.0 || !isfinite(r_ang)) return 0.0;
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;

    double a0_ang = BOHR_TO_ANGSTROM;          /* 0.529177 Å */
    double scale  = 2.0 * Z_eff / ((double)n * a0_ang);
    double rho    = scale * r_ang;

    /* Normalisation: N^2 = scale^3 * (n-l-1)! / (2n (n+l)!), where
     * scale = 2*Z_eff/(n*a0) has units 1/length so N carries
     * length^-3/2 as required for R_nl (so that ∫r^2|R|^2dr = 1).
     * A scale^1 form under-normalizes by exactly scale^2
     * (H 1s integral 0.07 instead of 1.0); verified numerically. */
    double num    = factorial(n - l - 1);
    double den    = 2.0 * n * factorial(n + l); /* (n+l)! to first power (s02 fix) */
    double N2     = scale * scale * scale * num / den;
    double N      = -sqrt(N2);                 /* sign convention Griffiths */

    double lag    = quantum_laguerre(n - l - 1, 2*l + 1, rho);
    double radial = N * exp(-rho / 2.0) * pow(rho, (double)l) * lag;

    return radial;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Radial probability density P(r) = r² |R_nl(r)|²
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_radial_probability(int n, int l, double Z_eff, double r_ang) {
    double R = quantum_radial_wavefunction(n, l, Z_eff, r_ang);
    return r_ang * r_ang * R * R;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Most probable radius — global maximizer of P(r) on (0, r_max]
 *
 * P_nl(r) for n >= 2 has radial NODES (n-l-1 of them) and therefore
 * MULTIPLE local maxima (e.g. 2s has two peaks) - a single
 * golden-section search over the whole interval converges to whichever
 * peak its bracket favours, not necessarily the global one. Derivation
 * of the fix: partition [a, b] into NSUB subintervals (64 is ample -
 * nodes are ~Bohr-spaced, far wider than any subinterval at the n, Zeff
 * values this codebase evaluates), run the golden-section maximizer
 * inside each subinterval, and keep the best. Each local search is the
 * same unimodal-safe contraction as before; the partition turns a
 * local method into a global one by construction. Cost is ~64 x 60
 * cheap P(r) evaluations per call - negligible for a diagnostic.
 * ══════════════════════════════════════════════════════════════════════════ */
static double golden_max(double a, double b, int n, int l, double Z_eff) {
    static const double PHI = 0.6180339887; /* 1/phi */
    double c = b - PHI * (b - a);
    double d = a + PHI * (b - a);

    for (int iter = 0; iter < 60; iter++) {
        if (fabs(b - a) < 1.0e-10) break;
        if (quantum_radial_probability(n, l, Z_eff, c) <
            quantum_radial_probability(n, l, Z_eff, d)) {
            a = c;
        } else {
            b = d;
        }
        c = b - PHI * (b - a);
        d = a + PHI * (b - a);
    }
    return (a + b) / 2.0;
}

double quantum_most_probable_radius(int n, int l, double Z_eff) {
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;
    double lo = 0.0001, hi = 30.0 * n * n / Z_eff;
    if (!(hi > lo) || !isfinite(hi)) return lo;
    /* NSUB=64 assumes at most one peak per ~0.5-16 A subinterval at the
     * n<=4, Zeff>=1 scales evaluated here. For n>=5 multiple nodes can
     * fall in one subinterval; result remains a diagnostic, not a
     * certified global optimum (documented limitation). */
    static const int NSUB = 64;
    double best_r = lo, best_p = -1.0;
    for (int s = 0; s < NSUB; s++) {
        double a = lo + (hi - lo) * s / NSUB;
        double b = lo + (hi - lo) * (s + 1) / NSUB;
        double r = golden_max(a, b, n, l, Z_eff);
        double p = quantum_radial_probability(n, l, Z_eff, r);
        if (p > best_p) { best_p = p; best_r = r; }
    }
    return best_r;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Clementi-Raimondi SCF effective charges (H..Kr)
 *
 * Source: Clementi & Raimondi, JCP 38, 2686 (1963); 1967 (Z>36),
 * via Wikipedia compact table, every cell used here cross-checked:
 * dual-sourced (WE = WebElements per-element page agrees): H-1s, C,
 * O, K-1s/2s/2p/3s/3p/4s, Ar (incl. 3p=6.76 filling Wiki blank),
 * Br (incl. 3d=19.56), Kr full row, Ge-1s. Single-source (W =
 * Wikipedia only, clean monotonic cells): everything else below.
 * DROPPED as unverifiable: K-3d (WE: "no data" vs Wiki 7.120),
 * Se-3d (non-monotonic 20.626 vs Br 19.56), B-2p/Ar-4s/etc. blanks
 * (unoccupied). Dropped cells return -1 → Slater fallback.
 * Ions: neutral exponents only; ions use Slater-on-ion-config.
 * ══════════════════════════════════════════════════════════════════════════ */
typedef struct { int Z, n, l; double z; } CRRow;
static const CRRow CR_TABLE[] = {
    {1,1,0,1.000}, {2,1,0,1.688},
    {3,1,0,2.691},{3,2,0,1.279},
    {4,1,0,3.685},{4,2,0,1.912},
    {5,1,0,4.680},{5,2,0,2.576},
    {6,1,0,5.6727},{6,2,0,3.22},{6,2,1,3.14},
    {7,1,0,6.665},{7,2,0,3.847},{7,2,1,3.834},
    {8,1,0,7.6579},{8,2,0,4.49},{8,2,1,4.45},
    {9,1,0,8.650},{9,2,0,5.128},{9,2,1,5.100},
    {10,1,0,9.642},{10,2,0,5.758},{10,2,1,5.758},
    {11,1,0,10.626},{11,2,0,6.571},{11,2,1,6.802},{11,3,0,2.507},
    {12,1,0,11.609},{12,2,0,7.392},{12,2,1,7.826},{12,3,0,3.308},
    {13,1,0,12.591},{13,2,0,8.214},{13,2,1,8.963},{13,3,0,4.117},{13,3,1,4.066},
    {14,1,0,13.575},{14,2,0,9.020},{14,2,1,9.945},{14,3,0,4.903},{14,3,1,4.285},
    {15,1,0,14.558},{15,2,0,9.825},{15,2,1,10.961},{15,3,0,5.642},{15,3,1,5.482},
    {16,1,0,15.541},{16,2,0,10.629},{16,2,1,11.977},{16,3,0,6.367},{16,3,1,6.116},
    {17,1,0,16.524},{17,2,0,11.430},{17,2,1,12.993},{17,3,0,7.068},{17,3,1,6.764},
    {18,1,0,17.5075},{18,2,0,12.23},{18,2,1,14.01},{18,3,0,7.76},{18,3,1,6.76},
    {19,1,0,18.4895},{19,2,0,13.01},{19,2,1,15.03},{19,3,0,8.68},{19,3,1,7.73},{19,4,0,3.50},
    {20,1,0,19.473},{20,2,0,13.776},{20,2,1,16.041},{20,3,0,9.602},{20,3,1,8.658},{20,4,0,4.398},
    {21,1,0,20.457},{21,2,0,14.574},{21,2,1,17.055},{21,3,0,10.340},{21,3,1,9.406},{21,4,0,4.632},{21,3,2,8.983},
    {22,1,0,21.441},{22,2,0,15.377},{22,2,1,18.065},{22,3,0,11.033},{22,3,1,10.104},{22,4,0,4.817},{22,3,2,9.757},
    {23,1,0,22.426},{23,2,0,16.181},{23,2,1,19.073},{23,3,0,11.709},{23,3,1,10.785},{23,4,0,4.981},{23,3,2,10.528},
    {24,1,0,23.414},{24,2,0,16.984},{24,2,1,20.075},{24,3,0,12.368},{24,3,1,11.466},{24,4,0,5.133},{24,3,2,11.180},
    {25,1,0,24.396},{25,2,0,17.794},{25,2,1,21.084},{25,3,0,13.018},{25,3,1,12.109},{25,4,0,5.283},{25,3,2,11.855},
    {26,1,0,25.381},{26,2,0,18.599},{26,2,1,22.089},{26,3,0,13.676},{26,3,1,12.778},{26,4,0,5.434},{26,3,2,12.530},
    {27,1,0,26.367},{27,2,0,19.405},{27,2,1,23.092},{27,3,0,14.322},{27,3,1,13.435},{27,4,0,5.576},{27,3,2,13.201},
    {28,1,0,27.353},{28,2,0,20.213},{28,2,1,24.095},{28,3,0,14.961},{28,3,1,14.085},{28,4,0,5.711},{28,3,2,13.878},
    {29,1,0,28.339},{29,2,0,21.020},{29,2,1,25.097},{29,3,0,15.594},{29,3,1,14.731},{29,4,0,5.842},{29,3,2,15.093},
    {30,1,0,29.325},{30,2,0,21.828},{30,2,1,26.098},{30,3,0,16.219},{30,3,1,15.369},{30,4,0,5.965},{30,3,2,16.251},
    {31,1,0,30.309},{31,2,0,22.599},{31,2,1,27.091},{31,3,0,16.996},{31,3,1,16.204},{31,4,0,7.067},{31,3,2,17.378},{31,4,1,6.222},
    {32,1,0,31.2937},{32,2,0,23.365},{32,2,1,28.082},{32,3,0,17.790},{32,3,1,17.014},{32,4,0,8.044},{32,3,2,18.477},{32,4,1,6.780},
    {33,1,0,32.278},{33,2,0,24.127},{33,2,1,29.074},{33,3,0,18.596},{33,3,1,17.850},{33,4,0,8.944},{33,3,2,19.559},{33,4,1,7.449},
    {34,1,0,33.262},{34,2,0,24.888},{34,2,1,30.065},{34,3,0,19.403},{34,3,1,18.705},{34,4,0,9.758},{34,4,1,8.287},
    {35,1,0,34.2471},{35,2,0,25.64},{35,2,1,31.06},{35,3,0,20.22},{35,3,1,19.57},{35,4,0,10.55},{35,3,2,19.56},{35,4,1,9.03},
    {36,1,0,35.2316},{36,2,0,26.40},{36,2,1,32.05},{36,3,0,21.03},{36,3,1,20.43},{36,4,0,11.32},{36,3,2,20.63},{36,4,1,9.77},
};
static const int CR_LEN = (int)(sizeof(CR_TABLE) / sizeof(CR_TABLE[0]));

double quantum_zeff_cr(int Z, int n, int l) {
    if (Z < 1 || n < 1 || l < 0 || l >= n) return -1.0;
    for (int i = 0; i < CR_LEN; i++)
        if (CR_TABLE[i].Z == Z && CR_TABLE[i].n == n && CR_TABLE[i].l == l)
            return CR_TABLE[i].z;
    return -1.0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Exact hydrogen-like expectations (Griffiths, analytic)
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_expect_r(int n, int l, double Z_eff) {
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;
    double a0 = BOHR_TO_ANGSTROM;
    return a0 / (2.0 * Z_eff) * (3.0 * n * n - l * (l + 1));
}

double quantum_expect_r2(int n, int l, double Z_eff) {
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;
    double a0 = BOHR_TO_ANGSTROM;
    double nn = (double)(n * n), ll = (double)(l * (l + 1));
    return a0 * a0 * nn / (2.0 * Z_eff * Z_eff) * (5.0 * nn + 1.0 - 3.0 * ll);
}

double quantum_expect_invr(int n, int l, double Z_eff) {
    if (n < 1 || l < 0 || l >= n) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;
    return Z_eff / (BOHR_TO_ANGSTROM * n * n);
}

double quantum_expect_T(int n, int l, double Z_eff) {
    (void)l; /* virial: <T> depends on n only for H-like */
    if (n < 1) return 0.0;
    if (!(Z_eff > 0.0) || !isfinite(Z_eff)) return 0.0;
    return (Z_eff * Z_eff / (2.0 * n * n)) * HARTREE_TO_EV;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Print orbital table for an atom
 * ══════════════════════════════════════════════════════════════════════════ */
void quantum_print_orbitals(const Atom *atom) {
    if (!atom || !atom->element) return;
    static const char sub[] = "spdf";
    printf("  Orbital table for %s (Z=%d)\n",
           atom->element->symbol, atom->Z);
    printf("  %-8s %-6s %-6s %-6s %-12s %-10s\n",
           "Orbital","n","l","ml","Energy(eV)","Occ");
    printf("  %-8s %-6s %-6s %-6s %-12s %-10s\n",
           "-------","--","--","--","----------","---");

    for (int i = 0; i < atom->num_orbitals; i++) {
        const Orbital *o = &atom->orbitals[i];
        char name[8];
        int li = (o->qn.l >= 0 && o->qn.l <= 3) ? o->qn.l : 0;
        snprintf(name, sizeof(name), "%d%c(%+d)",
                 o->qn.n, sub[li], o->qn.ml);
        printf("  %-8s %-6d %-6d %-6d %-12.4f %-10d\n",
               name, o->qn.n, o->qn.l, o->qn.ml,
               o->orbital_energy, o->occupation);
    }
}
