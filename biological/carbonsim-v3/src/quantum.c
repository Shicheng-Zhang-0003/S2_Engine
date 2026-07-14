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
 * Table from Slater (1930), Phys. Rev. 36, 57.
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_nstar(int n) {
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
 * Slater screening constant
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_zeff(int Z, int n, int l, const ElectronConfig *cfg) {
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

        /* Determine contribution to S based on Slater's rules */
        if (l <= 1) {
            /* Target is s or p electron */
            if (other_group == target_group) {
                /* Same group: 0.35 (0.30 for 1s) */
                double contrib = (n == 1) ? 0.30 : 0.35;
                S += contrib * electrons;
            } else if (other_group == target_group - 1) {
                /* One group inner */
                S += 0.85 * electrons;
            } else if (other_group < target_group - 1) {
                /* Two or more groups inner */
                S += 1.00 * electrons;
            }
            /* Higher groups do not screen (impossible in ground state) */
        } else {
            /* Target is d or f electron */
            if (other_group == target_group) {
                S += 0.35 * electrons;
            } else if (other_group < target_group) {
                S += 1.00 * electrons;
            }
        }
    }

    double Zeff = (double)Z - S;
    return (Zeff < 1.0) ? 1.0 : Zeff; /* physical minimum */
}

/* ══════════════════════════════════════════════════════════════════════════
 * Orbital energy
 * E_nl = -13.6058 eV × (Z_eff / n*)²
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_orbital_energy(int Z, int n, int l, const ElectronConfig *cfg) {
    double Zeff = quantum_zeff(Z, n, l, cfg);
    double nstar = quantum_nstar(n);
    return -13.6058 * (Zeff / nstar) * (Zeff / nstar);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Fill atom->orbitals[] from its electron configuration
 * ══════════════════════════════════════════════════════════════════════════ */
void quantum_fill_orbitals(Atom *atom) {
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
            orb->qn.ms   = 0.5;        /* representative spin */
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
    if (k < 0)  return 1.0;
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
 *  N  = -sqrt[(2Z_eff/na₀)³ × (n-l-1)! / (2n [(n+l)!]³)]
 *  R_nl(r) = N × e^(-ρ/2) × ρ^l × L_{n-l-1}^{2l+1}(ρ)
 *
 * Return units: Å^(-3/2)
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_radial_wavefunction(int n, int l, double Z_eff, double r_ang) {
    if (r_ang < 0.0) return 0.0;

    double a0_ang = BOHR_TO_ANGSTROM;          /* 0.529177 Å */
    double scale  = 2.0 * Z_eff / ((double)n * a0_ang);
    double rho    = scale * r_ang;

    /* Normalisation constant squared */
    double num    = factorial(n - l - 1);
    double den    = 2.0 * n * pow(factorial(n + l), 3.0);
    double N2     = pow(scale, 3.0) * num / den;
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
 * Most probable radius — golden-section search on P(r) in [0, r_max]
 * ══════════════════════════════════════════════════════════════════════════ */
double quantum_most_probable_radius(int n, int l, double Z_eff) {
    double a = 0.0001, b = 30.0 * n * n / Z_eff;
    static const double PHI = 0.6180339887; /* 1/φ */
    double c = b - PHI * (b - a);
    double d = a + PHI * (b - a);

    for (int iter = 0; iter < 200; iter++) {
        if (fabs(b - a) < 1.0e-8) break;
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

/* ══════════════════════════════════════════════════════════════════════════
 * Print orbital table for an atom
 * ══════════════════════════════════════════════════════════════════════════ */
void quantum_print_orbitals(const Atom *atom) {
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
        snprintf(name, sizeof(name), "%d%c(%+d)",
                 o->qn.n, sub[o->qn.l], o->qn.ml);
        printf("  %-8s %-6d %-6d %-6d %-12.4f %-10d\n",
               name, o->qn.n, o->qn.l, o->qn.ml,
               o->orbital_energy, o->occupation);
    }
}
