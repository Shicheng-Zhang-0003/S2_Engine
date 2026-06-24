#ifndef QUANTUM_H
#define QUANTUM_H

#include "types.h"

/*
 * quantum.h
 * Quantum mechanical calculations used to initialise and describe atoms.
 *
 * Strategy: hydrogen-like orbital model with Slater effective nuclear charge.
 * This is the standard approximation used in Hartree-Fock starting guesses
 * and gives chemically accurate orbital energies for main-group elements.
 *
 * For a multi-electron atom with Z protons, each electron experiences an
 * effective nuclear charge Z_eff = Z - S, where S is the Slater screening
 * constant calculated from all other electrons by their proximity.
 *
 * Orbital energy: E_nl = -13.6058 eV × (Z_eff / n*)²
 * where n* is the effective principal quantum number (Slater 1930).
 */

/* ── Effective principal quantum number (Slater 1930 table) ──────────────── */
double quantum_nstar(int n);

/* ── Slater screening constant ───────────────────────────────────────────── */
/*
 * Computes S for an electron in orbital (n, l) given the full electron
 * configuration. Returns Z_eff = Z - S.
 *
 * Rules:
 *  Group electrons: [1s][2s2p][3s3p][3d][4s4p][4d][4f][5s5p]…
 *  For s/p electron:
 *    same group   : +0.35 per electron (1s: +0.30)
 *    next inner   : +0.85 per electron
 *    deeper inner : +1.00 per electron
 *  For d/f electron:
 *    same group   : +0.35 per electron
 *    all inner    : +1.00 per electron
 */
double quantum_zeff(int Z, int n, int l, const ElectronConfig *cfg);

/* ── Orbital energy (eV, negative = bound) ───────────────────────────────── */
double quantum_orbital_energy(int Z, int n, int l, const ElectronConfig *cfg);

/* ── Populate orbital array for an atom from its ElectronConfig ──────────── */
/*
 * Fills atom->orbitals[] and sets atom->num_orbitals.
 * Each orbital gets its energy set via quantum_orbital_energy().
 * Orbital ml values are assigned in order: -l, -(l-1), …, +l.
 * Spins are assigned: first electron in each orbital gets ms=+0.5,
 * second gets ms=-0.5 (Hund's rule applied per subshell).
 */
void quantum_fill_orbitals(Atom *atom);

/* ── Hydrogen-like radial wave function ──────────────────────────────────── */
/*
 * Returns R_nl(r) for a hydrogen-like atom with effective charge Z_eff.
 * r is in Ångström. Uses the exact analytic formula:
 *
 *   R_nl(r) = N × exp(-ρ/2) × ρ^l × L_(n-l-1)^(2l+1)(ρ)
 *
 * where ρ = 2 Z_eff r / (n a₀), a₀ = 0.529177 Å,
 * N is the normalisation constant, and L_p^q are associated Laguerre
 * polynomials computed by three-term recurrence.
 *
 * Units of return value: Å^(-3/2) (so |R|²r²dr is dimensionless probability).
 */
double quantum_radial_wavefunction(int n, int l, double Z_eff, double r_angstrom);

/* ── Radial probability density P(r) = r² |R_nl(r)|² ───────────────────── */
double quantum_radial_probability(int n, int l, double Z_eff, double r_angstrom);

/* ── Most probable radius for orbital (n,l) with Z_eff ──────────────────── */
/*
 * Numerically finds r_max of P(r) by golden-section search in [0, 20 Å].
 */
double quantum_most_probable_radius(int n, int l, double Z_eff);

/* ── Print orbital energy table for atom ────────────────────────────────── */
void quantum_print_orbitals(const Atom *atom);

/* ── Associated Laguerre polynomial L_p^q(x) ───────────────────────────── */
/*    (exposed for unit testing)                                             */
double quantum_laguerre(int p, int q, double x);

#endif /* QUANTUM_H */
