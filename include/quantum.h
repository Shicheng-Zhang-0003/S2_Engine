#ifndef QUANTUM_H
#define QUANTUM_H

#include "types.h"

/*
 * quantum.h
 * Quantum mechanical calculations used to initialise and describe atoms.
 *
 * Strategy: hydrogen-like orbital model with Slater effective nuclear charge.
 * This is the standard approximation used in Hartree-Fock starting guesses
 * and gives qualitatively useful orbital shapes and orderings for
 * main-group elements (a diagnostic starting guess, not a substitute for
 * HF/DFT with exchange/correlation; do not cite energies beyond order of
 * magnitude. Slater-Hydrogen valence energies err by factors of ~2-5x,
 * e.g. O 2p computed -70.4 eV vs first IE 13.6 eV, C 2p -35.9 eV vs
 * 11.3 eV: useful for Zeff trends and radial scales, not spectroscopy).
 *
 * For a multi-electron atom with Z protons, each electron experiences an
 * effective nuclear charge Z_eff = Z - S, where S is the Slater screening
 * constant calculated from all other electrons by their proximity.
 *
 * Orbital energy: E_nl = -0.5*HARTREE_TO_EV eV × (Z_eff / n*)² (≈-13.6057; derived in-line, never hand-typed)
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

/* Raw Slater screening (no floor, no stderr): for finite-difference
 * probes like charge-response slopes that legitimately evaluate
 * over-screened censuses. Prefer quantum_zeff() elsewhere. */
double quantum_zeff_raw(int Z, int n, int l, const ElectronConfig *cfg);

/* ── Orbital energy (eV, negative = bound) ───────────────────────────────── */
double quantum_orbital_energy(int Z, int n, int l, const ElectronConfig *cfg);

/* ── Populate orbital array for an atom from its ElectronConfig ──────────── */
/*
 * Fills atom->orbitals[] and sets atom->num_orbitals.
 * Each orbital gets its energy set via quantum_orbital_energy().
 * Orbital ml values are assigned in order: -l, -(l-1), …, +l.
 * One table entry is stored per ml slot (not per electron): ms is always
 * a valid single-electron spin (+0.5 representative); a doubly-occupied
 * slot (occupation==2) holds a +0.5/-0.5 pair whose net is zero - the
 * pair count lives in occupation, not in ms. No energy/force path reads ms.
 */
void quantum_fill_orbitals(Atom *atom);

/*
 * As quantum_fill_orbitals, but reports whether the fixed-size orbital
 * table TRUNCATED (audit fix F14). The Atom carries MAX_ORBITALS=32
 * ml-slots, which covers the whole tabulated element range (Kr needs 18)
 * but not a fully occupied shell sequence up to 7p (60 slots), so the
 * lanthanides onward would silently lose electrons. The old code broke
 * out of the loop with no signal at all. *truncated is set to 1 and a
 * warning is emitted on first truncation.
 */
void quantum_fill_orbitals_checked(Atom *atom, int *truncated);

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
 * Global maximizer of P(r) on (0, 30n^2/Z_eff]: the interval is
 * partitioned into 64 subintervals, each searched by golden section,
 * and the best local maximum is kept - P_nl for n>=2 has n-l-1 nodes
 * and multiple peaks, so a single whole-interval search can land on a
 * secondary maximum instead of the global one.
 */
double quantum_most_probable_radius(int n, int l, double Z_eff);

/* ── Clementi-Raimondi SCF effective charges ───────────────────────────── */
/*
 * SCF-fitted Slater exponents (Clementi & Raimondi, JCP 1963; 1967 for
 * Z>36), dual-sourced: Wikipedia compact table cross-checked cell by
 * cell against WebElements per-element pages (C, O, K, Ar, Br, Kr,
 * Ge-1s all match; WebElements wins on conflict: K-3d dropped as
 * "no data", Se-3d dropped as non-monotonic, Kr row + Ar-3p taken
 * from WebElements). Returns CR Zeff for (Z,n,l), or -1.0 if that
 * orbital has no tabulated exponent (unoccupied/virtual or blank) —
 * caller falls back to Slater. Neutral ground-state atoms only; ions
 * use Slater rules on ion configs (documented limitation).
 *
 * SCOPE HONESTY: CR exponents reproduce SCF radial SHAPES (sizes,
 * ranges), not orbital energies: E=-13.6(Z/n*)^2 with CR charges is
 * still factors off measured IEs. Use CR for spatial quantities
 * (r_mp, overlap zeta, dispersion range), Slater for energy labels.
 */
double quantum_zeff_cr(int Z, int n, int l);

/* ── Exact hydrogen-like expectations (analytic, Griffiths) ─────────────── */
/*
 * For R_nl with nuclear charge Z_eff (a0 in Angstrom internally):
 *   <r>   = a0/(2Z)[3n^2 - l(l+1)]
 *   <r^2> = a0^2 n^2/(2Z^2)[5n^2 + 1 - 3l(l+1)]
 *   <1/r> = Z/(a0 n^2)
 *   <T>   = (Z^2/(2n^2)) Ha in eV (virial theorem, exact for H-like)
 * Return 0 on invalid input. Exact given Zeff — the approximation is
 * Zeff itself, not these formulas.
 */
double quantum_expect_r(int n, int l, double Z_eff);
double quantum_expect_r2(int n, int l, double Z_eff);
double quantum_expect_invr(int n, int l, double Z_eff);
double quantum_expect_T(int n, int l, double Z_eff);

/* ── Print orbital energy table for atom ────────────────────────────────── */
void quantum_print_orbitals(const Atom *atom);

/* ── Associated Laguerre polynomial L_p^q(x) ───────────────────────────── */
/*    (exposed for unit testing)                                             */
double quantum_laguerre(int p, int q, double x);

#endif /* QUANTUM_H */
