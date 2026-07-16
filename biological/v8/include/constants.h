#ifndef CONSTANTS_H
#define CONSTANTS_H

/*
 * constants.h
 * Physical constants for the carbon chemistry simulator.
 * All SI values are exact (2019 CODATA redefinition where applicable).
 * MD unit system: Length=Angstrom, Time=femtosecond, Energy=eV, Mass=AMU
 */

/* ── Fundamental constants (SI) ─────────────────────────────────────────── */
#define PLANCK_H            6.62607015e-34      /* J·s                       */
#define PLANCK_HBAR         1.054571817e-34     /* J·s  (h / 2π)             */
#define SPEED_OF_LIGHT      2.99792458e8        /* m/s  (exact)              */
#define ELEM_CHARGE         1.602176634e-19     /* C    (exact)              */
#define ELECTRON_MASS       9.1093837015e-31    /* kg                        */
#define PROTON_MASS         1.67262192369e-27   /* kg                        */
#define NEUTRON_MASS        1.67492749804e-27   /* kg                        */
#define BOLTZMANN_K         1.380649e-23        /* J/K  (exact)              */
#define AVOGADRO_N          6.02214076e23       /* mol⁻¹(exact)              */
#define VACUUM_PERMITTIVITY 8.8541878128e-12    /* F/m                       */
#define COULOMB_K           8.9875517923e9      /* N·m²/C²                   */

/* ── Atomic units ────────────────────────────────────────────────────────── */
#define BOHR_RADIUS         5.29177210903e-11   /* m  (a₀)                   */
#define HARTREE_ENERGY      4.3597447222071e-18 /* J  (Eₕ)                   */
#define AMU                 1.66053906660e-27   /* kg (atomic mass unit)     */

/* ── Conversion factors ──────────────────────────────────────────────────── */
#define ANGSTROM_TO_M       1.0e-10
#define M_TO_ANGSTROM       1.0e10
#define EV_TO_J             1.602176634e-19
#define J_TO_EV             6.241509074e18
#define EV_TO_HARTREE       0.0367493224
#define HARTREE_TO_EV       27.211386245988
#define KCAL_MOL_TO_EV      0.043363           /* 1 kcal/mol in eV           */
#define EV_TO_KCAL_MOL      23.060548          /* 1 eV in kcal/mol           */
#define BOHR_TO_ANGSTROM    0.529177210903
#define ANGSTROM_TO_BOHR    1.889726124626

/*
 * MD force-unit conversion factor:
 *   a [Å/fs²] = F [eV/Å] / m [AMU]  × MD_FORCE_CONV
 *
 * Derivation:
 *   1 eV/Å = 1.602176634e-19 J / 1e-10 m = 1.602176634e-9 N
 *   1 AMU  = 1.66053906660e-27 kg
 *   → a [m/s²] = (1.602176634e-9 / 1.66053906660e-27) × (F/m)
 *              = 9.64853322e17 × (F/m)
 *   1 Å/fs² = 1e-10 m / (1e-15 s)² = 1e20 m/s²
 *   → a [Å/fs²] = 9.64853322e17 / 1e20 × (F/m) = 9.64853322e-3 × (F/m)
 */
#define MD_FORCE_CONV       9.64853322e-3

/*
 * Coulomb prefactor in MD units:
 *   V [eV] = COULOMB_MD × q1 [e] × q2 [e] / r [Å]
 *
 * Derivation:
 *   k_e × e² / Å → eV
 *   = 8.9875517923e9 × (1.602176634e-19)² / 1e-10 / 1.602176634e-19
 *   = 14.3996 eV·Å/e²
 */
#define COULOMB_MD          14.3996             /* eV·Å / e²                 */

/* ── Simulation limits ───────────────────────────────────────────────────── */
#define MAX_ELECTRONS       128
#define MAX_ATOMS           100000
#define MAX_BONDS           200000
#define MAX_BONDS_PER_ATOM  8
#define MAX_SHELLS          7
#define MAX_ELEMENTS        118
#define MAX_ORBITALS        32   /* max orbital entries per atom             */

/* ── Quantum number helpers ──────────────────────────────────────────────── */
#define QN_L_S              0
#define QN_L_P              1
#define QN_L_D              2
#define QN_L_F              3
static const char QN_L_NAMES[] = "spdf";

/* orbital capacity = 2*(2l+1) */
#define ORBITAL_CAPACITY(l) (2*(2*(l)+1))

#endif /* CONSTANTS_H */
