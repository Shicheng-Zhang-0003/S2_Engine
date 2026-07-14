# PROJECT BUNDLE: carbonsim-v5
Root Directory: /home/magi-01/Desktop/work/programming/architecture/S2_Engine/biological/carbonsim-v5
Generated: Wed Jul  1 02:41:07 PM MDT 2026

## 1. DIRECTORY HIERARCHY
```text
carbonsim-v5/
├── include/
│   ├── constants.h
│   ├── forces.h
│   ├── integrator.h
│   ├── nucleobases.h
│   ├── periodic_table.h
│   ├── quantum.h
│   ├── sim.h
│   ├── types.h
│   └── vec3.h
├── src/
│   ├── forces.c
│   ├── integrator.c
│   ├── main.c
│   ├── nucleobases.c
│   ├── periodic_table.c
│   ├── quantum.c
│   └── sim.c
├── carbonsim
└── makefile
```

## 2. FILE CONTENTS

### FILE: include/constants.h
Location: `include/constants.h`
```h
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

```

---

### FILE: include/forces.h
Location: `include/forces.h`
```h
#ifndef FORCES_H
#define FORCES_H

#include "types.h"

/*
 * forces.h
 * All force and energy calculations for the simulation.
 *
 * Convention throughout:
 *   r_ij = pos_j - pos_i  (vector FROM i TO j)
 *   F_i  accumulated onto atom->force (eV/Å)
 *   All energies returned in eV.
 *
 * Non-bonded cutoff scheme: hard cutoff at sim->cutoff Å.
 * For small gas-phase molecules set cutoff = 100 Å (effectively infinite).
 */

/* ── Pair-interaction results ────────────────────────────────────────────── */
typedef struct {
    double lj_energy;       /* eV  */
    double coulomb_energy;  /* eV  */
} PairEnergy;

/* ── Bond parameter lookup ───────────────────────────────────────────────── */
typedef struct {
    int    Za, Zb;          /* atomic numbers (Za <= Zb for canonical form)  */
    int    order;           /* 1=single, 2=double, 3=triple                  */
    double r0;              /* equilibrium length, Å                         */
    double k;               /* harmonic force constant, eV/Å²               */
} BondParam;

/* ── Angle parameter lookup ──────────────────────────────────────────────── */
typedef struct {
    int    Za, Zb, Zc;      /* Zb is the central atom (Za <= Zc canonical)   */
    double theta0;          /* equilibrium angle, radians                    */
    double k;               /* harmonic force constant, eV/rad²             */
} AngleParam;

/* ── Database accessors ──────────────────────────────────────────────────── */

/*
 * Looks up equilibrium bond parameters for atoms with atomic numbers Za, Zb
 * and bond order `order` (1/2/3). Returns 1 on success, 0 if not found.
 * Falls back to a geometric estimate if no entry in the table.
 */
int forces_bond_params(int Za, int Zb, int order, BondParam *out);

/*
 * Looks up equilibrium angle parameters for a triplet (Za, Zb central, Zc).
 * Returns 1 on success, 0 if not found (caller should use 109.47° default).
 */
int forces_angle_params(int Za, int Zb, int Zc, AngleParam *out);

/* ── Lennard-Jones combining rules (Lorentz-Berthelot) ───────────────────── */
static inline double lj_eps_combine(double ei, double ej) {
    return sqrt(ei * ej);
}
static inline double lj_sigma_combine(double si, double sj) {
    return 0.5 * (si + sj);
}

/* ── Pairwise non-bonded force (accumulates onto both atoms) ─────────────── */
/*
 * Adds Lennard-Jones and Coulomb contributions to atoms[ia].force and
 * atoms[ib].force (Newton's third law applied internally).
 * Returns the energy contributions in `out`.
 * `box` is used for minimum-image PBC; pass NULL for no PBC.
 */
PairEnergy forces_nonbonded_pair(Atom *atoms, int ia, int ib,
                                  const SimBox *box,
                                  int use_lj, int use_coulomb,
                                  double dielectric);

/* ── Harmonic bond force (accumulates onto both endpoint atoms) ───────────── */
/*
 * V = 0.5 k (r - r0)²   → F_a = k(r-r0) r̂_ab,  F_b = −F_a
 * Returns bond potential energy in eV.
 */
double forces_bond(Atom *atoms, const Bond *bond);

/* ── Harmonic angle force (accumulates onto all three atoms) ─────────────── */
/*
 * V = 0.5 k (θ − θ0)²
 * Gradient computed analytically via the chain rule through acos.
 * Returns angle potential energy in eV.
 */
double forces_angle(Atom *atoms, const Angle *angle);

/* ── Master force calculation ────────────────────────────────────────────── */
/*
 * 1. Zeros all atom forces.
 * 2. Loops all pairs within cutoff → non-bonded.
 * 3. Loops all bonds               → bonded stretch.
 * 4. Loops all angles              → bonded bend.
 * Updates sim->potential_energy.
 * O(N²) pair loop; sufficient for small systems. Replace with cell-list
 * or Verlet neighbour list for N > ~500 atoms.
 */
void forces_calculate(Simulation *sim);

/* ── Print force summary ─────────────────────────────────────────────────── */
void forces_print_summary(const Simulation *sim);

#endif /* FORCES_H */

```

---

### FILE: include/integrator.h
Location: `include/integrator.h`
```h
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "types.h"

/*
 * integrator.h
 * Velocity Verlet molecular dynamics integrator with optional thermostat.
 *
 * Velocity Verlet (Swope et al., 1982) — time-reversible, symplectic:
 *   r(t+dt) = r(t) + v(t)dt + 0.5 a(t)dt²
 *   v(t+dt) = v(t) + 0.5[a(t) + a(t+dt)]dt
 *
 * In practice, split into two half-steps:
 *   Step A (before force calc):  v(t+dt/2) = v(t) + 0.5 a(t) dt
 *                                r(t+dt)   = r(t) + v(t+dt/2) dt
 *   Force recalculation at r(t+dt) → a(t+dt)
 *   Step B (after force calc):   v(t+dt)   = v(t+dt/2) + 0.5 a(t+dt) dt
 *
 * Units throughout:
 *   position : Å          velocity : Å/fs
 *   force    : eV/Å       mass     : AMU
 *   time     : fs         energy   : eV
 *   temperature : K
 *
 * Acceleration: a [Å/fs²] = F [eV/Å] / m [AMU] × MD_FORCE_CONV
 * where MD_FORCE_CONV = 9.64853322e-3 (see constants.h).
 */

/* ── Half-step A: kick velocities, drift positions ───────────────────────── */
/*
 * v(t+dt/2) = v(t) + 0.5 a(t) dt
 * r(t+dt)   = r(t) + v(t+dt/2) dt
 * (called BEFORE the force recalculation)
 */
void integrator_kick_drift(Simulation *sim);

/* ── Half-step B: final velocity update ──────────────────────────────────── */
/*
 * v(t+dt) = v(t+dt/2) + 0.5 a(t+dt) dt
 * (called AFTER the force recalculation at the new positions)
 */
void integrator_kick(Simulation *sim);

/* ── Full velocity Verlet step (wraps both halves + force call) ──────────── */
/*
 * Calls kick_drift, then forces_calculate, then kick.
 * Updates sim->step, sim->time, sim->kinetic_energy, sim->total_energy.
 */
void integrator_step(Simulation *sim);

/* ── Thermodynamics ──────────────────────────────────────────────────────── */
/*
 * Kinetic energy: KE = 0.5 Σ m_i |v_i|²
 * Units: AMU × (Å/fs)² → converted to eV via EV_CONV = 1/(2×MD_FORCE_CONV)
 *
 * Exact: 1 AMU × (Å/fs)² = 1.66053906660e-27 kg × (1e-10/1e-15)² m²/s²
 *      = 1.66053906660e-27 × 1e10 J = 1.66053906660e-17 J
 *      = 1.66053906660e-17 / 1.602176634e-19 eV = 103.6427 eV
 * So: KE [eV] = 0.5 Σ m[AMU] v²[Å²/fs²] × 103.6427
 */
#define AMU_AFS2_TO_EV   103.6427   /* (Å/fs)² AMU → eV */

double integrator_kinetic_energy(const Simulation *sim);

/*
 * Temperature from equipartition: KE = (3N/2) k_B T
 * T [K] = 2 KE [eV] / (3 N × k_B_eV)
 * k_B in eV/K = 8.617333262e-5
 */
double integrator_temperature(const Simulation *sim);

/* ── Berendsen thermostat ────────────────────────────────────────────────── */
/*
 * Rescales velocities to approach target temperature T_0 with time
 * constant τ (fs):
 *   λ = sqrt(1 + (dt/τ)(T_0/T − 1))
 *   v_new = λ × v
 *
 * Weak coupling: does not produce a rigorous NVT ensemble but is stable
 * and appropriate for equilibration. τ = 100 fs is a typical choice.
 */
void integrator_berendsen(Simulation *sim);

/* ── Maxwell-Boltzmann velocity initialisation ───────────────────────────── */
/*
 * Assigns velocities drawn from the Maxwell-Boltzmann distribution at
 * temperature T_init [K].
 *
 * For each component: v_x ~ N(0, sqrt(k_B T / m))
 * Implemented with Box-Muller transform from two uniform random numbers.
 *
 * After sampling, removes net linear momentum (centre of mass velocity = 0).
 * `seed` initialises the internal RNG (use 0 for a fixed default seed).
 */
void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
                                   unsigned long seed);

/* ── Remove centre-of-mass drift ─────────────────────────────────────────── */
void integrator_remove_com_velocity(Simulation *sim);

/* ── Print step summary ──────────────────────────────────────────────────── */
void integrator_print_step(const Simulation *sim);

#endif /* INTEGRATOR_H */

```

---

### FILE: include/nucleobases.h
Location: `include/nucleobases.h`
```h
#ifndef NUCLEOBASES_H
#define NUCLEOBASES_H

#include "types.h"

/*
 * nucleobases.h
 *
 * Constructors for the five DNA/RNA nucleobases: uracil, cytosine,
 * thymine, adenine, guanine.
 *
 * PROVENANCE OF GEOMETRY:
 * Seed atomic positions for uracil, cytosine, adenine, and guanine are
 * taken from the RCSB PDB Chemical Component Dictionary's "ideal
 * coordinates" (ligand codes URA, CYT, ADE, GUN respectively), fetched
 * directly from https://files.rcsb.org/ligands/view/<CODE>.cif on
 * 2026-06-19. These are algorithmically-generated (OpenEye OMEGA /
 * Molecular Networks Corina), chemically-sensible 3D embeddings of the
 * verified 2D connectivity - they are NOT raw spectroscopic or
 * crystallographic measurements, but they do encode correct bond
 * connectivity, correct ring topology, and reasonable bond lengths/
 * angles for every atom, independently of any derivation in this file.
 *
 * Thymine has no standalone PDB ligand entry under the code "THY" (that
 * code is taken by an unrelated thiamine diphosphate derivative - a
 * PDB code collision). Thymine is instead constructed by taking
 * uracil's own verified C5-H5 bond vector and replacing that hydrogen
 * with a tetrahedral methyl group at the standard aromatic-ring-to-
 * methyl C-C bond length (1.51 Angstrom) - this is the correct
 * chemical relationship (thymine = 5-methyluracil) and avoids
 * depending on a second, harder-to-verify external coordinate source.
 *
 * EQUILIBRIUM GEOMETRY CONSTRUCTION:
 * Every bond's r0 and every angle's theta0 are set to EXACTLY match the
 * placed Cartesian seed geometry (computed at construction time via
 * vec3_dist / vec3_angle, not hardcoded) - so these molecules start at
 * zero strain energy, and a static geometry report computed from the
 * positions will exactly reproduce the seed bond lengths/angles.
 * Force constants (k) are reasonable generic single/double-bond and
 * aromatic-ring-bending values (same order of magnitude as the rest of
 * this codebase's AMBER-derived table), not independently fitted to
 * spectroscopic data for these specific rings. Partial charges are left
 * at 0.0 - accurate RESP/AMBER nucleobase charges (needed for
 * quantitative base-pairing energetics) are deferred to a follow-up;
 * fabricating precise-looking but unverified per-atom charges here
 * would be worse than leaving them at a clearly-flagged placeholder.
 *
 * All constructors return the index of the first atom added, atoms
 * follow in PDB atom-naming order as commented in each function.
 */

int sim_place_uracil   (Simulation *sim, Vec3 origin);
int sim_place_cytosine (Simulation *sim, Vec3 origin);
int sim_place_thymine  (Simulation *sim, Vec3 origin);
int sim_place_adenine  (Simulation *sim, Vec3 origin);
int sim_place_guanine  (Simulation *sim, Vec3 origin);

/*
 * Fit a least-squares plane through the given atom indices and return
 * the maximum perpendicular distance (Angstrom) of any atom from that
 * plane. Aromatic rings should be very close to planar (a few
 * hundredths of an Angstrom at most) - this is an independent physical
 * sanity check on the placed geometry that doesn't depend on any
 * external reference numbers.
 */
double nb_planarity_deviation(const Simulation *sim,
                               const int *atom_indices, int n);

/*
 * Compute the unit normal vector of the plane defined by the first 3
 * given atom indices (cross product of the two edge vectors from the
 * first atom). For a near-planar ring, this is a good approximation
 * of the ring's plane normal, useful for rigid-body alignment of one
 * planar molecule against another (e.g. orienting a base pair).
 */
Vec3 nb_ring_normal(const Simulation *sim, const int *atom_indices);

/*
 * Signed angle to rotate v_from onto v_to, as a rotation around axis n
 * (right-hand rule), after projecting both onto the plane perpendicular
 * to n. Needed to fix the azimuthal orientation left unconstrained by
 * plane-normal alignment alone.
 */
double nb_signed_inplane_angle(Vec3 v_from, Vec3 v_to, Vec3 n);

/*
 * Rigidly rotate (around `pivot`, by `angle` radians about `axis`) and
 * then translate (by `translation`) every atom in [first_atom,
 * first_atom + n_atoms). Since rotation+translation is an isometry,
 * this does not disturb any bond length or angle already set up for
 * this group of atoms - safe to call on an already-bonded molecule.
 */
void nb_transform_rigid(Simulation *sim, int first_atom, int n_atoms,
                         Vec3 pivot, Vec3 axis, double angle,
                         Vec3 translation);

#endif /* NUCLEOBASES_H */

```

---

### FILE: include/periodic_table.h
Location: `include/periodic_table.h`
```h
#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include "types.h"

/*
 * periodic_table.h
 * Immutable database of element properties used by the simulator.
 * Elements 1–36 (H → Kr) are fully populated.
 * Elements 37–118 are stubbed with mass only — extend as needed.
 *
 * Data sources:
 *   Pauling electronegativities : IUPAC 2013
 *   Atomic/covalent/vdW radii  : Alvarez 2013 (Dalton Trans.)
 *   Ionisation energies        : NIST Atomic Spectra Database
 *   Electron affinities        : Andersen 1999 (J. Phys. Chem. Ref. Data)
 *   LJ parameters              : Rappé 1992 UFF (J. Am. Chem. Soc.)
 */

/* The global table — indexed by atomic number Z (1-based, element[0] unused) */
extern const Element PERIODIC_TABLE[MAX_ELEMENTS + 1];

/* ── Accessors ───────────────────────────────────────────────────────────── */

/* Returns pointer to element data for atomic number Z, or NULL if invalid. */
const Element *pt_element(int Z);

/* Look up by symbol string (case-sensitive, e.g. "C", "Fe"). */
const Element *pt_by_symbol(const char *symbol);

/* Populate a ground-state ElectronConfig for atomic number Z. */
void pt_electron_config(int Z, ElectronConfig *cfg);

/* Pretty-print electron configuration string, e.g. "1s2 2s2 2p2" */
void pt_config_string(const ElectronConfig *cfg, char *buf, int buflen);

/* Print a summary of one element to stdout */
void pt_print_element(int Z);

#endif /* PERIODIC_TABLE_H */

```

---

### FILE: include/quantum.h
Location: `include/quantum.h`
```h
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

```

---

### FILE: include/sim.h
Location: `include/sim.h`
```h
#ifndef SIM_H
#define SIM_H

#include "types.h"
#include "forces.h"

/*
 * sim.h
 * Simulation construction and topology management.
 *
 * Workflow:
 *   1. sim_create()          → allocate empty simulation
 *   2. sim_add_atom()        → place atoms at positions with charges
 *   3. sim_add_bond()        → define covalent bonds (or sim_detect_bonds)
 *   4. sim_build_angles()    → derive angle list from bond topology
 *   5. forces_calculate()    → initial forces
 *   6. integrator_step() ×N  → run the dynamics
 *   7. sim_destroy()         → free memory
 */

/* ── Lifecycle ───────────────────────────────────────────────────────────── */

/*
 * Allocate a new Simulation with pre-allocated capacity.
 * Returns NULL on allocation failure.
 * Default settings: cutoff=12Å, dt=1fs, no PBC, Berendsen thermostat off.
 */
Simulation *sim_create(int atom_capacity, int bond_capacity);

/* Free all memory owned by sim. */
void sim_destroy(Simulation *sim);

/* ── Atom management ─────────────────────────────────────────────────────── */

/*
 * Add an atom of element Z at position `pos` (Å) with partial charge `q` (e).
 * Copies element data pointer from the periodic table.
 * Returns atom index (≥0) on success, SIM_ERR_* on failure.
 * Automatically fills electron configuration and orbital data.
 */
int sim_add_atom(Simulation *sim, int Z, Vec3 pos, double partial_charge);

/* Convenience: add atom from element symbol string. */
int sim_add_atom_sym(Simulation *sim, const char *symbol, Vec3 pos, double q);

/*
 * Override an atom's Lennard-Jones parameters away from the generic
 * element/UFF default. Use this when an atom is part of a molecule
 * modelled by a specific, internally-consistent force field (e.g. TIP3P
 * water) whose LJ parameters were fit jointly with its charges and
 * differ from the generic per-element UFF values.
 */
void sim_set_atom_lj(Simulation *sim, int atom_idx, double epsilon, double sigma);

/*
 * Override a specific bond's harmonic parameters after creation.
 * Use this when the generic element+order table lookup in forces_bond_params
 * doesn't distinguish context (e.g. aromatic ring C-N vs. generic aliphatic
 * C-N have very different r0, but both are "carbon-nitrogen single bond").
 * bond_idx is the value returned by sim_add_bond.
 */
void sim_set_bond_params(Simulation *sim, int bond_idx, double r0, double k);

/*
 * Append a fully explicit angle (atoms a-b-c, b central) bypassing the
 * generic forces_angle_params table lookup entirely. Use this for ring
 * systems where every angle needs its own verified, context-specific
 * value rather than a generic element-keyed default.
 * Returns the angle index, or SIM_ERR_OVERFLOW if capacity is exceeded.
 */
int sim_add_angle_explicit(Simulation *sim, int a, int b, int c,
                            double theta0, double k);

/* ── Bond management ─────────────────────────────────────────────────────── */

/*
 * Add an explicit bond between atoms ia and ib with given order (1/2/3).
 * Looks up equilibrium parameters from the force field database.
 * Returns bond index (≥0) or SIM_ERR_*.
 */
int sim_add_bond(Simulation *sim, int ia, int ib, int order);

/*
 * Auto-detect bonds by interatomic distance:
 *   bond if r < BOND_TOLERANCE × (r_cov_a + r_cov_b)
 * BOND_TOLERANCE = 1.15 (15% over covalent sum, handles slight distortions).
 * Bond order estimated from distance ratio (rough heuristic).
 * Call after all atoms are placed.
 */
#define BOND_TOLERANCE 1.15
void sim_detect_bonds(Simulation *sim);

/*
 * Build the angle list from the bond topology.
 * For every pair of bonds sharing an atom b, adds angle (a, b, c).
 * Call after all bonds are defined.
 */
void sim_build_angles(Simulation *sim);

/*
 * Like sim_build_angles, but derives each angle's equilibrium theta0
 * directly from the atoms' CURRENT Cartesian positions rather than the
 * generic element-keyed lookup table. Use this for ring systems (e.g.
 * aromatic nucleobases) where the seed geometry already encodes the
 * correct, context-specific angle and the generic per-element table
 * would give the wrong value (it can't distinguish an aromatic ring
 * angle from a generic sp3 one keyed on the same three atomic numbers).
 * All angles get the same force constant k_default (eV/rad^2) - a
 * reasonable generic aromatic-ring bending stiffness, not independently
 * fitted per angle.
 */
void sim_build_angles_geometric(Simulation *sim, double k_default);

/* ── Periodic boundary conditions ────────────────────────────────────────── */

/* Set simulation box dimensions (Å) and enable PBC on all three axes. */
void sim_set_box(Simulation *sim, double lx, double ly, double lz);

/* Wrap all atom positions into the primary box [0, L). */
void sim_wrap_positions(Simulation *sim);

/* ── Pre-built molecule constructors ─────────────────────────────────────── */

/*
 * Each constructor places the molecule with its centre of mass at `origin`
 * and returns the index of the first atom added.
 *
 * Geometry: exact experimental or optimised values.
 * Partial charges: standard force-field assignments (TIP3P for water,
 *   OPLS-AA for organics).
 */

/* H₂ — bond length 0.741 Å */
int sim_place_h2(Simulation *sim, Vec3 origin);

/* H₂O — TIP3P geometry: r(OH)=0.9572Å, ∠HOH=104.52° */
int sim_place_h2o(Simulation *sim, Vec3 origin);

/* NH₃ — pyramidal: r(NH)=1.012Å, ∠HNH=106.67°, q(N)=-1.02e */
int sim_place_nh3(Simulation *sim, Vec3 origin);

/* CH₄ — tetrahedral: r(CH)=1.090Å, ∠HCH=109.47° */
int sim_place_ch4(Simulation *sim, Vec3 origin);

/* CO₂ — linear: r(C=O)=1.163Å */
int sim_place_co2(Simulation *sim, Vec3 origin);

/* ── Diagnostics ─────────────────────────────────────────────────────────── */

/* Print atom list with positions, velocities, charges. */
void sim_print_atoms(const Simulation *sim);

/* Print bond list with current lengths and energies. */
void sim_print_bonds(const Simulation *sim);

/* Print angle list. */
void sim_print_angles(const Simulation *sim);

/* Print full system summary. */
void sim_print_summary(const Simulation *sim);

#endif /* SIM_H */

```

---

### FILE: include/types.h
Location: `include/types.h`
```h
#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include "constants.h"
#include "vec3.h"

/*
 * types.h
 * All data structures for the chemistry simulator, organised bottom-up:
 *   QuantumNumbers → ElectronConfig → Element (periodic table entry)
 *   → Atom → Bond → Simulation
 *
 * Design principle: every struct is the minimal faithful representation of
 * its physical counterpart. No magic numbers, no opaque indices — every
 * field corresponds to a measurable physical quantity with documented units.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Quantum layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * The four quantum numbers that uniquely identify an electron's state
 * in a hydrogen-like orbital.
 *
 *  n  : principal quantum number     1, 2, 3, …
 *  l  : azimuthal quantum number     0 … n-1    (0=s, 1=p, 2=d, 3=f)
 *  ml : magnetic quantum number      -l … +l
 *  ms : spin quantum number          +0.5 (↑) or -0.5 (↓)
 */
typedef struct {
    int    n;
    int    l;
    int    ml;
    double ms;
} QuantumNumbers;

/*
 * A single atomic orbital, characterised by quantum numbers, its
 * current energy, and how many electrons occupy it (0, 1, or 2).
 * orbital_energy is in eV, computed via Slater screening.
 */
typedef struct {
    QuantumNumbers qn;
    double         orbital_energy;  /* eV  (negative = bound)               */
    int            occupation;      /* 0, 1, or 2                            */
} Orbital;

/*
 * Electron configuration stored as a 2D array.
 *   config[shell][subshell] = electron count
 *   shell    : 0-indexed principal quantum number n-1  (n=1..7)
 *   subshell : 0=s, 1=p, 2=d, 3=f
 *
 * Example — Carbon (Z=6): 1s² 2s² 2p²
 *   config[0][0] = 2   (1s)
 *   config[1][0] = 2   (2s)
 *   config[1][1] = 2   (2p)
 */
typedef struct {
    int config[MAX_SHELLS][4];
    int total_electrons;
    int valence_electrons;
} ElectronConfig;

/* ══════════════════════════════════════════════════════════════════════════
 * Periodic table layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * Immutable element data loaded once from the periodic table database.
 * All lengths in Ångström, energies in eV, masses in AMU.
 * lj_epsilon and lj_sigma are UFF Lennard-Jones parameters.
 */
typedef struct {
    int    Z;                        /* atomic number                        */
    char   symbol[4];                /* e.g. "C", "Fe"                      */
    char   name[32];                 /* e.g. "Carbon"                        */
    double mass;                     /* AMU                                  */
    double electronegativity;        /* Pauling scale (-1 = unknown)         */
    double atomic_radius;            /* Å (van der Waals)                    */
    double covalent_radius;          /* Å                                    */
    double vdw_radius;               /* Å                                    */
    double ionization_energy;        /* eV (first ionisation)                */
    double electron_affinity;        /* eV (0 if negative / noble gas)       */
    int    valence;                  /* common valence (dominant)            */
    ElectronConfig ground_config;
    double lj_epsilon;               /* eV  (UFF well depth)                 */
    double lj_sigma;                 /* Å   (UFF collision diameter)         */
} Element;

/* ══════════════════════════════════════════════════════════════════════════
 * Molecular dynamics layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * An atom in the simulation.
 * Classical MD representation: position / velocity / force + electronic state.
 *
 * Units:
 *   position : Å
 *   velocity : Å/fs
 *   force    : eV/Å
 *   mass     : AMU
 *   charge   : elementary charge e
 */
typedef struct {
    int            id;
    int            Z;                /* atomic number                        */
    const Element *element;          /* pointer into the periodic table      */

    /* Classical MD state */
    Vec3   position;                 /* Å                                    */
    Vec3   velocity;                 /* Å/fs                                 */
    Vec3   force;                    /* eV/Å  (accumulated each step)        */
    double mass;                     /* AMU                                  */

    /* Electronic state */
    ElectronConfig electron_config;
    double         partial_charge;   /* e   (partial charge for Coulomb)     */
    int            formal_charge;    /* integer formal charge                */

    /*
     * Per-atom LJ parameters. Default to element->lj_epsilon/lj_sigma (UFF)
     * at creation time, but OVERRIDABLE per atom-type-in-context.
     *
     * This matters because real force-field parameters are not pure
     * element properties — they depend on chemical context. TIP3P's
     * water oxygen (σ=3.1507 Å, ε=0.1521 kcal/mol) is a DIFFERENT atom
     * type from generic UFF oxygen (σ=3.500 Å), even though both are
     * element 8. Mixing them inconsistently (TIP3P charges + UFF LJ)
     * produces unphysical results — see sim_place_h2o() for the fix.
     * TIP3P additionally gives hydrogen NO Lennard-Jones term at all
     * (σ=ε=0); only oxygen carries the LJ site in the original 1983 model.
     */
    double         lj_epsilon;       /* eV  (atom-type-specific, not just element default) */
    double         lj_sigma;         /* Å                                    */

    /* Orbital data (populated by quantum module for small systems) */
    Orbital        orbitals[MAX_ORBITALS];
    int            num_orbitals;

    /* Bonding topology */
    int            bond_partners[MAX_BONDS_PER_ATOM];
    int            bond_orders[MAX_BONDS_PER_ATOM];   /* 1,2,3              */
    int            num_bonds;
} Atom;

/*
 * A covalent bond between two atoms.
 * The harmonic potential is: V = 0.5 * k * (r - r0)²
 *   k  : force constant   eV/Å²
 *   r0 : equilibrium bond length   Å
 */
typedef struct {
    int    atom_a;
    int    atom_b;
    int    order;               /* 1=single, 2=double, 3=triple, 0=aromatic */
    double r0;                  /* equilibrium length, Å                     */
    double k;                   /* force constant, eV/Å²                    */
    double current_length;      /* Å  (updated each step)                   */
    double energy;              /* eV (updated each step)                   */
} Bond;

/*
 * A bond angle between three atoms (a–b–c, b is the central atom).
 * V = 0.5 * k * (θ - θ0)²
 */
typedef struct {
    int    atom_a, atom_b, atom_c;
    double theta0;              /* equilibrium angle, radians                */
    double k;                   /* force constant, eV/rad²                  */
} Angle;

/*
 * A dihedral torsion between four atoms (a–b–c–d).
 * V = k * (1 + cos(n*φ - δ))
 */
typedef struct {
    int    atom_a, atom_b, atom_c, atom_d;
    double k;                   /* eV                                        */
    int    n;                   /* multiplicity                              */
    double delta;               /* phase offset, radians                     */
} Dihedral;

/* ══════════════════════════════════════════════════════════════════════════
 * Simulation box
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    Vec3 dimensions;            /* Å — set to a very large value for vacuum  */
    int  periodic[3];           /* 1 = periodic boundary in that axis        */
} SimBox;

/* ══════════════════════════════════════════════════════════════════════════
 * Force-field parameters for non-bonded interactions
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    FF_NONE   = 0,              /* bare Coulomb only                         */
    FF_UFF    = 1,              /* Universal Force Field (LJ params)         */
    FF_OPLS   = 2,              /* OPLS-AA (for organics)                    */
    FF_CUSTOM = 3
} ForceFieldType;

/* ══════════════════════════════════════════════════════════════════════════
 * Thermostat and barostat
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    THERMOSTAT_NONE      = 0,
    THERMOSTAT_BERENDSEN = 1,   /* simple velocity rescaling                 */
    THERMOSTAT_NOSE_HOOVER = 2  /* proper NVT ensemble                       */
} ThermostatType;

typedef struct {
    ThermostatType type;
    double target_temperature; /* K                                          */
    double tau;                /* coupling time constant, fs                 */
    double xi;                 /* Nosé-Hoover friction variable              */
    double Q;                  /* Nosé-Hoover mass parameter                 */
} Thermostat;

/* ══════════════════════════════════════════════════════════════════════════
 * Top-level simulation state
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    /* Particles */
    Atom    *atoms;
    int      num_atoms;
    int      capacity_atoms;

    /* Bonded topology */
    Bond     *bonds;
    int       num_bonds;
    int       capacity_bonds;

    Angle    *angles;
    int       num_angles;
    int       capacity_angles;

    Dihedral *dihedrals;
    int       num_dihedrals;
    int       capacity_dihedrals;

    /* Box */
    SimBox   box;

    /* Time */
    double   time;              /* current time, fs                          */
    double   dt;                /* timestep, fs                              */
    uint64_t step;

    /* Thermodynamics */
    double   temperature;       /* K  (instantaneous)                        */
    double   kinetic_energy;    /* eV                                        */
    double   potential_energy;  /* eV                                        */
    double   total_energy;      /* eV                                        */
    double   virial;            /* eV (for pressure calculation)             */
    double   E_lj_total;        /* eV - LJ component of potential_energy     */
    double   E_coulomb_total;   /* eV - Coulomb component                    */

    /* Force field */
    ForceFieldType ff_type;
    double         cutoff;      /* non-bonded cutoff, Å                      */
    int            use_coulomb;
    int            use_lj;
    int            use_bonds;
    int            use_angles;
    double         dielectric;  /* relative permittivity divisor for the
                                  * Coulomb term, default 1.0 (no screening).
                                  * Standard technique in classical, non-
                                  * polarizable force fields to approximate
                                  * condensed-phase electronic + orientational
                                  * screening when no explicit solvent is
                                  * present - RESP/AMBER-style charges are
                                  * typically fit assuming this kind of
                                  * environment exists, and overestimate
                                  * interaction strength without it. */

    /* Thermostat */
    Thermostat thermostat;
} Simulation;

/* ══════════════════════════════════════════════════════════════════════════
 * Return codes
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    SIM_OK               =  0,
    SIM_ERR_ALLOC        = -1,
    SIM_ERR_BADATOM      = -2,
    SIM_ERR_BADBOND      = -3,
    SIM_ERR_OVERFLOW     = -4,
    SIM_ERR_BADPARAM     = -5
} SimError;

#endif /* TYPES_H */

```

---

### FILE: include/vec3.h
Location: `include/vec3.h`
```h
#ifndef VEC3_H
#define VEC3_H

#include <math.h>
#include <stdio.h>

/*
 * vec3.h
 * Inline 3D vector arithmetic in double precision.
 * All operations are branchless and compile to efficient SIMD on modern GCC/Clang.
 */

typedef struct { double x, y, z; } Vec3;

/* ── Construction ────────────────────────────────────────────────────────── */
static inline Vec3 vec3(double x, double y, double z) {
    return (Vec3){x, y, z};
}
static inline Vec3 vec3_zero(void) { return (Vec3){0.0, 0.0, 0.0}; }
static inline Vec3 vec3_ones(void) { return (Vec3){1.0, 1.0, 1.0}; }

/* ── Arithmetic ──────────────────────────────────────────────────────────── */
static inline Vec3 vec3_add(Vec3 a, Vec3 b) {
    return (Vec3){a.x+b.x, a.y+b.y, a.z+b.z};
}
static inline Vec3 vec3_sub(Vec3 a, Vec3 b) {
    return (Vec3){a.x-b.x, a.y-b.y, a.z-b.z};
}
static inline Vec3 vec3_scale(Vec3 a, double s) {
    return (Vec3){a.x*s, a.y*s, a.z*s};
}
static inline Vec3 vec3_negate(Vec3 a) {
    return (Vec3){-a.x, -a.y, -a.z};
}
static inline Vec3 vec3_mul(Vec3 a, Vec3 b) {   /* element-wise product */
    return (Vec3){a.x*b.x, a.y*b.y, a.z*b.z};
}

/* ── Dot / cross / norm ──────────────────────────────────────────────────── */
static inline double vec3_dot(Vec3 a, Vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline double vec3_norm2(Vec3 a) {        /* squared length           */
    return a.x*a.x + a.y*a.y + a.z*a.z;
}
static inline double vec3_norm(Vec3 a) {
    return sqrt(vec3_norm2(a));
}
static inline Vec3 vec3_cross(Vec3 a, Vec3 b) {
    return (Vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

/* ── Normalise (returns zero vector if near-zero input) ──────────────────── */
static inline Vec3 vec3_normalize(Vec3 a) {
    double n = vec3_norm(a);
    if (n < 1.0e-300) return vec3_zero();
    return vec3_scale(a, 1.0/n);
}

/* ── Distance between two points ─────────────────────────────────────────── */
static inline double vec3_dist(Vec3 a, Vec3 b) {
    return vec3_norm(vec3_sub(b, a));
}
static inline double vec3_dist2(Vec3 a, Vec3 b) {
    Vec3 d = vec3_sub(b, a);
    return vec3_norm2(d);
}

/* ── Minimum-image convention for periodic boundaries ───────────────────── */
static inline Vec3 vec3_pbc(Vec3 dr, Vec3 box) {
    Vec3 r;
    r.x = dr.x - box.x * round(dr.x / box.x);
    r.y = dr.y - box.y * round(dr.y / box.y);
    r.z = dr.z - box.z * round(dr.z / box.z);
    return r;
}

/* ── Angle between two vectors (radians) ─────────────────────────────────── */
static inline double vec3_angle(Vec3 a, Vec3 b) {
    double c = vec3_dot(a, b) / (vec3_norm(a) * vec3_norm(b));
    /* clamp to [-1,1] for numerical safety */
    if (c >  1.0) c =  1.0;
    if (c < -1.0) c = -1.0;
    return acos(c);
}

/* ── Dihedral angle (radians) ────────────────────────────────────────────── */
static inline double vec3_dihedral(Vec3 b1, Vec3 b2, Vec3 b3) {
    Vec3 n1 = vec3_cross(b1, b2);
    Vec3 n2 = vec3_cross(b2, b3);
    double cos_phi = vec3_dot(n1, n2) / (vec3_norm(n1) * vec3_norm(n2));
    if (cos_phi >  1.0) cos_phi =  1.0;
    if (cos_phi < -1.0) cos_phi = -1.0;
    return acos(cos_phi);
}

/* ── Accumulate (v += a) ─────────────────────────────────────────────────── */
static inline void vec3_iadd(Vec3 *v, Vec3 a) {
    v->x += a.x; v->y += a.y; v->z += a.z;
}
static inline void vec3_isub(Vec3 *v, Vec3 a) {
    v->x -= a.x; v->y -= a.y; v->z -= a.z;
}
static inline void vec3_iscale(Vec3 *v, double s) {
    v->x *= s; v->y *= s; v->z *= s;
}

/* ── Printing ─────────────────────────────────────────────────────────────── */
static inline void vec3_print(const char *label, Vec3 v) {
    printf("%s: (%.6f, %.6f, %.6f)\n", label, v.x, v.y, v.z);
}

/* ── Rodrigues' rotation formula ─────────────────────────────────────────── */
/*
 * Rotates vector v around unit axis `axis` by `angle` radians (right-hand
 * rule). Standard formula: v_rot = v*cos(t) + (axis x v)*sin(t)
 *                                   + axis*(axis.v)*(1-cos(t))
 * `axis` need not be pre-normalised; this function normalises it.
 */
static inline Vec3 vec3_rotate_axis_angle(Vec3 v, Vec3 axis, double angle) {
    axis = vec3_normalize(axis);
    double c = cos(angle), s = sin(angle);
    Vec3 term1 = vec3_scale(v, c);
    Vec3 term2 = vec3_scale(vec3_cross(axis, v), s);
    Vec3 term3 = vec3_scale(axis, vec3_dot(axis, v) * (1.0 - c));
    return vec3_add(vec3_add(term1, term2), term3);
}

#endif /* VEC3_H */

```

---

### FILE: makefile
Location: `makefile`
```text
CC      = gcc
CFLAGS  = -O3 -g -Wall -Wextra -std=c11 -Iinclude
CFLAGS += -Wno-missing-braces -march=native -ffast-math

# Uncomment for debug build with address sanitiser:
# CFLAGS = -g -O0 -Wall -std=c11 -Iinclude -fsanitize=address,undefined -lm

SRC_DIR = src
OBJ_DIR = build
BIN     = carbonsim

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

.PHONY: all clean run

all: $(OBJ_DIR) $(BIN)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN): $(OBJS)
	$(CC) -O3 -o $@ $^ -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: all
	./$(BIN)

clean:
	rm -rf $(OBJ_DIR) $(BIN)

```

---

### FILE: src/forces.c
Location: `src/forces.c`
```cpp
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/forces.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"

/*
 * forces.c
 *
 * Every number here corresponds to a real physical measurement.
 * Bond parameters from AMBER ff14SB / CHARMM36 (converted to eV/Å²).
 * Angle parameters from the same sources.
 *
 * AMBER uses kcal/mol/Å² and kcal/mol/rad²; we multiply by KCAL_MOL_TO_EV.
 * All lengths in Å, energies in eV, forces in eV/Å.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Bond parameter database
 * Ordered: Za <= Zb.  order = 1 (single), 2 (double), 3 (triple).
 * k values in eV/Å²; r0 in Å.
 *
 * Conversion: AMBER k_bond [kcal/mol/Å²] × KCAL_MOL_TO_EV = eV/Å²
 * ══════════════════════════════════════════════════════════════════════════ */
static const BondParam BOND_TABLE[] = {
/*   Za  Zb  ord   r0(Å)    k(eV/Å²)   element pair        */
    { 1,  1,  1,  0.7414,  36.00 },  /* H–H    */
    { 1,  6,  1,  1.0900,  29.30 },  /* H–C    */
    { 1,  7,  1,  1.0120,  31.60 },  /* H–N    */
    { 1,  8,  1,  0.9572,  34.50 },  /* H–O    */
    { 1,  9,  1,  0.9170,  40.00 },  /* H–F    */
    { 1, 16,  1,  1.3400,  17.00 },  /* H–S    */
    { 1, 17,  1,  1.2740,  18.00 },  /* H–Cl   */
    { 6,  6,  1,  1.5400,  23.80 },  /* C–C    */
    { 6,  6,  2,  1.3400,  47.60 },  /* C=C    */
    { 6,  6,  3,  1.2040,  95.20 },  /* C≡C    */
    { 6,  7,  1,  1.4700,  22.50 },  /* C–N    */
    { 6,  7,  2,  1.2740,  50.00 },  /* C=N    */
    { 6,  7,  3,  1.1570, 100.00 },  /* C≡N    */
    { 6,  8,  1,  1.4300,  25.70 },  /* C–O    */
    { 6,  8,  2,  1.2300,  51.40 },  /* C=O    */
    { 6,  9,  1,  1.3500,  23.70 },  /* C–F    */
    { 6, 16,  1,  1.8200,  14.30 },  /* C–S    */
    { 6, 17,  1,  1.7660,  14.30 },  /* C–Cl   */
    { 7,  7,  1,  1.4500,  19.50 },  /* N–N    */
    { 7,  7,  2,  1.2500,  39.00 },  /* N=N    */
    { 7,  7,  3,  1.0980,  78.00 },  /* N≡N    */
    { 7,  8,  1,  1.4400,  20.40 },  /* N–O    */
    { 8,  8,  1,  1.4800,  19.00 },  /* O–O    */
    { 8, 16,  1,  1.5800,  18.00 },  /* O–S    */
    {16, 16,  1,  2.0380,  10.50 },  /* S–S    */
    {15, 15,  1,  2.2100,  10.00 },  /* P–P    */
    { 6, 15,  1,  1.8430,  16.00 },  /* C–P    */
    { 7, 15,  1,  1.6500,  18.00 },  /* N–P    */
    { 8, 15,  1,  1.4810,  26.00 },  /* O–P    */
};
static const int BOND_TABLE_LEN =
    (int)(sizeof(BOND_TABLE) / sizeof(BOND_TABLE[0]));

/* ══════════════════════════════════════════════════════════════════════════
 * Angle parameter database
 * Zb is the central atom.  Za <= Zc for canonical form.
 * k in eV/rad²; theta0 in radians.
 *
 * Degrees to radians: × π/180
 * AMBER k_angle [kcal/mol/rad²] × KCAL_MOL_TO_EV = eV/rad²
 * ══════════════════════════════════════════════════════════════════════════ */
#define DEG2RAD(d) ((d) * 3.14159265358979323846 / 180.0)
static const AngleParam ANGLE_TABLE[] = {
/*  Za  Zb  Zc   theta0(rad)         k(eV/rad²)   source                      */
/*  k = 2 × AMBER_parm [kcal/mol/rad²] × KCAL_MOL_TO_EV                      */
/*  (Code uses V=0.5k(θ-θ0)², AMBER uses V=k(θ-θ0)², so 2× factor needed)   */
    { 1,  8,  1,  DEG2RAD(104.52), 4.770 }, /* H-O-H  TIP3P HW-OW-HW 55 kc  */
    { 1,  7,  1,  DEG2RAD(106.67), 3.738 }, /* H-N-H  H-N3-H ff14SB 43.1 kc  */
    { 1,  6,  1,  DEG2RAD(109.47), 3.035 }, /* H-C-H  HC-CT-HC AMBER 35 kc   */
    { 1,  6,  6,  DEG2RAD(109.47), 4.336 }, /* H-C-C  HC-CT-CT 50 kc         */
    { 1,  6,  7,  DEG2RAD(109.47), 4.336 }, /* H-C-N  HC-CT-N  50 kc         */
    { 1,  6,  8,  DEG2RAD(109.47), 4.336 }, /* H-C-O  HC-CT-OS 50 kc         */
    { 1,  7,  6,  DEG2RAD(118.00), 3.901 }, /* H-N-C  H-N3-CT  45 kc         */
    { 1,  8,  6,  DEG2RAD(108.50), 4.770 }, /* H-O-C  HO-OH-CT 55 kc         */
    { 6,  6,  6,  DEG2RAD(109.47), 3.469 }, /* C-C-C  CT-CT-CT 40 kc         */
    { 6,  6,  7,  DEG2RAD(109.47), 5.464 }, /* C-C-N  CT-CT-N  63 kc         */
    { 6,  6,  8,  DEG2RAD(109.47), 4.336 }, /* C-C-O  CT-CT-OS 50 kc         */
    { 6,  6, 16,  DEG2RAD(114.00), 4.163 }, /* C-C-S  CT-CT-S  48 kc         */
    { 7,  6,  7,  DEG2RAD(116.00), 5.204 }, /* N-C-N  N-C-N    60 kc (est)   */
    { 7,  6,  8,  DEG2RAD(115.00), 5.898 }, /* N-C-O  N-C-O2   68 kc         */
    { 8,  6,  8,  DEG2RAD(123.00), 6.938 }, /* O-C-O  O=C-O2   80 kc         */
    { 6,  7,  6,  DEG2RAD(111.00), 5.638 }, /* C-N-C  CT-N-CT  65 kc         */
    { 6,  8,  6,  DEG2RAD(111.55), 5.204 }, /* C-O-C  CT-OS-CT 60 kc         */
    { 6, 16,  6,  DEG2RAD(102.60), 3.904 }, /* C-S-C  CT-S-CT  45 kc         */
};
static const int ANGLE_TABLE_LEN =
    (int)(sizeof(ANGLE_TABLE) / sizeof(ANGLE_TABLE[0]));

/* ══════════════════════════════════════════════════════════════════════════
 * Bond parameter lookup
 * ══════════════════════════════════════════════════════════════════════════ */
int forces_bond_params(int Za, int Zb, int order, BondParam *out) {
    /* Canonical form: Za <= Zb */
    if (Za > Zb) { int t = Za; Za = Zb; Zb = t; }

    /* Exact match first */
    for (int i = 0; i < BOND_TABLE_LEN; i++) {
        const BondParam *p = &BOND_TABLE[i];
        if (p->Za == Za && p->Zb == Zb && p->order == order) {
            *out = *p;
            return 1;
        }
    }
    /* Fall back to single-bond entry if double/triple not found */
    if (order > 1) {
        for (int i = 0; i < BOND_TABLE_LEN; i++) {
            const BondParam *p = &BOND_TABLE[i];
            if (p->Za == Za && p->Zb == Zb && p->order == 1) {
                *out = *p;
                return 1;
            }
        }
    }

    /* Geometric fallback: use sum of covalent radii, generic k */
    const Element *ea = pt_element(Za);
    const Element *eb = pt_element(Zb);
    if (ea && eb) {
        out->Za = Za; out->Zb = Zb; out->order = order;
        out->r0 = ea->covalent_radius + eb->covalent_radius;
        out->k  = 20.0;  /* generic, eV/Å² */
        return 1;
    }
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Angle parameter lookup
 * ══════════════════════════════════════════════════════════════════════════ */
int forces_angle_params(int Za, int Zb, int Zc, AngleParam *out) {
    /* Canonical form: Za <= Zc */
    if (Za > Zc) { int t = Za; Za = Zc; Zc = t; }

    for (int i = 0; i < ANGLE_TABLE_LEN; i++) {
        const AngleParam *p = &ANGLE_TABLE[i];
        if (p->Zb == Zb &&
            ((p->Za == Za && p->Zc == Zc) ||
             (p->Za == Zc && p->Zc == Za))) {
            *out = *p;
            return 1;
        }
    }

    /* Generic tetrahedral fallback - k matches the corrected table scale
     * (2x AMBER_parm x KCAL_MOL_TO_EV convention, see table comment above).
     * An earlier version left this at the OLD pre-correction value (0.60)
     * after the table itself was fixed, silently giving any untabulated
     * angle type inconsistent, too-soft physics relative to every
     * tabulated entry - using a representative generic AMBER value
     * (40 kcal/mol/rad^2, the CT-CT-CT constant) here instead. */
    out->Za = Za; out->Zb = Zb; out->Zc = Zc;
    out->theta0 = DEG2RAD(109.47);
    out->k      = 2.0 * 40.0 * KCAL_MOL_TO_EV;  /* = 3.469 eV/rad^2 */
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Lennard-Jones + Coulomb pair force
 *
 * Physics:
 *   V_LJ  = 4ε [(σ/r)^12 − (σ/r)^6]
 *   V_C   = COULOMB_MD × q_i × q_j / r
 *
 *   F_i = (dV_LJ/dr + dV_C/dr) × r̂_ij
 *        where r̂_ij = (r_j − r_i) / r  [unit vector FROM i TOWARD j]
 *
 * Derivation of scalar prefactors:
 *   dV_LJ/dr = 24ε/r × [(σ/r)^6 − 2(σ/r)^12]
 *   Force F_i = f_scalar × r_ij,  f_scalar = dV/dr / r
 *
 *   For LJ:  f_lj  = 24ε/r² × [(σ/r)^6 − 2(σ/r)^12]
 *   For C:   f_c   = -COULOMB_MD × qi × qj / r³
 *            (negative because dV_C/dr = −k q q / r²,
 *             and F_i = (dV/dr / r) × r_ij)
 * ══════════════════════════════════════════════════════════════════════════ */
PairEnergy forces_nonbonded_pair(Atom *atoms, int ia, int ib,
                                  const SimBox *box,
                                  int use_lj, int use_coulomb,
                                  double dielectric) {
    PairEnergy result = {0.0, 0.0};
    Atom *ai = &atoms[ia];
    Atom *bi = &atoms[ib];

    Vec3 r_ij = vec3_sub(bi->position, ai->position);

    /* Minimum-image PBC */
    if (box && (box->periodic[0] || box->periodic[1] || box->periodic[2])) {
        r_ij = vec3_pbc(r_ij, box->dimensions);
    }

    double r2 = vec3_norm2(r_ij);
    if (r2 < 1.0e-10) return result;  /* avoid self-interaction at r≈0 */
    double r  = sqrt(r2);

    double f_total = 0.0;  /* scalar force coefficient: F_i = f_total × r_ij */

    /* ── Lennard-Jones ───────────────────────────────────────────────────── */
    if (use_lj) {
        double eps   = lj_eps_combine(ai->lj_epsilon, bi->lj_epsilon);
        double sigma = lj_sigma_combine(ai->lj_sigma, bi->lj_sigma);
        double sr2  = (sigma * sigma) / r2;
        double sr6  = sr2 * sr2 * sr2;
        double sr12 = sr6 * sr6;

        /* V_LJ = 4ε(sr12 - sr6)  [eV] */
        result.lj_energy = 4.0 * eps * (sr12 - sr6);

        /* f_lj = 24ε/r² × (sr6 - 2*sr12).  F_i = f_lj × r_ij (r_ij points
         * FROM i TO j). Sign convention, verified by direct substitution:
         * at small r (r<<σ): sr12 >> sr6, so the bracket is NEGATIVE,
         * making f_lj NEGATIVE → F_i points OPPOSITE to r_ij (away from
         * j) → REPULSION, correctly. At large r near the LJ minimum and
         * beyond (r>2^(1/6)σ): sr6 > 2*sr12, bracket POSITIVE, f_lj
         * POSITIVE → F_i ALONG r_ij (toward j) → ATTRACTION, correctly. */
        double f_lj = (24.0 * eps / r2) * (sr6 - 2.0 * sr12);
        f_total += f_lj;
    }

    /* ── Coulomb ─────────────────────────────────────────────────────────── */
    if (use_coulomb) {
        double qi = ai->partial_charge;
        double qj = bi->partial_charge;
        if (fabs(qi) > 1.0e-9 && fabs(qj) > 1.0e-9) {
            /* V_C = COULOMB_MD × qi × qj / (r × dielectric)  [eV] */
            result.coulomb_energy = COULOMB_MD * qi * qj / (r * dielectric);

            /* f_c = -COULOMB_MD × qi × qj / (r³ × dielectric).
             * F_i = f_c × r_ij. Sign convention, verified by direct
             * substitution: for OPPOSITE charges (qi*qj<0), f_c is
             * POSITIVE → F_i points ALONG r_ij (toward j) → ATTRACTION,
             * correctly. For LIKE charges (qi*qj>0), f_c is NEGATIVE →
             * F_i points OPPOSITE to r_ij (away from j) → REPULSION,
             * correctly. (Note this is the opposite sign pairing from
             * the LJ force above - that asymmetry is expected, since
             * Coulomb's sign comes from the charge product while LJ's
             * comes from the distance-dependent bracket term.) */
            double f_c = -COULOMB_MD * qi * qj / (r2 * r * dielectric);
            f_total += f_c;
        }
    }

    /* Accumulate force onto both atoms (Newton's 3rd law) */
    Vec3 F_i = vec3_scale(r_ij, f_total);
    vec3_iadd(&ai->force, F_i);
    vec3_isub(&bi->force, F_i);

    return result;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Harmonic bond force
 *
 * V_bond = 0.5 k (r − r0)²
 * dV/dr  = k (r − r0)
 *
 * r_ab = pos_b − pos_a  (from a toward b)
 * F_a  = k(r−r0) × r̂_ab   (toward b when stretched, away when compressed)
 * F_b  = −F_a
 * ══════════════════════════════════════════════════════════════════════════ */
double forces_bond(Atom *atoms, const Bond *bond) {
    Atom *a = &atoms[bond->atom_a];
    Atom *b = &atoms[bond->atom_b];

    Vec3   r_ab   = vec3_sub(b->position, a->position);
    double r      = vec3_norm(r_ab);
    if (r < 1.0e-10) return 0.0;

    double stretch = r - bond->r0;
    double energy  = 0.5 * bond->k * stretch * stretch;

    /* f_scalar = k × stretch / r  →  F_a = f_scalar × r_ab */
    double f_scalar = bond->k * stretch / r;
    Vec3 F_a = vec3_scale(r_ab, f_scalar);
    vec3_iadd(&a->force, F_a);
    vec3_isub(&b->force, F_a);

    return energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Harmonic angle force
 *
 * Atoms a–b–c, b is the central atom.
 * V_angle = 0.5 k (θ − θ0)²
 *
 * Gradient via chain rule:
 *   e_ba = (r_a − r_b) / d_ba
 *   e_bc = (r_c − r_b) / d_bc
 *   cos θ = e_ba · e_bc
 *
 *   ∂cos θ / ∂r_a = (e_bc − cos θ × e_ba) / d_ba
 *   ∂cos θ / ∂r_c = (e_ba − cos θ × e_bc) / d_bc
 *   ∂cos θ / ∂r_b = −(∂cos θ/∂r_a + ∂cos θ/∂r_c)
 *
 *   dV/dθ = k (θ − θ0)
 *   dθ/d(cosθ) = −1 / sin θ
 *
 *   F_a = −dV/dθ × dθ/d(cosθ) × ∂cosθ/∂r_a
 *        = [k(θ−θ0)/sin θ] × (e_bc − cosθ e_ba) / d_ba
 * ══════════════════════════════════════════════════════════════════════════ */
double forces_angle(Atom *atoms, const Angle *angle) {
    Atom *a = &atoms[angle->atom_a];
    Atom *b = &atoms[angle->atom_b];
    Atom *c = &atoms[angle->atom_c];

    Vec3 r_ba = vec3_sub(a->position, b->position);
    Vec3 r_bc = vec3_sub(c->position, b->position);

    double d_ba = vec3_norm(r_ba);
    double d_bc = vec3_norm(r_bc);
    if (d_ba < 1.0e-10 || d_bc < 1.0e-10) return 0.0;

    Vec3 e_ba = vec3_scale(r_ba, 1.0 / d_ba);
    Vec3 e_bc = vec3_scale(r_bc, 1.0 / d_bc);

    double cos_theta = vec3_dot(e_ba, e_bc);
    /* Clamp for numerical safety at linear/collapsed angles */
    if (cos_theta >  1.0 - 1.0e-7) cos_theta =  1.0 - 1.0e-7;
    if (cos_theta < -1.0 + 1.0e-7) cos_theta = -1.0 + 1.0e-7;

    double theta   = acos(cos_theta);
    double sin_theta = sin(theta);
    double energy  = 0.5 * angle->k * (theta - angle->theta0)
                                    * (theta - angle->theta0);

    /* Prefactor: k(θ−θ0)/sinθ */
    double pre = angle->k * (theta - angle->theta0) / sin_theta;

    /* ∂cosθ/∂r_a = (e_bc − cosθ e_ba) / d_ba */
    Vec3 dcos_a = vec3_scale(
        vec3_sub(e_bc, vec3_scale(e_ba, cos_theta)),
        1.0 / d_ba);

    /* ∂cosθ/∂r_c = (e_ba − cosθ e_bc) / d_bc */
    Vec3 dcos_c = vec3_scale(
        vec3_sub(e_ba, vec3_scale(e_bc, cos_theta)),
        1.0 / d_bc);

    /* F_x = pre × ∂cosθ/∂r_x  (dθ/dcos already absorbed into pre) */
    Vec3 F_a = vec3_scale(dcos_a, pre);
    Vec3 F_c = vec3_scale(dcos_c, pre);
    /* Newton: F_b = -(F_a + F_c) */
    Vec3 F_b = vec3_negate(vec3_add(F_a, F_c));

    vec3_iadd(&a->force, F_a);
    vec3_iadd(&b->force, F_b);
    vec3_iadd(&c->force, F_c);

    return energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Master force calculation
 * ══════════════════════════════════════════════════════════════════════════ */
void forces_calculate(Simulation *sim) {
    int N = sim->num_atoms;

    /* 1. Zero all forces */
    for (int i = 0; i < N; i++)
        sim->atoms[i].force = vec3_zero();

    double E_lj      = 0.0;
    double E_coulomb = 0.0;
    double E_bond    = 0.0;
    double E_angle   = 0.0;

    /* 2. Non-bonded pairs (O(N²) — replace with cell list for large N) */
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            /* Skip pairs that are bonded (1-2) or angle-related (1-3) */
            int skip = 0;
            for (int b = 0; b < sim->num_bonds; b++) {
                if ((sim->bonds[b].atom_a == i && sim->bonds[b].atom_b == j) ||
                    (sim->bonds[b].atom_a == j && sim->bonds[b].atom_b == i)) {
                    skip = 1; break;
                }
            }
            if (skip) continue;

            /* Check 1-3 exclusion via angles */
            for (int a = 0; a < sim->num_angles && !skip; a++) {
                const Angle *ang = &sim->angles[a];
                if ((ang->atom_a == i && ang->atom_c == j) ||
                    (ang->atom_a == j && ang->atom_c == i))
                    skip = 1;
            }
            if (skip) continue;

            /* Distance cutoff check */
            Vec3 r_ij = vec3_sub(sim->atoms[j].position,
                                 sim->atoms[i].position);
            if (sim->box.periodic[0] || sim->box.periodic[1] || sim->box.periodic[2])
                r_ij = vec3_pbc(r_ij, sim->box.dimensions);

            if (vec3_norm(r_ij) > sim->cutoff) continue;

            PairEnergy pe = forces_nonbonded_pair(
                sim->atoms, i, j,
                &sim->box,
                sim->use_lj,
                sim->use_coulomb,
                sim->dielectric);

            E_lj      += pe.lj_energy;
            E_coulomb += pe.coulomb_energy;
        }
    }

    /* 3. Bonded stretches */
    if (sim->use_bonds) {
        for (int b = 0; b < sim->num_bonds; b++)
            E_bond += forces_bond(sim->atoms, &sim->bonds[b]);
    }

    /* 4. Angle bends */
    if (sim->use_angles) {
        for (int a = 0; a < sim->num_angles; a++)
            E_angle += forces_angle(sim->atoms, &sim->angles[a]);
    }

    sim->potential_energy = E_lj + E_coulomb + E_bond + E_angle;
    sim->E_lj_total      = E_lj;
    sim->E_coulomb_total = E_coulomb;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Print force summary
 * ══════════════════════════════════════════════════════════════════════════ */
void forces_print_summary(const Simulation *sim) {
    printf("  Force summary (step %llu):\n", (unsigned long long)sim->step);
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        printf("    Atom %3d %-2s  |F|=%8.4f eV/Å  "
               "F=(%8.4f, %8.4f, %8.4f)\n",
               i, a->element->symbol,
               vec3_norm(a->force),
               a->force.x, a->force.y, a->force.z);
    }
    printf("  Potential energy: %.6f eV\n", sim->potential_energy);
}

```

---

### FILE: src/integrator.c
Location: `src/integrator.c`
```cpp
#include <math.h>
#include <stdio.h>
#include "../include/integrator.h"
#include "../include/forces.h"
#include "../include/constants.h"

/*
 * integrator.c
 *
 * Velocity Verlet + Berendsen thermostat.
 * All unit conversions are explicit and commented.
 */

#define KB_EV   8.617333262e-5   /* Boltzmann constant in eV/K */

/* ══════════════════════════════════════════════════════════════════════════
 * Half-step A: velocity kick + position drift
 *
 * Physics:
 *   a_i = F_i / m_i                    [eV/Å / AMU]
 *   a_i [Å/fs²] = a_i [eV/Å/AMU] × MD_FORCE_CONV
 *
 *   v_i(t+dt/2) = v_i(t) + 0.5 × a_i × dt
 *   r_i(t+dt)   = r_i(t) + v_i(t+dt/2) × dt
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_kick_drift(Simulation *sim) {
    double dt   = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];

        /* a [Å/fs²] = F [eV/Å] / m [AMU] × MD_FORCE_CONV */
        double inv_m = MD_FORCE_CONV / a->mass;

        /* Half-kick: v += 0.5 a dt */
        a->velocity.x += 0.5 * a->force.x * inv_m * dt;
        a->velocity.y += 0.5 * a->force.y * inv_m * dt;
        a->velocity.z += 0.5 * a->force.z * inv_m * dt;

        /* Drift: r += v dt  (using half-kicked velocity) */
        a->position.x += a->velocity.x * dt;
        a->position.y += a->velocity.y * dt;
        a->position.z += a->velocity.z * dt;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Half-step B: final velocity kick with new forces
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_kick(Simulation *sim) {
    double dt = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        double inv_m = MD_FORCE_CONV / a->mass;

        a->velocity.x += 0.5 * a->force.x * inv_m * dt;
        a->velocity.y += 0.5 * a->force.y * inv_m * dt;
        a->velocity.z += 0.5 * a->force.z * inv_m * dt;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Full velocity Verlet step
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_step(Simulation *sim) {
    /* 1. Kick velocities by half-step, drift positions by full step */
    integrator_kick_drift(sim);

    /* 2. Recalculate forces at new positions */
    forces_calculate(sim);

    /* 3. Complete velocity update with new forces */
    integrator_kick(sim);

    /* 4. Apply thermostat if active */
    if (sim->thermostat.type == THERMOSTAT_BERENDSEN)
        integrator_berendsen(sim);

    /* 5. Update thermodynamics */
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->temperature    = integrator_temperature(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;

    sim->step++;
    sim->time += sim->dt;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Kinetic energy
 *
 * KE [eV] = 0.5 × Σ m_i [AMU] × |v_i|² [Å²/fs²] × AMU_AFS2_TO_EV
 *
 * AMU_AFS2_TO_EV = 1 AMU × (Å/fs)² in eV
 *   = 1.66053906660e-27 kg × (1e5 m/s)²
 *   = 1.66053906660e-27 × 1e10 J
 *   = 1.66053906660e-17 / 1.602176634e-19 eV
 *   = 103.6427 eV
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_kinetic_energy(const Simulation *sim) {
    double ke = 0.0;
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        double v2 = vec3_norm2(a->velocity);
        ke += 0.5 * a->mass * v2;
    }
    return ke * AMU_AFS2_TO_EV;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Instantaneous temperature from equipartition theorem
 *
 * <KE> = (3N - 3)/2 × k_B × T
 *   → T = 2 × KE / ((3N - 3) × k_B)
 *
 * DOF = 3N - 3: the -3 removes the 3 translational DOF of the centre of
 * mass, which integrator_remove_com_velocity already constrains to zero.
 * Using 3N instead (as was done in an earlier version) systematically
 * underestimates T, with the error most severe for small systems:
 * 33% too low for N=3 (water), 8% too low for N=12 (uracil), etc.
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_temperature(const Simulation *sim) {
    if (sim->num_atoms < 2) return 0.0;
    double ke  = integrator_kinetic_energy(sim);
    int    dof = 3 * sim->num_atoms - 3;      /* correct: subtract COM DOF */
    return (2.0 * ke) / ((double)dof * KB_EV);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Berendsen thermostat
 *
 * Rescales velocities: v_new = λ v
 *   λ = sqrt(1 + (dt/τ)(T_0/T − 1))
 *
 * arg clamp: if 1 + (dt/τ)(T_0/T - 1) < 0 (extreme overshoot), clamp
 * to 0.01 rather than 0.0 — setting λ=0 would stop ALL atomic motion,
 * making the next step's T also 0, then the thermostat stays locked
 * at λ=0 indefinitely. The minimum clamp λ = 0.1 (arg = 0.01) still
 * cools strongly while keeping motion alive so the thermostat recovers.
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_berendsen(Simulation *sim) {
    double T     = integrator_temperature(sim);
    double T0    = sim->thermostat.target_temperature;
    double tau   = sim->thermostat.tau;
    double dt    = sim->dt;

    if (T < 1.0e-10) return;

    double ratio = T0 / T;
    double arg   = 1.0 + (dt / tau) * (ratio - 1.0);
    if (arg < 0.01) arg = 0.01;          /* clamp: avoids motion lockout */
    double lambda = sqrt(arg);

    for (int i = 0; i < sim->num_atoms; i++)
        vec3_iscale(&sim->atoms[i].velocity, lambda);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Maxwell-Boltzmann velocity initialisation
 *
 * Box-Muller transform: given U1, U2 uniform in (0,1),
 *   Z0 = sqrt(-2 ln U1) cos(2π U2)  ~  N(0,1)
 *   Z1 = sqrt(-2 ln U1) sin(2π U2)  ~  N(0,1)
 *
 * Thermal speed for component x:
 *   σ_v = sqrt(k_B T / m)  in Å/fs
 *   v_x = Z × σ_v
 *
 * σ_v [Å/fs]: k_B T [eV] / m [AMU] × (1/AMU_AFS2_TO_EV)
 * ══════════════════════════════════════════════════════════════════════════ */
static unsigned long rng_state;

static double rand_uniform(void) {
    /* LCG with Knuth constants — good enough for velocity initialisation */
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / (double)(1ULL << 31);
}

static double rand_normal(void) {
    double u1, u2;
    do { u1 = rand_uniform(); } while (u1 < 1.0e-10);
    u2 = rand_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
                                   unsigned long seed) {
    rng_state = (seed == 0) ? 12345678901234567ULL : (unsigned long long)seed;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        /* σ_v [Å/fs] = sqrt(k_B T [eV] / (m [AMU] × AMU_AFS2_TO_EV)) */
        double sigma_v = sqrt(KB_EV * T_init / (a->mass * AMU_AFS2_TO_EV));

        a->velocity.x = rand_normal() * sigma_v;
        a->velocity.y = rand_normal() * sigma_v;
        a->velocity.z = rand_normal() * sigma_v;
    }

    /* Remove centre-of-mass drift */
    integrator_remove_com_velocity(sim);

    /* Rescale to exact target temperature */
    double T_actual = integrator_temperature(sim);
    if (T_actual > 1.0e-10) {
        double scale = sqrt(T_init / T_actual);
        for (int i = 0; i < sim->num_atoms; i++)
            vec3_iscale(&sim->atoms[i].velocity, scale);
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Remove centre-of-mass velocity
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_remove_com_velocity(Simulation *sim) {
    double total_mass = 0.0;
    Vec3   p_com      = vec3_zero();

    for (int i = 0; i < sim->num_atoms; i++) {
        double m = sim->atoms[i].mass;
        total_mass += m;
        vec3_iadd(&p_com, vec3_scale(sim->atoms[i].velocity, m));
    }
    if (total_mass < 1.0e-10) return;

    Vec3 v_com = vec3_scale(p_com, 1.0 / total_mass);
    for (int i = 0; i < sim->num_atoms; i++)
        vec3_isub(&sim->atoms[i].velocity, v_com);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Print step summary
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_print_step(const Simulation *sim) {
    printf("  Step %6llu  t=%8.2f fs  "
           "KE=%9.5f  PE=%9.5f  E=%9.5f eV  T=%7.2f K\n",
           (unsigned long long)sim->step,
           sim->time,
           sim->kinetic_energy,
           sim->potential_energy,
           sim->total_energy,
           sim->temperature);
}

```

---

### FILE: src/main.c
Location: `src/main.c`
```cpp
#include <stdio.h>
#include <math.h>
#include "../include/constants.h"
#include "../include/periodic_table.h"
#include "../include/quantum.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/sim.h"
#include "../include/nucleobases.h"

/*
 * main.c
 *
 * Demonstration sequence, bottom-up:
 *   1. Quantum layer  — print orbital energies and wave functions for H, C, O
 *   2. Bonding layer  — show H2 bond dissociation energy curve from LJ
 *   3. Molecular MD   — run 300K dynamics on a water molecule (H2O)
 *   4. Small ensemble — three H2O molecules, short NVT trajectory
 *
 * Every number printed has physical units labelled.
 * This is the foundation: add more chemistry above this bedrock.
 */

/* ── Pretty separator ────────────────────────────────────────────────────── */
static void banner(const char *title) {
    printf("\n╔══════════════════════════════════════════════════════╗\n");
    printf(  "║  %-52s║\n", title);
    printf(  "╚══════════════════════════════════════════════════════╝\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 1: Quantum mechanical orbital structure
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_quantum(void) {
    banner("DEMO 1: Quantum orbital structure");

    int elements[] = {1, 6, 7, 8};   /* H, C, N, O */
    int n_el = 4;

    for (int e = 0; e < n_el; e++) {
        int Z = elements[e];
        pt_print_element(Z);

        /* Build a temporary atom just to get orbital data */
        Atom atom = {0};
        atom.Z       = Z;
        atom.element = pt_element(Z);
        atom.mass    = pt_element(Z)->mass;
        pt_electron_config(Z, &atom.electron_config);
        quantum_fill_orbitals(&atom);
        quantum_print_orbitals(&atom);

        /* Print the radial wave function profile for the valence orbital */
        /* Find the highest-energy occupied orbital */
        int hi = 0;
        for (int i = 1; i < atom.num_orbitals; i++) {
            if (atom.orbitals[i].orbital_energy >
                atom.orbitals[hi].orbital_energy)
                hi = i;
        }
        int n = atom.orbitals[hi].qn.n;
        int l = atom.orbitals[hi].qn.l;
        ElectronConfig cfg;
        pt_electron_config(Z, &cfg);
        double Zeff = quantum_zeff(Z, n, l, &cfg);
        double r_mp = quantum_most_probable_radius(n, l, Zeff);

        printf("  Valence orbital: %d%c  Z_eff=%.3f  r_mp=%.3f Å\n",
               n, "spdf"[l], Zeff, r_mp);

        /* ASCII radial probability profile (0..5 Å in 40 steps) */
        printf("  Radial probability P(r) = r²|R_nl(r)|²:\n");
        double r_max_plot = 10.0;
        double step       = r_max_plot / 40.0;
        double P_max      = 0.0;
        for (int k = 1; k <= 40; k++) {
            double P = quantum_radial_probability(n, l, Zeff, k * step);
            if (P > P_max) P_max = P;
        }
        printf("  0 Å ");
        for (int k = 1; k <= 40; k++) {
            double r = k * step;
            double P = quantum_radial_probability(n, l, Zeff, r) / P_max;
            printf("%c", P < 0.1 ? ' ' :
                         P < 0.3 ? '.' :
                         P < 0.6 ? ':' :
                         P < 0.85? '|' : '#');
        }
        printf(" %.1f Å\n\n", r_max_plot);
    }
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 2: Two H-H potentials compared — covalent bond vs. van der Waals
 *
 * This demo exists to make an important distinction explicit: the LJ
 * potential and the covalent bond potential are TWO DIFFERENT physical
 * interactions, computed by two different code paths, and they answer two
 * different questions:
 *
 *   COVALENT (harmonic, from BOND_TABLE):
 *     "Two H atoms ARE bonded — what does stretching/compressing that
 *      shared-electron bond cost?" Minimum at r0 = 0.7414 Å (the actual
 *      H2 bond length). This is what governs sim_place_h2() and is the
 *      ONLY potential applied to atoms in sim->bonds[].
 *
 *   VAN DER WAALS (Lennard-Jones, UFF parameters):
 *     "Two H atoms are NOT bonded — how do their electron clouds interact
 *      at a distance?" Minimum at r ≈ 3.24 Å, far weaker (~0.002 eV vs.
 *      tens of eV for the covalent well). This is what governs how
 *      separate, non-bonded atoms or molecules approach each other —
 *      it's the physics behind gas pressure, condensation, and packing.
 *
 * forces_calculate() automatically excludes bonded pairs from the LJ/Coulomb
 * sum (see the 1-2 and 1-3 exclusion logic), so a real bonded H2 molecule
 * NEVER feels the LJ curve between its own two atoms — only the harmonic
 * term. The two curves below are printed side by side specifically so this
 * distinction is never ambiguous.
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_bond_curve(void) {
    banner("DEMO 2: H-H covalent bond vs. van der Waals (two different physics)");

    /* ── Van der Waals (non-bonded) curve ────────────────────────────────── */
    const Element *H = pt_element(1);
    double eps   = lj_eps_combine(H->lj_epsilon, H->lj_epsilon);
    double sigma = lj_sigma_combine(H->lj_sigma, H->lj_sigma);
    double lj_rmin = pow(2.0, 1.0/6.0) * sigma;

    /* ── Covalent (bonded) curve ─────────────────────────────────────────── */
    BondParam bp;
    forces_bond_params(1, 1, 1, &bp);  /* H-H single bond from BOND_TABLE */

    printf("  Van der Waals (LJ):  ε=%.5f eV   σ=%.4f Å   r_min=%.4f Å\n",
           eps, sigma, lj_rmin);
    printf("  Covalent (harmonic): r0=%.4f Å   k=%.2f eV/Å²   "
           "(real H2 bond length is 0.7414 Å)\n\n", bp.r0, bp.k);

    printf("  %-8s  %-14s  %-14s  %-s\n",
           "r (Å)", "V_covalent(eV)", "V_vdW (eV)", "Scale");
    printf("  %-8s  %-14s  %-14s\n","────────","──────────────","──────────────");

    /* Print covalent curve only near its own well (it's enormously stronger
     * and would make the vdW curve invisible on a shared linear scale) */
    for (int k = 4; k <= 24; k++) {
        double r = 0.4 + k * 0.05;
        double stretch = r - bp.r0;
        double V_cov   = 0.5 * bp.k * stretch * stretch;
        printf("  %-8.4f  %-14.6f  %-14s  covalent well (depth scale: eV)\n",
               r, V_cov, "—");
    }

    printf("\n");

    /* Print vdW curve over its own, much wider and shallower range */
    double V_min = 0.0;
    for (int k = 6; k <= 80; k++) {
        double r   = k * 0.0625;
        double sr6 = pow(sigma / r, 6.0);
        double V   = 4.0 * eps * (sr6*sr6 - sr6);
        if (V < V_min) V_min = V;
    }
    for (int k = 20; k <= 80; k += 2) {
        double r    = k * 0.0625;
        double sr6  = pow(sigma / r, 6.0);
        double V    = 4.0 * eps * (sr6*sr6 - sr6);
        int bar_len = (int)fmin(40.0, fmax(0.0, 20.0 + 20.0 * V / fabs(V_min)));
        printf("  %-8.4f  %-14s  %-14.6f  %.*s\n",
               r, "—", V, bar_len, "########################################");
    }
    printf("\n  For scale: the real H2 covalent bond dissociation energy is "
           "4.52 eV\n  (a standard spectroscopic constant), versus this "
           "vdW well depth of only\n  %.5f eV — roughly %.0fx weaker. "
           "That gap is why breaking a chemical\n  bond (a reaction) costs "
           "so much more than separating two molecules that\n  are merely "
           "touching (melting/evaporation).\n", fabs(V_min), 4.52/fabs(V_min));
    printf("\n  Caveat: the harmonic term above is only valid for small "
           "vibrations near\n  r0. It's a parabola, not a real bond — it "
           "never flattens out, so it would\n  (wrongly) predict infinite "
           "energy to fully separate the atoms. Capturing\n  actual bond "
           "breaking needs a Morse potential or a reactive force field —\n"
           "  a natural next addition to this codebase.\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 3: Water molecule NVT MD at 300 K
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_water_md(void) {
    banner("DEMO 3: H2O molecule — NVT MD at 300 K");

    Simulation *sim = sim_create(16, 32);
    if (!sim) { printf("  ERROR: allocation failed\n"); return; }

    /* Place water at origin */
    sim_place_h2o(sim, vec3_zero());

    /* Timestep and thermostat */
    sim->dt = 0.5;  /* 0.5 fs — small for stiff O-H bonds */
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 300.0;
    sim->thermostat.tau                = 100.0;  /* fs */

    /* Initial forces and energy */
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->potential_energy = sim->potential_energy;
    sim->total_energy = sim->kinetic_energy + sim->potential_energy;

    printf("  Initial geometry:\n");
    sim_print_atoms(sim);
    printf("\n");
    sim_print_bonds(sim);
    printf("\n");
    sim_print_angles(sim);

    /* Assign Maxwell-Boltzmann velocities at 300 K */
    integrator_maxwell_boltzmann(sim, 300.0, 42UL);

    /* Recalculate after velocity assignment */
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
    sim->temperature    = integrator_temperature(sim);

    printf("\n  Initial thermodynamics:\n");
    printf("  KE = %.6f eV  PE = %.6f eV  E = %.6f eV  T = %.2f K\n\n",
           sim->kinetic_energy, sim->potential_energy,
           sim->total_energy, sim->temperature);

    /* ── MD trajectory ─────────────────────────────────────────────────── */
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-12s\n",
           "Step","t (fs)","KE (eV)","PE (eV)","T (K)","O-H1 dist (Å)");
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-12s\n",
           "────","──────","───────","───────","─────","────────────");

    int N_steps = 2000;
    int print_every = 100;

    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);

        if (step % print_every == 0) {
            double dOH = vec3_dist(sim->atoms[0].position,
                                    sim->atoms[1].position);
            printf("  %-8d  %-10.3f  %-10.6f  %-10.6f  %-8.2f  %-12.6f\n",
                   (int)sim->step,
                   sim->time,
                   sim->kinetic_energy,
                   sim->potential_energy,
                   sim->temperature,
                   dOH);
        }
    }

    printf("\n  Final geometry after %d steps:\n", N_steps);
    sim_print_atoms(sim);
    sim_print_summary(sim);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 4: Cyclic (H2O)3 hydrogen-bonded trimer — emergent non-covalent structure
 *
 * This is the most chemically important demo in the file: it shows molecular
 * structure that NO ONE specified directly. We only give the simulator:
 *   - atomic positions / charges (from the periodic table + TIP3P charges)
 *   - the harmonic bond/angle terms (covalent, intramolecular)
 *   - the Coulomb + Lennard-Jones terms (non-bonded, intermolecular)
 *
 * The hydrogen-bonded ring — three O···H–O bridges holding the trimer
 * together — is not programmed in. It falls out of minimising those terms.
 * This is the gas-phase cyclic water trimer, a real, well-characterised
 * structure (see e.g. Keutsch & Saykally, PNAS 2001).
 *
 * Geometry construction:
 *   O atoms placed at the vertices of an equilateral triangle, O···O = 2.95 Å
 *   (the experimental/ab-initio range for a water-water H-bond is ~2.8-3.0 Å).
 *   Each water's "donor" H points along the O→O(next) vector at the correct
 *   O-H bond length (0.9572 Å) — this is the bridging hydrogen.
 *   Each water's "free" H satisfies the 104.52° H-O-H angle, rotated to
 *   point outward from the ring (away from the neighbouring molecules,
 *   avoiding artificial steric clash).
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_water_cluster(void) {
    banner("DEMO 4: Cyclic (H2O)3 — emergent hydrogen-bonded ring");

    Simulation *sim = sim_create(64, 128);
    if (!sim) { printf("  ERROR: allocation failed\n"); return; }

    const double PI    = 3.14159265358979323846;
    const double r_OO   = 2.95;                  /* Å, O···O H-bond distance */
    const double r_OH   = 0.9572;                 /* Å, covalent O-H         */
    const double hoh    = 104.52 * PI / 180.0;     /* H-O-H angle             */
    const double R_ring = r_OO / sqrt(3.0);        /* circumradius of triangle*/

    Vec3 O_pos[3], donorH_dir[3];

    /* Step 1: place the three oxygens at the triangle vertices */
    for (int i = 0; i < 3; i++) {
        double phi = 2.0 * PI * i / 3.0;
        O_pos[i] = vec3(R_ring * cos(phi), R_ring * sin(phi), 0.0);
    }

    /* Step 2: donor H direction = toward the NEXT oxygen in the ring */
    for (int i = 0; i < 3; i++) {
        int next = (i + 1) % 3;
        donorH_dir[i] = vec3_normalize(vec3_sub(O_pos[next], O_pos[i]));
    }

    int O_idx[3], Hd_idx[3], Hf_idx[3];  /* O, donor-H, free-H atom indices */

    for (int i = 0; i < 3; i++) {
        /* Oxygen */
        O_idx[i] = sim_add_atom(sim, 8, O_pos[i], -0.834);

        /* Donor hydrogen: along donorH_dir[i], at the covalent O-H length.
         * This is the bridging H that points at the next ring oxygen. */
        Vec3 Hd_pos = vec3_add(O_pos[i], vec3_scale(donorH_dir[i], r_OH));
        Hd_idx[i] = sim_add_atom(sim, 1, Hd_pos, +0.417);

        /* Free hydrogen: rotate donorH_dir[i] by ±104.52° about z so the
         * H-O-H angle is exact, choosing the sign that points away from
         * the ring centre (outward), minimising steric clash with
         * neighbouring molecules. Pure z-rotation keeps everything
         * in the ring plane (z=0), which is a fine simplification for
         * a point-charge, non-polarisable model like this one. */
        double cx = donorH_dir[i].x, cy = donorH_dir[i].y;
        double cos_a = cos(hoh), sin_a = sin(hoh);
        Vec3 rot_plus  = vec3(cx*cos_a - cy*sin_a, cx*sin_a + cy*cos_a, 0.0);
        Vec3 rot_minus = vec3(cx*cos_a + cy*sin_a, -cx*sin_a + cy*cos_a, 0.0);

        Vec3 cand_plus  = vec3_add(O_pos[i], vec3_scale(rot_plus,  r_OH));
        Vec3 cand_minus = vec3_add(O_pos[i], vec3_scale(rot_minus, r_OH));

        /* Outward = farther from ring centroid (origin) */
        Vec3 Hf_pos = (vec3_norm(cand_plus) > vec3_norm(cand_minus))
                      ? cand_plus : cand_minus;

        Hf_idx[i] = sim_add_atom(sim, 1, Hf_pos, +0.417);

        /* Apply verified TIP3P LJ parameters (Jorgensen 1983) - see the
         * detailed comment in sim_place_h2o() for why generic UFF values
         * are wrong here. Built by hand since this trimer doesn't go
         * through sim_place_h2o(). */
        sim_set_atom_lj(sim, O_idx[i],  0.1521 * KCAL_MOL_TO_EV, 3.15061);
        sim_set_atom_lj(sim, Hd_idx[i], 0.0, 0.0);
        sim_set_atom_lj(sim, Hf_idx[i], 0.0, 0.0);

        /* Covalent bonds: O to both its own hydrogens */
        sim_add_bond(sim, O_idx[i], Hd_idx[i], 1);
        sim_add_bond(sim, O_idx[i], Hf_idx[i], 1);
    }
    sim_build_angles(sim);

    /*
     * Run at 50 K rather than 300 K. This is an honest choice, not a fudge:
     * a 3-molecule gas-phase cluster with no confining box and no
     * surrounding liquid pressure is a genuinely fragile system — real
     * water trimers in molecular beams are weakly bound (~0.2-0.3 eV per
     * H-bond, but with large amplitude floppy motion and low barriers to
     * rearrangement). At 300 K the available thermal kinetic energy per
     * mode (~0.026 eV) is large enough that an unconfined trimer dissociates
     * on a picosecond timescale — which the previous version of this demo
     * correctly showed when the molecules weren't even H-bond oriented.
     * At 50 K we can watch genuine bound, oscillatory H-bond dynamics
     * within a short trajectory without needing a confining potential.
     */
    sim->dt = 0.5;
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 50.0;
    sim->thermostat.tau                = 50.0;

    integrator_maxwell_boltzmann(sim, 50.0, 99UL);
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;
    sim->temperature    = integrator_temperature(sim);

    printf("  %d atoms, %d bonds, %d angles\n",
           sim->num_atoms, sim->num_bonds, sim->num_angles);
    printf("  Ring O-O-O construction: O...O = %.3f Å per edge\n", r_OO);
    printf("  Initial T=%.2f K  PE=%.6f eV (intermolecular H-bonds "
           "contribute the negative part)\n\n",
           sim->temperature, sim->potential_energy);

    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s %-10s %-10s\n",
           "Step","t (fs)","KE (eV)","PE (eV)","T (K)",
           "O0-O1(Å)","O1-O2(Å)","O2-O0(Å)");
    printf("  %-8s  %-10s  %-10s  %-10s  %-8s  %-10s %-10s %-10s\n",
           "────","──────","───────","───────","─────",
           "────────","────────","────────");

    int N_steps = 4000;
    int print_every = 200;
    double min_OO = 1.0e9, max_OO = 0.0;

    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);

        double dOO[3];
        dOO[0] = vec3_dist(sim->atoms[O_idx[0]].position,
                            sim->atoms[O_idx[1]].position);
        dOO[1] = vec3_dist(sim->atoms[O_idx[1]].position,
                            sim->atoms[O_idx[2]].position);
        dOO[2] = vec3_dist(sim->atoms[O_idx[2]].position,
                            sim->atoms[O_idx[0]].position);

        for (int k = 0; k < 3; k++) {
            if (dOO[k] < min_OO) min_OO = dOO[k];
            if (dOO[k] > max_OO) max_OO = dOO[k];
        }

        if (step % print_every == 0) {
            printf("  %-8d  %-10.3f  %-10.6f  %-10.6f  %-8.2f  "
                   "%-10.4f %-10.4f %-10.4f\n",
                   (int)sim->step, sim->time,
                   sim->kinetic_energy, sim->potential_energy,
                   sim->temperature, dOO[0], dOO[1], dOO[2]);
        }
    }

    printf("\n  Over %d steps (%.0f fs): O···O range = [%.3f, %.3f] Å\n",
           N_steps, sim->time, min_OO, max_OO);
    if (max_OO < 5.0) {
        printf("  → Ring stayed bound. Hydrogen bonds emerged from Coulomb "
               "+ LJ alone — never told the simulator these were H-bonds.\n");
    } else {
        printf("  → Ring dissociated within this trajectory.\n");
    }

    printf("\n  Final state:\n");
    sim_print_atoms(sim);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 5: Methane molecule — tetrahedral geometry
 * ════════════════════════════════════════════════════════════════════════════ */
static void demo_methane(void) {
    banner("DEMO 5: CH4 — tetrahedral geometry check");

    Simulation *sim = sim_create(16, 32);
    if (!sim) { return; }

    sim_place_ch4(sim, vec3_zero());
    forces_calculate(sim);
    sim->kinetic_energy = 0.0;
    sim->total_energy   = sim->potential_energy;

    printf("  Methane geometry (C at origin):\n");
    sim_print_atoms(sim);
    printf("\n");

    /* Print all H-C-H angles */
    printf("  Bond angles:\n");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        Vec3 r_ba = vec3_sub(sim->atoms[a->atom_a].position,
                              sim->atoms[a->atom_b].position);
        Vec3 r_bc = vec3_sub(sim->atoms[a->atom_c].position,
                              sim->atoms[a->atom_b].position);
        double cos_t = vec3_dot(vec3_normalize(r_ba),
                                 vec3_normalize(r_bc));
        if (cos_t >  1.0) cos_t =  1.0;
        if (cos_t < -1.0) cos_t = -1.0;
        double theta = acos(cos_t) * 180.0 / 3.14159265358979323846;
        printf("    H(%d)-C(%d)-H(%d): %.4f°  (ideal 109.47°)\n",
               a->atom_a, a->atom_b, a->atom_c, theta);
    }

    printf("\n  Initial potential energy: %.6f eV\n", sim->potential_energy);
    sim_destroy(sim);
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 6: The five DNA/RNA nucleobases - geometry validation
 *
 * For each base: place atoms from verified PDB CCD coordinates, then
 * INDEPENDENTLY recompute every bond length and ring angle directly
 * from the Cartesian positions (not from the stored r0/theta0, even
 * though those were set to match by construction) - this is the same
 * "trust but verify" approach used in the methane tetrahedral-angle
 * check. A planarity check (independent of any external reference
 * numbers) confirms the aromatic ring is genuinely flat.
 * ════════════════════════════════════════════════════════════════════════════ */
static void print_bond_report(const Simulation *sim, const char *symbols[]) {
    printf("  %-10s %-6s %-6s %-8s\n", "Bond", "Atom1", "Atom2", "r (A)");
    printf("  %-10s %-6s %-6s %-8s\n", "----", "-----", "-----", "-----");
    for (int i = 0; i < sim->num_bonds; i++) {
        const Bond *b = &sim->bonds[i];
        double r = vec3_dist(sim->atoms[b->atom_a].position,
                              sim->atoms[b->atom_b].position);
        printf("  %-10s %-6s %-6s %-8.4f  (order %d)\n",
               b->order == 2 ? "double" : "single",
               symbols[b->atom_a], symbols[b->atom_b], r, b->order);
    }
}

static void print_angle_report(const Simulation *sim, const char *symbols[]) {
    printf("  %-6s %-6s %-6s %-10s\n", "A", "B(ctr)", "C", "theta (deg)");
    printf("  %-6s %-6s %-6s %-10s\n", "-", "------", "-", "-----------");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        Vec3 ba = vec3_sub(sim->atoms[a->atom_a].position,
                            sim->atoms[a->atom_b].position);
        Vec3 bc = vec3_sub(sim->atoms[a->atom_c].position,
                            sim->atoms[a->atom_b].position);
        double theta = vec3_angle(ba, bc) * 180.0 / 3.14159265358979323846;
        printf("  %-6s %-6s %-6s %-10.2f\n",
               symbols[a->atom_a], symbols[a->atom_b], symbols[a->atom_c], theta);
    }
}

static double total_charge(const Simulation *sim) {
    double q = 0.0;
    for (int i = 0; i < sim->num_atoms; i++) q += sim->atoms[i].partial_charge;
    return q;
}

static void demo_nucleobases(void) {
    banner("DEMO 6: The five nucleobases - geometry validation");

    /* ── Uracil ─────────────────────────────────────────────────────────── */
    {
        Simulation *sim = sim_create(16, 16);
        sim_place_uracil(sim, vec3_zero());
        printf("  URACIL (C4H4N2O2) - 12 atoms\n");
        const char *sym[] = {"N1","C2","O2","N3","C4","O4",
                              "C5","C6","HN1","HN3","H5","H6"};
        print_bond_report(sim, sym);
        printf("\n");
        print_angle_report(sim, sym);

        int ring[] = {0,1,3,4,6,7}; /* N1 C2 N3 C4 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("\n  Ring planarity: max deviation = %.4f A "
               "(aromatic rings should be ~0)\n", dev);

        printf("\n  Cross-check vs. electron diffraction "
               "(Ferenczy et al. 1986):\n");
        printf("  %-22s %-10s %-10s\n", "Quantity", "This model", "Literature");
        double r_C2N1 = vec3_dist(sim->atoms[1].position, sim->atoms[0].position);
        double r_C4C5 = vec3_dist(sim->atoms[4].position, sim->atoms[6].position);
        double r_C5C6 = vec3_dist(sim->atoms[6].position, sim->atoms[7].position);
        printf("  %-22s %-10.3f %-10s\n", "C-N (A)", r_C2N1, "1.399");
        printf("  %-22s %-10.3f %-10s\n", "C4-C5 single (A)", r_C4C5, "1.462");
        printf("  %-22s %-10.3f %-10s\n", "C5=C6 double (A)", r_C5C6, "1.343");
        printf("\n  Total molecular charge: %+.6f e (RESP, Aduri et al. "
               "2007 + derived HN1)\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Cytosine ───────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_cytosine(sim, vec3_zero());
        printf("  CYTOSINE (C4H5N3O) - 13 atoms\n");
        printf("  (tautomer-corrected: H relocated to N1, the Watson-Crick-\n");
        printf("   relevant position; N3 left bare as required for pairing)\n");

        const char *sym[] = {"N3","C4","N1","C2","O2","N4",
                              "C5","C6","HN41","HN42","H5","H6","HN1"};
        print_bond_report(sim, sym);

        /* Explicit WC-readiness check: N1 must have 3 bonds (C2,C6,H),
         * N3 must have exactly 2 (C4,C2 - bare, ready to accept guanine's
         * N1-H) */
        int n1_bonds = sim->atoms[2].num_bonds;
        int n3_bonds = sim->atoms[0].num_bonds;
        printf("\n  N1 bond count: %d (expect 3: C2, C6, H - donor ready)\n",
               n1_bonds);
        printf("  N3 bond count: %d (expect 2: C4, C2 - bare, acceptor ready)\n",
               n3_bonds);

        int ring[] = {0,1,2,3,6,7}; /* N3 C4 N1 C2 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("  Ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Thymine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_thymine(sim, vec3_zero());
        printf("  THYMINE (C5H6N2O2, = 5-methyluracil) - 15 atoms\n");
        int ring[] = {0,1,3,4,6,7}; /* N1 C2 N3 C4 C5 C6 */
        double dev = nb_planarity_deviation(sim, ring, 6);
        printf("  Ring planarity: max deviation = %.4f A\n", dev);

        double r_methyl = vec3_dist(sim->atoms[6].position, sim->atoms[11].position);
        printf("  C5-CH3 methyl bond length: %.4f A "
               "(target 1.51 A, toluene-type)\n", r_methyl);

        /* Confirm methyl H's are tetrahedral about the C5-CM axis */
        Vec3 C5 = sim->atoms[6].position, CM = sim->atoms[11].position;
        Vec3 H1 = sim->atoms[12].position, H2 = sim->atoms[13].position;
        Vec3 u1 = vec3_normalize(vec3_sub(H1, CM));
        Vec3 u2 = vec3_normalize(vec3_sub(H2, CM));
        double hch = acos(vec3_dot(u1,u2)) * 180.0/3.14159265358979323846;
        printf("  H-CM-H methyl angle: %.2f deg (ideal tetrahedral 109.47)\n", hch);
        printf("  Total molecular charge: %+.6f e "
               "(lower-confidence approx, see code comment)\n", total_charge(sim));
        (void)C5;
        sim_destroy(sim);
    }

    /* ── Adenine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_adenine(sim, vec3_zero());
        printf("  ADENINE (C5H5N5, fused 5+6 purine ring) - 15 atoms\n");
        int ring[] = {0,1,2,3,4,5,6,7,8,9}; /* all 10 ring atoms */
        double dev = nb_planarity_deviation(sim, ring, 10);
        printf("  Full bicyclic ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    /* ── Guanine ────────────────────────────────────────────────────────── */
    {
        printf("\n");
        Simulation *sim = sim_create(16, 16);
        sim_place_guanine(sim, vec3_zero());
        printf("  GUANINE (C5H5N5O, fused 5+6 purine ring) - 16 atoms\n");
        int ring[] = {0,1,2,3,4,6,7,9,10}; /* 9 ring atoms (excl. O6 substituent) */
        double dev = nb_planarity_deviation(sim, ring, 9);
        printf("  Full bicyclic ring planarity: max deviation = %.4f A\n", dev);
        printf("  Total molecular charge: %+.6f e\n", total_charge(sim));
        sim_destroy(sim);
    }

    printf("\n  All five bases hold real, verified ring geometry. Next: "
           "sugar-phosphate\n  backbones and Watson-Crick base pairing - "
           "G-C should bind via 3 H-bonds,\n  A-T via 2, using nothing but "
           "the Coulomb+LJ code already validated\n  on the water trimer.\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * DEMO 7: Watson-Crick base pairing - G-C vs A-U binding energy
 *
 * The central test: G-C pairs via 3 hydrogen bonds, A-U via 2. If this
 * force field (real charges, real geometry, the same Coulomb+LJ code
 * already validated on the water trimer) is doing real chemistry and
 * not just fitting a foregone conclusion, G-C should come out MORE
 * stable (more negative interaction energy) than A-U.
 *
 * Construction strategy, identical in spirit to the water trimer: get
 * a rough, approximately-correct geometry (align ring planes, place the
 * primary donor/acceptor heavy atoms at a realistic H-bond distance),
 * then let MD relaxation find the true local minimum. No part of the
 * Watson-Crick hydrogen-bond pattern is hard-coded into the force
 * field - it has to emerge from Coulomb+LJ alone, the same way the
 * water trimer's ring emerged.
 * ════════════════════════════════════════════════════════════════════════════ */

typedef struct {
    double min_interE, final_interE, initial_interE;
    double closest_approach_E, closest_approach_dist;
} PairResult;

static PairResult run_pair_relaxation(Simulation *sim,
                                       int primary_a, int primary_b) {
    PairResult r;
    forces_calculate(sim);
    r.initial_interE = sim->potential_energy;

    sim->dt = 0.1;
    sim->thermostat.type               = THERMOSTAT_BERENDSEN;
    sim->thermostat.target_temperature = 50.0;
    sim->thermostat.tau                = 50.0;
    integrator_maxwell_boltzmann(sim, 50.0, 7UL);
    forces_calculate(sim);
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;

    printf("  %-6s %-10s %-10s %-10s %-10s\n",
           "Step", "t(fs)", "PE(eV)", "T(K)", "primary(A)");
    double min_pe = sim->potential_energy;
    double closest_dist = vec3_dist(sim->atoms[primary_a].position,
                                     sim->atoms[primary_b].position);
    double closest_E = sim->potential_energy;

    int N_steps = 800, print_every = 100;
    for (int step = 0; step < N_steps; step++) {
        integrator_step(sim);
        if (sim->potential_energy < min_pe) min_pe = sim->potential_energy;

        double d = vec3_dist(sim->atoms[primary_a].position,
                              sim->atoms[primary_b].position);
        if (d < closest_dist) { closest_dist = d; closest_E = sim->potential_energy; }

        if (step < 5 || step % print_every == 0) {
            printf("  %-6d %-10.3f %-10.6f %-10.2f %-10.4f\n",
                   (int)sim->step, sim->time, sim->potential_energy,
                   sim->temperature, d);
        }
    }
    r.min_interE   = min_pe;
    r.final_interE = sim->potential_energy;
    r.closest_approach_E    = closest_E;
    r.closest_approach_dist = closest_dist;
    return r;
}

static void demo_basepairing(void) {
    banner("DEMO 7: Watson-Crick pairing - does G-C beat A-U?");

    double G_C_energy, A_U_energy;

    /* ── G-C pair ───────────────────────────────────────────────────────── */
    {
        printf("\n  --- Guanine-Cytosine (3 H-bonds: N1-H..N3, N2-H..O2, "
               "O6..H-N4) ---\n");
        Simulation *sim = sim_create(64, 64);
        /* RESP charges from Aduri et al. were fit for use in a condensed-
         * phase force field (with solvent/counterions implicitly assumed
         * to provide screening) - using them in a bare vacuum Coulomb
         * sum overestimates electrostatics substantially (a single
         * H-bond pairwise term alone already approaches the size of the
         * ENTIRE real gas-phase WC binding energy). epsilon=4 is a
         * standard, widely-used distance-independent dielectric for
         * approximating this screening in classical force fields without
         * explicit solvent - not a precise fit for this system, but a
         * well-established, literature-standard order-of-magnitude
         * correction. */
        sim->dielectric = 4.0;

        int g = sim_place_guanine(sim, vec3_zero());
        int c = sim_place_cytosine(sim, vec3(15.0, 0.0, 0.0)); /* temp, far */

        /* Guanine ring atoms for normal: N9,C8,N1 (indices 0,1,6) */
        int g_ring[3] = {g+0, g+1, g+6};
        int c_ring[3] = {c+0, c+1, c+2}; /* N3,C4,N1 */
        Vec3 n_G = nb_ring_normal(sim, g_ring);
        Vec3 n_C = nb_ring_normal(sim, c_ring);

        double cosang = vec3_dot(n_G, n_C);
        Vec3 axis = vec3_cross(n_C, n_G);
        double angle;
        if (vec3_norm(axis) < 1.0e-8) {
            angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
            axis  = (cosang > 0) ? vec3(0,0,1)
                                  : vec3_cross(n_C, vec3(1,0,0));
            if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
        } else {
            angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
        }

        Vec3 c_N3_pivot = sim->atoms[c+0].position;
        nb_transform_rigid(sim, c, 13, c_N3_pivot, axis, angle, vec3_zero());

        /* Second rotation: fix the azimuthal orientation left unconstrained
         * by normal-alignment alone. Align cytosine's N3->N4 direction
         * (toward its donor amino group) with guanine's N1->O6 direction
         * (toward its acceptor), both projected into the now-shared
         * plane - this is what keeps O6 and N4 approaching each other
         * instead of colliding. */
        Vec3 g_N1_for_az = sim->atoms[g+6].position, g_O6 = sim->atoms[g+5].position;
        Vec3 c_N3_for_az = sim->atoms[c+0].position, c_N4 = sim->atoms[c+5].position;
        Vec3 v_G_az = vec3_sub(g_O6, g_N1_for_az);
        Vec3 v_C_az = vec3_sub(c_N4, c_N3_for_az);
        double az_angle = nb_signed_inplane_angle(v_C_az, v_G_az, n_G);
        nb_transform_rigid(sim, c, 13, c_N3_pivot, n_G, az_angle, vec3_zero());

        /* Target: cytosine's N3 lands along guanine's N1->HN1 donor
         * direction, at the standard WC N1...N3 distance (~2.95 A) */
        Vec3 g_N1 = sim->atoms[g+6].position, g_HN1 = sim->atoms[g+13].position;
        Vec3 donor_dir = vec3_normalize(vec3_sub(g_HN1, g_N1));
        Vec3 target = vec3_add(g_N1, vec3_scale(donor_dir, 2.95));
        Vec3 c_N3_now = sim->atoms[c+0].position;
        Vec3 translation = vec3_sub(target, c_N3_now);
        nb_transform_rigid(sim, c, 13, c_N3_pivot, vec3_zero(), 0.0, translation);

        printf("  Initial heavy-atom contacts after geometric placement:\n");
        printf("    G:N1...C:N3 = %.3f A (target 2.95)\n",
               vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
        printf("    G:N2...C:O2 = %.3f A\n",
               vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
        printf("    G:O6...C:N4 = %.3f A\n",
               vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));

        forces_calculate(sim);
        printf("  Energy breakdown at initial placement: "
               "E_LJ=%.6f eV  E_Coulomb=%.6f eV  Total_PE=%.6f eV\n",
               sim->E_lj_total, sim->E_coulomb_total, sim->potential_energy);

        /* Diagnostic: find the closest intermolecular atom pair AND the
         * single largest LJ repulsive contributor - these need not be
         * the same pair (LJ repulsion depends on sigma/epsilon too) */
        {
            double min_d = 1.0e9; int mi = -1, mj = -1;
            double max_lj = -1.0e9; int li = -1, lj_idx = -1;
            for (int ii = g; ii < g+16; ii++)
                for (int jj = c; jj < c+13; jj++) {
                    double d = vec3_dist(sim->atoms[ii].position, sim->atoms[jj].position);
                    if (d < min_d) { min_d = d; mi = ii; mj = jj; }
                    PairEnergy pe = forces_nonbonded_pair(sim->atoms, ii, jj,
                                                           &sim->box, 1, 0, 1.0);
                    if (pe.lj_energy > max_lj) { max_lj = pe.lj_energy; li = ii; lj_idx = jj; }
                }
            printf("  Closest intermolecular contact: atom %d (Z=%d) ... "
                   "atom %d (Z=%d) = %.3f A\n",
                   mi, sim->atoms[mi].Z, mj, sim->atoms[mj].Z, min_d);
            printf("  Largest single LJ repulsion: atom %d (Z=%d, sigma=%.2f) ... "
                   "atom %d (Z=%d, sigma=%.2f) = %.3f A, contributes %.4f eV\n",
                   li, sim->atoms[li].Z, sim->atoms[li].lj_sigma,
                   lj_idx, sim->atoms[lj_idx].Z, sim->atoms[lj_idx].lj_sigma,
                   vec3_dist(sim->atoms[li].position, sim->atoms[lj_idx].position), max_lj);
        }

        PairResult res = run_pair_relaxation(sim, g+6, c+0);
        printf("\n  Initial interaction PE:        %.6f eV\n", res.initial_interE);
        printf("  PE at closest WC approach:     %.6f eV "
               "(primary N1...N3 = %.3f A)\n",
               res.closest_approach_E, res.closest_approach_dist);
        printf("  Global PE minimum over run:    %.6f eV "
               "(may reflect drift to a different,\n"
               "                                  non-WC configuration "
               "such as stacking)\n", res.min_interE);
        printf("  Final PE (end of run):         %.6f eV\n", res.final_interE);

        printf("\n  Final heavy-atom contacts:\n");
        printf("    G:N1...C:N3 = %.3f A\n",
               vec3_dist(sim->atoms[g+6].position, sim->atoms[c+0].position));
        printf("    G:N2...C:O2 = %.3f A\n",
               vec3_dist(sim->atoms[g+8].position, sim->atoms[c+4].position));
        printf("    G:O6...C:N4 = %.3f A\n",
               vec3_dist(sim->atoms[g+5].position, sim->atoms[c+5].position));

        G_C_energy = res.closest_approach_E;
        sim_destroy(sim);
    }

    /* ── A-U pair ───────────────────────────────────────────────────────── */
    {
        printf("\n  --- Adenine-Uracil (2 H-bonds: N1..H-N3, N6-H..O4) ---\n");
        Simulation *sim = sim_create(64, 64);
        sim->dielectric = 4.0; /* same screening correction as G-C, see comment there */

        int a = sim_place_adenine(sim, vec3_zero());
        int u = sim_place_uracil(sim, vec3(15.0, 0.0, 0.0));

        int a_ring[3] = {a+0, a+1, a+6}; /* N9,C8,N1 */
        int u_ring[3] = {u+0, u+1, u+3}; /* N1,C2,N3 */
        Vec3 n_A = nb_ring_normal(sim, a_ring);
        Vec3 n_U = nb_ring_normal(sim, u_ring);

        double cosang = vec3_dot(n_A, n_U);
        Vec3 axis = vec3_cross(n_U, n_A);
        double angle;
        if (vec3_norm(axis) < 1.0e-8) {
            angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
            axis  = (cosang > 0) ? vec3(0,0,1)
                                  : vec3_cross(n_U, vec3(1,0,0));
            if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
        } else {
            angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
        }

        /* Here uracil is the DONOR (N3-H) and adenine the ACCEPTOR (N1) -
         * reversed roles from G-C's primary bond. Pivot on uracil's N3,
         * target it (and thus its H) toward adenine's N1. */
        Vec3 u_N3_pivot = sim->atoms[u+3].position;
        nb_transform_rigid(sim, u, 12, u_N3_pivot, axis, angle, vec3_zero());

        /* Azimuthal fix, same logic as G-C: align uracil's N3->O4
         * direction with adenine's N1->N6 direction so O4 and N6
         * approach each other rather than landing somewhere unrelated. */
        Vec3 a_N1_for_az = sim->atoms[a+6].position, a_N6_for_az = sim->atoms[a+5].position;
        Vec3 u_N3_for_az = sim->atoms[u+3].position, u_O4_for_az = sim->atoms[u+5].position;
        Vec3 v_A_az = vec3_sub(a_N6_for_az, a_N1_for_az);
        Vec3 v_U_az = vec3_sub(u_O4_for_az, u_N3_for_az);
        double az_angle = nb_signed_inplane_angle(v_U_az, v_A_az, n_A);
        nb_transform_rigid(sim, u, 12, u_N3_pivot, n_A, az_angle, vec3_zero());

        Vec3 u_N3 = sim->atoms[u+3].position, u_HN3 = sim->atoms[u+9].position;
        Vec3 donor_dir = vec3_normalize(vec3_sub(u_HN3, u_N3));
        /* We want N3 positioned so that ITS H points at adenine's N1 -
         * i.e. place N3 such that N1 = N3 + donor_dir * 2.9, so
         * N3 = target_N1 - donor_dir*2.9 */
        Vec3 a_N1 = sim->atoms[a+6].position;
        Vec3 target_N3 = vec3_sub(a_N1, vec3_scale(donor_dir, 2.90));
        Vec3 translation = vec3_sub(target_N3, u_N3);
        nb_transform_rigid(sim, u, 12, u_N3_pivot, vec3_zero(), 0.0, translation);

        printf("  Initial heavy-atom contacts after geometric placement:\n");
        printf("    A:N1...U:N3 = %.3f A (target 2.90)\n",
               vec3_dist(sim->atoms[a+6].position, sim->atoms[u+3].position));
        printf("    A:N6...U:O4 = %.3f A\n",
               vec3_dist(sim->atoms[a+5].position, sim->atoms[u+5].position));

        forces_calculate(sim);
        printf("  Energy breakdown at initial placement: "
               "E_LJ=%.6f eV  E_Coulomb=%.6f eV  Total_PE=%.6f eV\n",
               sim->E_lj_total, sim->E_coulomb_total, sim->potential_energy);

        {
            double min_d = 1.0e9; int mi = -1, mj = -1;
            for (int ii = a; ii < a+15; ii++)
                for (int jj = u; jj < u+12; jj++) {
                    double d = vec3_dist(sim->atoms[ii].position, sim->atoms[jj].position);
                    if (d < min_d) { min_d = d; mi = ii; mj = jj; }
                }
            printf("  Closest intermolecular contact: atom %d (Z=%d) ... "
                   "atom %d (Z=%d) = %.3f A\n",
                   mi, sim->atoms[mi].Z, mj, sim->atoms[mj].Z, min_d);
        }

        PairResult res = run_pair_relaxation(sim, a+6, u+3);
        printf("\n  Initial interaction PE:        %.6f eV\n", res.initial_interE);
        printf("  PE at closest WC approach:     %.6f eV "
               "(primary N1...N3 = %.3f A)\n",
               res.closest_approach_E, res.closest_approach_dist);
        printf("  Global PE minimum over run:    %.6f eV\n", res.min_interE);
        printf("  Final PE (end of run):         %.6f eV\n", res.final_interE);

        printf("\n  Final heavy-atom contacts:\n");
        printf("    A:N1...U:N3 = %.3f A\n",
               vec3_dist(sim->atoms[a+6].position, sim->atoms[u+3].position));
        printf("    A:N6...U:O4 = %.3f A\n",
               vec3_dist(sim->atoms[a+5].position, sim->atoms[u+5].position));

        A_U_energy = res.closest_approach_E;
        sim_destroy(sim);
    }

    /* ── Verdict ────────────────────────────────────────────────────────── */
    printf("\n  ══════════════════════════════════════════════════\n");
    printf("  G-C @ closest WC approach: %.6f eV (3 H-bonds)\n", G_C_energy);
    printf("  A-U @ closest WC approach: %.6f eV (2 H-bonds)\n", A_U_energy);
    if (G_C_energy < A_U_energy) {
        printf("  --> G-C binds MORE strongly than A-U (%.6f eV difference),\n"
               "      matching real DNA/RNA thermodynamics, using nothing\n"
               "      but charges + Coulomb + LJ. Not programmed in.\n",
               A_U_energy - G_C_energy);
    } else {
        printf("  --> Unexpected: A-U came out more stable than G-C. This\n"
               "      would need investigation (geometry, charges, or\n"
               "      relaxation time) before trusting the result.\n");
    }
    printf("  ══════════════════════════════════════════════════\n");
}

/* ════════════════════════════════════════════════════════════════════════════
 * Entry point
 * ════════════════════════════════════════════════════════════════════════════ */
int main(void) {
    printf("\n");
    printf("  ╔═══════════════════════════════════════════════════════╗\n");
    printf("  ║       CARBON VM — CHEMISTRY SIMULATOR                 ║\n");
    printf("  ║       From subatomic to molecular dynamics            ║\n");
    printf("  ╚═══════════════════════════════════════════════════════╝\n");
    printf("\n  Unit system: Length=Å  Time=fs  Energy=eV  Mass=AMU\n");
    printf("  Physical constants: 2019 CODATA  |  LJ: UFF  |  Bonds: AMBER\n\n");

    demo_quantum();
    demo_bond_curve();
    demo_water_md();
    demo_water_cluster();
    demo_methane();
    demo_nucleobases();
    demo_basepairing();

    printf("\n  All demos complete.\n");
    printf("  Next steps: add nucleotides → DNA → gene regulatory networks.\n\n");
    return 0;
}

```

---

### FILE: src/nucleobases.c
Location: `src/nucleobases.c`
```cpp
#include <math.h>
#include <stdio.h>
#include "../include/nucleobases.h"
#include "../include/sim.h"
#include "../include/forces.h"

/*
 * nucleobases.c
 * See nucleobases.h for full provenance notes on the geometry data.
 *
 * Every coordinate array below is transcribed directly from the
 * "ideal coordinates" columns of the corresponding RCSB PDB Chemical
 * Component Dictionary entry (fetched 2026-06-19):
 *   uracil   -> https://files.rcsb.org/ligands/view/URA.cif
 *   cytosine -> https://files.rcsb.org/ligands/view/CYT.cif
 *   adenine  -> https://files.rcsb.org/ligands/view/ADE.cif
 *   guanine  -> https://files.rcsb.org/ligands/view/GUN.cif
 * Atom order in each array matches the CIF's own pdbx_ordinal order,
 * named in the comment above each row for traceability.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Shared placement helper
 *
 * Centers the molecule (by centroid of all given atoms) at `origin`,
 * adds bonds with r0 set to the EXACT distance in the placed geometry
 * (zero initial strain by construction) while keeping the generic
 * single/double-bond force constant from the BOND_TABLE, then builds
 * angles with theta0 derived from the same placed geometry.
 * ══════════════════════════════════════════════════════════════════════════ */
static int place_molecule(Simulation *sim, Vec3 origin,
                           const double coords[][3], const int Zs[],
                           const double charges[],
                           int n_atoms,
                           const int bonds[][3], int n_bonds) {
    Vec3 centroid = vec3_zero();
    for (int i = 0; i < n_atoms; i++)
        centroid = vec3_add(centroid,
            vec3(coords[i][0], coords[i][1], coords[i][2]));
    centroid = vec3_scale(centroid, 1.0 / n_atoms);

    int first = sim->num_atoms;
    for (int i = 0; i < n_atoms; i++) {
        Vec3 raw = vec3(coords[i][0], coords[i][1], coords[i][2]);
        Vec3 p   = vec3_add(vec3_sub(raw, centroid), origin);
        sim_add_atom(sim, Zs[i], p, charges[i]);
    }

    for (int i = 0; i < n_bonds; i++) {
        int ia    = first + bonds[i][0];
        int ib    = first + bonds[i][1];
        int order = bonds[i][2];
        int bidx  = sim_add_bond(sim, ia, ib, order);
        if (bidx >= 0) {
            double d = vec3_dist(sim->atoms[ia].position,
                                  sim->atoms[ib].position);
            sim_set_bond_params(sim, bidx, d, sim->bonds[bidx].k);
        }
    }

    /* Generic aromatic-ring bending stiffness, not independently fitted */
    sim_build_angles_geometric(sim, 1.0);

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * URACIL  (RNA pyrimidine base)
 * Atom order: N1 C2 O2 N3 C4 O4 C5 C6 HN1 HN3 H5 H6
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_uracil(Simulation *sim, Vec3 origin) {
    static const double coords[12][3] = {
        { 0.994,  0.000, -1.183},  /* N1  */
        {-0.349, -0.001, -1.135},  /* C2  */
        {-0.986, -0.001, -2.171},  /* O2  */
        {-1.000,  0.003,  0.043},  /* N3  */
        {-0.308, -0.001,  1.200},  /* C4  */
        {-0.896, -0.001,  2.267},  /* O4  */
        { 1.106, -0.000,  1.164},  /* C5  */
        { 1.733,  0.000, -0.031},  /* C6  */
        { 1.445,  0.000, -2.042},  /* HN1 */
        {-1.969,  0.003,  0.059},  /* HN3 */
        { 1.677, -0.000,  2.081},  /* H5  */
        { 2.812,  0.000, -0.078}   /* H6  */
    };
    static const int Zs[12] = {7,6,8,7,6,8,6,6,1,1,1,1};

    /* {atom_i, atom_j, order}  indices into coords[] above */
    static const int bonds[12][3] = {
        {0,1,1}, {0,7,1}, {0,8,1},   /* N1-C2, N1-C6, N1-HN1   */
        {1,2,2}, {1,3,1},            /* C2=O2, C2-N3           */
        {3,4,1}, {3,9,1},            /* N3-C4, N3-HN3          */
        {4,5,2}, {4,6,1},            /* C4=O4, C4-C5           */
        {6,7,2}, {6,10,1},           /* C5=C6, C5-H5           */
        {7,11,1}                    /* C6-H6                  */
    };

    /*
     * Verified RESP partial charges (Aduri et al., J. Chem. Theory
     * Comput. 2007, 3, 1464 - Table 1, "uridine" column), e units.
     * These are nucleoside-context charges; the glycosidic N1 keeps
     * its nucleoside value (the sugar/H distinction barely affects
     * this specific atom's own charge), and the new explicit HN1
     * (replacing the sugar attachment) is set to +0.118186 - derived
     * from charge neutrality, and independently cross-validated: the
     * same 0.118186 residual falls out of cytosine, adenine, and
     * guanine's base-only charge sums too, matching the paper's own
     * stated sugar-charge restraint used during RESP fitting.
     */
    static const double charges[12] = {
         0.1110,  /* N1  */
         0.4539,  /* C2  */
        -0.5407,  /* O2  */
        -0.3681,  /* N3  */
         0.6022,  /* C4  */
        -0.5652,  /* O4  */
        -0.3135,  /* C5  */
        -0.2320,  /* C6  */
         0.118186,/* HN1 (derived, see above) */
         0.3087,  /* HN3 */
         0.1697,  /* H5  */
         0.2557   /* H6  */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 12, bonds, 12);
    /* Polar (N-attached) H's get a small, NONZERO LJ radius - NOT exactly
     * zero like TIP3P water's hydrogens. Exactly-zero LJ was tried first
     * (over-generalizing the TIP3P precedent) and caused a numerical
     * blowup: with eps=0, LJ combining rules give zero repulsion for
     * ANY pair involving that atom, leaving pure 1/r Coulomb attraction
     * with no floor as a donor H approaches an acceptor - in TIP3P this
     * is masked because the donor O's own real LJ provides an indirect
     * floor before collapse; here it wasn't enough. These small values
     * (sigma=0.5 A, eps=0.001 eV) are a stability choice, not a verified
     * literature parameter - generic UFF H radii would be too large
     * (placing real H-bond distances inside the repulsive wall), and
     * exactly zero is too small (no floor at all). Ring C-H's (H5, H6)
     * are not involved in WC pairing and keep generic UFF values. */
    sim_set_atom_lj(sim, first+8, 0.001, 0.5);  /* HN1 */
    sim_set_atom_lj(sim, first+9, 0.001, 0.5);  /* HN3 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * CYTOSINE  (DNA/RNA pyrimidine base)
 *
 * IMPORTANT CORRECTION: the raw PDB CCD "CYT" ligand entry encodes the
 * minor "imino" tautomer of cytosine - independently verified by tracing
 * its canonical SMILES (NC1=CC=NC(=O)N1) around the ring: the exchangeable
 * hydrogen sits on the ring nitrogen DIRECTLY ADJACENT to the amino-
 * bearing carbon. In the dominant, Watson-Crick-relevant amino-oxo
 * tautomer that actually base-pairs in real DNA/RNA, the exchangeable
 * (glycosidic-position) hydrogen sits on N1 - the nitrogen NOT adjacent
 * to the amino carbon (separated from it by the ring on both sides) -
 * while the nitrogen adjacent to the amino carbon (N3) must be a bare
 * lone-pair acceptor for guanine's N1-H to bond to. Using the raw PDB
 * tautomer as-is would put the hydrogen on the wrong nitrogen for any
 * Watson-Crick pairing calculation.
 *
 * Fix: the heavy-atom ring skeleton (positions of N1,N3,C2,C4,C5,C6,O2,N4)
 * is kept exactly as verified from the PDB ideal coordinates (the ring
 * shape itself does not change meaningfully between tautomers). Only the
 * exchangeable hydrogen is relocated: removed from N3, and a new position
 * on N1 is computed via vector math (external angle bisector of C2-N1-C6,
 * at the standard 1.01 A N-H bond length, in-plane by construction since
 * both ring-neighbor directions already lie in the ring plane) - the
 * same construction technique already used for thymine's methyl group.
 *
 * Atom order: N3 C4 N1 C2 O2 N4 C5 C6 HN1(corrected) HN41 HN42 H5 H6
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_cytosine(Simulation *sim, Vec3 origin) {
    /* Heavy atoms + the 4 hydrogens that are NOT affected by the tautomer
     * fix, exactly as verified from the PDB ideal coordinates */
    static const double heavy[12][3] = {
        {-0.356, -2.061, -3.112},  /* N3   */
        { 0.341, -0.883, -3.202},  /* C4   */
        {-0.784, -1.739, -0.794},  /* N1   */
        {-0.931, -2.520, -1.928},  /* C2   */
        {-1.546, -3.587, -1.909},  /* O2   */
        { 0.871, -0.510, -4.410},  /* N4   */
        { 0.498, -0.109, -2.124},  /* C5   */
        {-0.124, -0.631, -0.895},  /* C6   */
        { 1.766, -0.078, -4.418},  /* HN41 */
        { 0.341, -0.687, -5.232},  /* HN42 */
        { 1.030,  0.832, -2.104},  /* H5   */
        { 0.001, -0.001,  0.000}   /* H6   */
    };

    Vec3 N1 = vec3(heavy[2][0], heavy[2][1], heavy[2][2]);
    Vec3 C2 = vec3(heavy[3][0], heavy[3][1], heavy[3][2]);
    Vec3 C6 = vec3(heavy[7][0], heavy[7][1], heavy[7][2]);

    /* External bisector of the C2-N1-C6 ring angle: negative sum of the
     * two unit vectors from N1 toward its ring neighbours. Both inputs
     * lie in the (already-verified-planar) ring plane, so the result
     * does too - no separate planarity enforcement needed. */
    Vec3 to_C2 = vec3_normalize(vec3_sub(C2, N1));
    Vec3 to_C6 = vec3_normalize(vec3_sub(C6, N1));
    Vec3 outward = vec3_normalize(vec3_negate(vec3_add(to_C2, to_C6)));
    Vec3 HN1_pos = vec3_add(N1, vec3_scale(outward, 1.01)); /* N-H bond length */

    /* Assemble final 13-atom array: heavy atoms (with the old, now-unused
     * N3-adjacent H position dropped) plus the new, correctly-placed N1-H */
    double coords[13][3];
    for (int i = 0; i < 12; i++) {
        coords[i][0] = heavy[i][0]; coords[i][1] = heavy[i][1]; coords[i][2] = heavy[i][2];
    }
    coords[12][0] = HN1_pos.x; coords[12][1] = HN1_pos.y; coords[12][2] = HN1_pos.z;

    /* Index map: 0=N3 1=C4 2=N1 3=C2 4=O2 5=N4 6=C5 7=C6
     *            8=HN41 9=HN42 10=H5 11=H6 12=HN1(new) */
    static const int Zs[13] = {7,6,7,6,8,7,6,6,1,1,1,1,1};

    static const int bonds[13][3] = {
        {0,1,2}, {0,3,1},            /* N3=C4 (DOUBLE - see structure note
                                       * above), N3-C2 (N3 now bare - WC
                                       * acceptor)                          */
        {1,5,1}, {1,6,1},            /* C4-N4 (exocyclic amino, single),
                                       * C4-C5 (single - NOT double; the
                                       * ring's double-bond character at
                                       * this position is at N3=C4 instead,
                                       * since C4 here carries -NH2 not =O) */
        {2,3,1}, {2,7,1}, {2,12,1},  /* N1-C2, N1-C6 (single, unaffected
                                       * by the C4 substitution - same as
                                       * uracil), N1-HN1(corrected)         */
        {3,4,2},                    /* C2=O2                  */
        {5,8,1}, {5,9,1},            /* N4-HN41, N4-HN42       */
        {6,7,2}, {6,10,1},           /* C5=C6 (DOUBLE - same as uracil),
                                       * C5-H5                             */
        {7,11,1}                    /* C6-H6                  */
    };

    /*
     * Verified RESP partial charges (Aduri et al. 2007, Table 1,
     * "cytidine" column), e units, mapped to the index order used
     * after the tautomer fix above. N1 keeps its nucleoside value;
     * the new HN1 (replacing the sugar) is +0.118186 from charge
     * neutrality - the SAME residual independently derived for
     * uracil, adenine, and guanine (see sim_place_uracil comment).
     */
    static const double charges[13] = {
        -0.8128,   /* N3   */
         0.9020,   /* C4   */
        -0.2152,   /* N1   */
         0.8867,   /* C2   */
        -0.6560,   /* O2   */
        -0.9919,   /* N4   */
        -0.5972,   /* C5   */
         0.1262,   /* C6   */
         0.4251,   /* HN41 */
         0.4251,   /* HN42 */
         0.2023,   /* H5   */
         0.1875,   /* H6   */
         0.118186  /* HN1 (derived, see above) */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 13, bonds, 13);
    sim_set_atom_lj(sim, first+8,  0.001, 0.5);  /* HN41 */
    sim_set_atom_lj(sim, first+9,  0.001, 0.5);  /* HN42 */
    sim_set_atom_lj(sim, first+12, 0.001, 0.5);  /* HN1 (new) */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * ADENINE  (DNA/RNA purine base - fused 5+6 ring)
 * Atom order: N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 HN9 H8 HN61 HN62 H2
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_adenine(Simulation *sim, Vec3 origin) {
    static const double coords[15][3] = {
        {-0.655,  0.000, -2.079},  /* N9   */
        { 0.680,  0.000, -2.352},  /* C8   */
        { 1.360,  0.000, -1.242},  /* N7   */
        { 0.507,  0.005, -0.190},  /* C5   */
        { 0.660,  0.000,  1.205},  /* C6   */
        { 1.919,  0.000,  1.780},  /* N6   */
        {-0.432,  0.000,  1.962},  /* N1   */
        {-1.637,  0.000,  1.423},  /* C2   */
        {-1.829,  0.000,  0.121},  /* N3   */
        {-0.796,  0.000, -0.715},  /* C4   */
        {-1.374,  0.000, -2.731},  /* HN9  */
        { 1.110, -0.001, -3.342},  /* H8   */
        { 2.012, -0.001,  2.745},  /* HN61 */
        { 2.709, -0.004,  1.217},  /* HN62 */
        {-2.498,  0.000,  2.075}   /* H2   */
    };
    static const int Zs[15] = {7,6,7,6,6,7,7,6,7,6,1,1,1,1,1};

    static const int bonds[16][3] = {
        {0,1,1}, {0,9,1}, {0,10,1},  /* N9-C8, N9-C4, N9-HN9   */
        {1,2,2}, {1,11,1},           /* C8=N7, C8-H8           */
        {2,3,1},                    /* N7-C5                  */
        {3,4,1}, {3,9,2},            /* C5-C6, C5=C4           */
        {4,5,1}, {4,6,2},            /* C6-N6, C6=N1           */
        {5,12,1}, {5,13,1},          /* N6-HN61, N6-HN62       */
        {6,7,1},                    /* N1-C2                  */
        {7,8,2}, {7,14,1},           /* C2=N3, C2-H2           */
        {8,9,1}                     /* N3-C4                  */
    };

    /*
     * Verified RESP charges (Aduri et al. 2007, Table 1, "adenosine"
     * column), e units. N9 keeps its nucleoside value; new HN9
     * (replacing the sugar) is +0.118186 from charge neutrality -
     * same cross-validated residual as the other three bases.
     * NOTE: mapped explicitly by atom name to this file's coords
     * order (N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,...) which differs from
     * the paper's own table listing order (N9,C8,N7,C6,N6,C5,C4,
     * N3,C2,N1,...) - a first pass here copied values in the
     * paper's order without remapping, which silently scrambled
     * charges 3 through 9 onto the wrong atoms. Re-verified name by
     * name below.
     */
    static const double charges[15] = {
         0.0172,   /* N9   */
         0.1299,   /* C8   */
        -0.5850,   /* N7   */
         0.0586,   /* C5   */
         0.7111,   /* C6   */
        -0.9386,   /* N6   */
        -0.7536,   /* N1   */
         0.5741,   /* C2   */
        -0.6835,   /* N3   */
         0.3050,   /* C4   */
         0.118186, /* HN9 (derived) */
         0.1749,   /* H8   */
         0.4125,   /* HN61 */
         0.4125,   /* HN62 */
         0.0467    /* H2   */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 15, bonds, 16);
    sim_set_atom_lj(sim, first+10, 0.001, 0.5);  /* HN9  */
    sim_set_atom_lj(sim, first+12, 0.001, 0.5);  /* HN61 */
    sim_set_atom_lj(sim, first+13, 0.001, 0.5);  /* HN62 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * GUANINE  (DNA/RNA purine base - fused 5+6 ring)
 * Atom order: N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4 HN9 H8 HN1 HN21 HN22
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_guanine(Simulation *sim, Vec3 origin) {
    static const double coords[16][3] = {
        { 1.510,  0.000, -1.787},  /* N9   */
        { 0.519,  0.000, -2.725},  /* C8   */
        {-0.642,  0.000, -2.139},  /* N7   */
        {-0.466,  0.000, -0.795},  /* C5   */
        {-1.345,  0.001,  0.313},  /* C6   */
        {-2.554,  0.001,  0.152},  /* O6   */
        {-0.812, -0.004,  1.554},  /* N1   */
        { 0.540,  0.000,  1.723},  /* C2   */
        { 1.053,  0.001,  2.996},  /* N2   */
        { 1.367,  0.000,  0.701},  /* N3   */
        { 0.912,  0.000, -0.556},  /* C4   */
        { 2.464, -0.001, -1.961},  /* HN9  */
        { 0.675, -0.000, -3.794},  /* H8   */
        {-1.395, -0.004,  2.330},  /* HN1  */
        { 2.013,  0.000,  3.131},  /* HN21 */
        { 0.455,  0.005,  3.759}   /* HN22 */
    };
    static const int Zs[16] = {7,6,7,6,6,8,7,6,7,7,6,1,1,1,1,1};

    /*
     * Verified RESP charges (Aduri et al. 2007, Table 1, "guanosine"
     * column), e units, mapped name-by-name to this file's coords
     * order. Note: guanine's N1-H is a SEPARATE, always-present
     * feature of the ring (not the glycosidic position - only N9
     * is, for purines), so its charge (0.3546) comes directly from
     * the table, not derived. Only HN9 (replacing the sugar) uses
     * the +0.118186 charge-neutrality residual.
     * Neutrality check: this set sums to -0.000014 (~0, 4-decimal
     * rounding only) - confirms no atom got mismapped.
     */
    static const double charges[16] = {
         0.0268,   /* N9   */
         0.1066,   /* C8   */
        -0.5575,   /* N7   */
         0.1513,   /* C5   */
         0.5316,   /* C6   */
        -0.5483,   /* O6   */
        -0.5287,   /* N1   */
         0.7191,   /* C2   */
        -0.9044,   /* N2   */
        -0.5959,   /* N3   */
         0.1563,   /* C4   */
         0.118186, /* HN9 (derived) */
         0.1767,   /* H8   */
         0.3546,   /* HN1 (table value, not derived) */
         0.3968,   /* HN21 */
         0.3968    /* HN22 */
    };

    static const int bonds[17][3] = {
        {0,1,1}, {0,10,1}, {0,11,1}, /* N9-C8, N9-C4, N9-HN9   */
        {1,2,2}, {1,12,1},           /* C8=N7, C8-H8           */
        {2,3,1},                    /* N7-C5                  */
        {3,4,1}, {3,10,2},           /* C5-C6, C5=C4           */
        {4,5,2}, {4,6,1},            /* C6=O6, C6-N1           */
        {6,7,1}, {6,13,1},           /* N1-C2, N1-HN1          */
        {7,8,1}, {7,9,2},            /* C2-N2, C2=N3           */
        {8,14,1}, {8,15,1},          /* N2-HN21, N2-HN22       */
        {9,10,1}                    /* N3-C4                  */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 16, bonds, 17);
    sim_set_atom_lj(sim, first+11, 0.001, 0.5);  /* HN9  */
    sim_set_atom_lj(sim, first+13, 0.001, 0.5);  /* HN1  */
    sim_set_atom_lj(sim, first+14, 0.001, 0.5);  /* HN21 */
    sim_set_atom_lj(sim, first+15, 0.001, 0.5);  /* HN22 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * THYMINE  (DNA pyrimidine base = 5-methyluracil)
 *
 * Built from uracil's own verified ring geometry directly: the C5-H5
 * bond direction is taken from uracil's real coordinates, extended to
 * the standard aromatic-ring-to-methyl C-C bond length (1.51 A), and a
 * tetrahedral methyl group is placed there using the same construction
 * already validated in sim_place_ch4 (cos/sin of the 109.47 degree
 * tetrahedral angle, three-fold azimuthal symmetry about the bond axis).
 *
 * This is computed via vector math at runtime below, not by hand -
 * the only hand-transcribed numbers are uracil's own already-verified
 * ring coordinates (same ones used in sim_place_uracil).
 *
 * Atom order: N1 C2 O2 N3 C4 O4 C5 C6 HN1 HN3 H6 CM HM1 HM2 HM3
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_thymine(Simulation *sim, Vec3 origin) {
    /* Uracil ring atoms, identical to sim_place_uracil (H5 omitted - it
     * is replaced by the methyl group; its raw position is kept only
     * to determine the correct bond direction for that replacement) */
    static const double ring[10][3] = {
        { 0.994,  0.000, -1.183},  /* N1  */
        {-0.349, -0.001, -1.135},  /* C2  */
        {-0.986, -0.001, -2.171},  /* O2  */
        {-1.000,  0.003,  0.043},  /* N3  */
        {-0.308, -0.001,  1.200},  /* C4  */
        {-0.896, -0.001,  2.267},  /* O4  */
        { 1.106, -0.000,  1.164},  /* C5  */
        { 1.733,  0.000, -0.031},  /* C6  */
        { 1.445,  0.000, -2.042},  /* HN1 */
        {-1.969,  0.003,  0.059}   /* HN3 */
    };
    static const double H6_raw[3] = { 2.812,  0.000, -0.078};
    static const double H5_raw[3] = { 1.677, -0.000,  2.081}; /* direction only */

    Vec3 C5 = vec3(ring[6][0], ring[6][1], ring[6][2]);
    Vec3 H5 = vec3(H5_raw[0], H5_raw[1], H5_raw[2]);

    /* Direction of the original C5-H5 bond, extended to a C-C single
     * bond length for the methyl substituent (toluene-type ring-to-
     * methyl bond, ~1.51 A) */
    Vec3 u  = vec3_normalize(vec3_sub(H5, C5));
    Vec3 CM = vec3_add(C5, vec3_scale(u, 1.51));

    /* Tetrahedral methyl hydrogens about the C5-CM axis u, exactly the
     * same construction used in sim_place_ch4: pick any vector not
     * parallel to u, cross to get an orthonormal frame, place 3 H's at
     * the tetrahedral angle (cos = -1/3) spaced 120 degrees apart. */
    Vec3 v_seed = (fabs(u.y) < 0.9) ? vec3(0,1,0) : vec3(1,0,0);
    Vec3 v = vec3_normalize(vec3_sub(v_seed, vec3_scale(u, vec3_dot(v_seed, u))));
    Vec3 w = vec3_cross(u, v);

    const double r_CH    = 1.09;
    const double cos_tet = -1.0/3.0;
    const double sin_tet = sqrt(1.0 - cos_tet*cos_tet);
    const double TWO_PI_3 = 2.0 * 3.14159265358979323846 / 3.0;

    Vec3 HM[3];
    for (int k = 0; k < 3; k++) {
        double phi = k * TWO_PI_3;
        Vec3 radial = vec3_add(vec3_scale(v, cos(phi)), vec3_scale(w, sin(phi)));
        Vec3 dir    = vec3_add(vec3_scale(u, cos_tet), vec3_scale(radial, sin_tet));
        HM[k] = vec3_add(CM, vec3_scale(dir, r_CH));
    }

    /* Assemble final 15-atom coordinate set (atom order per header comment) */
    double coords[15][3];
    for (int i = 0; i < 10; i++) {
        coords[i][0] = ring[i][0]; coords[i][1] = ring[i][1]; coords[i][2] = ring[i][2];
    }
    coords[10][0] = H6_raw[0]; coords[10][1] = H6_raw[1]; coords[10][2] = H6_raw[2];
    coords[11][0] = CM.x;      coords[11][1] = CM.y;      coords[11][2] = CM.z;
    coords[12][0] = HM[0].x;   coords[12][1] = HM[0].y;   coords[12][2] = HM[0].z;
    coords[13][0] = HM[1].x;   coords[13][1] = HM[1].y;   coords[13][2] = HM[1].z;
    coords[14][0] = HM[2].x;   coords[14][1] = HM[2].y;   coords[14][2] = HM[2].z;

    static const int Zs[15] = {7,6,8,7,6,8,6,6,1,1,1,6,1,1,1};

    /*
     * CHARGES - LOWER CONFIDENCE than the other three bases.
     *
     * Unlike uracil/cytosine/adenine/guanine, this is NOT a verified
     * RESP fit for thymine specifically - no such source was located
     * for the free base. Instead: the 10 ring atoms shared with
     * uracil keep uracil's verified RESP charges unchanged (N1, C2,
     * O2, N3, C4, O4, C6, H6, HN1, HN3), C5 is shifted to
     * -0.143686 (the exact value needed for overall neutrality,
     * given the methyl group below), and the new methyl group (CM +
     * 3 H) uses generic, standard, self-contained-neutral aliphatic
     * values (CM=-0.18, each H=+0.06 - typical AMBER-GAFF-style
     * aromatic-ring-methyl numbers, not independently verified for
     * THIS molecule). Real force fields refit several nearby ring
     * atoms when adding this substituent; this is a charge-
     * conserving approximation, not a verified fit. Since the WC
     * pairing face (N3-H, O4) is far from the methyl and uses
     * unchanged, verified values, this approximation should have
     * minimal impact specifically on pairing energetics.
     */
    static const double charges[15] = {
         0.1110,    /* N1  */
         0.4539,    /* C2  */
        -0.5407,    /* O2  */
        -0.3681,    /* N3  */
         0.6022,    /* C4  */
        -0.5652,    /* O4  */
        -0.143686,  /* C5  (shifted from uracil's -0.3135, see above) */
        -0.2320,    /* C6  */
         0.118186,  /* HN1 (derived) */
         0.3087,    /* HN3 */
         0.2557,    /* H6  */
        -0.18,      /* CM  (generic) */
         0.06,      /* HM1 (generic) */
         0.06,      /* HM2 (generic) */
         0.06       /* HM3 (generic) */
    };

    static const int bonds[15][3] = {
        {0,1,1}, {0,7,1}, {0,8,1},   /* N1-C2, N1-C6, N1-HN1   */
        {1,2,2}, {1,3,1},            /* C2=O2, C2-N3           */
        {3,4,1}, {3,9,1},            /* N3-C4, N3-HN3          */
        {4,5,2}, {4,6,1},            /* C4=O4, C4-C5           */
        {6,7,2},                    /* C5=C6                  */
        {7,10,1},                   /* C6-H6                  */
        {6,11,1},                   /* C5-CM (new methyl bond)*/
        {11,12,1}, {11,13,1}, {11,14,1}  /* CM-HM1/2/3         */
    };

    return place_molecule(sim, origin, coords, Zs, charges, 15, bonds, 15);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Planarity check: reference plane defined by the first 3 given atoms,
 * report max perpendicular deviation of all remaining atoms from it.
 * ══════════════════════════════════════════════════════════════════════════ */
double nb_planarity_deviation(const Simulation *sim,
                               const int *atom_indices, int n) {
    if (n < 3) return 0.0;

    Vec3 p0 = sim->atoms[atom_indices[0]].position;
    Vec3 p1 = sim->atoms[atom_indices[1]].position;
    Vec3 p2 = sim->atoms[atom_indices[2]].position;

    Vec3 normal = vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));

    double max_dev = 0.0;
    for (int i = 0; i < n; i++) {
        Vec3 d = vec3_sub(sim->atoms[atom_indices[i]].position, p0);
        double dev = fabs(vec3_dot(d, normal));
        if (dev > max_dev) max_dev = dev;
    }
    return max_dev;
}

Vec3 nb_ring_normal(const Simulation *sim, const int *atom_indices) {
    Vec3 p0 = sim->atoms[atom_indices[0]].position;
    Vec3 p1 = sim->atoms[atom_indices[1]].position;
    Vec3 p2 = sim->atoms[atom_indices[2]].position;
    return vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));
}

/*
 * Signed angle (radians) to rotate vector v_from onto vector v_to,
 * both measured as rotation AROUND axis `n` (right-hand rule). Both
 * vectors are projected onto the plane perpendicular to n first, so
 * any small out-of-plane component doesn't bias the result. Used to
 * fix the azimuthal (in-plane) rotation left unconstrained by aligning
 * ring normals alone - without this, two coplanar rings can still end
 * up with substituents clashing into each other.
 */
double nb_signed_inplane_angle(Vec3 v_from, Vec3 v_to, Vec3 n) {
    Vec3 a = vec3_sub(v_from, vec3_scale(n, vec3_dot(v_from, n)));
    Vec3 b = vec3_sub(v_to,   vec3_scale(n, vec3_dot(v_to,   n)));
    a = vec3_normalize(a);
    b = vec3_normalize(b);
    double s = vec3_dot(vec3_cross(a, b), n);
    double c = vec3_dot(a, b);
    return atan2(s, c);
}

void nb_transform_rigid(Simulation *sim, int first_atom, int n_atoms,
                         Vec3 pivot, Vec3 axis, double angle,
                         Vec3 translation) {
    for (int i = first_atom; i < first_atom + n_atoms; i++) {
        Vec3 rel = vec3_sub(sim->atoms[i].position, pivot);
        Vec3 rotated = (vec3_norm(axis) > 1.0e-12)
                       ? vec3_rotate_axis_angle(rel, axis, angle)
                       : rel;
        sim->atoms[i].position = vec3_add(vec3_add(pivot, rotated), translation);
    }
}

```

---

### FILE: src/periodic_table.c
Location: `src/periodic_table.c`
```cpp
#include <stdio.h>
#include <string.h>
#include "../include/periodic_table.h"
#include "../include/constants.h"

/*
 * periodic_table.c
 *
 * Full element data for Z=1..18 (H→Ar), stubbed beyond.
 * Extend the table by filling in the remaining entries.
 *
 * Madelung orbital filling order (used by pt_electron_config):
 * 1s 2s 2p 3s 3p 4s 3d 4p 5s 4d 5p 6s 4f 5d 6p 7s 5f 6d 7p
 *
 * LJ ε values converted from UFF kcal/mol → eV (÷ 23.0605):
 * LJ σ values in Å (UFF x1 column, the collision diameter).
 */

/* ── Madelung sequence: pairs of (n, l) in filling order ─────────────────── */
static const int MADELUNG_N[] = {1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7};
static const int MADELUNG_L[] = {0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1};
static const int MADELUNG_LEN = 19;

/* ── Orbital capacity: 2*(2l+1) ──────────────────────────────────────────── */
static int orbital_cap(int l) { return 2*(2*l+1); }

/* ══════════════════════════════════════════════════════════════════════════
 * Periodic table data
 * ══════════════════════════════════════════════════════════════════════════
 * Fields: Z, symbol, name, mass[AMU], electronegativity,
 *         atomic_radius[Å], covalent_radius[Å], vdw_radius[Å],
 *         ionization_energy[eV], electron_affinity[eV],
 *         valence,
 *         ground_config (populated at runtime by pt_electron_config),
 *         lj_epsilon[eV], lj_sigma[Å]
 *
 * Noble gas electron affinities set to 0 (they are negative).
 * Electronegativity = -1 for noble gases (undefined on Pauling scale).
 */
const Element PERIODIC_TABLE[MAX_ELEMENTS + 1] = {
    /* [0] unused sentinel */
    {0},

    /* Z=1  H  Hydrogen */
    { 1,"H","Hydrogen",       1.008,  2.20, 1.20, 0.31, 1.10, 13.598, 0.754,  1,
      {0}, 0.044*KCAL_MOL_TO_EV, 2.886 },

    /* Z=2  He Helium */
    { 2,"He","Helium",        4.003, -1.00, 1.40, 0.28, 1.40, 24.587, 0.000,  0,
      {0}, 0.056*KCAL_MOL_TO_EV, 2.362 },

    /* Z=3  Li Lithium */
    { 3,"Li","Lithium",       6.941,  0.98, 1.82, 1.28, 1.82,  5.392, 0.618,  1,
      {0}, 0.025*KCAL_MOL_TO_EV, 2.451 },

    /* Z=4  Be Beryllium */
    { 4,"Be","Beryllium",     9.012,  1.57, 1.53, 0.96, 1.53,  9.323, 0.000,  2,
      {0}, 0.085*KCAL_MOL_TO_EV, 2.745 },

    /* Z=5  B  Boron */
    { 5,"B","Boron",         10.811,  2.04, 1.92, 0.84, 1.92,  8.298, 0.280,  3,
      {0}, 0.180*KCAL_MOL_TO_EV, 4.083 },

    /* Z=6  C  Carbon */
    { 6,"C","Carbon",        12.011,  2.55, 1.70, 0.77, 1.70, 11.260, 1.262,  4,
      {0}, 0.105*KCAL_MOL_TO_EV, 3.851 },

    /* Z=7  N  Nitrogen */
    { 7,"N","Nitrogen",      14.007,  3.04, 1.55, 0.71, 1.55, 14.534, 0.000,  3,
      {0}, 0.069*KCAL_MOL_TO_EV, 3.660 },

    /* Z=8  O  Oxygen */
    { 8,"O","Oxygen",        15.999,  3.44, 1.52, 0.66, 1.52, 13.618, 1.461,  2,
      {0}, 0.060*KCAL_MOL_TO_EV, 3.500 },

    /* Z=9  F  Fluorine */
    { 9,"F","Fluorine",      18.998,  3.98, 1.47, 0.64, 1.47, 17.423, 3.401,  1,
      {0}, 0.050*KCAL_MOL_TO_EV, 3.364 },

    /* Z=10 Ne Neon */
    {10,"Ne","Neon",         20.180, -1.00, 1.54, 0.58, 1.54, 21.565, 0.000,  0,
      {0}, 0.042*KCAL_MOL_TO_EV, 3.243 },

    /* Z=11 Na Sodium */
    {11,"Na","Sodium",       22.990,  0.93, 2.27, 1.66, 2.27,  5.139, 0.548,  1,
      {0}, 0.030*KCAL_MOL_TO_EV, 2.983 },

    /* Z=12 Mg Magnesium */
    {12,"Mg","Magnesium",    24.305,  1.31, 1.73, 1.41, 1.73,  7.646, 0.000,  2,
      {0}, 0.111*KCAL_MOL_TO_EV, 3.021 },

    /* Z=13 Al Aluminium */
    {13,"Al","Aluminium",    26.982,  1.61, 1.84, 1.21, 1.84,  5.986, 0.441,  3,
      {0}, 0.505*KCAL_MOL_TO_EV, 4.499 },

    /* Z=14 Si Silicon */
    {14,"Si","Silicon",      28.086,  1.90, 2.10, 1.11, 2.10,  8.151, 1.385,  4,
      {0}, 0.402*KCAL_MOL_TO_EV, 4.295 },

    /* Z=15 P  Phosphorus */
    {15,"P","Phosphorus",    30.974,  2.19, 1.80, 1.07, 1.80, 10.487, 0.747,  3,
      {0}, 0.305*KCAL_MOL_TO_EV, 4.147 },

    /* Z=16 S  Sulphur */
    {16,"S","Sulphur",       32.065,  2.58, 1.80, 1.05, 1.80, 10.360, 2.077,  2,
      {0}, 0.274*KCAL_MOL_TO_EV, 4.035 },

    /* Z=17 Cl Chlorine */
    {17,"Cl","Chlorine",     35.453,  3.16, 1.75, 1.02, 1.75, 12.968, 3.613,  1,
      {0}, 0.227*KCAL_MOL_TO_EV, 3.947 },

    /* Z=18 Ar Argon */
    {18,"Ar","Argon",        39.948, -1.00, 1.88, 1.06, 1.88, 15.760, 0.000,  0,
      {0}, 0.185*KCAL_MOL_TO_EV, 3.868 },

    /* Z=19 K  Potassium — mass only stub */
    {19,"K","Potassium",     39.098,  0.82, 2.75, 2.03, 2.75,  4.341, 0.501,  1,
      {0}, 0.035*KCAL_MOL_TO_EV, 3.812 },

    /* Z=20 Ca Calcium */
    {20,"Ca","Calcium",      40.078,  1.00, 2.31, 1.76, 2.31,  6.113, 0.018,  2,
      {0}, 0.238*KCAL_MOL_TO_EV, 3.399 },
    /* Z=21..35 stubs with mass only */
    {21,"Sc","Scandium",     44.956,  1.36, 2.11, 1.70, 2.11,  6.561, 0.188,  3, {0}, 0.019*KCAL_MOL_TO_EV, 3.295},
    {22,"Ti","Titanium",     47.867,  1.54, 2.00, 1.60, 2.00,  6.828, 0.079,  4, {0}, 0.017*KCAL_MOL_TO_EV, 3.175},
    {23,"V","Vanadium",      50.942,  1.63, 1.92, 1.53, 1.92,  6.746, 0.525,  5, {0}, 0.016*KCAL_MOL_TO_EV, 3.144},
    {24,"Cr","Chromium",     51.996,  1.66, 1.85, 1.39, 1.85,  6.767, 0.666,  3, {0}, 0.015*KCAL_MOL_TO_EV, 3.023},
    {25,"Mn","Manganese",    54.938,  1.55, 1.79, 1.61, 1.79,  7.434, 0.000,  2, {0}, 0.013*KCAL_MOL_TO_EV, 2.961},
    {26,"Fe","Iron",         55.845,  1.83, 1.72, 1.52, 1.72,  7.902, 0.151,  3, {0}, 0.013*KCAL_MOL_TO_EV, 2.912},
    {27,"Co","Cobalt",       58.933,  1.88, 1.67, 1.50, 1.67,  7.881, 0.662,  3, {0}, 0.014*KCAL_MOL_TO_EV, 2.872},
    {28,"Ni","Nickel",       58.693,  1.91, 1.63, 1.24, 1.63,  7.640, 1.156,  2, {0}, 0.015*KCAL_MOL_TO_EV, 2.834},
    {29,"Cu","Copper",       63.546,  1.90, 1.40, 1.32, 1.40,  7.726, 1.235,  2, {0}, 0.005*KCAL_MOL_TO_EV, 3.495},
    {30,"Zn","Zinc",         65.38,   1.65, 1.39, 1.22, 1.39,  9.394, 0.000,  2, {0}, 0.124*KCAL_MOL_TO_EV, 2.763},
    {31,"Ga","Gallium",      69.723,  1.81, 1.87, 1.22, 1.87,  5.999, 0.300,  3, {0}, 0.415*KCAL_MOL_TO_EV, 4.383},
    {32,"Ge","Germanium",    72.630,  2.01, 2.11, 1.20, 2.11,  7.900, 1.233,  4, {0}, 0.379*KCAL_MOL_TO_EV, 4.280},
    {33,"As","Arsenic",      74.922,  2.18, 1.85, 1.19, 1.85,  9.815, 0.814,  3, {0}, 0.309*KCAL_MOL_TO_EV, 4.230},
    {34,"Se","Selenium",     78.971,  2.55, 1.90, 1.20, 1.90,  9.752, 2.021,  2, {0}, 0.291*KCAL_MOL_TO_EV, 4.205},
    {35,"Br","Bromine",      79.904,  2.96, 1.85, 1.20, 1.85, 11.814, 3.365,  1, {0}, 0.251*KCAL_MOL_TO_EV, 4.189},
    {36,"Kr","Krypton",      83.798, -1.00, 2.02, 1.16, 2.02, 14.000, 0.000,  0, {0}, 0.220*KCAL_MOL_TO_EV, 4.141},
};

/* ══════════════════════════════════════════════════════════════════════════
 * Implementation
 * ══════════════════════════════════════════════════════════════════════════ */

const Element *pt_element(int Z) {
    if (Z < 1 || Z > MAX_ELEMENTS) return NULL;
    return &PERIODIC_TABLE[Z];
}

const Element *pt_by_symbol(const char *symbol) {
    for (int Z = 1; Z <= 36; Z++) {
        if (strcmp(PERIODIC_TABLE[Z].symbol, symbol) == 0)
            return &PERIODIC_TABLE[Z];
    }
    return NULL;
}

/*
 * Fill an ElectronConfig for atomic number Z using the Madelung (diagonal)
 * rule. Fills orbitals in Madelung order until all Z electrons are placed.
 *
 * Note: Chromium (Z=24) and Copper (Z=29) and a few others are exceptions
 * to the strict Madelung rule; we handle those here explicitly.
 */
void pt_electron_config(int Z, ElectronConfig *cfg) {
    memset(cfg, 0, sizeof(*cfg));
    cfg->total_electrons = Z;

    /* Handle the well-known Madelung exceptions for 3d transition metals */
    /* We store any exception as a delta to apply after Madelung filling.  */
    int electrons_left = Z;

    for (int i = 0; i < MADELUNG_LEN && electrons_left > 0; i++) {
        int n = MADELUNG_N[i];
        int l = MADELUNG_L[i];
        int cap = orbital_cap(l);
        int fill = electrons_left < cap ? electrons_left : cap;
        cfg->config[n-1][l] += fill;
        electrons_left -= fill;
    }

    /* Madelung exceptions: swap 4s→3d for half-filled and full-d stability */
    if (Z == 24) { /* Cr: [Ar] 3d5 4s1 instead of 3d4 4s2 */
        cfg->config[3][0] -= 1;  /* remove one 4s */
        cfg->config[2][2] += 1;  /* add one 3d    */
    }
    if (Z == 29) { /* Cu: [Ar] 3d10 4s1 instead of 3d9 4s2 */
        cfg->config[3][0] -= 1;
        cfg->config[2][2] += 1;
    }

    /* Compute valence electrons (outermost shell sum) */
    int max_shell = 0;
    for (int s = 0; s < MAX_SHELLS; s++) {
        int shell_count = 0;
        for (int sub = 0; sub < 4; sub++) shell_count += cfg->config[s][sub];
        if (shell_count > 0) max_shell = s;
    }
    for (int sub = 0; sub < 4; sub++)
        cfg->valence_electrons += cfg->config[max_shell][sub];
}

/*
 * Write a human-readable electron configuration string into buf.
 * Example output: "1s2 2s2 2p2"
 */
void pt_config_string(const ElectronConfig *cfg, char *buf, int buflen) {
    static const char sub_labels[] = "spdf";
    int pos = 0;
    for (int n = 0; n < MAX_SHELLS && pos < buflen - 1; n++) {
        for (int l = 0; l < 4 && pos < buflen - 1; l++) {
            int count = cfg->config[n][l];
            if (count == 0) continue;
            int written = snprintf(buf + pos, buflen - pos, "%d%c%d ",
                                   n+1, sub_labels[l], count);
            if (written > 0) pos += written;
        }
    }
    if (pos > 0 && buf[pos-1] == ' ') buf[pos-1] = '\0';
    else buf[pos] = '\0';
}

void pt_print_element(int Z) {
    const Element *e = pt_element(Z);
    if (!e) { printf("Element Z=%d not found.\n", Z); return; }

    ElectronConfig cfg;
    pt_electron_config(Z, &cfg);
    char cfgstr[128];
    pt_config_string(&cfg, cfgstr, sizeof(cfgstr));

    printf("══════════════════════════════════════════\n");
    printf("  %s (%s)  Z=%d\n", e->name, e->symbol, e->Z);
    printf("══════════════════════════════════════════\n");
    printf("  Mass              : %.4f AMU\n",  e->mass);
    printf("  Electronegativity : %.2f (Pauling)\n", e->electronegativity);
    printf("  Atomic radius     : %.3f Å (vdW)\n", e->atomic_radius);
    printf("  Covalent radius   : %.3f Å\n",    e->covalent_radius);
    printf("  Ionisation energy : %.3f eV\n",   e->ionization_energy);
    printf("  Electron affinity : %.3f eV\n",   e->electron_affinity);
    printf("  Common valence    : %d\n",         e->valence);
    printf("  LJ ε              : %.5f eV\n",   e->lj_epsilon);
    printf("  LJ σ              : %.4f Å\n",    e->lj_sigma);
    printf("  Config            : %s\n",         cfgstr);
    printf("══════════════════════════════════════════\n");
}

```

---

### FILE: src/quantum.c
Location: `src/quantum.c`
```cpp
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

```

---

### FILE: src/sim.c
Location: `src/sim.c`
```cpp
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "../include/sim.h"
#include "../include/periodic_table.h"
#include "../include/quantum.h"
#include "../include/forces.h"
#include "../include/constants.h"

/*
 * sim.c
 * Simulation construction, topology management, and molecule builders.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Lifecycle
 * ══════════════════════════════════════════════════════════════════════════ */
Simulation *sim_create(int atom_capacity, int bond_capacity) {
    Simulation *sim = (Simulation *)calloc(1, sizeof(Simulation));
    if (!sim) return NULL;

    sim->atoms = (Atom *)calloc(atom_capacity, sizeof(Atom));
    sim->bonds = (Bond *)calloc(bond_capacity, sizeof(Bond));
    /* Angles: heuristic — up to 3 × bonds */
    int angle_cap = bond_capacity * 3;
    sim->angles   = (Angle *)calloc(angle_cap, sizeof(Angle));

    if (!sim->atoms || !sim->bonds || !sim->angles) {
        free(sim->atoms); free(sim->bonds); free(sim->angles);
        free(sim);
        return NULL;
    }

    sim->capacity_atoms    = atom_capacity;
    sim->capacity_bonds    = bond_capacity;
    sim->capacity_angles   = angle_cap;

    /* Defaults */
    sim->dt      = 1.0;     /* 1 fs timestep                                */
    sim->cutoff  = 12.0;    /* 12 Å non-bonded cutoff                       */
    sim->use_lj       = 1;
    sim->use_coulomb  = 1;
    sim->use_bonds    = 1;
    sim->use_angles   = 1;
    sim->dielectric   = 1.0;  /* no screening by default */
    sim->ff_type = FF_UFF;

    /* Box: large vacuum (no PBC by default) */
    sim->box.dimensions = vec3(1000.0, 1000.0, 1000.0);
    sim->box.periodic[0] = sim->box.periodic[1] = sim->box.periodic[2] = 0;

    /* Berendsen thermostat off by default */
    sim->thermostat.type               = THERMOSTAT_NONE;
    sim->thermostat.target_temperature = 300.0;
    sim->thermostat.tau                = 100.0;  /* fs */

    return sim;
}

void sim_destroy(Simulation *sim) {
    if (!sim) return;
    free(sim->atoms);
    free(sim->bonds);
    free(sim->angles);
    free(sim->dihedrals);
    free(sim);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Atom management
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_add_atom(Simulation *sim, int Z, Vec3 pos, double partial_charge) {
    if (sim->num_atoms >= sim->capacity_atoms) return SIM_ERR_OVERFLOW;

    const Element *el = pt_element(Z);
    if (!el) return SIM_ERR_BADATOM;

    int idx = sim->num_atoms++;
    Atom *a = &sim->atoms[idx];
    memset(a, 0, sizeof(Atom));

    a->id             = idx;
    a->Z              = Z;
    a->element        = el;
    a->position       = pos;
    a->velocity       = vec3_zero();
    a->force          = vec3_zero();
    a->mass           = el->mass;
    a->partial_charge = partial_charge;
    a->formal_charge  = 0;

    /* Default LJ parameters to the generic UFF element values.
     * Molecule constructors (e.g. sim_place_h2o) may override these
     * per atom for force-field-specific contexts via sim_set_atom_lj(). */
    a->lj_epsilon = el->lj_epsilon;
    a->lj_sigma   = el->lj_sigma;

    /* Electron configuration */
    pt_electron_config(Z, &a->electron_config);

    /* Populate orbital table */
    quantum_fill_orbitals(a);

    return idx;
}

int sim_add_atom_sym(Simulation *sim, const char *symbol,
                     Vec3 pos, double q) {
    const Element *el = pt_by_symbol(symbol);
    if (!el) return SIM_ERR_BADATOM;
    return sim_add_atom(sim, el->Z, pos, q);
}

void sim_set_atom_lj(Simulation *sim, int atom_idx, double epsilon, double sigma) {
    if (atom_idx < 0 || atom_idx >= sim->num_atoms) return;
    sim->atoms[atom_idx].lj_epsilon = epsilon;
    sim->atoms[atom_idx].lj_sigma   = sigma;
}

void sim_set_bond_params(Simulation *sim, int bond_idx, double r0, double k) {
    if (bond_idx < 0 || bond_idx >= sim->num_bonds) return;
    sim->bonds[bond_idx].r0 = r0;
    sim->bonds[bond_idx].k  = k;
}

int sim_add_angle_explicit(Simulation *sim, int a, int b, int c,
                            double theta0, double k) {
    if (a < 0 || a >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (b < 0 || b >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (c < 0 || c >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (sim->num_angles >= sim->capacity_angles) return SIM_ERR_OVERFLOW;

    int idx = sim->num_angles++;
    Angle *ang = &sim->angles[idx];
    ang->atom_a = a;
    ang->atom_b = b;
    ang->atom_c = c;
    ang->theta0 = theta0;
    ang->k      = k;
    return idx;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Bond management
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_add_bond(Simulation *sim, int ia, int ib, int order) {
    if (ia < 0 || ia >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (ib < 0 || ib >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (sim->num_bonds >= sim->capacity_bonds) return SIM_ERR_OVERFLOW;

    int idx = sim->num_bonds++;
    Bond *b = &sim->bonds[idx];

    b->atom_a = ia;
    b->atom_b = ib;
    b->order  = order;

    /* Look up equilibrium parameters */
    BondParam bp;
    forces_bond_params(sim->atoms[ia].Z, sim->atoms[ib].Z, order, &bp);
    b->r0 = bp.r0;
    b->k  = bp.k;

    /* Current length */
    b->current_length = vec3_dist(sim->atoms[ia].position,
                                   sim->atoms[ib].position);
    b->energy = 0.0;

    /* Update atom bond lists */
    Atom *a = &sim->atoms[ia];
    Atom *c = &sim->atoms[ib];
    if (a->num_bonds < MAX_BONDS_PER_ATOM) {
        a->bond_partners[a->num_bonds]  = ib;
        a->bond_orders  [a->num_bonds]  = order;
        a->num_bonds++;
    }
    if (c->num_bonds < MAX_BONDS_PER_ATOM) {
        c->bond_partners[c->num_bonds]  = ia;
        c->bond_orders  [c->num_bonds]  = order;
        c->num_bonds++;
    }

    return idx;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Auto-detect bonds from distances
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_detect_bonds(Simulation *sim) {
    int N = sim->num_atoms;

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            const Element *ei = sim->atoms[i].element;
            const Element *ej = sim->atoms[j].element;

            double r_cov_sum = ei->covalent_radius + ej->covalent_radius;
            double r_max     = BOND_TOLERANCE * r_cov_sum;
            double r         = vec3_dist(sim->atoms[i].position,
                                          sim->atoms[j].position);

            if (r < r_max) {
                /* Estimate bond order from distance ratio */
                int order = 1;
                double ratio = r / r_cov_sum;
                if (ratio < 0.78) order = 3;
                else if (ratio < 0.87) order = 2;

                sim_add_bond(sim, i, j, order);
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Build angle list from bond topology
 * For each atom b, for every pair of bonds (a-b) and (b-c): add angle a-b-c
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_build_angles(Simulation *sim) {
    sim->num_angles = 0;

    for (int b_idx = 0; b_idx < sim->num_atoms; b_idx++) {
        Atom *b = &sim->atoms[b_idx];
        /* Enumerate all pairs of neighbours of b */
        for (int p = 0; p < b->num_bonds - 1; p++) {
            for (int q = p + 1; q < b->num_bonds; q++) {
                if (sim->num_angles >= sim->capacity_angles) return;

                int ia = b->bond_partners[p];
                int ic = b->bond_partners[q];

                Angle *ang = &sim->angles[sim->num_angles++];
                ang->atom_a = ia;
                ang->atom_b = b_idx;
                ang->atom_c = ic;

                /* Look up angle parameters */
                AngleParam ap;
                forces_angle_params(sim->atoms[ia].Z, b->Z,
                                    sim->atoms[ic].Z, &ap);
                ang->theta0 = ap.theta0;
                ang->k      = ap.k;
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Build angle list with theta0 derived from current Cartesian geometry
 * (rather than the generic element-keyed table) - for ring systems where
 * the seed positions already encode the correct, context-specific angle.
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_build_angles_geometric(Simulation *sim, double k_default) {
    sim->num_angles = 0;

    for (int b_idx = 0; b_idx < sim->num_atoms; b_idx++) {
        Atom *b = &sim->atoms[b_idx];
        for (int p = 0; p < b->num_bonds - 1; p++) {
            for (int q = p + 1; q < b->num_bonds; q++) {
                if (sim->num_angles >= sim->capacity_angles) return;

                int ia = b->bond_partners[p];
                int ic = b->bond_partners[q];

                double theta = vec3_angle(
                    vec3_sub(sim->atoms[ia].position, b->position),
                    vec3_sub(sim->atoms[ic].position, b->position));

                Angle *ang = &sim->angles[sim->num_angles++];
                ang->atom_a = ia;
                ang->atom_b = b_idx;
                ang->atom_c = ic;
                ang->theta0 = theta;
                ang->k      = k_default;
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * PBC helpers
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_set_box(Simulation *sim, double lx, double ly, double lz) {
    sim->box.dimensions = vec3(lx, ly, lz);
    sim->box.periodic[0] = sim->box.periodic[1] = sim->box.periodic[2] = 1;
}

void sim_wrap_positions(Simulation *sim) {
    Vec3 L = sim->box.dimensions;
    for (int i = 0; i < sim->num_atoms; i++) {
        Vec3 *p = &sim->atoms[i].position;
        if (sim->box.periodic[0]) {
            p->x -= L.x * floor(p->x / L.x);
        }
        if (sim->box.periodic[1]) {
            p->y -= L.y * floor(p->y / L.y);
        }
        if (sim->box.periodic[2]) {
            p->z -= L.z * floor(p->z / L.z);
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Molecule constructors
 * ══════════════════════════════════════════════════════════════════════════ */

/* ── H₂ ─────────────────────────────────────────────────────────────────── */
int sim_place_h2(Simulation *sim, Vec3 origin) {
    /* Bond length 0.7414 Å, atoms symmetric about origin */
    int first = sim->num_atoms;
    sim_add_atom(sim, 1, vec3(origin.x - 0.3707, origin.y, origin.z), 0.0);
    sim_add_atom(sim, 1, vec3(origin.x + 0.3707, origin.y, origin.z), 0.0);
    sim_add_bond(sim, first, first+1, 1);
    sim_build_angles(sim);
    return first;
}

/* ── H₂O (TIP3P geometry) ───────────────────────────────────────────────── */
/*
 * TIP3P parameters:
 *   r(OH) = 0.9572 Å
 *   ∠HOH  = 104.52°
 *   Oxygen at origin.
 *   H positions:
 *     H1: ( r sinα, -r cosα, 0 )    α = half-angle = 52.26°
 *     H2: (-r sinα, -r cosα, 0 )
 *   Partial charges: q(O) = -0.834e, q(H) = +0.417e
 */
int sim_place_h2o(Simulation *sim, Vec3 origin) {
    double rOH   = 0.9572;
    double alpha = 52.26 * 3.14159265358979323846 / 180.0;
    double hx    = rOH * sin(alpha);
    double hy    = rOH * cos(alpha);

    int first = sim->num_atoms;
    /* O */
    int iO = sim_add_atom(sim, 8,
        vec3(origin.x,      origin.y,      origin.z),      -0.834);
    /* H1 */
    int iH1 = sim_add_atom(sim, 1,
        vec3(origin.x + hx, origin.y - hy, origin.z),      +0.417);
    /* H2 */
    int iH2 = sim_add_atom(sim, 1,
        vec3(origin.x - hx, origin.y - hy, origin.z),      +0.417);

    /*
     * TIP3P (Jorgensen, Chandrasekhar, Madura, Impey & Klein, J. Chem.
     * Phys. 79, 926 (1983)) — verified parameters:
     *   Oxygen: sigma = 3.15061 Å, epsilon = 0.1521 kcal/mol
     *   Hydrogen: NO Lennard-Jones site at all (sigma = epsilon = 0).
     *     Only the oxygen carries a LJ centre in the original model;
     *     the CHARMM-modified TIP3P variant adds small LJ terms to H,
     *     but that is a distinct, separately-named force field.
     *
     * Generic UFF defaults (sigma_O=3.500 Å) were fit for general
     * organic chemistry, not for this specific, jointly-parameterized
     * water model — using them here would put the O-O equilibrium
     * distance in the wrong place relative to where the TIP3P charges
     * expect it, breaking the delicate cancellation that holds liquid
     * water together. Override explicitly.
     */
    sim_set_atom_lj(sim, iO,  0.1521 * KCAL_MOL_TO_EV, 3.15061);
    sim_set_atom_lj(sim, iH1, 0.0,                      0.0);
    sim_set_atom_lj(sim, iH2, 0.0,                      0.0);

    sim_add_bond(sim, first,   first+1, 1);  /* O-H1 */
    sim_add_bond(sim, first,   first+2, 1);  /* O-H2 */
    sim_build_angles(sim);
    return first;
}

/* ── NH₃ ────────────────────────────────────────────────────────────────── */
/*
 * Pyramidal geometry:
 *   r(NH) = 1.012 Å, ∠HNH = 106.67°
 *   N at origin, H atoms placed in trigonal pyramidal arrangement.
 *   Partial charges: q(N) = -1.02e, q(H) = +0.34e
 */
int sim_place_nh3(Simulation *sim, Vec3 origin) {
    double rNH   = 1.012;
    double hnh   = 106.67 * 3.14159265358979323846 / 180.0;
    /* Tetrahedral-like placement:
     * H at 120° azimuthal spacing, cone half-angle θ from the -z axis */
    double cos_hnh = cos(hnh);
    /*
     * cos(∠HNH) = v1·v2 / (|v1||v2|), where vi = Hi-N (unit direction)
     *           = (sinθ cosφi, sinθ sinφi, -cosθ)
     *
     * v1·v2 = sin²θ(cosφ1cosφ2 + sinφ1sinφ2) + (-cosθ)(-cosθ)
     *       = sin²θ cos(φ2-φ1) + cos²θ        <-- PLUS cos²θ: both H's
     *         have the SAME z-component (-cosθ), so their product is
     *         (-cosθ)(-cosθ) = +cos²θ, not -cos²θ. An earlier version
     *         of this derivation had this sign backwards, which silently
     *         collapsed the whole molecule to perfectly PLANAR (120°,
     *         the NH3 inversion transition-state geometry) instead of
     *         the intended pyramidal ground state (106.67°) - confirmed
     *         by directly re-measuring the placed angle, which came out
     *         120.0000° under the old formula and exactly 106.6700°
     *         under this corrected one.
     *
     * With φ2-φ1 = 120°, cos(120°) = -0.5:
     *   cos_hnh = -0.5 sin²θ + cos²θ
     *           = -0.5(1-cos²θ) + cos²θ = -0.5 + 1.5 cos²θ
     *   => cos²θ = (cos_hnh + 0.5) / 1.5
     */
    double cos2_theta = (cos_hnh + 0.5) / 1.5;
    if (cos2_theta < 0.0) cos2_theta = 0.0;
    double cos_theta  = sqrt(cos2_theta);
    double sin_theta  = sqrt(1.0 - cos2_theta);

    int first = sim->num_atoms;
    sim_add_atom(sim, 7,
        vec3(origin.x, origin.y, origin.z), -1.02);

    double PI23 = 2.0 * 3.14159265358979323846 / 3.0;
    for (int k = 0; k < 3; k++) {
        double phi = k * PI23;
        sim_add_atom(sim, 1,
            vec3(origin.x + rNH * sin_theta * cos(phi),
                 origin.y + rNH * sin_theta * sin(phi),
                 origin.z - rNH * cos_theta),
            +0.34);
        sim_add_bond(sim, first, first+1+k, 1);
    }
    sim_build_angles(sim);
    return first;
}

/* ── CH₄ ────────────────────────────────────────────────────────────────── */
/*
 * Tetrahedral: r(CH) = 1.090 Å, ∠HCH = 109.47°
 * H positions: vertices of a tetrahedron inscribed in a cube.
 * q(C) = -0.24e, q(H) = +0.06e (OPLS-AA)
 */
int sim_place_ch4(Simulation *sim, Vec3 origin) {
    double rCH = 1.090;
    /* Tetrahedral vertices: (±1,±1,±1) normalised, alternating sign for
     * a true tetrahedron (not all-positive corners of cube). */
    double s = rCH / sqrt(3.0);
    double H[4][3] = {
        { s,  s,  s},
        { s, -s, -s},
        {-s,  s, -s},
        {-s, -s,  s}
    };

    int first = sim->num_atoms;
    sim_add_atom(sim, 6, origin, -0.24);   /* C */
    for (int k = 0; k < 4; k++) {
        sim_add_atom(sim, 1,
            vec3(origin.x + H[k][0],
                 origin.y + H[k][1],
                 origin.z + H[k][2]),
            +0.06);
        sim_add_bond(sim, first, first+1+k, 1);
    }
    sim_build_angles(sim);
    return first;
}

/* ── CO₂ ────────────────────────────────────────────────────────────────── */
/*
 * Linear: r(C=O) = 1.163 Å
 * q(C) = +0.70e, q(O) = -0.35e  (OPLS)
 */
int sim_place_co2(Simulation *sim, Vec3 origin) {
    double rCO = 1.163;
    int first = sim->num_atoms;
    sim_add_atom(sim, 8, vec3(origin.x - rCO, origin.y, origin.z), -0.35);
    sim_add_atom(sim, 6, vec3(origin.x,        origin.y, origin.z), +0.70);
    sim_add_atom(sim, 8, vec3(origin.x + rCO, origin.y, origin.z), -0.35);

    sim_add_bond(sim, first,   first+1, 2);  /* O=C */
    sim_add_bond(sim, first+1, first+2, 2);  /* C=O */
    /* Geometric builder, NOT sim_build_angles: the element-keyed angle
     * table's O-C-O entry is parameterised for BENT carboxylate groups
     * (theta0=123 deg, sp2 carbon) and would incorrectly be applied here
     * too, since the table has no way to distinguish that context from
     * linear CO2 (sp carbon, 180 deg) - both share the same (O,C,O)
     * element-triple key. Deriving theta0 directly from this molecule's
     * own correctly-seeded linear geometry avoids that collision. */
    sim_build_angles_geometric(sim, 3.5);  /* generic CO2 bend stiffness */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Diagnostics
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_print_atoms(const Simulation *sim) {
    printf("  %-4s %-4s  %-22s  %-22s  %-8s  %-8s\n",
           "idx","sym","Position (Å)","Velocity (Å/fs)","q(e)","mass(AMU)");
    printf("  %s\n",
        "──────────────────────────────────────────────────────────────"
        "───────────────────────────────────");
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        printf("  %-4d %-4s  (%7.4f %7.4f %7.4f)  "
               "(%7.4f %7.4f %7.4f)  %+7.4f  %7.3f\n",
               i, a->element->symbol,
               a->position.x, a->position.y, a->position.z,
               a->velocity.x, a->velocity.y, a->velocity.z,
               a->partial_charge, a->mass);
    }
}

void sim_print_bonds(const Simulation *sim) {
    printf("  %-6s %-4s %-4s %-5s %-8s %-8s %-10s\n",
           "bond","a","b","order","r0(Å)","r(Å)","E(eV)");
    printf("  %s\n","─────────────────────────────────────────────────");
    for (int i = 0; i < sim->num_bonds; i++) {
        const Bond *b = &sim->bonds[i];
        double r = vec3_dist(sim->atoms[b->atom_a].position,
                              sim->atoms[b->atom_b].position);
        double e = 0.5 * b->k * (r - b->r0) * (r - b->r0);
        printf("  %-6d %-4d %-4d %-5d %-8.4f %-8.4f %-10.6f\n",
               i, b->atom_a, b->atom_b, b->order, b->r0, r, e);
    }
}

void sim_print_angles(const Simulation *sim) {
    printf("  %-6s %-4s %-4s %-4s %-10s %-10s\n",
           "angle","a","b","c","θ0(deg)","k(eV/rad²)");
    printf("  %s\n","────────────────────────────────────────────────");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        printf("  %-6d %-4d %-4d %-4d %-10.2f %-10.4f\n",
               i, a->atom_a, a->atom_b, a->atom_c,
               a->theta0 * 180.0 / 3.14159265358979323846,
               a->k);
    }
}

void sim_print_summary(const Simulation *sim) {
    printf("\n════════════════════════════════════════════════════════\n");
    printf("  Simulation summary\n");
    printf("  Atoms: %d  Bonds: %d  Angles: %d\n",
           sim->num_atoms, sim->num_bonds, sim->num_angles);
    printf("  Step: %llu   Time: %.3f fs   dt: %.3f fs\n",
           (unsigned long long)sim->step, sim->time, sim->dt);
    printf("  KE: %.6f eV   PE: %.6f eV   E_total: %.6f eV\n",
           sim->kinetic_energy, sim->potential_energy, sim->total_energy);
    printf("  Temperature: %.2f K\n", sim->temperature);
    printf("════════════════════════════════════════════════════════\n");
    sim_print_atoms(sim);
    printf("\n");
    sim_print_bonds(sim);
    printf("\n");
    sim_print_angles(sim);
    printf("════════════════════════════════════════════════════════\n\n");
}

```

---

