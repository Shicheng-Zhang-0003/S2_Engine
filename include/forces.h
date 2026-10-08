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
 * Coulomb 1/r is long-ranged: a hard cutoff truncates its tail and
 * creates a discontinuity (energies are cutoff- and size-dependent).
 * Acceptable here only because all current demos are small gas-phase
 * clusters whose pair distances largely fall inside 12 A; do not cite
 * magnitudes as converged electrostatics. Use Ewald/PME or
 * reaction-field before any condensed-phase quantitative claim.
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
/* AMBER combines R_min-half arithmetically and eps geometrically; with
 * sigma = 2 Rstar over 2^(1/6) this is exactly sigma-arithmetic plus
 * eps-geometric below. Requires correct (halved-otherwise) sigma. */
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
PairEnergy forces_nonbonded_pair(Atom *atoms, int num_atoms, int ia, int ib,
                                   const SimBox *box,
                                   int use_lj, int use_coulomb,
                                   double dielectric);

/* Side-effect-free pair energy (no force accumulation).
 * Use for diagnostics that must not perturb dynamics state.
 * Same potential, same cutoff handling via box; forces untouched. */
PairEnergy forces_nonbonded_energy(const Atom *atoms, int num_atoms, int ia, int ib,
                                   const SimBox *box,
                                   int use_lj, int use_coulomb,
                                   double dielectric);

/* ── Harmonic bond force (accumulates onto both endpoint atoms) ───────────── */
/*
 * V = 0.5 k (r - r0)²   → F_a = k(r-r0) r̂_ab,  F_b = −F_a
 * Returns bond potential energy in eV.
 */
double forces_bond(Atom *atoms, int num_atoms, const Bond *bond);

/* ── Harmonic angle force (accumulates onto all three atoms) ─────────────── */
/*
 * V = 0.5 k (θ − θ0)²
 * Gradient computed analytically via the chain rule through acos.
 * Returns angle potential energy in eV.
 */
double forces_angle(Atom *atoms, int num_atoms, const Angle *angle);

/* ── Dihedral (torsion) force (accumulates onto all four atoms) ──────────── */
/*
 * V = k * (1 + cos(n*phi - delta)), the standard AMBER/CHARMM torsion
 * functional form - general enough to express any periodic torsion
 * preference (n=1 for a single minimum, n=2/3 for the multiple minima
 * typical of real sp3-sp3 or partial-double-bond rotations).
 *
 * The FORCE is the exact analytic chain-rule gradient through
 * phi = atan2(y, x) (derivation in-source at forces_dihedral; the
 * translation-invariance identity dx_a+dx_b+dx_c+dx_d = 0 holds
 * exactly by construction of the formulas). Audit P2 closed: the
 * former central-finite-difference default is kept as
 * forces_dihedral_fd() (the validation oracle) and
 * tests/test_forces.c asserts analytic-vs-FD agreement to 1e-6 plus
 * net-zero force on generic, helical, and near-planar geometries.
 * Collinear plane-normals (|n| < 1e-12) return energy with zero
 * forces rather than FD's garbage-scale values.
 *
 * Returns the dihedral potential energy in eV.
 */
double forces_dihedral(Atom *atoms, int num_atoms, const Dihedral *dihedral);

/* Finite-difference torsion force (former default, validation oracle).
 * Same energy; forces via central differences (h=1e-5 A, 24 evals).
 * See tests/test_forces.c for the agreement contract. */
double forces_dihedral_fd(Atom *atoms, int num_atoms, const Dihedral *dihedral);

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
