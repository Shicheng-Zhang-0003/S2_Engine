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
