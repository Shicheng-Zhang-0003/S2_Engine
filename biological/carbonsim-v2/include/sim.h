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
