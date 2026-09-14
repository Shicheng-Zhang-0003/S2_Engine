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
 * where MD_FORCE_CONV is derived in-line in constants.h from
 * EV_TO_J/ANGSTROM_TO_M/AMU (≈9.64853322e-3).
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
 * Units: AMU × (Å/fs)² → converted to eV via AMU_AFS2_TO_EV.
 *
 * Exact: 1 AMU × (Å/fs)² = 1.66053906660e-27 kg × (1e-10/1e-15)² m²/s²
 *      = 1.66053906660e-27 × 1e10 J = 1.66053906660e-17 J
 *      = 1.66053906660e-17 / 1.602176634e-19 eV = 103.6427 eV
 * So: KE [eV] = 0.5 Σ m[AMU] v²[Å²/fs²] × 103.6427
 * (the 0.5 lives in the summation loop in integrator.c, NOT in the
 * conversion factor - AMU_AFS2_TO_EV is the full single-unit value.)
 */
#define AMU_AFS2_TO_EV   (AMU * 1.0e10 / EV_TO_J)   /* (Å/fs)² AMU → eV, derived in-line */

double integrator_kinetic_energy(const Simulation *sim);

/*
 * Temperature from equipartition theorem: KE = (dof/2) k_B T
 * T [K] = 2 KE [eV] / (dof x k_B_eV)
 * k_B (eV/K) = BOLTZMANN_K/EV_TO_J (derived in-line; ≈8.617333262e-5)
 *
 * dof = 3N - 3 - constrained: the -3 removes the 3 centre-of-mass
 * translational degrees of freedom, already constrained to zero by
 * integrator_remove_com_velocity(). `constrained` counts FURTHER
 * removed degrees of freedom the caller knows about:
 *   - 3 per frozen/immobilised atom (see integrator_minimize_frozen);
 *   - 2 more for a strictly linear molecule (rotation about its own
 *     axis carries no energy: 3N-5, not 3N-3).
 * Harmonic positional restraints do NOT remove degrees of freedom (the
 * atom still moves in 3D about its anchor); leave constrained=0 for
 * restrained systems. sim->num_constrained_dof carries that count
 * (0 by default, i.e. the long-standing 3N-3 behaviour). Getting this
 * wrong does not crash - it silently biases every reported temperature,
 * so count honestly.
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

/*
 * Steepest-descent energy minimization: repeatedly moves every atom a
 * small distance along its own force vector (NOT full Newtonian
 * dynamics), shrinking the step size whenever a step would increase
 * the energy, until forces are small or max_iterations is reached.
 * Standard MD practice - resolves initial-condition clashes (like a
 * freshly-built chain with unavoidable local steric overlaps) before
 * running velocity-based dynamics, which has no safe way to absorb a
 * severe initial clash. Returns the final potential energy.
 */
double integrator_minimize(Simulation *sim, int max_iterations,
                            double initial_step, double force_tolerance);

/*
 * Same as integrator_minimize, but atoms with frozen[i] != 0 are
 * excluded from all position updates. `frozen` must be an array of
 * length sim->num_atoms.
 */
double integrator_minimize_frozen(Simulation *sim, const int *frozen,
                                   int max_iterations, double initial_step,
                                   double force_tolerance);

#endif /* INTEGRATOR_H */
