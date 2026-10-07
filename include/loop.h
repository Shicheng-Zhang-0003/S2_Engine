#ifndef LOOP_H
#define LOOP_H

/*
 * loop.h — bio/QC/QM loop-closure bridge.
 *
 * The three tracks (nucleic acids, proteins, electrophysiology) run
 * standalone in the dynamics: gating is not derived from protein
 * structure and the filter is rigid. That is a scope separation, not a
 * physics inconsistency — all three consume the SAME primaries
 * (CODATA h/e/k/NA, COULOMB_MD, KCAL_MOL_TO_EV, Shannon VIII radii,
 * Marcus 1991 TATB free energies) through the SAME derive-in-line
 * rule as constants.h.
 *
 * This module makes that claim checkable: every scale-bridging number
 * (Nernst reversal, thermal voltage, 1000:1 selectivity scale, LJ
 * minima, Coulomb prefactor) is recomputed here from primaries, and
 * the irreducible reduced-model choices are enumerated so no
 * abstraction stays hidden.
 *
 * Record impact: none. Pure functions over constants; not called from
 * the Demo record path. Provenance for the values below lives in the
 * cited literature, retrieval notes in tests/test_loop.c.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* Thermal voltage RT/F in mV at T_K (R=N_A*k, F=N_A*e). NaN on bad T. */
double loop_thermal_mV(double T_K);

/* Thermal energy kT in eV at T_K. NaN on bad T. */
double loop_thermal_eV(double T_K);

/* Nernst reversal in mV: (RT/zF) ln(c_out/c_in). NaN on bad input. */
double loop_nernst_mV(double z, double c_out, double c_in, double T_K);

/* Selectivity free-energy scale in eV: -kT ln(ratio). NaN on bad input.
 * ratio=1000 at 300K must give ~-0.179 eV (the Demo 12 expt reference). */
double loop_selectivity_scale_eV(double ratio, double T_K);

/* LJ 12-6 minimum position from collision diameter: 2^(1/6)*sigma. */
double loop_lj_minimum(double sigma);

/* Irreducible reduced-model choices. Count + names are asserted in
 * tests/test_loop.c so documentation cannot drift from code. */
int loop_abstraction_count(void);
const char *loop_abstraction_name(int i);

#ifdef __cplusplus
}
#endif

#endif /* LOOP_H */
