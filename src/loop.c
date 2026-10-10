#include <math.h>
#include "../include/loop.h"
#include "../include/constants.h"
#include "../include/amber_lj.h"

/*
 * loop.c — single-source recomputation of every cross-track number.
 *
 * No magic numbers: R=N_A*k, F=N_A*e, kT and RT/F derived in-line.
 * Physiological squid concentrations below are the Hodgkin-Huxley
 * 1952 external/internal values (440/50 Na, 20/400 K mM at 6.3 C =
 * 279.45 K); they are inputs to the check, not outputs, and the
 * 2-5 mV residual against HH_E_NA/HH_E_K is activity + leak mixing,
 * documented in the test, not hidden.
 */

double loop_thermal_mV(double T_K) {
    if (!(T_K > 0.0) || !isfinite(T_K)) return NAN;
    /* RT/F = N_A*k*T / (N_A*e) = k*T/e, in V; x1000 for mV. */
    double rt_over_f = (BOLTZMANN_K * T_K) / ELEM_CHARGE;
    return rt_over_f * 1000.0;
}

double loop_thermal_eV(double T_K) {
    if (!(T_K > 0.0) || !isfinite(T_K)) return NAN;
    /* kT in J / (J per eV). EV_TO_J is exact by SI redefinition. */
    return (BOLTZMANN_K * T_K) / EV_TO_J;
}

double loop_nernst_mV(double z, double c_out, double c_in, double T_K) {
    if (!(z != 0.0) || !isfinite(z)) return NAN;
    if (!(c_out > 0.0) || !(c_in > 0.0)) return NAN;
    if (!isfinite(c_out) || !isfinite(c_in)) return NAN;
    if (!(T_K > 0.0) || !isfinite(T_K)) return NAN;
    double vt = loop_thermal_mV(T_K);
    if (!isfinite(vt)) return NAN;
    return (vt / z) * log(c_out / c_in);
}

double loop_selectivity_scale_eV(double ratio, double T_K) {
    if (!(ratio > 0.0) || !isfinite(ratio)) return NAN;
    if (!(T_K > 0.0) || !isfinite(T_K)) return NAN;
    double kt = loop_thermal_eV(T_K);
    if (!isfinite(kt)) return NAN;
    return -kt * log(ratio);
}

double loop_lj_minimum(double sigma) {
    if (!(sigma > 0.0) || !isfinite(sigma)) return NAN;
    return TWOPOW_SIXTH * sigma;
}

/* Irreducible choices — mirrors readme §7. Adding a new reduced-model
 * assumption without adding it here breaks test_loop, by design.
 * FULL-AUDIT Q1/Q2 (third pass): two assumptions the earlier inventory did
 * not name, both surfaced by the independent-oracle sweep. */
static const char *const LOOP_ABSTRACTIONS[] = {
    "no-bulk-solvent: vacuum + dielectric divisor, 6-water clusters only",
    "rigid-filter: deposited 1K4C TVGYG, no flexible-filter free energy",
    "tracks-bridged-at-constants: dynamics unconnected, primaries shared",
    "base-pair-magnitudes-qualitative: ordering validated, ~4x overshoot",
    "harmonic-bonds-cannot-break: Morse/reactive is future work",
    "non-nucleobase-charges-approximate: charge-balanced, not RESP",
    "qm-overlap-heuristic: Slater prefactor dropped, m maximised per atom",
    "no-1-4-scaling-by-design: AMBER scaling needs co-fitted torsions",
    "cutoff-switch-off-by-default: gas-phase correct, opt-in for condensed",
    "o-n2-no-neighbour-list: fine at current sizes, exclusion O(degree)",
    "no-membrane-potential-in-dynamics: Nernst supplied here from same kT",
    "qeq-unscreened-coulomb: bare 1/r, no Rappe screening; refused when the"
    " ordering inverts (Q1) — screened variant is future work",
    "minimum-image-one-convention: every pair term, periodic or not, uses the"
    " same separation (Q2); no long-range QEq correction",
    "diagnostic-velocity-cap-on-by-default: 2000 K rescales velocities, so a"
    " hand-seeded NVE run above it is NOT NVE (Q6); disable with maxtemp 0"
};

int loop_abstraction_count(void) {
    return (int)(sizeof(LOOP_ABSTRACTIONS) / sizeof(LOOP_ABSTRACTIONS[0]));
}

const char *loop_abstraction_name(int i) {
    if (i < 0 || i >= loop_abstraction_count()) return 0;
    return LOOP_ABSTRACTIONS[i];
}
