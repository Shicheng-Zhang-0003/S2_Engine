#ifndef NEURON_H
#define NEURON_H

/*
 * neuron.h
 * Hodgkin-Huxley (1952) point-neuron model of the squid giant axon.
 *
 * This is a genuinely independent track from the chemistry/MD code in
 * the rest of this codebase - it does not depend on atoms, bonds, or
 * force fields. It models the SAME underlying physical idea one level
 * up: instead of simulating every ion explicitly, it uses the
 * empirically-measured, voltage-dependent conductance model that
 * Hodgkin and Huxley derived from real voltage-clamp experiments on
 * the squid giant axon (Nobel Prize, 1963).
 *
 * PARAMETER SOURCING (verified against 2+ independent, mutually
 * consistent literature sources before implementation, not from
 * memory alone):
 *   - Conductances/capacitance/reversal potentials: cross-checked
 *     against Ostojic (arXiv:0801.3963) and Kirst/Timme/Battaglia
 *     (arXiv:1506.05828), both giving identical values.
 *   - Gating rate functions (alpha/beta for m,h,n): cross-checked
 *     against Tulane University's HH course reference and the same
 *     Kirst/Timme/Battaglia review, both giving algebraically
 *     identical formulas (confirmed by direct substitution, e.g.
 *     "-0.0556" in one source equals "1/18" in the other).
 *   - Resting steady-state gating values (m0,h0,n0): independently
 *     recomputed from the rate functions themselves via the
 *     equilibrium condition m_inf = alpha_m/(alpha_m+beta_m) at
 *     V=-65 mV, and cross-checked against a published worked example
 *     (Proteus Analytics HH walkthrough) - matched to 4 decimal
 *     places (0.0529, 0.5961, 0.3177), confirming the sourced
 *     formulas are correct before any C implementation was written.
 *
 * Convention: "modern" voltage convention (V = actual membrane
 * potential in mV, resting near -65 mV, depolarization is POSITIVE),
 * not the original 1952 paper's own internal convention (which used
 * V=0 at rest with depolarization negative) - these are mathematically
 * equivalent under V_1952 = -(V_modern + 65), but mixing formulas from
 * sources using different conventions without converting would be a
 * real, silent, order-of-magnitude-scale bug. Every source used here
 * was confirmed to already use the modern convention before use.
 *
 * Units: V in mV, t in ms, C in uF/cm^2, g in mS/cm^2, currents in
 * uA/cm^2 (all standard HH-model units, consistent throughout).
 */

/* ── Fixed model parameters (squid giant axon, 6.3 C) ────────────────────── */
#define HH_C_M       1.0     /* membrane capacitance, uF/cm^2            */
#define HH_G_NA      120.0   /* max Na+ conductance, mS/cm^2             */
#define HH_G_K       36.0    /* max K+ conductance, mS/cm^2              */
#define HH_G_L       0.3     /* leak conductance, mS/cm^2                */
#define HH_E_NA      50.0    /* Na+ reversal potential, mV               */
#define HH_E_K       -77.0   /* K+ reversal potential, mV                */
#define HH_E_L       -54.4   /* leak reversal potential, mV              */
#define HH_V_REST    -65.0   /* resting potential, mV (initial V)        */

/* ── Gating rate functions (modern convention, V in mV) ──────────────────── */
double hh_alpha_m(double V);
double hh_beta_m (double V);
double hh_alpha_h(double V);
double hh_beta_h (double V);
double hh_alpha_n(double V);
double hh_beta_n (double V);

/* ── Single HH neuron state ───────────────────────────────────────────────── */
typedef struct {
    double V;        /* membrane potential, mV                          */
    double m, h, n;   /* gating variables (dimensionless, 0-1)           */
    double t;         /* simulation time, ms                            */
    double I_ext;     /* externally applied current, uA/cm^2 (settable) */

    /* Diagnostics, updated each step */
    double I_Na, I_K, I_L; /* instantaneous ionic currents, uA/cm^2      */
} HHNeuron;

/* Initialise a neuron at its analytically-correct resting steady state
 * (V=-65mV, m/h/n at their equilibrium values for that voltage - NOT
 * hardcoded copies of the reference numbers, but computed live from
 * the same rate functions used throughout, so the neuron starts in a
 * genuine, self-consistent equilibrium rather than an approximation). */
void hh_init(HHNeuron *n);

/* Advance the neuron by one Runge-Kutta 4th-order step of size dt (ms).
 * RK4 is used rather than simple Euler because the HH equations are
 * numerically stiff (the m gating variable can change on a much
 * faster timescale than h or n during the sodium spike upstroke) -
 * Euler would require an impractically small dt to remain stable and
 * accurate; RK4's higher-order error term tolerates dt=0.01ms
 * comfortably, which is standard practice for this model. */
void hh_step(HHNeuron *n, double dt);

/* Detect whether a spike (action potential threshold crossing) is in
 * progress: returns 1 if V > threshold_mV, 0 otherwise. A simple
 * level-crossing check, not a sophisticated spike-detection algorithm -
 * sufficient for identifying whether the model produced real spikes. */
int hh_is_spiking(const HHNeuron *n, double threshold_mV);

#endif /* NEURON_H */
