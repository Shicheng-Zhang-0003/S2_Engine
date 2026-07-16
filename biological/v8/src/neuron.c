#include <math.h>
#include "../include/neuron.h"

/*
 * neuron.c
 * See neuron.h for full parameter/formula sourcing and verification.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Gating rate functions (modern convention, V in mV)
 *
 * Each has a removable singularity at its own specific voltage (where
 * the denominator (1 - exp(...)) -> 0): alpha_m at V=-40, alpha_n at
 * V=-55. L'Hopital's rule gives the correct limiting value there;
 * implemented via a small-argument Taylor expansion of the
 * x/(1-exp(-x)) form (which -> 1 as x->0) to avoid a literal 0/0
 * division at those exact voltages.
 * ══════════════════════════════════════════════════════════════════════════ */

/* Numerically stable x / (1 - exp(-x)), correct in the x->0 limit
 * (where the true value is 1, by L'Hopital / Taylor expansion:
 * 1 - exp(-x) ~ x - x^2/2 for small x, so x/(1-exp(-x)) ~ 1/(1-x/2) ~ 1). */
static double stable_ratio(double x) {
    if (fabs(x) < 1.0e-6) {
        /* Taylor expansion near the singularity: limit is 1, with a
         * small linear correction for the next-order term. */
        return 1.0 - x / 2.0;
    }
    return x / (1.0 - exp(-x));
}

double hh_alpha_m(double V) {
    /* True form: 0.1*(V+40) / (1 - exp(-(V+40)/10)), with a removable
     * singularity at V=-40. Rewritten as 0.1 * 10 * stable_ratio(x)
     * where x=(V+40)/10, since (V+40) = 10x and
     * (1-exp(-(V+40)/10)) = (1-exp(-x)), so the ratio (V+40)/(1-exp(-x))
     * = 10x/(1-exp(-x)) = 10*stable_ratio(x). */
    double x = (V + 40.0) / 10.0;
    return 0.1 * 10.0 * stable_ratio(x);
}

double hh_beta_m(double V) {
    return 4.0 * exp(-(V + 65.0) / 18.0);
}

double hh_alpha_h(double V) {
    return 0.07 * exp(-(V + 65.0) / 20.0);
}

double hh_beta_h(double V) {
    return 1.0 / (1.0 + exp(-(V + 35.0) / 10.0));
}

double hh_alpha_n(double V) {
    /* Same pattern as alpha_m: true form 0.01*(V+55)/(1-exp(-(V+55)/10)),
     * removable singularity at V=-55, rewritten via stable_ratio. */
    double x = (V + 55.0) / 10.0;
    return 0.01 * 10.0 * stable_ratio(x);
}

double hh_beta_n(double V) {
    return 0.125 * exp(-(V + 65.0) / 80.0);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Derivative bundle for the 4-variable HH system, used by the RK4
 * integrator. Autonomous system (no explicit t dependence beyond the
 * fixed I_ext parameter), so this only needs the current state.
 * ══════════════════════════════════════════════════════════════════════════ */
typedef struct { double dV, dm, dh, dn; } HHDeriv;

static HHDeriv hh_derivatives(double V, double m, double h, double n,
                               double I_ext) {
    HHDeriv d;

    double I_Na = HH_G_NA * m*m*m * h * (V - HH_E_NA);
    double I_K  = HH_G_K  * n*n*n*n * (V - HH_E_K);
    double I_L  = HH_G_L  * (V - HH_E_L);

    d.dV = (I_ext - I_Na - I_K - I_L) / HH_C_M;
    d.dm = hh_alpha_m(V) * (1.0 - m) - hh_beta_m(V) * m;
    d.dh = hh_alpha_h(V) * (1.0 - h) - hh_beta_h(V) * h;
    d.dn = hh_alpha_n(V) * (1.0 - n) - hh_beta_n(V) * n;

    return d;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Initialise at the analytically correct resting equilibrium
 * ══════════════════════════════════════════════════════════════════════════ */
void hh_init(HHNeuron *neuron) {
    neuron->V = HH_V_REST;
    neuron->t = 0.0;
    neuron->I_ext = 0.0;

    double am = hh_alpha_m(HH_V_REST), bm = hh_beta_m(HH_V_REST);
    double ah = hh_alpha_h(HH_V_REST), bh = hh_beta_h(HH_V_REST);
    double an = hh_alpha_n(HH_V_REST), bn = hh_beta_n(HH_V_REST);

    neuron->m = am / (am + bm);
    neuron->h = ah / (ah + bh);
    neuron->n = an / (an + bn);

    neuron->I_Na = HH_G_NA * neuron->m*neuron->m*neuron->m * neuron->h
                   * (neuron->V - HH_E_NA);
    neuron->I_K  = HH_G_K * neuron->n*neuron->n*neuron->n*neuron->n
                   * (neuron->V - HH_E_K);
    neuron->I_L  = HH_G_L * (neuron->V - HH_E_L);
}

/* ══════════════════════════════════════════════════════════════════════════
 * RK4 integration step
 * ══════════════════════════════════════════════════════════════════════════ */
void hh_step(HHNeuron *neuron, double dt) {
    double V = neuron->V, m = neuron->m, h = neuron->h, n = neuron->n;
    double I = neuron->I_ext;

    HHDeriv k1 = hh_derivatives(V, m, h, n, I);

    HHDeriv k2 = hh_derivatives(V + 0.5*dt*k1.dV, m + 0.5*dt*k1.dm,
                                 h + 0.5*dt*k1.dh, n + 0.5*dt*k1.dn, I);

    HHDeriv k3 = hh_derivatives(V + 0.5*dt*k2.dV, m + 0.5*dt*k2.dm,
                                 h + 0.5*dt*k2.dh, n + 0.5*dt*k2.dn, I);

    HHDeriv k4 = hh_derivatives(V + dt*k3.dV, m + dt*k3.dm,
                                 h + dt*k3.dh, n + dt*k3.dn, I);

    neuron->V += (dt/6.0) * (k1.dV + 2.0*k2.dV + 2.0*k3.dV + k4.dV);
    neuron->m += (dt/6.0) * (k1.dm + 2.0*k2.dm + 2.0*k3.dm + k4.dm);
    neuron->h += (dt/6.0) * (k1.dh + 2.0*k2.dh + 2.0*k3.dh + k4.dh);
    neuron->n += (dt/6.0) * (k1.dn + 2.0*k2.dn + 2.0*k3.dn + k4.dn);

    neuron->t += dt;

    /* Update diagnostic currents at the new state */
    neuron->I_Na = HH_G_NA * neuron->m*neuron->m*neuron->m * neuron->h
                   * (neuron->V - HH_E_NA);
    neuron->I_K  = HH_G_K * neuron->n*neuron->n*neuron->n*neuron->n
                   * (neuron->V - HH_E_K);
    neuron->I_L  = HH_G_L * (neuron->V - HH_E_L);
}

int hh_is_spiking(const HHNeuron *neuron, double threshold_mV) {
    return neuron->V > threshold_mV;
}
