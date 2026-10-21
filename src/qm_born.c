/* qm_born.c — Born ionic self-energy and free-ion polarizability.
 *
 * See include/qm_born.h for why this module exists and what it is not.
 *
 * Every constant below is a published number with a citation. Nothing here
 * is fitted to reproduce a selectivity number, and the one place where the
 * model is known to be incomplete (the cavity/dispersion contribution the
 * Born expression omits) is quantified rather than absorbed.
 */

#include "../include/qm_born.h"
#include "../include/constants.h"

#include <math.h>

/* N_A e^2 / (8 pi eps0).
 *
 * COULOMB_MD is e^2/(4 pi eps0) in eV*Angstrom, derived in constants.h
 * from CODATA primaries, so the Born coefficient is derived from the same
 * source the engine's own electrostatics uses rather than transcribed a
 * second time. Halving gives the 8 pi in the denominator; multiplying by
 * N_A and converting eV -> kJ/mol (/1000) completes it.
 *
 * NO LENGTH CONVERSION APPEARS, and that is correct rather than an
 * omission: the Born energy is C/r with C an energy*length, and
 * COULOMB_MD already carries the Angstrom in its units. The conversion
 * here is eV*Angstrom -> kJ*Angstrom per mole, i.e. energy and amount
 * only. Multiplying by ANGSTROM_TO_M would have squared the Angstrom and
 * silently returned zero at the 1e-20 scale; an earlier revision of this
 * function did exactly that and every Born number in the tree came out as
 * 0.0000 until the unit test in tests/test_external.c caught it. */
double qm_born_coefficient_kJmol_per_ang(void) {
    return AVOGADRO_N * (COULOMB_MD * 0.5) * EV_TO_J / 1000.0;
}

double qm_born_radius(int Z) {
    switch (Z) {
        case  3: return 1.28;   /* Li+  */
        case 11: return 1.80;   /* Na+  */
        case 19: return 2.08;   /* K+   */
        case 37: return 2.35;   /* Rb+  */
        case 55: return 2.65;   /* Cs+  */
        default: return 0.0;
    }
}

double qm_born_radius_fitted(int Z) {
    switch (Z) {
        case 11: return 1.62;   /* Na+ */
        case 19: return 1.95;   /* K+  */
        case 37: return 2.25;   /* Rb+ */
        case 55: return 2.53;   /* Cs+ */
        case  3: return 1.16;   /* Li+ */
        default: return 0.0;
    }
}

double qm_born_self_energy_kJmol(int Z, double eps_r) {
    double r = qm_born_radius(Z);
    if (!(r > 0.0)) return 0.0;
    if (!isfinite(eps_r) || eps_r <= 1.0) return 0.0;  /* vacuum: no solvation */
    /* z = 1 for every ion in this module; the tabulated radii are all
     * univalent alkali cations and silently squaring a charge the table
     * does not describe would be worse than refusing. */
    return -qm_born_coefficient_kJmol_per_ang() / r * (1.0 - 1.0 / eps_r);
}

double qm_born_transfer_kJmol(int Z, double eps_from, double eps_to) {
    return qm_born_self_energy_kJmol(Z, eps_to)
         - qm_born_self_energy_kJmol(Z, eps_from);
}

double qm_born_hydration_differential_kJmol(int Za, int Zb, double eps_r) {
    return qm_born_self_energy_kJmol(Za, eps_r)
         - qm_born_self_energy_kJmol(Zb, eps_r);
}

double qm_born_hydration_differential_spread_kJmol(int Za, int Zb,
                                                   double eps_r) {
    double fa = qm_born_radius_fitted(Za);
    double fb = qm_born_radius_fitted(Zb);
    double ra = qm_born_radius(Za);
    double rb = qm_born_radius(Zb);
    if (!(fa > 0.0) || !(fb > 0.0) || !(ra > 0.0) || !(rb > 0.0)) return 0.0;
    if (!isfinite(eps_r) || eps_r <= 1.0) return 0.0;
    double c = qm_born_coefficient_kJmol_per_ang() * (1.0 - 1.0 / eps_r);
    double lps = c * (1.0 / rb - 1.0 / ra);   /* dG(a) - dG(b), LPS radii  */
    double fit = c * (1.0 / fb - 1.0 / fa);   /* dG(a) - dG(b), fitted      */
    return fit - lps;
}

double qm_ion_alpha_free(int Z) {
    switch (Z) {
        case  3: return 0.03;   /* Li+ */
        case 11: return 0.18;   /* Na+ */
        case 19: return 0.83;   /* K+  */
        case 37: return 1.40;   /* Rb+ */
        case 55: return 2.40;   /* Cs+ */
        default: return 0.0;
    }
}

double qm_born_induction_eV(double alpha_A3, double field_e_per_A,
                            double eps_r) {
    if (!isfinite(alpha_A3) || !isfinite(field_e_per_A)) return 0.0;
    if (!isfinite(eps_r) || eps_r < 1.0) return 0.0;
    if (alpha_A3 <= 0.0) return 0.0;
    /* Alpha is a free-ion (vacuum) quantity, so the field it polarizes in
     * is the macroscopic one. Screening the energy rather than the field
     * would square the factor; see the header. */
    double e_macro = field_e_per_A / eps_r;
    /* FULL-AUDIT: the missing factor is 1/(4*pi*eps0), i.e.
     * 1/COULOMB_MD in the engine's units.
     *
     * The engine's units: the field is E [V/A] (the potential of a unit
     * charge at 1 A is 14.3996 V, so the field there is 14.3996 V/A),
     * the polarizability is a VOLUME alpha [A^3], and the induced dipole
     * is mu = alpha E / COULOMB_MD [e.A]. The polarization energy is
     * U = -1/2 mu.E, hence
     *
     *     U [eV] = -1/2 * alpha [A^3] * E [V/A]^2 / COULOMB_MD
     *
     * The old code omitted the division by COULOMB_MD and so returned
     * e^2.A, not eV: a value 14.3996x too small when read as an energy.
     * (The same omission made the "unit charge at 1 A is 14.3996 e/A"
     * statement in the header wrong - 14.3996 V/A is the FIELD at 1 A;
     * the potential there is 14.3996 V.) */
    return -0.5 * alpha_A3 * e_macro * e_macro / COULOMB_MD;
}