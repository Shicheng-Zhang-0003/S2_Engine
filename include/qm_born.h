#ifndef QM_BORN_H
#define QM_BORN_H

#include "types.h"

/* ══════════════════════════════════════════════════════════════════════════
 * qm_born — the ion's own electrostatic self-energy, and the ion's own
 * electronic polarizability.
 * ══════════════════════════════════════════════════════════════════════════
 *
 * WHY THIS EXISTS AT ALL
 *
 * The KcsA selectivity question is not a question about electrostatics. It
 * is a question about how much energy it costs to take a bare ion out of
 * water, and a vacuum force field has no water. Until this module existed,
 * the repository's only handle on that quantity was a transcribed table
 * (kcsa_hydration_free_energy_kJmol), which is right to use and wrong to
 * call physics: nothing in the code could check it, vary it, or explain
 * why its K/Na difference is what it is.
 *
 * This module supplies the missing internal model:
 *
 *   1. A Born ionic self-energy — the continuum limit of the ion's own
 *      field energy, -z^2 e^2 (1 - 1/eps) / (8 pi eps0 r_Born), with the
 *      effective radius from Latimer, Pitzer & Slansky (1936). This is
 *      what makes the desolvation leg derivable rather than transcribed,
 *      and it makes the ion radius enter the thermodynamics on its own.
 *
 *   2. Free-ion electronic polarizabilities (Pauling's table) for the
 *      induction term. The engine already has an atom polarizability
 *      model (qm_alpha), but it is a NEUTRAL-ATOM estimate used for
 *      dispersion; a free cation's polarizability is a different physical
 *      quantity with a different value (K+ 0.83 A^3 against neutral K
 *      2.73 A^3), and the K/Na ratio is what the induction leg turns on.
 *
 * WHAT THIS MODULE IS NOT
 *
 *   - It is not a force field and adds no energy to any trajectory. It is
 *     a set of closed-form reference functions plus tabulated constants.
 *   - The Born self-energy is a linear-response continuum result. It omits
 *     the cavity/dispersion contribution that the real ion–water
 *     interaction has, and it therefore UNDER-predicts the K/Na hydration
 *     differential. That deficit is quantified by
 *     qm_born_hydration_differential_kJmol() against the tabulated value
 *     rather than hidden: see the module's own caveat block and the
 *     Born-vs-tabulated row the KcsA leg prints.
 *   - Nothing here is fitted to a selectivity result. Every constant is a
 *     published number with a citation in qm_born.c.
 *
 * SCOPE OF THE BORRADII: there is no such thing as a unique Born radius.
 * The values below are the classic Latimer–Pitzer–Slansky set, which is
 * the one most commonly tabulated and the one Hummer, Pratt & Garcia
 * (arXiv:physics/9709042) explicitly confirm for Na+ (1.80 A there, 1.88 A
 * in Marcus's later compilation). A fitted alternative set exists
 * (Valisko: Na 1.62, K 1.95) and gives a different differential; the
 * spread between them is reported by
 * qm_born_hydration_differential_spread_kJmol() so the sensitivity of the
 * hydration leg to this single unsolvable choice is a number, not a
 * footnote.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * N_A e^2 / (8 pi eps0), the coefficient the Born expression divides by
 * the Born radius. kJ/mol per Angstrom for a univalent ion, derived from
 * CODATA primaries already in constants.h (N_A * e^2/(8 pi eps0) with the
 * Angstrom<->metre conversion) rather than transcribed, so it cannot drift
 * from the Coulomb constant the engine's electrostatics already use.
 */
double qm_born_coefficient_kJmol_per_ang(void);

/*
 * Effective Born radius in Angstrom — Latimer, Pitzer & Slansky (1936),
 * the standard tabulation:
 *
 *     Li+ 1.28    Na+ 1.80    K+ 2.08    Rb+ 2.35    Cs+ 2.65
 *
 * These are NOT the Shannon crystal radii the LJ sites use
 * (kcsa_cation_radius: K 1.51, Na 1.18 for CN=8) and the two must not be
 * confused: the crystal radius describes the ion in a lattice, the Born
 * radius describes the effective cavity that reproduces the measured
 * hydration free energy. Confirmed for Na+ as 0.180 nm by Hummer, Pratt &
 * Garcia, arXiv:physics/9709042, who recover it from explicit-solvent
 * finite-size analysis of the Hummer et al. single-ion energies.
 * Returns 0.0 for an ion with no tabulated value.
 */
double qm_born_radius(int Z);

/* Valisko's alternative Born radii obtained by fitting measured hydration
 * free energies (reported in Antonijevic, Electrochem. Sci. 15 (2010)
 * 1511, Table 1): Na+ 1.62, K+ 1.95. Kept as the sensitivity bound, not as
 * a second opinion to average with. Returns 0.0 when unlisted. */
double qm_born_radius_fitted(int Z);

/*
 * Born ionic self-energy (solvation free energy) of a charged sphere in a
 * dielectric continuum of relative permittivity eps_r:
 *
 *     dG = -N_A z^2 e^2 / (8 pi eps0 r_Born) * (1 - 1/eps_r)
 *
 * in kJ/mol, negative for eps_r > 1. Returns 0.0 for an unlisted ion or a
 * non-physical eps_r (<= 1 is vacuum, where the expression is identically
 * zero, so it is returned as 0.0 rather than rejected).
 *
 * Sanity value: Na+ in water (eps_r = 78.5, r_Born = 1.80 A) gives
 * -381 kJ/mol against a measured -365.3 kJ/mol — the linear-response
 * continuum over-binding the larger hydration terms fix.
 */
double qm_born_self_energy_kJmol(int Z, double eps_r);

/*
 * Born self-energy TRANSFER between two dielectric media, i.e. the cost of
 * moving an ion from eps_from to eps_to:
 *
 *     dG_transfer = dG(eps_to) - dG(eps_from)
 *
 * Positive to bury a cation in a low-dielectric medium. With the standard
 * protein-interior eps = 4 and water eps = 78.5 this is the term that says
 * why a bare charge prefers the pore over the protein matrix.
 */
double qm_born_transfer_kJmol(int Z, double eps_from, double eps_to);

/*
 * Differential Born self-energy dG(Za) - dG(Zb) for two ions in the same
 * medium — the shape of the Born hydration leg as a function of the two
 * Born radii alone. For K+/Na+ in water this is about -51 kJ/mol, against
 * a measured -70 kJ/mol on the TATB absolute scale: the missing 19 kJ/mol
 * is the cavity-plus-dispersion term the Born expression has no room for.
 * That gap is a property of the model, not an error to be tuned away.
 */
double qm_born_hydration_differential_kJmol(int Za, int Zb, double eps_r);

/*
 * The same differential evaluated on the fitted (Valisko) Born radii,
 * subtracted from the Latimer–Pitzer–Slansky value: a signed sensitivity
 * bar on the hydration leg attributable solely to the Born-radius choice.
 */
double qm_born_hydration_differential_spread_kJmol(int Za, int Zb,
                                                   double eps_r);

/*
 * Free-ion electronic polarizability in Angstrom^3 — Pauling's table,
 * which is the standard compilation for ISOLATED ions (as opposed to the
 * in-crystal values of Fowler & Madden or the fitted solution values of
 * Li et al., J. Phys. Chem. B 121 (2017) 6416, which differ by tens of
 * percent and are a different question):
 *
 *     Li+ 0.03   Na+ 0.18   K+ 0.83   Rb+ 1.40   Cs+ 2.40
 *
 * The K+/Na+ ratio 4.6 is the physical handle the induction leg uses:
 * the same field polarizes K+ 4.6x harder than Na+. Returns 0.0 when
 * unlisted.
 */
double qm_ion_alpha_free(int Z);

/*
 * Ionic polarization (induced-dipole) energy -1/2 alpha E^2 in eV, from a
 * free-ion polarizability in A^3 and an electric field magnitude in e/A
 * (the engine's native Coulomb units: 14.3996 e/A is the field of a unit
 * charge at 1 A). The 1/2 is the self-consistency factor — the dipole is
 * built by the field it responds to, so the work done is half of alpha E^2.
 *
 * `eps_r` SCREENS THE FIELD and is not optional. Alpha is tabulated for a
 * free ion in vacuum, so the field it responds to is the macroscopic one,
 * E_macro = E_bare / eps_r. Leaving it out is not a small approximation
 * error: at the ~20 e/A a cation sees inside an eight-carbonyl cage, the
 * unscreened induction energy is a couple of hundred eV, which is four
 * orders of magnitude larger than any physically meaningful ion-ligand
 * binding energy and would silently dominate every total it was added to.
 * Use eps_r = 4 for the protein interior (Warshel et al., PNAS 105
 * (2008) 11845, and Gong, Hocky & Freed there) and eps_r = 1 only for a
 * genuinely gas-phase ion.
 *
 * This is the isolated-ion half; the ligand response is the induction
 * solver's job (qm_induction_energy / qm_solve_dipoles), which divides by
 * the same dielectric internally.
 */
double qm_born_induction_eV(double alpha_A3, double field_e_per_A,
                            double eps_r);

#ifdef __cplusplus
}
#endif

#endif /* QM_BORN_H */