#ifndef QM_FLUCTUATION_H
#define QM_FLUCTUATION_H

#include "types.h"

/* ══════════════════════════════════════════════════════════════════════════
 * qm_fluctuation — the entropic half of ion binding, and the strain that
 * produces it.
 * ══════════════════════════════════════════════════════════════════════════
 *
 * WHY THIS MODULE EXISTS
 *
 * Until now the repository's KcsA story was a mean-energy story: a filter
 * binding energy, a tabulated hydration free energy, and a sum. That sum
 * is not KcsA, and Dixit & Asthagiri (J. Chem. Theory Comput. 5 (2009)
 * 1037, doi:10.1021/jct800426g, with Asthagiri, Roux & Cross, Biophys. J.
 * 95 (2008) 1269, doi:10.1529/biophysj.108.138638) showed why: in the real
 * filter the MEAN binding energy slightly FAVOURS Na+, not K+. Na+ is the
 * better fit for an eight-carbonyl oxygen cage. Selectivity does not come
 * from the mean; it comes from the VARIANCE,
 *
 *     dG_bind = d(beta*<E_bind>)  +  d(beta*sigma^2 / 2)
 *
 * and the second term is large enough to overturn the first. The physical
 * origin they identified is strain: Na+ only reaches its favourable
 * contact distances (2.55 A, versus 2.93 A for K+) by pulling the filter
 * oxygens inward, which compresses the carbonyl-carbonyl distances and
 * costs repulsion strain that shows up as an increased fluctuation.
 *
 * A rigid-cage model cannot produce that term. It has no variance to speak
 * of, by construction. This module supplies it by running the filter as
 * what it is — a breathing peptide with the ion inside — and measuring
 * the variance of the ion's binding energy over an explicit NVT
 * trajectory.
 *
 * WHAT IS MEASURED HERE
 *
 *   E_bind(t) = E(totals, ion present) - E(totals, ion absent, same
 *               coordinates)
 *
 * The second term is computed by zeroing the ion's charge AND its
 * Lennard-Jones epsilon in place and re-evaluating at the same
 * coordinates. It is an exact counterfactual at that configuration, not an
 * average over a different trajectory, so the difference isolates the
 * ion's contribution to the energy of that conformation and nothing else.
 *
 *   O-O shell separation and its Lennard-Jones energy
 *
 * are the strain readout proper: the direct measure of what the ion did to
 * the carbonyl-carbonyl distances, which is the term the mean total energy
 * hides behind the tethers and the bonded terms. The literature value to
 * compare against is a Gly77 diagonal O-O held at roughly 4.5-5.5 A by the
 * DFT/MM models even with Na+ present.
 *
 *   E_filter(t) = E(totals, ion absent)
 *
 * is the same quantity on its own. Its mean IS the ligand-ligand strain
 * energy: comparing it between a run with K+ and a run with Na+ measures
 * directly how much the smaller cation compresses the carbonyl network.
 * That is the term Dixit et al. call dU_COCO.
 *
 * CAUTION ON WHAT THIS IS NOT
 *
 *   - A short trajectory of a tethered filter is not a converged ensemble.
 *     The standard error is reported with every variance and the caller is
 *     expected to print it. A variance whose error bar is comparable to the
 *     difference between the two ions is not a result.
 *   - The Langevin thermostat is exactly canonical but the trajectory is
 *     short; these are canonical-ensemble samples, not equilibrium
 *     guarantees. Independent re-runs with different seeds are the check,
 *     and the caller can run more.
 *   - The ion is bound at a FIXED axial site and the cage is tethered.
 *     This measures the site thermodynamics of the deposited geometry
 *     family, not the permeation path.
 *
 * DETERMINISM
 *
 *   Every number here is a function of the initial coordinates, the step
 *   counts and the seed alone. The engine's PCG64 stream is seeded through
 *   integrator_maxwell_boltzmann(), so two runs with the same seed replay
 *   bit-for-bit. Using the SAME seed for the K+ and Na+ runs is deliberate:
 *   it makes them common random numbers, so the difference of the two
 *   variances is far better determined than either variance alone.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Descriptive statistics for an energy series, with the error bars the
 * caller needs in order to decide whether a difference is a result.
 * Returns 0 on success, -1 on bad input.
 */
typedef struct {
    int    n;
    double mean_eV;
    double var_eV2;      /* unbiased sample variance, ddof = 1           */
    double sigma_eV;     /* sqrt(var)                                     */
    double sem_eV;       /* standard error of the mean                    */
    double sem_var_eV2;  /* standard error of the variance estimate       */
    double min_eV;
    double max_eV;
} QmSeries;

int qm_series_stats(const double *series_eV, int n, QmSeries *out);

/*
 * The fluctuation (entropic) contribution to a binding free energy,
 *
 *     d(beta sigma^2 / 2)  ->  half the difference of variances
 *
 * in eV, given the two variances in eV^2. Sign convention, used
 * consistently everywhere in the KcsA code: the first argument is the ion
 * whose preference the result expresses, so
 *
 *     qm_fluctuation_term_eV(sigma2_a, sigma2_b) < 0  means  A is favoured.
 *
 * For the K+/Na+ comparison this is 0.5 * (sigma^2_K - sigma^2_Na) in
 * eV, and a NEGATIVE value means K+ is favoured by the fluctuation term,
 * because Na+ compresses the filter harder and its binding energy swings
 * wider.
 */
/* Returns (var_a - var_b)/2, or NaN if either variance is not finite
 * (FULL-AUDIT: a non-finite input used to be silently converted to a
 * clean 0.0 eV, which is the one answer this term must never give by
 * accident). Callers MUST check isfinite() before using the result. */
double qm_fluctuation_term_eV(double var_a_eV2, double var_b_eV2);

/*
 * A sampled KcsA site: one ion in the real 1K4C filter under Langevin NVT.
 * All energies in eV, all distances in Angstrom.
 */
typedef struct {
    int    ok;
    int    n_samples;
    double mean_bind_eV;      /* <E_bind>  : the MEAN binding term          */
    double var_bind_eV2;      /* sigma^2   : the FLUCTUATION term's input   */
    double sem_bind_eV;
    double sem_var_bind_eV2;
    double mean_filter_eV;    /* <E_filter>, ion absent: strain reference   */
    double var_filter_eV2;
    double mean_iono_A;       /* <ion-O> to the coordinating oxygens       */
    double cn_mean;           /* mean coordination number at the site       */
    double mean_oo_A;         /* <O-O> across the coordinating shell        */
    double oo_lj_eV;          /* shell carbonyl-carbonyl LJ (dU_COCO)        */
    int    site_maintained;   /* 1 iff the filter held its 8-fold site      */
} QmKcsaSeries;

/*
 * The deposited 1K4C filter is 8-coordinate at every one of its four ion
 * sites; that is the statement kcsa_coord_stats verifies and the one the
 * literature's O-O diagonals (roughly 4.5-5.5 A) are measured against.
 * A run whose mean coordination has fallen below
 * QM_KCSA_SITE_MIN_CN was not a measurement of the site — it was a
 * measurement of a cage that came apart under its own thermal bath, and
 * its binding energies describe that.
 *
 * FULL-AUDIT: the gate used to be a flag with no consumer, so a run in
 * which the cage came apart still handed back a binding-energy variance
 * that a caller could quote as the fluctuation term of a held site. It
 * now ENFORCES: qm_kcsa_binding_series() returns -1 when the site was
 * not maintained, and writes no usable numbers (out->ok is cleared and
 * the statistics in `out` describe the collapse, so a caller can report
 * it).
 *
 * This fires readily, and that is measured rather than asserted: the
 * mean coordination falls from the deposited 8 to 4.4-5.6 and the O-O
 * diagonals close from the crystallographic 4.5-5.5 A to 3.7-4.1 A, at
 * every temperature (50-300 K), equilibration length (100-2000 steps)
 * and tether stiffness (0.5-8 eV/A^2) tested, with a filter strain of
 * 6-9 eV. The backbone and side-chain torsions are built and the ion is
 * held on the pore axis, but neither stops the ion's Coulomb pull (5.2 eV
 * per charge pair at 2.78 A) from closing the shell against a generic
 * angle stiffness and a generic dihedral force constant. That is a
 * property of this engine's reduced protein model, and it is reported,
 * not hidden.
 */
#define QM_KCSA_SITE_MIN_CN 7

/*
 * Run one ion's binding series in the real filter and fill `out`.
 *
 * Builds the deposited 1K4C TVGYG tetramer (via kcsa_filter.h), completes
 * the side-chain hydrogens, tethers the C-alpha trace at protein-like
 * stiffness, places the ion at `site`, minimises, then runs Langevin NVT
 * at `temp_K` for `n_equil` equilibration steps followed by `n_sample`
 * sampled steps.
 *
 * If `bind_series` is non-NULL it receives the per-sample E_bind values
 * (needs series_cap >= n_sample) for plotting or a second statistic; the
 * `filter_series` buffer does the same for E_filter. Either may be NULL.
 *
 * Returns 0 on success, -1 on bad input or allocation failure.
 */
int qm_kcsa_binding_series(int ion_Z, Vec3 site, double temp_K,
                           int n_equil, int n_sample, unsigned long seed,
                           double tethers_k,
                           double *bind_series, double *filter_series,
                           int series_cap, QmKcsaSeries *out);

/*
 * Ion-induced strain of the coordinating oxygen shell, as a linear elastic
 * response. This is the quantity Dixit & Asthagiri call dU_COCO and the
 * one their whole selectivity argument rests on.
 *
 * WHY IT IS DONE THIS WAY RATHER THAN BY DYNAMICS
 *
 * qm_kcsa_binding_series is the direct measurement, and on this engine it
 * does not converge: the filter holds its 8-fold site through
 * minimisation (CN 8 at the deposited 2.777 A mean ion-O distance) and
 * loses it within a few hundred steps of Langevin dynamics at every
 * temperature from 50 K to 300 K and every C-alpha restraint stiffness
 * from 0.5 to 16 eV/A^2.
 *
 * FULL-AUDIT: the paragraph below blamed the failure on the ABSENCE of
 * side-chain torsion terms. The torsions exist (kcsa_add_dihedrals) and
 * they are built here, so that reason was stale - and the conclusion was
 * still right, which is the dangerous kind of wrong. Measured with the
 * torsions in place: the site holds for the first ~50 fs (CN 8.000,
 * O-O 4.17-4.28 A) and then the ion's Coulomb pull closes the shell -
 * CN 6.2 at 100 fs, CN 4.0 with O-O at 3.28 A at 200 fs - at every
 * temperature, equilibration length and tether stiffness tested, with a
 * filter strain of 6-9 eV. Stiffening the angle term makes it WORSE
 * (CN 3.2 at 1.0 eV/rad^2, 2.7 at 2.2), so this is the limit of the
 * reduced protein model rather than a missing-stiffness parameter. The
 * two endpoints of that table are pinned by tests/test_audit4.c, and the
 * site-integrity gate below REFUSES to return a variance for a site that
 * was not maintained.
 *
 * The reason the torsions were added in the first place:
 * this engine's reduced protein model has no side-chain torsion terms, and
 * four of the eight filter oxygens belong to Thr75's side chain. Without
 * torsions, the hydroxyl rotates out of the shell and the coordination
 * number falls to 3-5 while the O-O diagonals close from the
 * crystallographic 4.5-5.5 A to ~3.5 A. No thermostat stiffness and no
 * restraint rescues it, because the missing term is the physics.
 *
 * Rather than report a "strain" number that is really the filter coming
 * apart, the strain is computed directly from the ion's forces on the
 * coordinating oxygens under a stated linear-elastic model of the
 * protein's response:
 *
 *     d_i     = F_ion,i / k            static displacement of oxygen i
 *     U_strain = sum_i |F_ion,i|^2 / (2k)
 *     dOO     = change in the shell's O-O diagonals at those displacements
 *
 * `k` is the protein's effective stiffness against a local oxygen
 * displacement and is THE parameter of the model. It is not knowable from
 * first principles here, so the caller is expected to SCAN it and report
 * the range; a single k presented as an answer would be a fitted number
 * wearing a derived number's clothes.
 *
 * Valid in the linear regime only, which the scan reveals: at small k the
 * oxygens displace further than the harmonic response is entitled to and
 * the estimate saturates, which shows up as the O-O diagonal falling
 * below the DFT/MM values.
 *
 * EACH ION IS PLACED AT ITS OWN RADIAL EQUILIBRIUM FIRST, AND THE FILTER
 * IS NEVER RELAXED. The deposited ion coordinates are K+'s; dropping Na+
 * into them and reading the forces compares two ions at one geometry,
 * where the Coulomb pulls are identical by construction and only the
 * Lennard-Jones sites differ. That produces strains differing by ~1% and
 * answers nothing, because the physical claim being tested is precisely
 * that Na+ sits closer to the same eight oxygens and so pulls harder.
 *
 * The ion's radius is relaxed by golden section with its axial coordinate
 * held at the deposited site, as kcsa_site_binding does and for the reason
 * it documents: the C4 filter makes azimuth immaterial, the problem is
 * one-dimensional, and a free ion ranges to whichever site is the global
 * minimum, which is not a site-binding measurement. The filter atoms stay
 * exactly where the structure refined them — the only state in which the
 * 8-fold shell is known to exist. Relaxing the filter instead, with no
 * bonded terms live, was measured here as CN 4.
 */
typedef struct {
    int    ok;
    int    shell_n;          /* coordinating oxygens used                  */
    double f_max_eV_A;       /* largest single ion-O force magnitude        */
    double f_rms_eV_A;       /* rms over the shell                         */
    double delta_A;          /* static oxygen displacement at stiffness k   */
    double strain_eV;        /* sum |F|^2 / 2k over the shell               */
    double oo_before_A;      /* mean shell O-O diagonal, undisplaced       */
    double oo_after_A;       /* mean shell O-O diagonal, displaced         */
    double cn;               /* coordination number at the placed site     */
    double mean_iono_A;      /* mean ion-O distance there                   */
    double target_contact_A; /* its preferred contact, Shannon + 1.40 A     */
    double contact_mismatch_A; /* achieved mean ion-O minus preferred      */
    double shell_mean_iono_A;  /* that mean, over the frozen 8-ligand set  */
    double mean_bind_eV;     /* filter binding energy there                 */
} QmKcsaElastic;

int qm_kcsa_elastic_strain(int ion_Z, Vec3 site, double stiffness,
                           QmKcsaElastic *out);

#ifdef __cplusplus
}
#endif

#endif /* QM_FLUCTUATION_H */