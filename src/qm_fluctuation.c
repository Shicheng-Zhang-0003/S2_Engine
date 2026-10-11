/* qm_fluctuation.c — the entropic half of ion binding, and the strain that
 * produces it. See include/qm_fluctuation.h for the physics and the
 * caveats.
 *
 * The measurement protocol is deliberately the simplest thing that can be
 * defended:
 *
 *   1. Build the deposited 1K4C TVGYG tetramer and complete its side-chain
 *      hydrogens (the same construction the record's flexible-filter legs
 *      already use, so this is not a new model, it is the old model with
 *      the ion left free to strain it).
 *   2. Tether the C-alpha trace so the peptide breathes around its
 *      crystallographic geometry instead of being frozen. Protein-like
 *      stiffness, from the same guide the duplex glycosidic proxy uses.
 *   3. Minimise with the ion present, so we start from the ion's own
 *      minimum rather than from a strained deposit.
 *   4. Equilibrate and sample under Langevin NVT.
 *   5. At every sample compute E_totals twice: once with the ion live,
 *      once with the ion's charge and LJ switched off at the SAME
 *      coordinates. The difference is that conformation's ion binding
 *      energy, exactly; the second value is the same-conformation strain
 *      reference.
 *
 * Step 5 is what makes the variance meaningful. Averaging the ion-free
 * energy over its own trajectory instead would conflate the strain with
 * the ordinary thermal wander of the peptide and wash out the signal;
 * the paired, same-coordinate difference removes the shared wander.
 */

#include "../include/qm_fluctuation.h"
#include "../include/kcsa_filter.h"
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/integrator.h"
#include "../include/constants.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/*───────────────────────────────────────────────────────────────────────────
 * Descriptive statistics
 * ───────────────────────────────────────────────────────────────────────────*/

/*
 * The variance error bar is the normal-approximation one,
 *
 *     SE(s^2) ~ s^2 * sqrt(2/n)
 *
 * which is exact in the limit for Gaussian samples and is the honest
 * thing to attach to a variance estimated from a trajectory. It is NOT a
 * block average, and for a correlated time series a block average would
 * give a LARGER and more correct bar. The correlation time is not
 * estimated here, so this number is a lower bound on the true error. It
 * is reported as such rather than presented as the error on the variance.
 */
static double series_sem_var(double var, int n) {
    if (n < 2 || !(var >= 0.0)) return INFINITY;
    return var * sqrt(2.0 / (double)n);
}

int qm_series_stats(const double *series_eV, int n, QmSeries *out) {
    if (!series_eV || !out || n < 1) return -1;
    double sum = 0.0, lo = INFINITY, hi = -INFINITY;
    for (int i = 0; i < n; i++) {
        double e = series_eV[i];
        if (!isfinite(e)) return -1;
        sum += e;
        if (e < lo) lo = e;
        if (e > hi) hi = e;
    }
    double mean = sum / (double)n;
    double ss = 0.0;
    for (int i = 0; i < n; i++) {
        double d = series_eV[i] - mean;
        ss += d * d;
    }
    /* Unbiased (ddof=1) for n>1. With a single sample the variance is not
     * estimable and is reported as 0 with an infinite error bar, which is
     * the honest statement rather than a fabricated number. */
    double var = (n > 1) ? ss / (double)(n - 1) : 0.0;
    out->n          = n;
    out->mean_eV    = mean;
    out->var_eV2    = var;
    out->sigma_eV   = sqrt(var);
    out->sem_eV     = (n > 1) ? sqrt(var / (double)n) : INFINITY;
    out->sem_var_eV2 = series_sem_var(var, n);
    out->min_eV     = lo;
    out->max_eV     = hi;
    return 0;
}

double qm_fluctuation_term_eV(double var_a_eV2, double var_b_eV2) {
    if (!isfinite(var_a_eV2) || !isfinite(var_b_eV2)) return 0.0;
    return 0.5 * (var_a_eV2 - var_b_eV2);
}

/*───────────────────────────────────────────────────────────────────────────
 * The KcsA binding series
 * ────────────────────────────────────────────────────────────────────────────*/

/*
 * The scan cutoff used for the coordination number and the mean ion-O
 * distance. 3.4 A is the filter's own deposited range (K-O spread across
 * all four sites is 2.688-3.082 A, header of kcsa_filter.h), so this
 * counts exactly the first shell the filter actually forms and does not
 * reach the next shell up.
 */
#define QM_KCSA_CN_CUTOFF 3.4

/* Langevin friction, 1/fs. 1/ps: slow compared with the 1 fs bond
 * vibrations (an OH stretch is 10 fs) and fast compared with the 100 ps
 * breathing motions being sampled, so it damps atomic detail without
 * freezing the conformational mode that carries the signal. */
#define QM_KCSA_GAMMA (1.0 / 1000.0)

/* Carbonyl oxygen radius, Angstrom — the same 1.40 A the KcsA contact
 * distances are built from (kcsa_filter.c's KCSA_O_RADIUS, Shannon O2- VI).
 * Mirrored here as a named constant rather than reaching into the .c
 * file, because the contact placement below has to state its own inputs;
 * tests/test_external.c asserts the two agree so they cannot drift. */
#define QM_KCSA_O_RADIUS 1.40

/* Azimuth samples for the in-plane site search. Eight is the symmetry of
 * the site: the filter is C4, so the in-plane energy is 90-degree
 * periodic, and eight samples resolve it with margin while keeping the
 * search deterministic and cheap. */
#define QM_KCSA_AZIMUTHS 8

/* Time step. 1 fs is the engine default and is comfortably inside the
 * stability limit for these bonded terms (types.h documents ~34 fs for the
 * stiffness in play). */
#define QM_KCSA_DT 1.0

/* Generic angle stiffness, eV/rad^2, applied to every angle by
 * sim_rebuild_angles_geometric. ~1 kcal/mol/rad^2: soft enough to stay
 * well inside the 1 fs stability limit for the light atoms in the filter,
 * firm enough that the amide plane does not fold. The engine cannot
 * assign a real amide stiffness here because its angle table is keyed on
 * element triple alone; this is a stated model limitation rather than a
 * parameter fitted to a selectivity number. */
#define QM_KCSA_ANGLE_K 0.35

/* Mean ion-oxygen distance over the coordinating oxygens.
 *
 * OXYGENS ONLY, matching kcsa_coord_stats' convention. This matters once
 * the side-chain hydrogens are in place: an unrestricted count inside the
 * same 3.4 A cutoff returns CN = 10.5, because the Thr75 OG1-H and Tyr78
 * OH hydrogen sit at ion-proximal distances and are not ligands. The
 * filter's defining statement is "eight oxygens", so the diagnostic has
 * to count the same species the filter statement counts, or the two
 * numbers cannot be compared. */
static double qm_kcsa_mean_iono(const Simulation *sim, int ion) {
    double sum = 0.0;
    int n = 0;
    for (int i = 0; i < sim->num_atoms; i++) {
        if (i == ion || sim->atoms[i].Z != 8) continue;
        double r = vec3_dist(sim->atoms[ion].position, sim->atoms[i].position);
        if (r <= QM_KCSA_CN_CUTOFF) { sum += r; n++; }
    }
    return (n > 0) ? sum / (double)n : 0.0;
}

static int qm_kcsa_cn(const Simulation *sim, int ion) {
    int n = 0;
    for (int i = 0; i < sim->num_atoms; i++) {
        if (i == ion || sim->atoms[i].Z != 8) continue;
        double r = vec3_dist(sim->atoms[ion].position, sim->atoms[i].position);
        if (r <= QM_KCSA_CN_CUTOFF) n++;
    }
    return n;
}

/*
 * Carbonyl-carbonyl strain, measured directly.
 *
 * The engine's TOTAL ion-free potential is a poor strain probe: it carries
 * the tethers, the bonded terms against the deposited r0, and every
 * non-coordinating atom's LJ, and those terms can move either way
 * independently of what the ion did to the filter. The term Dixit,
 * Asthagiri & Roux call dU_COCO is specifically the extra repulsion
 * between the FILTER's own coordinating oxygens, so that is what gets
 * measured: the mean O-O separation across the coordinating shell and the
 * Lennard-Jones energy between those oxygens, using the engine's own
 * combining rules and the engine's own 4*eps*((sigma/r)^12 - (sigma/r)^6).
 *
 * A shell-relative mean: the nearest-coordination-shell oxygens only.
 * Comparing a K+ shell against a Na+ shell that has a different
 * coordination number would otherwise compare populations, not strains.
 */
typedef struct {
    double mean_oo_A;   /* mean O-O distance across the shell (the diagonal) */
    double oo_lj_eV;    /* sum of LJ energy between shell oxygens             */
    int    n_shell;
} QmShellStrain;

static QmShellStrain qm_kcsa_shell_strain(const Simulation *sim, int ion) {
    QmShellStrain out = {0.0, 0.0, 0};
    int idx[64];
    int n = 0;
    for (int i = 0; i < sim->num_atoms && n < 64; i++) {
        if (i == ion || sim->atoms[i].Z != 8) continue;
        double r = vec3_dist(sim->atoms[ion].position, sim->atoms[i].position);
        if (r <= QM_KCSA_CN_CUTOFF) idx[n++] = i;
    }
    out.n_shell = n;
    if (n < 2) return out;
    double sum_r = 0.0, sum_e = 0.0;
    int npair = 0;
    for (int a = 0; a < n; a++) {
        for (int b = a + 1; b < n; b++) {
            double r = vec3_dist(sim->atoms[idx[a]].position,
                                 sim->atoms[idx[b]].position);
            sum_r += r;
            npair++;
            double eps   = lj_eps_combine(sim->atoms[idx[a]].lj_epsilon,
                                          sim->atoms[idx[b]].lj_epsilon);
            double sigma = lj_sigma_combine(sim->atoms[idx[a]].lj_sigma,
                                             sim->atoms[idx[b]].lj_sigma);
            if (r > 1e-9) {
                double sr2 = (sigma * sigma) / (r * r);
                double sr6 = sr2 * sr2 * sr2;
                sum_e += 4.0 * eps * (sr6 * sr6 - sr6);
            }
        }
    }
    if (npair > 0) {
        out.mean_oo_A = sum_r / (double)npair;
        out.oo_lj_eV  = sum_e;
    }
    return out;
}

/*
 * Run one ion's series. Shared by the K+ and Na+ legs so the protocol
 * cannot drift between them — the single most important thing this
 * function does, because a comparison of two runs made by two code paths
 * is a comparison of two protocols.
 */
static int kcsa_series_run(int ion_Z, Vec3 site, double temp_K,
                           int n_equil, int n_sample, unsigned long seed,
                           double tethers_k,
                           double *bind_series, double *filter_series,
                           int series_cap, QmKcsaSeries *out,
                           double *ref_filter_mean) {
    if (!out || n_sample < 2 || n_equil < 0) return -1;
    if (!(temp_K > 0.0) || !isfinite(temp_K)) return -1;
    if (bind_series && series_cap < n_sample) return -1;
    if (filter_series && series_cap < n_sample) return -1;
    if (!isfinite(tethers_k) || tethers_k <= 0.0) return -1;

    /* Capacity. Two counters, and both have to be right:
 *
 *   atoms: 4 segments x 41 heavy + 4 x 26 side-chain hydrogens + the ion.
 *   bonds: 4 segments x 41 + 4 x 26 (one per appended hydrogen).
 *
 * The BOND capacity is the one that bites. sim_create takes it as a
 * second argument and kcsa_add_hydrogens appends a bond per hydrogen; ask
 * for KCSA_FILTER_BONDS*4 + 16, as the record's other filter builders do,
 * and the builder fills the filter's own 164 bonds and then fails on the
 * 18th hydrogen with -1 — a hydrogen-less Thr/Val/Gly and a TYR with one
 * missing aromatic hydrogen, returned as a failure rather than silently
 * half-built, but an easy trap for the next caller. kcsa_add_hydrogens's
 * header now states the requirement. */
    const int cap_atoms  = KCSA_FILTER_ATOMS * 4 + 26 * 4 + 8;
    const int cap_bonds  = KCSA_FILTER_BONDS * 4 + 26 * 4 + 16;
    Simulation *sim = sim_create(cap_atoms, cap_bonds);
    if (!sim) return -1;

    int f = kcsa_build_filter(sim, vec3(0, 0, 0), 4);
    if (f < 0) { sim_destroy(sim); return -1; }

    /* The flexible-filter leg of the record already establishes that this
     * configuration (side-chain hydrogens completed, C-alpha tethers,
     * Langevin) HOLDS the filter around K+ — CN 8 -> 8 over 1 ps — and
     * that the same protocol collapses the empty filter, which is the
     * 1K4D inactivated state and not an artifact. That result is the
     * licence to measure a variance here: without it the "strain" being
     * measured would be the collapse of the model. Reused, not
     * rediscovered. */
    if (kcsa_add_hydrogens(sim, f, 4) < 0) { sim_destroy(sim); return -1; }

    sim->use_lj = 1; sim->use_coulomb = 1; sim->cutoff = 30.0;
    sim->use_bonds = 1; sim->use_angles = 1; sim->use_dihedrals = 0;

    /* ZERO-STRAIN the bonded terms against the deposited geometry BEFORE
     * anything is switched on. kcsa_build_filter deliberately does not
     * rebuild angles, because the engine's angle table falls back to a
     * generic 109.47 deg tetrahedral for an untabulated triple and a
     * peptide is planar — so as built, the bonded terms want to crush it.
     * Enabling them unzeroed collapsed the 8-fold site to CN 4 with the
     * shell diagonals closing from ~4.5-5.5 A to ~3.7 A, and every
     * "strain" number that came out of that was the deposit being
     * deformed by the force field rather than by the ion.
     *
     * kcsa_site_binding solves this by disabling the bonded terms
     * entirely, which is right for a rigid cage and useless here: the
     * strain mechanism REQUIRES a filter that can be strained. Zero-strain
     * is the third option and the one this measurement needs: the bonded
     * terms are live, and zero at the deposited geometry, so what they
     * report afterwards is deviation from experiment.
     *
     * The angle stiffness is a single generic value (QM_KCSA_ANGLE_K,
     * eV/rad^2) because the engine's per-angle table cannot distinguish a
     * planar amide from a tetrahedral centre on element alone. That is a
     * stated limitation of the reduced protein model, not a fitted value. */
    if (sim_zero_strain_bonds(sim) < 1) { sim_destroy(sim); return -1; }
    if (sim_rebuild_angles_geometric(sim, QM_KCSA_ANGLE_K) < 1) {
        sim_destroy(sim); return -1;
    }

    if (kcsa_add_calpha_tethers(sim, f, 4, tethers_k) < 0) {
        sim_destroy(sim); return -1;
    }

    int ion = (ion_Z > 0) ? sim_add_ion(sim, ion_Z, 1, site, 1.0) : -1;
    if (ion < 0) { sim_destroy(sim); return -1; }
    kcsa_set_ion_radius(sim, ion, ion_Z);

    /* Minimise with the ion present. This settles the deposited geometry
     * into the ion's own minimum before any sampling, so the series
     * measures the site rather than the transient from a strained start.
     * 600 steepest-descent iterations is enough: the potential is
     * dominated by LJ contacts and settles long before the budget, and
     * the Langevin equilibration that follows mops up the rest. This
     * engine has no neighbour list, so a force call on the full 269-atom
     * hydrogenated tetramer costs ~16 ms and the minimisation budget is
     * the single largest term in the leg's wall clock. */
    (void)integrator_minimize(sim, 600, 0.01, 1e-4);
    forces_calculate(sim);

    /* Seeded stream. Same seed for both ions => common random numbers, so
     * the variance DIFFERENCE is far better determined than either
     * variance. Both runs have identical atom counts and identical force
     * evaluation order, so the streams line up step for step. */
    integrator_maxwell_boltzmann(sim, temp_K, seed);
    integrator_remove_com_velocity(sim);
    sim->dt = QM_KCSA_DT;
    sim->thermostat.type = THERMOSTAT_LANGEVIN;
    sim->thermostat.target_temperature = temp_K;
    sim->thermostat.gamma = QM_KCSA_GAMMA;

    for (int s = 0; s < n_equil; s++) integrator_step(sim);
    forces_calculate(sim);

    /* Welford's online algorithm for both moments.
     *
     * This is NOT a convenience. An earlier revision computed the moments
     * in a second pass over the caller's optional series buffers and fell
     * back to the running mean when no buffer was supplied — which made
     * every variance exactly 0.0000 for any caller that did not ask for
     * the series, silently turning the module's entire reason for
     * existing into a constant. Moments must not depend on whether
     * anybody asked to see the samples. Welford is also the numerically
     * better form here: the binding energies run to tens of eV with
     * fluctuations of order 0.1, so a sum-of-squares about a separately
     * accumulated mean loses precision to cancellation. */
    double mean_b = 0.0, m2_b = 0.0;
    double mean_f = 0.0, m2_f = 0.0;
    double sum_oo = 0.0, sum_oolj = 0.0;
    int sum_cn = 0;
    for (int s = 0; s < n_sample; s++) {
        /* integrator_step's internal forces_calculate leaves
         * sim->potential_energy valid for the post-drift coordinates, and
         * the Langevin bath only touches velocities — so e_live is the
         * converged energy of this conformation with no extra call. The
         * counterfactual below is the only other one needed: the ion is
         * present but has no charge and no Lennard-Jones, so the field it
         * exerted is gone and the conformation's strain energy remains.
         * Turning the parameters off in place, rather than re-running a
         * separate ion-free trajectory, is what keeps the two samples
         * paired conformation for conformation. */
        integrator_step(sim);
        double e_live = sim->potential_energy;

        Atom *a = &sim->atoms[ion];
        double q_save = a->partial_charge;
        double e_save = a->lj_epsilon;
        double s_save = a->lj_sigma;
        a->partial_charge = 0.0;
        a->lj_epsilon = 0.0;
        forces_calculate(sim);
        double e_bare = sim->potential_energy;
        a->partial_charge = q_save;
        a->lj_epsilon = e_save;
        a->lj_sigma = s_save;

        double eb = e_live - e_bare;
        if (bind_series)   bind_series[s]   = eb;
        if (filter_series) filter_series[s] = e_bare;

        double nb = (double)(s + 1);
        double db = eb - mean_b;
        mean_b += db / nb;
        m2_b   += db * (eb - mean_b);
        double df = e_bare - mean_f;
        mean_f += df / nb;
        m2_f   += df * (e_bare - mean_f);

        sum_cn += qm_kcsa_cn(sim, ion);
        QmShellStrain sh = qm_kcsa_shell_strain(sim, ion);
        sum_oo   += sh.mean_oo_A;
        sum_oolj += sh.oo_lj_eV;
    }

    double mean_bind   = mean_b;
    double mean_filter = mean_f;
    double var_b = m2_b / (double)(n_sample - 1);
    double var_f = m2_f / (double)(n_sample - 1);

    out->ok = 1;
    out->n_samples = n_sample;
    out->mean_bind_eV = mean_bind;
    out->var_bind_eV2 = var_b;
    out->sem_bind_eV = sqrt(var_b / (double)n_sample);
    out->sem_var_bind_eV2 = series_sem_var(var_b, n_sample);
    out->mean_filter_eV = mean_filter;
    out->var_filter_eV2 = var_f;
    out->mean_iono_A = qm_kcsa_mean_iono(sim, ion);
    out->cn_mean = (double)sum_cn / (double)n_sample;
    out->mean_oo_A = sum_oo / (double)n_sample;
    out->oo_lj_eV  = sum_oolj / (double)n_sample;
    out->site_maintained = (out->cn_mean >= (double)QM_KCSA_SITE_MIN_CN) ? 1 : 0;
    if (ref_filter_mean) *ref_filter_mean = mean_filter;

    sim_destroy(sim);
    return 0;
}

/*
 * Run one ion's binding series. See the header for the protocol and the
 * caveats. `filter_series` (if given) must have room for n_sample.
 */
int qm_kcsa_binding_series(int ion_Z, Vec3 site, double temp_K,
                           int n_equil, int n_sample, unsigned long seed,
                           double tethers_k,
                           double *bind_series, double *filter_series,
                           int series_cap, QmKcsaSeries *out) {
    if (ion_Z <= 0) return -1;
    return kcsa_series_run(ion_Z, site, temp_K, n_equil, n_sample, seed,
                           tethers_k, bind_series, filter_series,
                           series_cap, out, NULL);
}

/*───────────────────────────────────────────────────────────────────────────
 * Ion-induced shell strain, linear elastic
 * ────────────────────────────────────────────────────────────────────────────*/

/* One shell oxygen's equilibrium position, saved before displacement so
 * the displaced and undisplaced O-O diagonals come from the same
 * reference geometry. */
#define QM_ELASTIC_MAX 32

static int qm_kcsa_shell_indices(const Simulation *sim, int ion, int *idx,
                                  int max_idx) {
    int n = 0;
    for (int i = 0; i < sim->num_atoms && n < max_idx; i++) {
        if (i == ion || sim->atoms[i].Z != 8) continue;
        double r = vec3_dist(sim->atoms[ion].position, sim->atoms[i].position);
        if (r <= QM_KCSA_CN_CUTOFF) idx[n++] = i;
    }
    return n;
}

/* Mean O-O separation across the given shell, at the CURRENT positions. */
static double qm_shell_mean_oo(const Simulation *sim, const int *idx, int n) {
    double sum = 0.0;
    int np = 0;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            sum += vec3_dist(sim->atoms[idx[a]].position,
                             sim->atoms[idx[b]].position);
            np++;
        }
    return (np > 0) ? sum / (double)np : 0.0;
}

int qm_kcsa_elastic_strain(int ion_Z, Vec3 site, double stiffness,
                           QmKcsaElastic *out) {
    if (!out || ion_Z <= 0) return -1;
    if (!isfinite(stiffness) || stiffness <= 0.0) return -1;
    memset(out, 0, sizeof *out);
    out->oo_before_A = 0.0;

    /* The rigid-filter build used by every binding leg: no bonded terms,
     * no hydrogens, no tethers. This measures the ion's pull on the
     * deposited geometry, which is the reference state, not a disturbed
     * one. */
    Simulation *sim = sim_create(KCSA_FILTER_ATOMS * 4 + 8,
                                 KCSA_FILTER_BONDS * 4 + 16);
    if (!sim) return -1;
    int f = kcsa_build_filter(sim, vec3(0, 0, 0), 4);
    if (f < 0) { sim_destroy(sim); return -1; }
    sim->use_lj = 1; sim->use_coulomb = 1; sim->cutoff = 30.0;
    sim->use_bonds = 0; sim->use_angles = 0; sim->use_dihedrals = 0;
    /* The overlap (exchange) repulsion and the dispersion term are ON.
     *
     * Without them the two ions are indistinguishable and the strain is
     * meaningless: with the filter held at its crystallographic geometry
     * every ion sits 2.777 A from the eight oxygens whatever its radius,
     * so the Coulomb pull on the shell is identical to 1% and the
     * Lennard-Jones wall (0.05 kcal/mol) is four orders of magnitude too
     * weak to register against it. Measured with LJ + Coulomb alone, the
     * K and Na strain energies differ by 1% because they are the same
     * geometry twice.
     *
     * That is not a modelling shortcut, it is the point. The filter's
     * GEOMETRY sets the ion-oxygen distance — the oxygens are where the
     * protein put them — so the only thing that can distinguish a K+ from
     * a Na+ at the same distance is how hard each one's electron cloud
     * is pushed back by the oxygen's. That is the exchange term, it is
     * size-sensitive by construction, and it is a quantum-mechanical
     * quantity rather than a fitted Lennard-Jones depth. qm_pauli
     * computes it from the Slater overlap of the two real atoms, so this
     * is where the engine's quantum layer earns its place in the KcsA
     * story at all. */
    sim->use_pauli = 1;
    sim->use_disp  = 1;

    int ion = sim_add_ion(sim, ion_Z, 1, site, 1.0);
    if (ion < 0) { sim_destroy(sim); return -1; }
    kcsa_set_ion_radius(sim, ion, ion_Z);

    int lig[QM_ELASTIC_MAX];
    int n_lig = 0;

    /* EACH ION FINDS ITS OWN SITE GEOMETRY, IN THE PLANE OF ITS SITE.
     *
     * Three placements were tried and each one is wrong in a way worth
     * recording.
     *
     * (1) Relax the ion's radius to the filter's potential minimum. That
     * returns r = 0 for BOTH ions, at the same on-axis point with the same
     * 2.777 A mean ion-O distance, because the ion-O interaction here is a
     * 0.05 kcal/mol Lennard-Jones wall with no minimum of its own and the
     * Coulomb term simply parks the cation on the axis. Measured that way
     * the two ions' strains agree to 1%, because it is one measurement
     * taken twice.
     *
     * (2) Impose each ion's Shannon contact distance (2.91 A for K, 2.58 A
     * for Na) by bisecting the ion's radius. This is where the real
     * content of KcsA shows up. On the axis, the mean ion-O distance over
     * the eight crystallographic oxygens is 2.777 A and is already the
     * SHORTEST value the geometry permits — r cannot go below zero. So
     * K+'s 2.91 A target is unreachable on the short side by 0.13 A,
     * while Na+'s 2.58 A target is unreachable by 0.20 A. Na+ cannot get
     * within 0.20 A of its own preferred contact anywhere on this filter's
     * axis. That is the whole exclusion mechanism, in one number, and it
     * is a property of the deposited geometry rather than of any force
     * field.
     *
     * (3) The bisection itself cannot satisfy an unreachable target, so it
     * pins both ions at r = 0 and the contact mismatch is silently lost.
     *
     * What is actually measured instead is each ion's own minimum in the
     * PLANE of its site — radius and azimuth free, axial coordinate held.
     * K+ is predicted to stay on the axis; Na+ is predicted to slide off
     * it, because being small in a cage sized for a larger ion is exactly
     * the situation in which a cation rattles, and a rattling ion is what
     * gives its binding energy a wide variance. The axial coordinate is
     * held because letting it range slides the cation to a neighbouring
     * site and stops measuring this one (kcsa_site_binding documents that
     * trap at length).
     *
     * The azimuth is scanned on a fixed grid and the radius is then golden
     * -sectioned at each azimuth: deterministic, no gradient projection,
     * and the reported minimum is the engine's own energy at its own
     * minimum rather than a reimplementation of it. */
    {
        sim->atoms[ion].position = vec3(0.0, 0.0, site.z);
        forces_calculate(sim);
        n_lig = qm_kcsa_shell_indices(sim, ion, lig, QM_ELASTIC_MAX);
        if (n_lig < 2) { sim_destroy(sim); return -1; }

        const double PHI = 0.6180339887498949;
        double best_e = 1e300, best_x = 0.0, best_y = 0.0;
        for (int t = 0; t < QM_KCSA_AZIMUTHS; t++) {
            double az = 2.0 * S2_PI * (double)t / (double)QM_KCSA_AZIMUTHS;
            double ca = cos(az), sa = sin(az);
            double a = 0.0, b = 4.5;
            for (int it = 0; it < 40 && (b - a) > 1e-5; it++) {
                double c1 = b - PHI * (b - a);
                double c2 = a + PHI * (b - a);
                sim->atoms[ion].position = vec3(c1 * ca, c1 * sa, site.z);
                forces_calculate(sim);
                double e1 = sim->potential_energy;
                sim->atoms[ion].position = vec3(c2 * ca, c2 * sa, site.z);
                forces_calculate(sim);
                double e2 = sim->potential_energy;
                if (e1 < e2) b = c2; else a = c1;
            }
            double r = 0.5 * (a + b);
            sim->atoms[ion].position = vec3(r * ca, r * sa, site.z);
            forces_calculate(sim);
            if (sim->potential_energy < best_e) {
                best_e = sim->potential_energy;
                best_x = r * ca;
                best_y = r * sa;
            }
        }
        sim->atoms[ion].position = vec3(best_x, best_y, site.z);
        out->target_contact_A = kcsa_cation_radius(ion_Z) + QM_KCSA_O_RADIUS;
    }

    /* NO RELAXATION OF THE FILTER, BY DESIGN. The strain being measured is the ion's
     * pull on the coordinating oxygens AT THE CRYSTALLOGRAPHIC
     * GEOMETRY, which is the state the structure was refined in and the
     * only state in which the 8-fold shell is known to exist. Minimising
     * here does not refine the measurement, it destroys it: with no bonded
     * terms the filter is a rigid cloud and the ION is the only free body,
     * so an unconstrained minimiser slides the cation off axis or down to
     * a neighbouring site and reports the coordination of whatever it
     * landed on (CN 4 at this site, measured). That is the trap
     * kcsa_site_binding documents at length, and this function walks
     * straight into it if it relaxes.
     *
     * The forces are well defined without relaxation — a force field is a
     * field — and at the deposited geometry they are the field of a real
     * site. */
    forces_calculate(sim);
    out->cn = qm_kcsa_cn(sim, ion);
    out->mean_iono_A = qm_kcsa_mean_iono(sim, ion);
    out->mean_bind_eV = sim->potential_energy;
    /* The contact mismatch: how far short of its own preferred ion-oxygen
     * distance this ion is forced to sit in a filter sized for K+. It is
     * computed over the FROZEN eight-ligand set, so it means the same
     * thing for both ions and does not drift with whichever oxygens happen
     * to be in range at each ion's own minimum. */
    {
        double sum = 0.0;
        for (int j = 0; j < n_lig; j++)
            sum += vec3_dist(sim->atoms[ion].position,
                             sim->atoms[lig[j]].position);
        out->shell_mean_iono_A = sum / (double)n_lig;
        out->contact_mismatch_A =
            out->shell_mean_iono_A - out->target_contact_A;
    }

    /* The SAME frozen ligand set the contact was placed against: the
     * strain is measured on the shell whose geometry was matched to the
     * ion's contact distance, not on a shell that re-selected itself
     * around wherever the ion ended up. */
    int *idx = lig;
    int n = n_lig;
    out->shell_n = n;

    /* Force on each atom WITH the ion. */
    double *live_all = (double *)malloc(sizeof(double) * 3 *
                                        (size_t)sim->num_atoms);
    if (!live_all) { sim_destroy(sim); return -1; }
    for (int i = 0; i < sim->num_atoms; i++) {
        live_all[3 * i + 0] = sim->atoms[i].force.x;
        live_all[3 * i + 1] = sim->atoms[i].force.y;
        live_all[3 * i + 2] = sim->atoms[i].force.z;
    }

    /* The SAME forces with the ion neutralised in place: what is left is
     * the filter's own force at these coordinates, and the difference is
     * the ion's pull on every atom. */
    Atom *a = &sim->atoms[ion];
    double q_save = a->partial_charge, e_save = a->lj_epsilon;
    double s_save = a->lj_sigma;
    a->partial_charge = 0.0; a->lj_epsilon = 0.0;
    forces_calculate(sim);
    a->partial_charge = q_save; a->lj_epsilon = e_save; a->lj_sigma = s_save;

    double fsum2 = 0.0, fmax2 = 0.0, dsum2 = 0.0;
    Vec3 ref[QM_ELASTIC_MAX];
    for (int j = 0; j < n; j++) {
        int i = idx[j];
        ref[j] = sim->atoms[i].position;
        Vec3 f = vec3(live_all[3 * i + 0] - sim->atoms[i].force.x,
                      live_all[3 * i + 1] - sim->atoms[i].force.y,
                      live_all[3 * i + 2] - sim->atoms[i].force.z);
        double f2 = vec3_norm2(f);
        if (f2 > fmax2) fmax2 = f2;
        fsum2 += f2;
        Vec3 d = vec3_scale(f, 1.0 / stiffness);   /* d = F / k */
        dsum2 += vec3_norm2(d);
        sim->atoms[i].position = vec3_add(ref[j], d);
    }
    out->f_max_eV_A = sqrt(fmax2);
    out->f_rms_eV_A = sqrt(fsum2 / (double)n);
    out->delta_A    = sqrt(dsum2 / (double)n);
    out->strain_eV  = fsum2 / (2.0 * stiffness);

    /* Geometry before and after, from the same reference. */
    for (int j = 0; j < n; j++) sim->atoms[idx[j]].position = ref[j];
    out->oo_before_A = qm_shell_mean_oo(sim, idx, n);
    for (int j = 0; j < n; j++) {
        int i = idx[j];
        Vec3 f = vec3(live_all[3 * i + 0] - sim->atoms[i].force.x,
                      live_all[3 * i + 1] - sim->atoms[i].force.y,
                      live_all[3 * i + 2] - sim->atoms[i].force.z);
        sim->atoms[i].position = vec3_add(ref[j], vec3_scale(f, 1.0 / stiffness));
    }
    out->oo_after_A = qm_shell_mean_oo(sim, idx, n);

    out->ok = 1;
    free(live_all);
    sim_destroy(sim);
    return 0;
}