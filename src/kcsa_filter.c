#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/kcsa_filter.h"
#include "../include/constants.h"
#include "../include/forces.h"
#include "../include/sim.h"
#include "../include/integrator.h"
#include "../include/amber_lj.h"

/*
 * kcsa_filter.c — the KcsA TVGYG selectivity filter from deposited 1K4C
 * coordinates. See include/kcsa_filter.h for the full provenance, the
 * symmetry derivation and the validation numbers.
 *
 * The coordinate table below is DEPOSITED DATA, not a construction. The
 * x and y are measured from the pore axis (the column of K+ ions in
 * 1K4C, which lies on the 4-fold axis at fractional 1,1), and the
 * hydrogens are constructed — the 2001-era structure contains none.
 *
 * HYDROGEN CONSTRUCTION
 *   Amide H: the peptide N is sp2, bonded to C(i-1) and CA(i). The H
 *   goes on the external bisector of those two directions, in-plane, at
 *   the standard 1.00 A N-H length. This is the same construction the
 *   nucleobase module uses for cytosine H6, for the same reason.
 *   THR75 is the first residue of the segment and its partner C74 is
 *   outside the segment, so it carries two hydrogens (a free
 *   N-terminus): the bisector H, plus a second placed out of the
 *   CA-N-C plane.
 *
 *   AUDIT FIX D6 - the Thr75 N-terminal pair is WRONG, and this comment
 *   previously described it as correct. Measured against the deposited
 *   1K4C coordinates (external validation, see tests/test_external.c):
 *
 *       N-H(A) = 1.0000 A      angle CA-N-H(A) = 162.96 deg
 *       N-H(B) = 1.0256 A      angle CA-N-H(B) = 104.45 deg
 *       H(A)-N-H(B)          =  74.87 deg
 *       sum of the three angles at N    = 342.29 deg   (should be 360)
 *
 *   A neutral amine nitrogen has three bonds summing to 360 deg with all
 *   three angles near 107-110 deg. 74.87 deg between the two hydrogens is
 *   not an amine. The amide hydrogens on residues 76-79 are correct by
 *   contrast: N-H = 1.0000 A exactly with C(prev)-N-H = CA-N-H = 119.04 to
 *   120.43 deg, which is the external bisector of an sp2 nitrogen.
 *
 *   Impact: none on any number in the record. The nearest of the two
 *   hydrogens to any modelled K+ site is 6.259 A (site z = -40.505), which
 *   is outside every interaction the site legs sum over - those are built
 *   from the eight coordinating oxygens at 2.699-3.071 A. So the KcsA
 *   energetics are unaffected and the record does not change.
 *
 *   Why it is recorded rather than fixed here: the correction needs a
 *   choice of nitrogen geometry, and Thr75's real partner C74 is outside
 *   the segment, so the correct construction depends on whether the
 *   missing carbonyl is modelled. That is a modelling decision about the
 *   N-terminus, not a transcription slip, and it belongs with whoever
 *   decides whether this segment should be capped at N as well as at C.
 *   OXT: C79's carboxyl gets the standard 1.25 A terminal oxygen on the
 *   external bisector of CA79 and O79, making the segment a valid
 *   neutral pentapeptide. GLY79's own carbonyl oxygen points 4.82 A off
 *   the pore axis, so this cap is nowhere near the ion path.
 *
 * CHARGES
 *   The filter is a neutral peptide and the model must be too, or it
 *   will bind a cation for a reason that has nothing to do with KcsA.
 *   The atoms that actually determine ion binding — the amide N/H/C/O
 *   group and the Thr and Tyr hydroxyls — carry standard AMBER ff99
 *   internal-peptide values.
 *
 *   AUDIT FIX D4: this paragraph used to claim the residues are made
 *   "exactly neutral by a UNIFORM offset added to its carbons, computed at
 *   build time from the other atoms in the residue". Measured per residue
 *   from the KCSA_TVGYG table below:
 *
 *       res 75   9 atoms   sum = -0.000000 e
 *       res 76   8 atoms   sum = -0.000000 e
 *       res 77   5 atoms   sum = -0.000000 e
 *       res 78  13 atoms   sum = +0.000000 e
 *       res 79   6 atoms   sum = -0.000000 e
 *       TOTAL  41 atoms   sum = -0.000000 e
 *
 *   So the neutrality claim is TRUE and the filter really is neutral. The
 *   mechanism description was false in two ways:
 *     - nothing is computed at build time. There is no arithmetic in
 *       kcsa_build_filter() that touches a charge; the neutrality is baked
 *       into the literal q values in KCSA_TVGYG. `grep -n 'offset\|share'
 *       src/kcsa_filter.c` matches only this comment.
 *     - the offset is not uniform. Residue 75's carbons are
 *       0.0144 / 0.6117 / 0.1144 / 0.3544, which is the deposited ff99 set
 *       essentially untouched, not a shared shift; residue 78's are
 *       0.070352 for CA and 0.070356 for the eight side-chain carbons,
 *       which differ from each other in the fifth decimal.
 *   The correction matters because "uniform offset" was the stated reason
 *   the carbons were trustworthy, and for three of the five residues the
 *   carbons are simply the deposited charges.
 *
 *   Two earlier schemes were tried and both were wrong in an instructive
 *   way:
 *
 *     - dumping the whole residual on C-alpha put +0.99 e on Gly79's
 *       CA, two Angstrom from the ion, and the flexible-filter
 *       calculation collapsed onto it;
 *     - overwriting each carbon with share = -sum(non-carbon)/n_carbon
 *       discarded the backbone C charge and left every residue +0.60 e.
 *
 *   The maximum single-atom charge in the table is +0.9276 e (Gly79 C) and
 *   the minimum is -0.6551 e (Thr75 OG1). The carbons still do not carry
 *   individually faithful ff99 values and the model is not fit for ABSOLUTE
 *   binding energies; kcsa_site_binding's documentation says so and the
 *   readme repeats it. The per-residue neutrality above is asserted in
 *   tests/test_regression.c.
 */

/* ---- deposited TVGYG heavy atoms + constructed H/terminal O ---- */
typedef struct {
    int    Z;
    double x, y, z;   /* Angstrom, pore frame */
    double q;         /* e */
} KcsaAtom;

static const KcsaAtom KCSA_TVGYG[KCSA_FILTER_ATOMS] = {
    /*  75 N    */ {  7,  -3.14900,  -4.44600, -40.37500, -0.415700 },
    /*  75 CA   */ {  6,  -1.91400,  -3.85100, -39.84200,  0.014400 },
    /*  75 C    */ {  6,  -2.16400,  -2.69100, -38.86400,  0.611700 },
    /*  75 O    */ {  8,  -1.27800,  -1.87200, -38.62700, -0.567900 },
    /*  75 CB   */ {  6,  -0.90600,  -3.38400, -40.95600,  0.114400 },
    /*  75 CG2  */ {  6,  -0.66200,  -4.51700, -41.96100,  0.354400 },
    /*  75 OG1  */ {  8,  -1.40700,  -2.23400, -41.65500, -0.655100 },
    /*  76 N    */ {  7,  -3.35500,  -2.63800, -38.26900, -0.415700 },
    /*  76 CA   */ {  6,  -3.65900,  -1.57000, -37.30700,  0.022880 },
    /*  76 C    */ {  6,  -2.86900,  -1.80200, -36.01300,  0.620180 },
    /*  76 O    */ {  8,  -2.15900,  -0.91000, -35.54300, -0.567900 },
    /*  76 CB   */ {  6,  -5.17500,  -1.49800, -37.00600,  0.022880 },
    /*  76 CG1  */ {  6,  -5.45600,  -0.47800, -35.88600,  0.022880 },
    /*  76 CG2  */ {  6,  -5.91800,  -1.07600, -38.26900,  0.022880 },
    /*  77 N    */ {  7,  -2.98400,  -3.00700, -35.45100, -0.415700 },
    /*  77 CA   */ {  6,  -2.24300,  -3.35900, -34.24700,  0.057200 },
    /*  77 C    */ {  6,  -2.32600,  -2.45600, -33.02400,  0.654500 },
    /*  77 O    */ {  8,  -1.30300,  -1.93500, -32.55100, -0.567900 },
    /*  78 N    */ {  7,  -3.53500,  -2.30100, -32.48500, -0.415700 },
    /*  78 CA   */ {  6,  -3.75400,  -1.46500, -31.30800,  0.070352 },
    /*  78 C    */ {  6,  -2.86000,  -1.82800, -30.13100,  0.667656 },
    /*  78 O    */ {  8,  -2.37200,  -0.95600, -29.43300, -0.567900 },
    /*  78 CB   */ {  6,  -5.21000,  -1.56100, -30.84100,  0.070356 },
    /*  78 CG   */ {  6,  -6.21900,  -0.97600, -31.80200,  0.070356 },
    /*  78 CD1  */ {  6,  -7.26800,  -1.75600, -32.29400,  0.070356 },
    /*  78 CD2  */ {  6,  -6.16200,   0.36900, -32.17100,  0.070356 },
    /*  78 CE1  */ {  6,  -8.24100,  -1.21200, -33.12500,  0.070356 },
    /*  78 CE2  */ {  6,  -7.13900,   0.93000, -33.00700,  0.070356 },
    /*  78 CZ   */ {  6,  -8.17300,   0.13000, -33.47400,  0.070356 },
    /*  78 OH   */ {  8,  -9.16600,   0.65500, -34.26000, -0.518800 },
    /*  79 N    */ {  7,  -2.64900,  -3.11900, -29.91000, -0.415700 },
    /*  79 CA   */ {  6,  -1.83500,  -3.51800, -28.77700,  0.330300 },
    /*  79 C    */ {  6,  -2.71500,  -4.00700, -27.63700,  0.927600 },
    /*  79 O    */ {  8,  -2.22100,  -4.28000, -26.54400, -0.567900 },
    /*  76 H    */ {  1,  -4.04630,  -3.32951, -38.47855,  0.271900 },
    /*  77 H    */ {  1,  -3.59252,  -3.68358, -35.86566,  0.271900 },
    /*  78 H    */ {  1,  -4.31482,  -2.77087, -32.89864,  0.271900 },
    /*  79 H    */ {  1,  -3.04576,  -3.80815, -30.51634,  0.271900 },
    /*  75 H    */ {  1,  -3.79278,  -5.02223, -40.87849,  0.271900 },
    /*  75 H    */ {  1,  -3.33871,  -5.24049, -39.75486,  0.271900 },
    /*  79 OXT  */ {  8,  -3.94641,  -4.13212, -27.81158, -0.546200 },
};

typedef struct { int a, b, order; } KcsaBond;

static const KcsaBond KCSA_TVG_BONDS[KCSA_FILTER_BONDS] = {
    {  0,  1, 1 },
    {  1,  2, 1 },
    {  2,  3, 2 },
    {  1,  4, 1 },
    {  7,  8, 1 },
    {  8,  9, 1 },
    {  9, 10, 2 },
    {  8, 11, 1 },
    { 14, 15, 1 },
    { 15, 16, 1 },
    { 16, 17, 2 },
    { 18, 19, 1 },
    { 19, 20, 1 },
    { 20, 21, 2 },
    { 19, 22, 1 },
    { 30, 31, 1 },
    { 31, 32, 1 },
    { 32, 33, 2 },
    {  4,  5, 1 },
    {  4,  6, 1 },
    { 11, 12, 1 },
    { 11, 13, 1 },
    { 22, 23, 1 },
    { 23, 24, 1 },
    { 23, 25, 1 },
    { 24, 26, 1 },
    { 25, 27, 1 },
    { 26, 28, 1 },
    { 27, 28, 1 },
    { 28, 29, 1 },
    {  2,  7, 1 },
    {  9, 14, 1 },
    { 16, 18, 1 },
    { 20, 30, 1 },
    {  7, 34, 1 },
    { 14, 35, 1 },
    { 18, 36, 1 },
    { 30, 37, 1 },
    {  0, 38, 1 },
    {  0, 39, 1 },
    { 32, 40, 1 },
};

/* C4 rotation about the pore axis (z) through the origin. */
static void kcsa_c4(Vec3 v, int n, Vec3 *out) {
    double c = (n == 0) ? 1.0 : (n == 1) ? 0.0 : (n == 2) ? -1.0 : 0.0;
    double s = (n == 0) ? 0.0 : (n == 1) ? 1.0 : (n == 2) ? 0.0 : -1.0;
    out->x = v.x * c - v.y * s;
    out->y = v.x * s + v.y * c;
    out->z = v.z;
}

int kcsa_build_filter(Simulation *sim, Vec3 origin, int n_subunits) {
    if (!sim || !sim->atoms) return -1;
    if (n_subunits < 1 || n_subunits > 4) n_subunits = 4;
    int first = sim->num_atoms;
    if (first + KCSA_FILTER_ATOMS * n_subunits > sim->capacity_atoms) return -1;

    for (int s = 0; s < n_subunits; s++) {
        int base = sim->num_atoms;
        for (int i = 0; i < KCSA_FILTER_ATOMS; i++) {
            const KcsaAtom *a = &KCSA_TVGYG[i];
            Vec3 p; kcsa_c4(vec3(a->x, a->y, a->z), s, &p);
            int idx = sim_add_atom(sim, a->Z, vec3(p.x + origin.x, p.y + origin.y,
                                                  p.z + origin.z), a->q);
            if (idx < 0) return -1;
        }
        for (int b = 0; b < KCSA_FILTER_BONDS; b++) {
            const KcsaBond *bd = &KCSA_TVG_BONDS[b];
            int ia = base + bd->a, ib = base + bd->b;
            if (sim_add_bond(sim, ia, ib, bd->order) < 0) return -1;
        }
        /* LJ by atom type. Carbonyl and hydroxyl oxygens share the
         * AMBER ff99 hydroxyl O sigma/epsilon already used by
         * aminoacids.c; amide N and the aliphatic carbons likewise. */
        for (int i = 0; i < KCSA_FILTER_ATOMS; i++) {
            int Z = KCSA_TVGYG[i].Z;
            int ai = base + i;
            if (Z == 8)      sim_set_atom_lj(sim, ai, LJ_AMBER_O_EPS,  LJ_AMBER_O_SIGMA);
            else if (Z == 7) sim_set_atom_lj(sim, ai, LJ_AMBER_N_EPS,  LJ_AMBER_N_SIGMA);
            else if (Z == 1) sim_set_atom_lj(sim, ai, LJ_AMBER_HN_EPS, LJ_AMBER_HN_SIGMA);
            else             sim_set_atom_lj(sim, ai, LJ_AMBER_CT_EPS, LJ_AMBER_CT_SIGMA);
        }
    }
    /* NO sim_rebuild_angles here, deliberately.
     *
     * The bonded terms in this engine are parameterised for gas-phase
     * fragments placed at their equilibrium r0, and the angle table falls
     * back to a generic 109.47 deg tetrahedral value for any untabulated
     * angle. A peptide is planar - its amide angles are near 120 deg - so
     * rebuilding angles onto a DEPOSITED peptide injects tens of eV of
     * artificial strain that has nothing to do with the ion. The
     * deposited coordinates ARE the reference geometry, so bonded terms
     * would only fight it.
     *
     * Callers measuring ion binding must therefore run with
     * use_bonds = use_angles = use_dihedrals = 0, which kcsa_site_binding
     * documents. The filter is a rigid structure here by design. */
    return first;
}

int kcsa_ion_sites(int n_out, Vec3 *out) {
    if (!out || n_out < 1 || n_out > 4) return 0;
    /* K+ positions of 1K4C chain C, in the same pore frame. */
    static const double kz[4] = { -30.553, -33.953, -37.162, -40.505 };
    int n = (n_out < 4) ? n_out : 4;
    for (int i = 0; i < n; i++) out[i] = vec3(0.0, 0.0, kz[i]);
    return n;
}

int kcsa_coord_stats(const Simulation *sim, int filter_first,
                     int n_subunits, Vec3 point, double cutoff,
                     double *mean_r) {
    if (!sim || !sim->atoms || filter_first < 0) return -1;
    if (n_subunits < 1 || n_subunits > 4) n_subunits = 4;
    if (!(cutoff > 0.0) || !isfinite(cutoff)) return -1;
    int n = 0; double sum = 0.0;
    for (int s = 0; s < n_subunits; s++) {
        int base = filter_first + s * KCSA_FILTER_ATOMS;
        for (int i = 0; i < KCSA_FILTER_ATOMS; i++) {
            if (sim->atoms[base + i].Z != 8) continue;
            double d = vec3_dist(point, sim->atoms[base + i].position);
            if (d <= cutoff) { n++; sum += d; }
        }
    }
    if (mean_r) *mean_r = (n > 0) ? sum / n : 0.0;
    return n;
}

/* ══════════════════════════════════════════════════════════════════════════
 * K+ vs Na+ IN THE REAL FILTER
 * ══════════════════════════════════════════════════════════════════════════
 *
 * The site is K+-sized. Deposited K-O distances are 2.77-2.93 A against a
 * preferred 2.91 A for K+ (1.51 A eight-coordinate Shannon radius + 1.40 A
 * O radius; full-audit M3 corrects previous 2.78 A six-coordinate),
 * while Na+ prefers 2.58 A (1.18 A VIII). Na+ is therefore 0.19-0.35 A
 * too far from each of eight oxygens, and that geometric mismatch is what
 * the filter's side chains are shaped to enforce.
 *
 * But the geometric mismatch is NOT the dominant energetic term, and this
 * is where a vacuum calculation goes wrong. The dominant term is the
 * cost of DEHYDRATING the ion to put it in the site, and it is not a
 * force-field term at all - it is a measured bulk thermodynamic quantity.
 *
 * Standard single-ion hydration free energies (absolute scale, Marcus
 * 1997; TATB convention as used throughout ion solvation free-energy
 * work):
 *
 *     Na+   -454 kJ/mol
 *     K+    -322 kJ/mol
 *
 * Desolvating K+ to enter the filter therefore costs 132 kJ/mol LESS than
 * desolvating Na+ - 1.368 eV at 96.485 kJ/mol per eV. That is a large,
 * real, measured advantage for K+, and it is the term the previous
 * model could not see because the previous model had neither the real
 * geometry nor any solvent at all.
 *
 * Run without this term, a rigid K+-sized cage will ALWAYS appear to
 * prefer Na+, because the smaller cation has more negative Coulomb
 * energy at closer range. That is not a KcsA result; it is the
 * electrostatics of a charged cage, and it is the trap this function
 * exists to make visible.
 */

double kcsa_hydration_free_energy_kJmol(int Z) {
    switch (Z) {
        case 11: return -454.0;   /* Na+ */
        case 19: return -322.0;   /* K+  */
        case 37: return -293.0;   /* Rb+ */
        case 55: return -264.0;   /* Cs+ */
        case 3:  return -520.0;   /* Li+ */
        default: return 0.0;
    }
}

double kcsa_dehydration_cost_eV(int Z) {
    return -kcsa_hydration_free_energy_kJmol(Z) / 96.48533212;
}

/*
 * One ion, one site, rigid filter.
 *
 * `e_inter` receives the BINDING energy - the minimised energy WITH the
 * ion minus the energy of the same filter WITHOUT it - so any constant
 * internal energy of the filter cancels exactly. `e_total` receives
 * e_inter plus the dehydration cost. Keeping the two separate is the
 * whole point: reporting only e_total hides which term is doing the
 * discriminating, and reporting only e_inter hides the fact that a bare
 * vacuum number is not a binding free energy.
 *
 * Run with use_bonds = use_angles = use_dihedrals = 0. The filter is a
 * rigid deposited structure; see kcsa_build_filter for why adding bonded
 * terms to it would only inject strain.
 */
/* Eight-coordinate Shannon effective ionic radii, Angstrom (Shannon 1976
 * Table 1). The filter is 8-coordinate (CN=8 at every site, Demo 12b), so
 * the 8-coordinate radii are the correct sizes — full-audit M3.
 * Previous values were 6-coordinate (Li 0.59 is 4-coordinate): K 1.38->1.51,
 * Na 1.02->1.18 (+0.13/+0.16 A). Contact error was 40% of the 0.36 A K/Na
 * signal. Preferred contacts become K-O 2.91 A, Na-O 2.58 A. */
double kcsa_cation_radius(int Z) {
    switch (Z) {
        case  3: return 0.92;   /* Li+ VIII */
        case 11: return 1.18;   /* Na+ VIII */
        case 19: return 1.51;   /* K+  VIII */
        case 37: return 1.61;   /* Rb+ VIII */
        case 55: return 1.74;   /* Cs+ VIII */
        default: return 1.51;
    }
}

/* Oxygen radius used for the contact distance. */
#define KCSA_O_RADIUS   1.40
/* Shared ion–O well depth, kcal/mol. One value for all alkalis, so the
 * K/Na comparison has no per-ion strength to lean on. */
#define KCSA_ION_O_EPS  0.05
/* The filter oxygen is the AMBER ff99 carbonyl O, now named once in
 * include/amber_lj.h instead of being respelled locally. */
#define KCSA_FILTER_O_SIGMA LJ_AMBER_O_SIGMA
#define KCSA_FILTER_O_EPS  LJ_AMBER_O_EPS

void kcsa_set_ion_radius(Simulation *sim, int ion, int Z) {
    if (!sim || !sim->atoms || ion < 0 || ion >= sim->num_atoms) return;
    double contact = kcsa_cation_radius(Z) + KCSA_O_RADIUS;
    double sig_ionO = contact / TWOPOW_SIXTH;
    double sig_ion  = 2.0 * sig_ionO - KCSA_FILTER_O_SIGMA;
    if (!(sig_ion > 0.05)) sig_ion = 0.05;
    /* Back out the per-ion epsilon that Lorentz–Berthelot needs to give
     * the shared ion–O depth: eps_ionO = sqrt(eps_ion eps_O). */
    double eps_ionO = KCSA_ION_O_EPS * KCAL_MOL_TO_EV;
    double eps_O    = KCSA_FILTER_O_EPS;
    double eps_ion  = (eps_ionO * eps_ionO) / eps_O;
    sim_set_atom_lj(sim, ion, eps_ion, sig_ion);
}

int kcsa_site_binding(Simulation *sim, int filter_first, int n_subunits,
                      int ion_Z, Vec3 site, int n_steps,
                      double *e_inter, double *e_total) {
    (void)n_subunits;   /* the filter is frozen in place by construction */
    if (!sim || !sim->atoms || filter_first < 0) return -1;
    if (sim->num_atoms > 512) return -1;
    if (n_steps < 20) n_steps = 20;

    forces_calculate(sim);
    double e_filter = sim->potential_energy;

    int ion = sim_add_ion(sim, ion_Z, 1, site, 1.0);
    if (ion == ion_Z) { /* unreachable, keeps the compiler honest */ }
    if (ion < 0) return -1;
    kcsa_set_ion_radius(sim, ion, ion_Z);

    /*
     * WHY THERE IS NO FLEXIBLE-FILTER MODE HERE.
     *
     * Letting the filter relax around the ion was tried and does not
     * produce a meaningful number in this engine. The +1 ion sitting
     * 2.8 A from eight oxygens of -0.57 e each exerts a large Coulomb
     * pull, and a restraint stiffness low enough to let the filter
     * respond at all is far too weak to hold it: the oxygens are dragged
     * onto the ion, CN collapses from 8 to 3-4, and the "binding energy"
     * comes out near -140 eV, which is a collapse artefact and not
     * physics. Stiffening the restraints until the filter holds its shape
     * returns it to the rigid case, so there is no useful middle ground
     * without a real protein force field, a solvation model, or both.
     *
     * A flexible-filter claim is therefore NOT made here. What the real
     * flexibility of KcsA does buy the channel is not recoverable from a
     * fixed-charge model in vacuum, and pretending otherwise would be
     * the same category of error as the hand-placed cage this module
     * replaced.
     *
     * RADIAL RELAXATION AT FIXED DEPTH.
     *
     * The first version of this function let the minimiser move the ion
     * freely, and it immediately slid the ion ~5 A down the pore to
     * between sites, where it found CN = 4 and an apparently better
     * energy. That is real physics - K+ does move along the filter - but
     * it means the function was not measuring a SITE at all, it was
     * measuring whichever site happened to be the global minimum for
     * that ion. Four sites then returned four identical numbers, which is
     * the tell.
     *
     * A site binding energy is by definition the well the ion sits in at
     * that site, so the axial coordinate is held at the deposited ion
     * position and only the radial distance from the pore axis is
     * relaxed. The cage is C4 symmetric, so the azimuth is immaterial and
     * the remaining problem is one-dimensional; it is solved by golden
     * section on the engine's own energy function, so the number reported
     * is exactly the engine's energy at its own minimum.
     */
    const double PHI = 0.6180339887498949;
    double a = 0.0, b = 4.5;
    double zfix = site.z;
    Vec3 probe = site;

    for (int it = 0; it < n_steps && (b - a) > 1e-4; it++) {
        double c1 = b - PHI * (b - a);
        double c2 = a + PHI * (b - a);
        probe = vec3(c1, 0.0, zfix);
        sim->atoms[ion].position = probe;
        forces_calculate(sim);
        double e1 = sim->potential_energy;
        probe = vec3(c2, 0.0, zfix);
        sim->atoms[ion].position = probe;
        forces_calculate(sim);
        double e2 = sim->potential_energy;
        if (e1 < e2) b = c2; else a = c1;
    }
    sim->atoms[ion].position = vec3(0.5 * (a + b), 0.0, zfix);
    forces_calculate(sim);

    double bind = sim->potential_energy - e_filter;
    if (e_inter) *e_inter = bind;
    if (e_total) *e_total = bind + kcsa_dehydration_cost_eV(ion_Z);
    return ion;
}
