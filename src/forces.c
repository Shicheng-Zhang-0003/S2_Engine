#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/forces.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"
#include "../include/qm.h"

/*
* forces.c
*
* UNITS: lengths in Å, energies in eV, forces in eV/Å. Where a tabulated
* source gives kcal/mol/Å² or kcal/mol/rad², conversion is by
* multiplication with KCAL_MOL_TO_EV.
*
* PROVENANCE (audit fix F2 - this block replaces an earlier header that
* claimed the bond parameters came from AMBER ff14SB / CHARMM36, which
* spot-checking does NOT support):
*
* BOND_TABLE k values are generic spectroscopic-order stiffness
* constants, not force-field-fitted parameters. Spot checks
* (eV/Å² -> kcal/mol/Å², dividing by KCAL_MOL_TO_EV):
*   H-H 36.00 -> 830 kcal/mol/Å²  the real spectroscopic H2 force
*          constant (~575 N/m). AMBER has no H-H bond term at all.
*   C-C 23.80 -> 549 kcal/mol/Å²  ff14SB CT-CT is 310, CHARMM36 ~222.
*          Matches neither.
*   H-O 34.50 -> 796 kcal/mol/Å²  TIP3P/AMBER HW-OW is 553,
*          spectroscopic O-H ~1100. Matches neither.
* Practical impact is limited: every molecule constructor overrides r0
* with the exact placed-geometry distance (zero initial strain), so
* these k's govern bond stiffness and vibration frequency only, never
* equilibrium structure. The claim is corrected anyway, because this
* codebase's standard is "every number sourced, every approximation
* flagged".
*
* ANGLE_TABLE rows are AMBER-sourced where the row's own comment cites
* the specific constant (e.g. TIP3P HW-OW-HW 55, ff14SB H-N3-H 43.1),
* converted with the 2x energy-convention factor documented in the
* table below; untabulated angle types fall back to the flagged
* generic tetrahedral value in forces_angle_params().
*/

/* ══════════════════════════════════════════════════════════════════════════
 * Bond parameter database
 * Ordered: Za <= Zb.  order = 1 (single), 2 (double), 3 (triple).
 * k values in eV/Å²; r0 in Å.
 *
 * PROVENANCE (audit fix F2): k values are generic spectroscopic-order
 * stiffness constants, NOT AMBER/CHARMM-fitted - the spot-check table
 * lives in the file header comment above. r0 entries here are reference
 * lengths only: every molecule constructor overrides r0 with the exact
 * placed-geometry distance (zero initial strain), so only k affects
 * dynamics (bond stiffness and vibration frequency, never equilibrium
 * structure).
 * ══════════════════════════════════════════════════════════════════════════ */
static const BondParam BOND_TABLE[] = {
/*   Za  Zb  ord   r0(Å)    k(eV/Å²)   element pair        */
    { 1,  1,  1,  0.7414,  36.00 },  /* H–H    */
    { 1,  6,  1,  1.0900,  29.30 },  /* H–C    */
    { 1,  7,  1,  1.0120,  31.60 },  /* H–N    */
    { 1,  8,  1,  0.9572,  34.50 },  /* H–O    */
    { 1,  9,  1,  0.9170,  40.00 },  /* H–F    */
    { 1, 16,  1,  1.3400,  17.00 },  /* H–S    */
    { 1, 17,  1,  1.2740,  18.00 },  /* H–Cl   */
    { 6,  6,  1,  1.5400,  23.80 },  /* C–C    */
    { 6,  6,  2,  1.3400,  47.60 },  /* C=C    */
    { 6,  6,  3,  1.2040,  95.20 },  /* C≡C    */
    { 6,  7,  1,  1.4700,  22.50 },  /* C–N    */
    { 6,  7,  2,  1.2740,  50.00 },  /* C=N    */
    { 6,  7,  3,  1.1570, 100.00 },  /* C≡N    */
    { 6,  8,  1,  1.4300,  25.70 },  /* C–O    */
    { 6,  8,  2,  1.2300,  51.40 },  /* C=O    */
    { 6,  9,  1,  1.3500,  23.70 },  /* C–F    */
    { 6, 16,  1,  1.8200,  14.30 },  /* C–S    */
    { 6, 17,  1,  1.7660,  14.30 },  /* C–Cl   */
    { 7,  7,  1,  1.4500,  19.50 },  /* N–N    */
    { 7,  7,  2,  1.2500,  39.00 },  /* N=N    */
    { 7,  7,  3,  1.0980,  78.00 },  /* N≡N    */
    { 7,  8,  1,  1.4400,  20.40 },  /* N–O    */
    { 8,  8,  1,  1.4800,  19.00 },  /* O–O    */
    { 8, 16,  1,  1.5800,  18.00 },  /* O–S    */
    {16, 16,  1,  2.0380,  10.50 },  /* S–S    */
    {15, 15,  1,  2.2100,  10.00 },  /* P–P    */
    { 6, 15,  1,  1.8430,  16.00 },  /* C–P    */
    { 7, 15,  1,  1.6500,  18.00 },  /* N–P    */
    { 8, 15,  1,  1.4810,  26.00 },  /* O–P    */
};
static const int BOND_TABLE_LEN =
    (int)(sizeof(BOND_TABLE) / sizeof(BOND_TABLE[0]));

/* ══════════════════════════════════════════════════════════════════════════
 * Angle parameter database
 * Zb is the central atom.  Za <= Zc for canonical form.
 * k in eV/rad²; theta0 in radians.
 *
 * Degrees to radians: × π/180
 * AMBER k_angle [kcal/mol/rad²] × KCAL_MOL_TO_EV = eV/rad²
 * ══════════════════════════════════════════════════════════════════════════ */
#define DEG2RAD(d) ((d) * 3.14159265358979323846 / 180.0)
static const AngleParam ANGLE_TABLE[] = {
/*  Za  Zb  Zc   theta0(rad)         k(eV/rad²)   source                      */
/*  k = 2 × AMBER_parm [kcal/mol/rad²] × KCAL_MOL_TO_EV                      */
/*  (Code uses V=0.5k(θ-θ0)², AMBER uses V=k(θ-θ0)², so 2× factor needed)   */
    { 1,  8,  1,  DEG2RAD(104.52), 4.770 }, /* H-O-H  TIP3P HW-OW-HW 55 kc  */
    { 1,  7,  1,  DEG2RAD(106.67), 3.738 }, /* H-N-H  H-N3-H ff14SB 43.1 kc  */
    { 1,  6,  1,  DEG2RAD(109.47), 3.035 }, /* H-C-H  HC-CT-HC AMBER 35 kc   */
    { 1,  6,  6,  DEG2RAD(109.47), 4.336 }, /* H-C-C  HC-CT-CT 50 kc         */
    { 1,  6,  7,  DEG2RAD(109.47), 4.336 }, /* H-C-N  HC-CT-N  50 kc         */
    { 1,  6,  8,  DEG2RAD(109.47), 4.336 }, /* H-C-O  HC-CT-OS 50 kc         */
    { 1,  7,  6,  DEG2RAD(118.00), 3.901 }, /* H-N-C  H-N3-CT  45 kc         */
    { 1,  8,  6,  DEG2RAD(108.50), 4.770 }, /* H-O-C  HO-OH-CT 55 kc         */
    { 6,  6,  6,  DEG2RAD(109.47), 3.469 }, /* C-C-C  CT-CT-CT 40 kc         */
    { 6,  6,  7,  DEG2RAD(109.47), 5.464 }, /* C-C-N  CT-CT-N  63 kc         */
    { 6,  6,  8,  DEG2RAD(109.47), 4.336 }, /* C-C-O  CT-CT-OS 50 kc         */
    { 6,  6, 16,  DEG2RAD(114.00), 4.163 }, /* C-C-S  CT-CT-S  48 kc         */
    { 7,  6,  7,  DEG2RAD(116.00), 5.204 }, /* N-C-N  N-C-N    60 kc (est)   */
    { 7,  6,  8,  DEG2RAD(115.00), 5.898 }, /* N-C-O  N-C-O2   68 kc         */
    { 8,  6,  8,  DEG2RAD(123.00), 6.938 }, /* O-C-O  O=C-O2   80 kc         */
    { 6,  7,  6,  DEG2RAD(111.00), 5.638 }, /* C-N-C  CT-N-CT  65 kc         */
    { 6,  8,  6,  DEG2RAD(111.55), 5.204 }, /* C-O-C  CT-OS-CT 60 kc         */
    { 6, 16,  6,  DEG2RAD(102.60), 3.904 }, /* C-S-C  CT-S-CT  45 kc         */
};
static const int ANGLE_TABLE_LEN =
    (int)(sizeof(ANGLE_TABLE) / sizeof(ANGLE_TABLE[0]));

/* ══════════════════════════════════════════════════════════════════════════
 * Bond parameter lookup
 * ══════════════════════════════════════════════════════════════════════════ */
int forces_bond_params(int Za, int Zb, int order, BondParam *out) {
    /* Canonical form: Za <= Zb */
    if (Za > Zb) { int t = Za; Za = Zb; Zb = t; }

    /* Exact match first */
    for (int i = 0; i < BOND_TABLE_LEN; i++) {
        const BondParam *p = &BOND_TABLE[i];
        if (p->Za == Za && p->Zb == Zb && p->order == order) {
            *out = *p;
            return 1;
        }
    }
    /* Fall back to single-bond entry if double/triple not found */
    if (order > 1) {
        for (int i = 0; i < BOND_TABLE_LEN; i++) {
            const BondParam *p = &BOND_TABLE[i];
            if (p->Za == Za && p->Zb == Zb && p->order == 1) {
                static volatile int bond_order_warned[118][118][4] = {{0}};
                if (!bond_order_warned[Za][Zb][order]) {
                    bond_order_warned[Za][Zb][order] = 1;
                    fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for Z%d-Z%d order %d - falling back to single-bond parameters instead of the requested bond order\n", Za, Zb, order);
                }
                *out = *p;
                return 1;
            }
        }
    }

    /* Geometric fallback: use sum of covalent radii, generic k */
    const Element *ea = pt_element(Za);
    const Element *eb = pt_element(Zb);
    if (ea && eb) {
        static volatile int bond_geom_warned[118][118][4] = {{0}};
        if (!bond_geom_warned[Za][Zb][order]) {
            bond_geom_warned[Za][Zb][order] = 1;
            fprintf(stderr, "forces: WARNING (audit fix F3): no BOND_TABLE entry for %s(Z%d)-%s(Z%d) order %d - geometric fallback: r0 from covalent-radii sum, generic k = 20 eV/A^2\n", ea->symbol, Za, eb->symbol, Zb, order);
        }
        out->Za = Za; out->Zb = Zb; out->order = order;
        out->r0 = ea->covalent_radius + eb->covalent_radius;
        out->k  = 20.0;  /* generic, eV/Å² */
        return 1;
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Angle parameter lookup
 * ══════════════════════════════════════════════════════════════════════════ */
int forces_angle_params(int Za, int Zb, int Zc, AngleParam *out) {
    /* Canonical form: Za <= Zc */
    if (Za > Zc) { int t = Za; Za = Zc; Zc = t; }

    for (int i = 0; i < ANGLE_TABLE_LEN; i++) {
        const AngleParam *p = &ANGLE_TABLE[i];
        if (p->Zb == Zb &&
            ((p->Za == Za && p->Zc == Zc) ||
             (p->Za == Zc && p->Zc == Za))) {
            *out = *p;
            return 1;
        }
    }

    /* Generic tetrahedral fallback - k matches the corrected table scale
     * (2x AMBER_parm x KCAL_MOL_TO_EV convention, see table comment above).
     * An earlier version left this at the OLD pre-correction value (0.60)
     * after the table itself was fixed, silently giving any untabulated
     * angle type inconsistent, too-soft physics relative to every
     * tabulated entry - using a representative generic AMBER value
     * (40 kcal/mol/rad^2, the CT-CT-CT constant) here instead. */
    static volatile int angle_warned[118][118][118] = {{0}};
    if (!angle_warned[Za][Zb][Zc]) {
        angle_warned[Za][Zb][Zc] = 1;
        fprintf(stderr, "forces: WARNING (audit fix F3): no ANGLE_TABLE entry for angle Z%d-Z%d-Z%d - generic tetrahedral fallback (109.47 deg, k = 3.469 eV/rad^2)\n", Za, Zb, Zc);
    }
    out->Za = Za; out->Zb = Zb; out->Zc = Zc;
    out->theta0 = DEG2RAD(109.47);
    out->k      = 2.0 * 40.0 * KCAL_MOL_TO_EV;  /* = 3.469 eV/rad^2 */
    return 0;
}

/* ======================================================================
 * Lennard-Jones + Coulomb pair force, with optional smooth switching
 * (audit F5).
 *
 * Physics:
 *   V_LJ  = 4 eps [(sig/r)^12 - (sig/r)^6]
 *   V_C   = COULOMB_MD * qi * qj / (r * dielectric)
 *   F_i = f_total * r_ij,  where r_ij points FROM i TOWARD j.
 *
 * Force sign convention (unchanged from the original): for LJ, f_lj is
 * negative at small r (repulsive, pushes i away from j) and positive
 * beyond the LJ minimum (attractive). For Coulomb, f_c is positive for
 * opposite charges (attractive, pulls i toward j) and negative for like
 * charges (repulsive).
 *
 * SWITCHING (audit F5): when do_switch is set, both the LJ and Coulomb
 * potential and force are multiplied by a smooth switching function S(r)
 * that is 1 for r <= r_switch, falls smoothly to 0 at r_cutoff, and has
 * a continuous first derivative. This removes the energy AND force
 * discontinuity of a hard cutoff. When do_switch is 0, S=1 and dS/dr=0
 * and the original hard-cutoff behaviour is recovered exactly.
 * ====================================================================== */

/* CHARMM-style switching function. Returns S(r) and, via *dSdr, dS/dr.
 *   S(r) = 1                                             for r <= r_switch
 *   S(r) = (roff^2-r^2)^2 (roff^2+2r^2-3ron^2)/(roff^2-ron^2)^3
 *                                                        for ron < r < roff
 *   S(r) = 0                                             for r >= r_cutoff
 * S and dS/dr are continuous everywhere, so both the switched potential
 * and force go smoothly to zero at r_cutoff. */
static double lj_switch_fn(double r, double r_switch, double r_cutoff,
                           double *dSdr) {
    if (r <= r_switch) { if (dSdr) *dSdr = 0.0; return 1.0; }
    if (r >= r_cutoff) { if (dSdr) *dSdr = 0.0; return 0.0; }
    double ron2  = r_switch * r_switch;
    double roff2 = r_cutoff * r_cutoff;
    double r2    = r * r;
    double d     = roff2 - ron2;
    double denom = d * d * d;
    double a     = roff2 - r2;
    double S     = (a * a) * (roff2 + 2.0 * r2 - 3.0 * ron2) / denom;
    if (dSdr) {
        /* dS/dr = -12 r (roff^2 - r^2)(r^2 - ron^2) / (roff^2 - ron^2)^3 */
        *dSdr = -12.0 * r * a * (r2 - ron2) / denom;
    }
    return S;
}

/* Core non-bonded pair interaction with optional switching. When
 * do_switch is 0 this is identical to the original hard-cutoff path.
 * When do_switch is set, the switched force coefficient for each term is
 * f*S + V*(dS/dr)/r, which reduces to f when S=1 and dS/dr=0. */
static PairEnergy pair_nonbonded_core(Atom *atoms, int ia, int ib,
                                      const SimBox *box,
                                      int use_lj, int use_coulomb,
                                      double dielectric,
                                      int do_switch,
                                      double r_switch, double r_cutoff) {
    PairEnergy result = {0.0, 0.0};
    if (!atoms || ia < 0 || ib < 0 || ia == ib) return result;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    Atom *ai = &atoms[ia];
    Atom *bi = &atoms[ib];
    Vec3 r_ij = vec3_sub(bi->position, ai->position);
    if (!isfinite(r_ij.x) || !isfinite(r_ij.y) || !isfinite(r_ij.z)) return result;
    if (box && (box->periodic[0] || box->periodic[1] || box->periodic[2]))
        r_ij = vec3_pbc_box(r_ij, box->dimensions, box->periodic);
    double r2 = vec3_norm2(r_ij);
    if (r2 < 1.0e-10) return result;
    double r = sqrt(r2);
    double f_total = 0.0;

    double S = 1.0, dSdr = 0.0;
    if (do_switch) S = lj_switch_fn(r, r_switch, r_cutoff, &dSdr);

    if (use_lj) {
        double eps   = lj_eps_combine(ai->lj_epsilon, bi->lj_epsilon);
        double sigma = lj_sigma_combine(ai->lj_sigma, bi->lj_sigma);
        double sr2  = (sigma * sigma) / r2;
        double sr6  = sr2 * sr2 * sr2;
        double sr12 = sr6 * sr6;
        double V_lj = 4.0 * eps * (sr12 - sr6);
        double f_lj = (24.0 * eps / r2) * (sr6 - 2.0 * sr12); /* (1/r) dV_lj/dr */
        result.lj_energy = V_lj * S;
        f_total += f_lj * S + V_lj * dSdr / r;
    }
    if (use_coulomb) {
        double qi = ai->partial_charge;
        double qj = bi->partial_charge;
        if (fabs(qi) > 1.0e-9 && fabs(qj) > 1.0e-9) {
            double V_c = COULOMB_MD * qi * qj / (r * dielectric);
            double f_c = -COULOMB_MD * qi * qj / (r2 * r * dielectric); /* (1/r) dV_c/dr */
            result.coulomb_energy = V_c * S;
            f_total += f_c * S + V_c * dSdr / r;
        }
    }

    Vec3 F_i = vec3_scale(r_ij, f_total);
    vec3_iadd(&ai->force, F_i);
    vec3_isub(&bi->force, F_i);
    return result;
}

/* Public, unswitched entry point (used by the base-pairing diagnostics,
 * which inspect specific close pairs and must not be affected by the
 * cutoff switch). Identical to the pre-F5 behaviour. */
PairEnergy forces_nonbonded_pair(Atom *atoms, int ia, int ib,
                                 const SimBox *box,
                                 int use_lj, int use_coulomb,
                                 double dielectric) {
    return pair_nonbonded_core(atoms, ia, ib, box, use_lj, use_coulomb,
                               dielectric, 0, 0.0, 0.0);
}

/* Side-effect-free variant: identical potential, no force writes.
 * Implemented directly (not via save/restore) so diagnostics cannot
 * leak residue even if future re-ordering forgets the zeroing step. */
PairEnergy forces_nonbonded_energy(const Atom *atoms, int ia, int ib,
                                   const SimBox *box,
                                   int use_lj, int use_coulomb,
                                   double dielectric) {
    PairEnergy out = {0.0, 0.0};
    if (!atoms || ia < 0 || ib < 0 || ia == ib) return out;
    if (!(dielectric > 1e-9) || !isfinite(dielectric)) dielectric = 1.0;
    Vec3 r_ij = vec3_sub(atoms[ib].position, atoms[ia].position);
    if (!isfinite(r_ij.x) || !isfinite(r_ij.y) || !isfinite(r_ij.z)) return out;
    if (box && (box->periodic[0] || box->periodic[1] || box->periodic[2]))
        r_ij = vec3_pbc_box(r_ij, box->dimensions, box->periodic);
    double r2 = vec3_norm2(r_ij);
    if (r2 < 1.0e-10) return out;
    double r = sqrt(r2);
    if (use_lj) {
        double eps = lj_eps_combine(atoms[ia].lj_epsilon, atoms[ib].lj_epsilon);
        double sigma = lj_sigma_combine(atoms[ia].lj_sigma, atoms[ib].lj_sigma);
        double sr2 = (sigma * sigma) / r2;
        double sr6 = sr2 * sr2 * sr2;
        out.lj_energy = 4.0 * eps * (sr6 * sr6 - sr6);
    }
    if (use_coulomb) {
        double qi = atoms[ia].partial_charge, qj = atoms[ib].partial_charge;
        if (fabs(qi) > 1.0e-9 && fabs(qj) > 1.0e-9)
            out.coulomb_energy = COULOMB_MD * qi * qj / (r * dielectric);
    }
    return out;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Harmonic bond force
 *
 * V_bond = 0.5 k (r − r0)²
 * dV/dr  = k (r − r0)
 *
 * r_ab = pos_b − pos_a  (from a toward b)
 * F_a  = k(r−r0) × r̂_ab   (toward b when stretched, away when compressed)
 * F_b  = −F_a
 * ══════════════════════════════════════════════════════════════════════════ */
double forces_bond(Atom *atoms, const Bond *bond) {
    if (!atoms || !bond) return 0.0;
    if (bond->atom_a < 0 || bond->atom_b < 0) return 0.0;
    Atom *a = &atoms[bond->atom_a];
    Atom *b = &atoms[bond->atom_b];

    Vec3   r_ab   = vec3_sub(b->position, a->position);
    double r      = vec3_norm(r_ab);
    if (r < 1.0e-10) return 0.0;

    double stretch = r - bond->r0;
    double energy  = 0.5 * bond->k * stretch * stretch;

    /* f_scalar = k × stretch / r  →  F_a = f_scalar × r_ab */
    double f_scalar = bond->k * stretch / r;
    Vec3 F_a = vec3_scale(r_ab, f_scalar);
    vec3_iadd(&a->force, F_a);
    vec3_isub(&b->force, F_a);

    return energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Harmonic angle force
 *
 * Atoms a–b–c, b is the central atom.
 * V_angle = 0.5 k (θ − θ0)²
 *
 * Gradient via chain rule:
 *   e_ba = (r_a − r_b) / d_ba
 *   e_bc = (r_c − r_b) / d_bc
 *   cos θ = e_ba · e_bc
 *
 *   ∂cos θ / ∂r_a = (e_bc − cos θ × e_ba) / d_ba
 *   ∂cos θ / ∂r_c = (e_ba − cos θ × e_bc) / d_bc
 *   ∂cos θ / ∂r_b = −(∂cos θ/∂r_a + ∂cos θ/∂r_c)
 *
 *   dV/dθ = k (θ − θ0)
 *   dθ/d(cosθ) = −1 / sin θ
 *
 *   F_a = −dV/dθ × dθ/d(cosθ) × ∂cosθ/∂r_a
 *        = [k(θ−θ0)/sin θ] × (e_bc − cosθ e_ba) / d_ba
 * ══════════════════════════════════════════════════════════════════════════ */
double forces_angle(Atom *atoms, const Angle *angle) {
    if (!atoms || !angle) return 0.0;
    if (angle->atom_a < 0 || angle->atom_b < 0 || angle->atom_c < 0) return 0.0;
    Atom *a = &atoms[angle->atom_a];
    Atom *b = &atoms[angle->atom_b];
    Atom *c = &atoms[angle->atom_c];

    Vec3 r_ba = vec3_sub(a->position, b->position);
    Vec3 r_bc = vec3_sub(c->position, b->position);

    double d_ba = vec3_norm(r_ba);
    double d_bc = vec3_norm(r_bc);
    if (d_ba < 1.0e-10 || d_bc < 1.0e-10) return 0.0;

    Vec3 e_ba = vec3_scale(r_ba, 1.0 / d_ba);
    Vec3 e_bc = vec3_scale(r_bc, 1.0 / d_bc);

    double cos_theta = vec3_dot(e_ba, e_bc);
    /* Clamp for numerical safety at linear/collapsed angles */
    if (cos_theta >  1.0 - 1.0e-7) cos_theta =  1.0 - 1.0e-7;
    if (cos_theta < -1.0 + 1.0e-7) cos_theta = -1.0 + 1.0e-7;

    double theta   = acos(cos_theta);
    double sin_theta = sin(theta);
    double energy  = 0.5 * angle->k * (theta - angle->theta0)
                                    * (theta - angle->theta0);

    /* Prefactor: k(θ−θ0)/sinθ */
    double pre = angle->k * (theta - angle->theta0) / sin_theta;

    /* ∂cosθ/∂r_a = (e_bc − cosθ e_ba) / d_ba */
    Vec3 dcos_a = vec3_scale(
        vec3_sub(e_bc, vec3_scale(e_ba, cos_theta)),
        1.0 / d_ba);

    /* ∂cosθ/∂r_c = (e_ba − cosθ e_bc) / d_bc */
    Vec3 dcos_c = vec3_scale(
        vec3_sub(e_ba, vec3_scale(e_bc, cos_theta)),
        1.0 / d_bc);

    /* F_x = pre × ∂cosθ/∂r_x  (dθ/dcos already absorbed into pre) */
    Vec3 F_a = vec3_scale(dcos_a, pre);
    Vec3 F_c = vec3_scale(dcos_c, pre);
    /* Newton: F_b = -(F_a + F_c) */
    Vec3 F_b = vec3_negate(vec3_add(F_a, F_c));

    vec3_iadd(&a->force, F_a);
    vec3_iadd(&b->force, F_b);
    vec3_iadd(&c->force, F_c);

    return energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Dihedral (torsion) energy only - used internally by the finite-
 * difference force calculation below, and exposed via forces_dihedral
 * for direct energy queries.
 *
 * V = k * (1 + cos(n*phi - delta))
 * ══════════════════════════════════════════════════════════════════════════ */
static double dihedral_energy_only(const Atom *atoms, const Dihedral *dh) {
    Vec3 pa = atoms[dh->atom_a].position;
    Vec3 pb = atoms[dh->atom_b].position;
    Vec3 pc = atoms[dh->atom_c].position;
    Vec3 pd = atoms[dh->atom_d].position;

    Vec3 b1 = vec3_sub(pb, pa);
    Vec3 b2 = vec3_sub(pc, pb);
    Vec3 b3 = vec3_sub(pd, pc);

    /* Guard against a degenerate (near-zero-length or collinear) case,
     * where the dihedral angle is undefined - contributes no energy
     * rather than risk a NaN from dividing by a near-zero norm inside
     * vec3_dihedral's normalize step. Collinearity checked via plane
     * normals: |b1×b2| or |b2×b3| near zero means undefined torsion. */
    if (vec3_norm(b1) < 1.0e-10 || vec3_norm(b2) < 1.0e-10 ||
        vec3_norm(b3) < 1.0e-10) return 0.0;
    if (vec3_norm(vec3_cross(b1, b2)) < 1.0e-12 ||
        vec3_norm(vec3_cross(b2, b3)) < 1.0e-12) return 0.0;
    if (!isfinite(b1.x) || !isfinite(b2.x) || !isfinite(b3.x)) return 0.0;

    double phi = vec3_dihedral(b1, b2, b3);
    if (!isfinite(phi)) return 0.0;
    return dh->k * (1.0 + cos(dh->n * phi - dh->delta));
}

double forces_dihedral(Atom *atoms, const Dihedral *dh) {
    if (!atoms || !dh) return 0.0;
    if (dh->atom_a < 0 || dh->atom_b < 0 || dh->atom_c < 0 || dh->atom_d < 0) return 0.0;
    /* Analytic torsion gradient (audit P2 closed): exact chain rule
     * through phi = atan2(y, x), x = n1.n2, y = m1.n2, with
     * n1 = b1xb2, n2 = b2xb3, m1 = n1 x b2hat.
     * dphi = (x*dy - y*dx)/(x^2+y^2); F = -dV/dphi * dphi/dr with
     * dV/dphi = -k*n*sin(n*phi-delta). Per-atom dx/dr, dy/dr below are
     * derived term-by-term from dn1 = db1xb2 + b1xdb2 etc. via the
     * cyclic identity (AxB).C = A.(BxC):
     *   a: dx = n2xb2,              dy = (bhat2xn2)xb2
     *   d: dx = n1xb2,              dy = m1xb2
     *   b: dx = (b2xn2)-(n2xb1)+(n1xb3),
     *      dy = (b2xW2)-(W2xb1)+(m1xb3)+(-W3+bhat2(bhat2.W3))/|b2|
     *   c: dx = (n2xb1)+(b3xn1)+(b2xn1),
     *      dy = (W2xb1)+(b3xm1)+(b2xm1)+(W3-bhat2(bhat2.W3))/|b2|
     * with W2 = bhat2xn2, W3 = n2xn1. The m1.dn2 terms (present only
     * for the middle atoms b, c) carry (b3xm1)/(b2xm1) order — the
     * cyclic identity A.(BxC) = B.(CxA), NOT B.(AxC): writing (m1xb3)
     * here was a real sign bug caught by tests/test_forces.c.
     * Translation invariance (sums = 0) holds exactly — asserted.
     * guard (|n| < 1e-12) returns energy with zero forces (FD gives
     * garbage-scale values there; zero is the honest choice). */
    Vec3 pa = atoms[dh->atom_a].position;
    Vec3 pb = atoms[dh->atom_b].position;
    Vec3 pc = atoms[dh->atom_c].position;
    Vec3 pd = atoms[dh->atom_d].position;
    Vec3 b1 = vec3_sub(pb, pa);
    Vec3 b2 = vec3_sub(pc, pb);
    Vec3 b3 = vec3_sub(pd, pc);
    double b2n = vec3_norm(b2);
    if (vec3_norm(b1) < 1.0e-10 || b2n < 1.0e-10 ||
        vec3_norm(b3) < 1.0e-10) return 0.0;
    Vec3 n1 = vec3_cross(b1, b2);
    Vec3 n2 = vec3_cross(b2, b3);
    double n1n = vec3_norm(n1), n2n = vec3_norm(n2);
    Vec3 bhat = vec3_scale(b2, 1.0 / b2n);
    Vec3 m1 = vec3_cross(n1, bhat);
    double x = vec3_dot(n1, n2);
    double y = vec3_dot(m1, n2);
    double denom = x * x + y * y;
    double phi = atan2(y, x);
    if (!isfinite(phi)) return 0.0;
    double sarg = dh->n * phi - dh->delta;
    double energy = dh->k * (1.0 + cos(sarg));
    if (!isfinite(energy)) return 0.0;
    if (n1n < 1.0e-12 || n2n < 1.0e-12) return energy;
    if (!(denom > 0.0) || !isfinite(denom)) return energy;
    double Vp = -dh->k * dh->n * sin(sarg); /* dV/dphi */
    if (!isfinite(Vp)) return energy;

    Vec3 W2 = vec3_cross(bhat, n2);
    Vec3 W3 = vec3_cross(n2, n1);
    double bhW3 = vec3_dot(bhat, W3);

    Vec3 dxa = vec3_cross(n2, b2);
    Vec3 dya = vec3_cross(W2, b2);
    Vec3 dxd = vec3_cross(n1, b2);
    Vec3 dyd = vec3_cross(m1, b2);
    Vec3 dxb = vec3_add(vec3_sub(vec3_cross(b2, n2), vec3_cross(n2, b1)),
                        vec3_cross(n1, b3));
    Vec3 dyb = vec3_add(vec3_add(vec3_sub(vec3_cross(b2, W2), vec3_cross(W2, b1)),
                                 vec3_cross(m1, b3)),
                        vec3_scale(vec3_add(vec3_negate(W3),
                                            vec3_scale(bhat, bhW3)),
                                   1.0 / b2n));
    Vec3 dxc = vec3_add(vec3_add(vec3_cross(n2, b1), vec3_cross(b3, n1)),
                        vec3_cross(b2, n1));
    Vec3 dyc = vec3_add(vec3_add(vec3_cross(W2, b1), vec3_cross(b3, m1)),
                        vec3_add(vec3_cross(b2, m1),
                                 vec3_scale(vec3_sub(W3, vec3_scale(bhat, bhW3)),
                                            1.0 / b2n)));

    Vec3 dpa = vec3_scale(vec3_sub(vec3_scale(dya, x), vec3_scale(dxa, y)),
                          -Vp / denom);
    Vec3 dpb = vec3_scale(vec3_sub(vec3_scale(dyb, x), vec3_scale(dxb, y)),
                          -Vp / denom);
    Vec3 dpc = vec3_scale(vec3_sub(vec3_scale(dyc, x), vec3_scale(dxc, y)),
                          -Vp / denom);
    Vec3 dpd = vec3_scale(vec3_sub(vec3_scale(dyd, x), vec3_scale(dxd, y)),
                          -Vp / denom);
    if (!isfinite(dpa.x + dpa.y + dpa.z + dpb.x + dpb.y + dpb.z +
                  dpc.x + dpc.y + dpc.z + dpd.x + dpd.y + dpd.z))
        return energy;
    vec3_iadd(&atoms[dh->atom_a].force, dpa);
    vec3_iadd(&atoms[dh->atom_b].force, dpb);
    vec3_iadd(&atoms[dh->atom_c].force, dpc);
    vec3_iadd(&atoms[dh->atom_d].force, dpd);
    return energy;
}

/* Finite-difference torsion force (former default, kept as the
 * validation oracle): perturbs each coordinate by +-h and
 * re-evaluates dihedral_energy_only. tests/test_forces.c asserts the
 * analytic forces_dihedral() above agrees with this to 1e-6 on
 * generic, helical, and near-planar geometries. */
double forces_dihedral_fd(Atom *atoms, const Dihedral *dh) {
    if (!atoms || !dh) return 0.0;
    if (dh->atom_a < 0 || dh->atom_b < 0 || dh->atom_c < 0 || dh->atom_d < 0) return 0.0;
    double energy = dihedral_energy_only(atoms, dh);
    if (!isfinite(energy)) return 0.0;

    const double h = 1.0e-5; /* Angstrom */
    int idx[4] = {dh->atom_a, dh->atom_b, dh->atom_c, dh->atom_d};

    for (int a = 0; a < 4; a++) {
        Atom *atom = &atoms[idx[a]];
        double *coords[3] = {&atom->position.x, &atom->position.y, &atom->position.z};

        for (int c = 0; c < 3; c++) {
            double original = *coords[c];

            *coords[c] = original + h;
            double E_plus = dihedral_energy_only(atoms, dh);

            *coords[c] = original - h;
            double E_minus = dihedral_energy_only(atoms, dh);

            *coords[c] = original; /* restore exactly */

            double dVdx = (E_plus - E_minus) / (2.0 * h);
            double force_component = -dVdx;

            if (c == 0) atom->force.x += force_component;
            else if (c == 1) atom->force.y += force_component;
            else atom->force.z += force_component;
        }
    }

    return energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Master force calculation
 * ══════════════════════════════════════════════════════════════════════════ */
void forces_calculate(Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    if (!sim->bonds && sim->num_bonds > 0) return;
    if (!sim->angles && sim->num_angles > 0) return;
    if (!sim->dihedrals && sim->num_dihedrals > 0) return;
    /* v3 SCF: converge charges (with dipole feedback) before forces. */
    if (sim->use_scf) {
        if (sim->scf_pinned_idx >= 0)
            qm_scf_charges(sim, sim->scf_total_q, sim->dielectric,
                           sim->scf_pinned_idx, sim->scf_pinned_q);
        else {
            double qq[128];
            int nn = sim->num_atoms < 128 ? sim->num_atoms : 128;
            if (qm_qeq(sim, sim->scf_total_q, sim->dielectric, qq) == 0)
                for (int i = 0; i < nn; i++)
                    if (isfinite(qq[i])) sim->atoms[i].partial_charge = qq[i];
        }
    }
    int N = sim->num_atoms;

    /* 1. Zero all forces */
    for (int i = 0; i < N; i++)
        sim->atoms[i].force = vec3_zero();

    double E_lj      = 0.0;
    double E_coulomb = 0.0;
    double E_bond    = 0.0;
    double E_angle   = 0.0;
    double E_dihedral = 0.0;
    double E_restraint = 0.0;

    /* ── 1-4 scaling: deliberately NOT applied (audit F1 resolution) ────────
     * Investigated and rejected; kept in-source so the reasoning survives.
     * AMBER scales 1-4 non-bonded pairs (LJ /SCNB with SCNB=2, Coulomb /SCEE
     * with SCEE=1.2) because its torsion parameters were co-fitted WITH that
     * scaling already in place. This force field has no fitted torsion
     * potentials - its dihedrals are harmonic restraints toward fixed
     * textbook angles and carry no 1-4 interaction physics. Scaling the 1-4
     * non-bonded terms here would remove real physics with nothing to
     * compensate, and when tested it destabilised the validated alpha-helix
     * i,i+4 backbone hydrogen bond (Demo 11 failed to form it). 1-4 pairs
     * are therefore deliberately kept at full strength. Re-introduce 1-4
     * scaling ONLY alongside properly fitted torsion parameters, and
     * re-validate Demo 11 when doing so. */

    /* Audit F5 follow-through: LJ/Coulomb cutoff switching is OPT-IN
     * via sim->use_switching (default OFF). The switching function
     * smoothly tapers the non-bonded potential and force to zero
     * between r_switch and r_cutoff, removing a REAL condensed-phase
     * cutoff's discontinuity. Every current demo is gas-phase/vacuum,
     * where the plain unswitched potential is the correct treatment.
     * Keying only on cutoff>2.0 was a bug: with the default cutoff=12,
     * r_switch=9.6 and the switch wrongly attenuated every pair in the
     * 9.6-12 A band in the gas-phase demos (Demos 7 and 11),
     * contradicting the old 'never activates for gas phase' claim and
     * perturbing validated results. Enable use_switching only when a
     * condensed-phase system with a real cutoff is added. */
    double r_cutoff = sim->cutoff;
    double r_switch = 0.8 * r_cutoff;
    int do_switch = (sim->use_switching && r_cutoff > 2.0) ? 1 : 0;
    /* 2. Non-bonded pairs.
     * Exclusion via per-atom bond_partners lists: O(degree) per pair
     * instead of O(B+A) scans. 1-2: j in i's partners. 1-3: share a
     * common bonded neighbor (angle center or either endpoint's
     * partner walk). Degree is bounded by MAX_BONDS_PER_ATOM (8), so
     * per-pair work is O(1). Falls back to exact angle-list check only
     * if topology and partner lists ever diverge (defensive). */
    double cutoff2 = sim->cutoff * sim->cutoff;
    int cutoff_ok = (sim->cutoff > 0.0) && isfinite(sim->cutoff);
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            /* Skip pairs that are bonded (1-2) or angle-related (1-3) */
            int skip = 0;
            const Atom *ai0 = &sim->atoms[i];
            const Atom *aj0 = &sim->atoms[j];
            for (int p = 0; p < ai0->num_bonds; p++) {
                if (ai0->bond_partners[p] == j) { skip = 1; break; }
            }
            if (!skip) {
                /* 1-3: any common bonded neighbor k of i and j. */
                for (int p = 0; p < ai0->num_bonds && !skip; p++) {
                    int k = ai0->bond_partners[p];
                    if (k < 0 || k >= N) continue;
                    const Atom *ak = &sim->atoms[k];
                    for (int q = 0; q < ak->num_bonds; q++) {
                        if (ak->bond_partners[q] == j) { skip = 1; break; }
                    }
                }
                /* Defensive: if partner lists are incomplete, consult
                 * the authoritative angle list (covers rebuilt topologies
                 * where partner sync may lag). */
                if (!skip) {
                    for (int a = 0; a < sim->num_angles; a++) {
                        const Angle *ang = &sim->angles[a];
                        if ((ang->atom_a == i && ang->atom_c == j) ||
                            (ang->atom_a == j && ang->atom_c == i)) {
                            skip = 1; break;
                        }
                    }
                }
                /* Silence unused-variable warning when partner path hits. */
                (void)aj0;
            }
            if (skip) continue;

            /* Distance cutoff check (r2 compare: identical branching
             * for valid cutoffs, skips sqrt for out-of-range pairs;
             * degenerate cutoff (<=0/NaN) disables pairs — garbage in,
             * documented out, never silent wrong physics). */
            Vec3 r_ij = vec3_sub(sim->atoms[j].position,
                                 sim->atoms[i].position);
            if (sim->box.periodic[0] || sim->box.periodic[1] || sim->box.periodic[2])
                r_ij = vec3_pbc_box(r_ij, sim->box.dimensions, sim->box.periodic);

            if (!cutoff_ok || vec3_norm2(r_ij) > cutoff2) continue;

            PairEnergy pe = pair_nonbonded_core(
                sim->atoms, i, j,
                &sim->box,
                sim->use_lj,
                sim->use_coulomb,
                sim->dielectric,
                do_switch, r_switch, r_cutoff);

            E_lj      += pe.lj_energy;
            E_coulomb += pe.coulomb_energy;
        }
    }

    /* 3. Bonded stretches */
    if (sim->use_bonds) {
        for (int b = 0; b < sim->num_bonds; b++)
            E_bond += forces_bond(sim->atoms, &sim->bonds[b]);
    }

    /* 4. Angle bends */
    if (sim->use_angles) {
        for (int a = 0; a < sim->num_angles; a++)
            E_angle += forces_angle(sim->atoms, &sim->angles[a]);
    }

    /* 5. Dihedral torsions */
    if (sim->use_dihedrals) {
        for (int d = 0; d < sim->num_dihedrals; d++)
            E_dihedral += forces_dihedral(sim->atoms, &sim->dihedrals[d]);
    }

    /* 6. Harmonic positional restraints: V = 0.5 k |r - anchor|^2,
     *    F = -k (r - anchor). Conservative, integrated by the same
     *    Verlet step as every term above. */
    for (int r = 0; r < sim->num_restraints; r++) {
        int ia = sim->restraint_atom[r];
        if (ia < 0 || ia >= N) continue; /* defensive: stale index */
        Vec3 disp = vec3_sub(sim->atoms[ia].position,
                             sim->restraint_anchor[r]);
        double k = sim->restraint_k[r];
        if (!isfinite(k) || !isfinite(disp.x)) continue;
        E_restraint += 0.5 * k * vec3_norm2(disp);
        vec3_isub(&sim->atoms[ia].force, vec3_scale(disp, k));
    }

    /* 7. qm v2/v4: induced-dipole polarization (analytic forces).
     * use_pol_scf (coupled dipoles) supersedes use_polar (first-order). */
    double E_polar = 0.0;
    if (sim->use_pol_scf) {
        E_polar = qm_induction_scf_forces(sim, sim->dielectric, NULL);
        if (!isfinite(E_polar)) E_polar = 0.0;
    } else if (sim->use_polar) {
        E_polar = qm_induction_forces(sim, sim->dielectric);
        if (!isfinite(E_polar)) E_polar = 0.0;
    }

    /* 8. qm v2: overlap-Pauli repulsion (FD forces). */
    double E_pauli = 0.0;
    if (sim->use_pauli) {
        E_pauli = qm_pauli_forces(sim);
        if (!isfinite(E_pauli)) E_pauli = 0.0;
    }

    /* 9. v3: QM Slater-Kirkwood damped dispersion (analytic forces). */
    double E_disp = 0.0;
    if (sim->use_disp) {
        E_disp = qm_dispersion_forces(sim);
        if (!isfinite(E_disp)) E_disp = 0.0;
    }

    sim->potential_energy = E_lj + E_coulomb + E_bond + E_angle + E_dihedral
                          + E_restraint + E_polar + E_pauli + E_disp;
    sim->E_lj_total      = E_lj;
    sim->E_coulomb_total = E_coulomb;
    sim->E_restraint_total = E_restraint;
    sim->E_polar_total   = E_polar;
    sim->E_pauli_total   = E_pauli;
    sim->E_disp_total    = E_disp;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Print force summary
 * ══════════════════════════════════════════════════════════════════════════ */
void forces_print_summary(const Simulation *sim) {
    printf("  Force summary (step %llu):\n", (unsigned long long)sim->step);
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        printf("    Atom %3d %-2s  |F|=%8.4f eV/Å  "
               "F=(%8.4f, %8.4f, %8.4f)\n",
               i, a->element->symbol,
               vec3_norm(a->force),
               a->force.x, a->force.y, a->force.z);
    }
    printf("  Potential energy: %.6f eV\n", sim->potential_energy);
}
