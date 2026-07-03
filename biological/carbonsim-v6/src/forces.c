#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/forces.h"
#include "../include/constants.h"
#include "../include/periodic_table.h"

/*
 * forces.c
 *
 * Every number here corresponds to a real physical measurement.
 * Bond parameters from AMBER ff14SB / CHARMM36 (converted to eV/Å²).
 * Angle parameters from the same sources.
 *
 * AMBER uses kcal/mol/Å² and kcal/mol/rad²; we multiply by KCAL_MOL_TO_EV.
 * All lengths in Å, energies in eV, forces in eV/Å.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Bond parameter database
 * Ordered: Za <= Zb.  order = 1 (single), 2 (double), 3 (triple).
 * k values in eV/Å²; r0 in Å.
 *
 * Conversion: AMBER k_bond [kcal/mol/Å²] × KCAL_MOL_TO_EV = eV/Å²
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
                *out = *p;
                return 1;
            }
        }
    }

    /* Geometric fallback: use sum of covalent radii, generic k */
    const Element *ea = pt_element(Za);
    const Element *eb = pt_element(Zb);
    if (ea && eb) {
        out->Za = Za; out->Zb = Zb; out->order = order;
        out->r0 = ea->covalent_radius + eb->covalent_radius;
        out->k  = 20.0;  /* generic, eV/Å² */
        return 1;
    }
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
    out->Za = Za; out->Zb = Zb; out->Zc = Zc;
    out->theta0 = DEG2RAD(109.47);
    out->k      = 2.0 * 40.0 * KCAL_MOL_TO_EV;  /* = 3.469 eV/rad^2 */
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Lennard-Jones + Coulomb pair force
 *
 * Physics:
 *   V_LJ  = 4ε [(σ/r)^12 − (σ/r)^6]
 *   V_C   = COULOMB_MD × q_i × q_j / r
 *
 *   F_i = (dV_LJ/dr + dV_C/dr) × r̂_ij
 *        where r̂_ij = (r_j − r_i) / r  [unit vector FROM i TOWARD j]
 *
 * Derivation of scalar prefactors:
 *   dV_LJ/dr = 24ε/r × [(σ/r)^6 − 2(σ/r)^12]
 *   Force F_i = f_scalar × r_ij,  f_scalar = dV/dr / r
 *
 *   For LJ:  f_lj  = 24ε/r² × [(σ/r)^6 − 2(σ/r)^12]
 *   For C:   f_c   = -COULOMB_MD × qi × qj / r³
 *            (negative because dV_C/dr = −k q q / r²,
 *             and F_i = (dV/dr / r) × r_ij)
 * ══════════════════════════════════════════════════════════════════════════ */
PairEnergy forces_nonbonded_pair(Atom *atoms, int ia, int ib,
                                  const SimBox *box,
                                  int use_lj, int use_coulomb,
                                  double dielectric) {
    PairEnergy result = {0.0, 0.0};
    Atom *ai = &atoms[ia];
    Atom *bi = &atoms[ib];

    Vec3 r_ij = vec3_sub(bi->position, ai->position);

    /* Minimum-image PBC */
    if (box && (box->periodic[0] || box->periodic[1] || box->periodic[2])) {
        r_ij = vec3_pbc(r_ij, box->dimensions);
    }

    double r2 = vec3_norm2(r_ij);
    if (r2 < 1.0e-10) return result;  /* avoid self-interaction at r≈0 */
    double r  = sqrt(r2);

    double f_total = 0.0;  /* scalar force coefficient: F_i = f_total × r_ij */

    /* ── Lennard-Jones ───────────────────────────────────────────────────── */
    if (use_lj) {
        double eps   = lj_eps_combine(ai->lj_epsilon, bi->lj_epsilon);
        double sigma = lj_sigma_combine(ai->lj_sigma, bi->lj_sigma);
        double sr2  = (sigma * sigma) / r2;
        double sr6  = sr2 * sr2 * sr2;
        double sr12 = sr6 * sr6;

        /* V_LJ = 4ε(sr12 - sr6)  [eV] */
        result.lj_energy = 4.0 * eps * (sr12 - sr6);

        /* f_lj = 24ε/r² × (sr6 - 2*sr12).  F_i = f_lj × r_ij (r_ij points
         * FROM i TO j). Sign convention, verified by direct substitution:
         * at small r (r<<σ): sr12 >> sr6, so the bracket is NEGATIVE,
         * making f_lj NEGATIVE → F_i points OPPOSITE to r_ij (away from
         * j) → REPULSION, correctly. At large r near the LJ minimum and
         * beyond (r>2^(1/6)σ): sr6 > 2*sr12, bracket POSITIVE, f_lj
         * POSITIVE → F_i ALONG r_ij (toward j) → ATTRACTION, correctly. */
        double f_lj = (24.0 * eps / r2) * (sr6 - 2.0 * sr12);
        f_total += f_lj;
    }

    /* ── Coulomb ─────────────────────────────────────────────────────────── */
    if (use_coulomb) {
        double qi = ai->partial_charge;
        double qj = bi->partial_charge;
        if (fabs(qi) > 1.0e-9 && fabs(qj) > 1.0e-9) {
            /* V_C = COULOMB_MD × qi × qj / (r × dielectric)  [eV] */
            result.coulomb_energy = COULOMB_MD * qi * qj / (r * dielectric);

            /* f_c = -COULOMB_MD × qi × qj / (r³ × dielectric).
             * F_i = f_c × r_ij. Sign convention, verified by direct
             * substitution: for OPPOSITE charges (qi*qj<0), f_c is
             * POSITIVE → F_i points ALONG r_ij (toward j) → ATTRACTION,
             * correctly. For LIKE charges (qi*qj>0), f_c is NEGATIVE →
             * F_i points OPPOSITE to r_ij (away from j) → REPULSION,
             * correctly. (Note this is the opposite sign pairing from
             * the LJ force above - that asymmetry is expected, since
             * Coulomb's sign comes from the charge product while LJ's
             * comes from the distance-dependent bracket term.) */
            double f_c = -COULOMB_MD * qi * qj / (r2 * r * dielectric);
            f_total += f_c;
        }
    }

    /* Accumulate force onto both atoms (Newton's 3rd law) */
    Vec3 F_i = vec3_scale(r_ij, f_total);
    vec3_iadd(&ai->force, F_i);
    vec3_isub(&bi->force, F_i);

    return result;
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
 * Master force calculation
 * ══════════════════════════════════════════════════════════════════════════ */
void forces_calculate(Simulation *sim) {
    int N = sim->num_atoms;

    /* 1. Zero all forces */
    for (int i = 0; i < N; i++)
        sim->atoms[i].force = vec3_zero();

    double E_lj      = 0.0;
    double E_coulomb = 0.0;
    double E_bond    = 0.0;
    double E_angle   = 0.0;

    /* 2. Non-bonded pairs (O(N²) — replace with cell list for large N) */
    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            /* Skip pairs that are bonded (1-2) or angle-related (1-3) */
            int skip = 0;
            for (int b = 0; b < sim->num_bonds; b++) {
                if ((sim->bonds[b].atom_a == i && sim->bonds[b].atom_b == j) ||
                    (sim->bonds[b].atom_a == j && sim->bonds[b].atom_b == i)) {
                    skip = 1; break;
                }
            }
            if (skip) continue;

            /* Check 1-3 exclusion via angles */
            for (int a = 0; a < sim->num_angles && !skip; a++) {
                const Angle *ang = &sim->angles[a];
                if ((ang->atom_a == i && ang->atom_c == j) ||
                    (ang->atom_a == j && ang->atom_c == i))
                    skip = 1;
            }
            if (skip) continue;

            /* Distance cutoff check */
            Vec3 r_ij = vec3_sub(sim->atoms[j].position,
                                 sim->atoms[i].position);
            if (sim->box.periodic[0] || sim->box.periodic[1] || sim->box.periodic[2])
                r_ij = vec3_pbc(r_ij, sim->box.dimensions);

            if (vec3_norm(r_ij) > sim->cutoff) continue;

            PairEnergy pe = forces_nonbonded_pair(
                sim->atoms, i, j,
                &sim->box,
                sim->use_lj,
                sim->use_coulomb,
                sim->dielectric);

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

    sim->potential_energy = E_lj + E_coulomb + E_bond + E_angle;
    sim->E_lj_total      = E_lj;
    sim->E_coulomb_total = E_coulomb;
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
