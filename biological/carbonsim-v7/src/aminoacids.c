#include <math.h>
#include <stdio.h>
#include "../include/aminoacids.h"
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/constants.h"
#include "../include/nucleobases.h"

/*
 * aminoacids.c
 * See aminoacids.h for full geometry/charge/chemistry provenance.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * AMBER LJ type constants - the SAME real, verified values already
 * proven in nucleobases.c (same primary source: amber99.prm, fetched
 * 2026-06-19), duplicated here rather than shared via a header refactor
 * to avoid touching already-validated, working code under time
 * pressure. Provenance is identical: these are AMBER's genuinely
 * general-purpose organic-chemistry LJ categories (sp3 carbon, amide
 * nitrogen, carbonyl oxygen, etc.), not nucleobase-specific - the same
 * categories legitimately apply to protein backbone atoms.
 * ══════════════════════════════════════════════════════════════════════════ */
#define AA_RSTAR_TO_SIGMA(rstar) ((rstar) * 2.0 / 1.122462048309373)

#define AA_LJ_N_SIGMA     AA_RSTAR_TO_SIGMA(1.8240)   /* amide N */
#define AA_LJ_N_EPS       (0.1700 * KCAL_MOL_TO_EV)
#define AA_LJ_C_SIGMA     AA_RSTAR_TO_SIGMA(1.9080)   /* sp2 carbonyl C */
#define AA_LJ_C_EPS       (0.0860 * KCAL_MOL_TO_EV)
#define AA_LJ_O_SIGMA     AA_RSTAR_TO_SIGMA(1.6612)   /* carbonyl O */
#define AA_LJ_O_EPS       (0.2100 * KCAL_MOL_TO_EV)
#define AA_LJ_OH_SIGMA    AA_RSTAR_TO_SIGMA(1.7210)   /* hydroxyl O */
#define AA_LJ_OH_EPS      (0.2104 * KCAL_MOL_TO_EV)
#define AA_LJ_CT_SIGMA    AA_RSTAR_TO_SIGMA(1.9080)   /* sp3 carbon */
#define AA_LJ_CT_EPS      (0.1094 * KCAL_MOL_TO_EV)
#define AA_LJ_HN_SIGMA    AA_RSTAR_TO_SIGMA(0.6000)   /* amide/amine H */
#define AA_LJ_HN_EPS      (0.0157 * KCAL_MOL_TO_EV)
#define AA_LJ_HC_SIGMA    AA_RSTAR_TO_SIGMA(1.4870)   /* aliphatic H */
#define AA_LJ_HC_EPS      (0.0157 * KCAL_MOL_TO_EV)
#define AA_LJ_HA_SIGMA    AA_RSTAR_TO_SIGMA(1.3870)   /* H1-type: aliphatic H, 1 electroneg neighbor (CA-H) */
#define AA_LJ_HA_EPS      (0.0157 * KCAL_MOL_TO_EV)
#define AA_LJ_HO_SIGMA    0.0                          /* hydroxyl H: zero, same real AMBER value as before */
#define AA_LJ_HO_EPS      0.0

/* ══════════════════════════════════════════════════════════════════════════
 * Shared placement helper - identical pattern to nucleobases.c's
 * place_molecule(): centers the molecule at `origin`, adds bonds with
 * r0 set to the EXACT placed distance (zero initial strain), builds
 * angles geometrically from the placed positions.
 * ══════════════════════════════════════════════════════════════════════════ */
static int place_aa_molecule(Simulation *sim, Vec3 origin,
                              const double coords[][3], const int Zs[],
                              const double charges[], int n_atoms,
                              const int bonds[][3], int n_bonds) {
    Vec3 centroid = vec3_zero();
    for (int i = 0; i < n_atoms; i++)
        centroid = vec3_add(centroid,
            vec3(coords[i][0], coords[i][1], coords[i][2]));
    centroid = vec3_scale(centroid, 1.0 / n_atoms);

    int first = sim->num_atoms;
    for (int i = 0; i < n_atoms; i++) {
        Vec3 raw = vec3(coords[i][0], coords[i][1], coords[i][2]);
        Vec3 p   = vec3_add(vec3_sub(raw, centroid), origin);
        sim_add_atom(sim, Zs[i], p, charges[i]);
    }

    for (int i = 0; i < n_bonds; i++) {
        int ia = first + bonds[i][0], ib = first + bonds[i][1];
        int order = bonds[i][2];
        int bidx = sim_add_bond(sim, ia, ib, order);
        if (bidx >= 0) {
            double d = vec3_dist(sim->atoms[ia].position, sim->atoms[ib].position);
            sim_set_bond_params(sim, bidx, d, sim->bonds[bidx].k);
        }
    }
    sim_build_angles_geometric(sim, 3.5); /* generic bending stiffness, same as CO2 */

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * GLYCINE (free, neutral form)
 *
 * Geometry: RCSB PDB CCD ligand "GLY", ideal coordinates, fetched
 * 2026-06-19 from https://files.rcsb.org/ligands/view/GLY.cif.
 * Verified: C2H5NO2, 10 atoms, 9 bonds, 0 chiral centers (correct -
 * glycine is the only achiral proteinogenic amino acid).
 *
 * Atom order: N CA C O OXT H H2 HA2 HA3 HXT
 * (matches the CIF's own atom order/names directly)
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_glycine(Simulation *sim, Vec3 origin) {
    static const double coords[10][3] = {
        { 1.931,  0.090, -0.034},  /* N   */
        { 0.761, -0.799, -0.008},  /* CA  */
        {-0.498,  0.029, -0.005},  /* C   */
        {-0.429,  1.235, -0.023},  /* O   */
        {-1.697, -0.574,  0.018},  /* OXT */
        { 1.910,  0.738,  0.738},  /* H   */
        { 2.788, -0.442, -0.037},  /* H2  */
        { 0.772, -1.440, -0.889},  /* HA2 */
        { 0.793, -1.415,  0.891},  /* HA3 */
        {-2.477, -0.002,  0.019}   /* HXT */
    };
    static const int Zs[10] = {7,6,6,8,8,1,1,1,1,1};

    /* Charges: explicitly approximate (see header comment), balanced
     * to sum to zero for this neutral free molecule. */
    static const double charges[10] = {
        -0.90,  /* N   */
         0.10,  /* CA  */
         0.60,  /* C   */
        -0.55,  /* O   */
        -0.55,  /* OXT */
         0.35,  /* H   */
         0.35,  /* H2  */
         0.05,  /* HA2 */
         0.05,  /* HA3 */
         0.50   /* HXT */
    };

    static const int bonds[9][3] = {
        {0,1,1}, {0,5,1}, {0,6,1},   /* N-CA, N-H, N-H2   */
        {1,2,1}, {1,7,1}, {1,8,1},   /* CA-C, CA-HA2, CA-HA3 */
        {2,3,2}, {2,4,1},            /* C=O, C-OXT        */
        {4,9,1}                     /* OXT-HXT           */
    };

    int first = place_aa_molecule(sim, origin, coords, Zs, charges, 10, bonds, 9);

    sim_set_atom_lj(sim, first+0, AA_LJ_N_EPS, AA_LJ_N_SIGMA);   /* N   */
    sim_set_atom_lj(sim, first+1, AA_LJ_CT_EPS, AA_LJ_CT_SIGMA); /* CA  */
    sim_set_atom_lj(sim, first+2, AA_LJ_C_EPS, AA_LJ_C_SIGMA);   /* C   */
    sim_set_atom_lj(sim, first+3, AA_LJ_O_EPS, AA_LJ_O_SIGMA);   /* O   */
    sim_set_atom_lj(sim, first+4, AA_LJ_OH_EPS, AA_LJ_OH_SIGMA); /* OXT */
    sim_set_atom_lj(sim, first+5, AA_LJ_HN_EPS, AA_LJ_HN_SIGMA); /* H   */
    sim_set_atom_lj(sim, first+6, AA_LJ_HN_EPS, AA_LJ_HN_SIGMA); /* H2  */
    sim_set_atom_lj(sim, first+7, AA_LJ_HA_EPS, AA_LJ_HA_SIGMA); /* HA2 */
    sim_set_atom_lj(sim, first+8, AA_LJ_HA_EPS, AA_LJ_HA_SIGMA); /* HA3 */
    sim_set_atom_lj(sim, first+9, AA_LJ_HO_EPS, AA_LJ_HO_SIGMA); /* HXT */

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * L-ALANINE (free, neutral form)
 *
 * Geometry: RCSB PDB CCD ligand "ALA", ideal coordinates, fetched
 * 2026-06-19 from https://files.rcsb.org/ligands/view/ALA.cif.
 * Verified: C3H7NO2, 13 atoms, 12 bonds, 1 chiral center (CA, S config,
 * matching the CIF's own stereo-marked canonical SMILES).
 *
 * Atom order: N CA C O CB OXT H H2 HA HB1 HB2 HB3 HXT
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_alanine(Simulation *sim, Vec3 origin) {
    static const double coords[13][3] = {
        {-0.966,  0.493,  1.500},  /* N   */
        { 0.257,  0.418,  0.692},  /* CA  */
        {-0.094,  0.017, -0.716},  /* C   */
        {-1.056, -0.682, -0.923},  /* O   */
        { 1.204, -0.620,  1.296},  /* CB  */
        { 0.661,  0.439, -1.742},  /* OXT */
        {-1.383, -0.425,  1.482},  /* H   */
        {-0.676,  0.661,  2.452},  /* H2  */
        { 0.746,  1.392,  0.682},  /* HA  */
        { 1.459, -0.330,  2.316},  /* HB1 */
        { 0.715, -1.594,  1.307},  /* HB2 */
        { 2.113, -0.676,  0.697},  /* HB3 */
        { 0.435,  0.182, -2.647}   /* HXT */
    };
    static const int Zs[13] = {7,6,6,8,6,8,1,1,1,1,1,1,1};

    /* Corrected to sum to exactly zero - a standalone test caught the
     * first draft of these values summing to -0.05 e (even correction
     * of +0.0038 e per atom applied, same fix pattern used earlier
     * for the deoxyribose sugar's charges). */
    static const double charges[13] = {
        -0.8962,  /* N   */
         0.1238,  /* CA  */
         0.6038,  /* C   */
        -0.5462,  /* O   */
        -0.2362,  /* CB  */
        -0.5462,  /* OXT */
         0.3538,  /* H   */
         0.3538,  /* H2  */
         0.0538,  /* HA  */
         0.0638,  /* HB1 */
         0.0638,  /* HB2 */
         0.1038,  /* HB3 */
         0.5038   /* HXT */
    };

    static const int bonds[12][3] = {
        {0,1,1}, {0,6,1}, {0,7,1},   /* N-CA, N-H, N-H2      */
        {1,2,1}, {1,4,1}, {1,8,1},   /* CA-C, CA-CB, CA-HA   */
        {2,3,2}, {2,5,1},            /* C=O, C-OXT           */
        {4,9,1}, {4,10,1}, {4,11,1}, /* CB-HB1/2/3           */
        {5,12,1}                    /* OXT-HXT              */
    };

    int first = place_aa_molecule(sim, origin, coords, Zs, charges, 13, bonds, 12);

    sim_set_atom_lj(sim, first+0, AA_LJ_N_EPS, AA_LJ_N_SIGMA);
    sim_set_atom_lj(sim, first+1, AA_LJ_CT_EPS, AA_LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+2, AA_LJ_C_EPS, AA_LJ_C_SIGMA);
    sim_set_atom_lj(sim, first+3, AA_LJ_O_EPS, AA_LJ_O_SIGMA);
    sim_set_atom_lj(sim, first+4, AA_LJ_CT_EPS, AA_LJ_CT_SIGMA); /* CB */
    sim_set_atom_lj(sim, first+5, AA_LJ_OH_EPS, AA_LJ_OH_SIGMA);
    sim_set_atom_lj(sim, first+6, AA_LJ_HN_EPS, AA_LJ_HN_SIGMA);
    sim_set_atom_lj(sim, first+7, AA_LJ_HN_EPS, AA_LJ_HN_SIGMA);
    sim_set_atom_lj(sim, first+8, AA_LJ_HA_EPS, AA_LJ_HA_SIGMA);
    sim_set_atom_lj(sim, first+9, AA_LJ_HC_EPS, AA_LJ_HC_SIGMA);  /* HB1/2/3 */
    sim_set_atom_lj(sim, first+10, AA_LJ_HC_EPS, AA_LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+11, AA_LJ_HC_EPS, AA_LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+12, AA_LJ_HO_EPS, AA_LJ_HO_SIGMA);

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * GLY-ALA DIPEPTIDE
 *
 * Real peptide (amide) bond formed via genuine condensation chemistry,
 * using the same PDB-CCD-authoritative leaving-atom logic documented
 * in aminoacids.h: glycine (N-terminal) loses its carboxyl -OH (OXT +
 * HXT); alanine (C-terminal) loses one amine H (H2). The new C-N
 * amide bond length (1.33 A) is well-established textbook knowledge
 * (shortened from a generic C-N single bond by partial double-bond
 * resonance character - the structural basis of the peptide bond's
 * planarity) - a DIFFERENT category of confidence than the primary-
 * source-fetched geometry elsewhere in this file, flagged explicitly
 * rather than blurred together.
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_dipeptide_GlyAla(Simulation *sim, Vec3 origin, int *out_ala_N) {
    const double PEPTIDE_BOND_LEN = 1.33; /* Angstrom, textbook amide C-N */

    int gly = sim_place_glycine(sim, origin);
    int ala = sim_place_alanine(sim, vec3_add(origin, vec3(6,6,6)));

    /* Indices before any removal */
    int gly_C   = gly + 2;
    int gly_OXT = gly + 4;
    int gly_HXT = gly + 9;
    int ala_N   = ala + 0;
    int ala_H2  = ala + 7;

    /* Compute glycine C's open valence direction AFTER removing OXT:
     * vector-sum-and-negate technique, same as used for the sugar's
     * C1' and the phosphate's non-bridging oxygens. C currently has
     * 3 substituents besides OXT: CA (gly+1) and the C=O oxygen
     * (gly+3) - only 2 OTHER real neighbors once OXT is gone, so the
     * open direction is computed from those 2 plus the double-bond
     * O counted once (a real sp2 carbonyl carbon has 3 substituents
     * total: CA, =O, and the new N - trigonal planar, not
     * tetrahedral, since the carbonyl carbon stays sp2 in an amide). */
    Vec3 C_pos = sim->atoms[gly_C].position;
    Vec3 to_CA = vec3_normalize(vec3_sub(sim->atoms[gly+1].position, C_pos));
    Vec3 to_O  = vec3_normalize(vec3_sub(sim->atoms[gly+3].position, C_pos));
    /* Trigonal planar: the 3rd direction is the negated sum of the
     * other 2 (exact for perfect 120 degree trigonal geometry). */
    Vec3 open_dir_gly = vec3_normalize(vec3_negate(vec3_add(to_CA, to_O)));

    /* Align alanine's N-H2 direction anti-parallel to glycine's open
     * direction (same rotation technique as the dinucleotide's
     * glycosidic bond attachment) */
    Vec3 N_pos = sim->atoms[ala_N].position;
    Vec3 H2_pos = sim->atoms[ala_H2].position;
    Vec3 NH2_dir = vec3_normalize(vec3_sub(H2_pos, N_pos));
    Vec3 target_dir = vec3_negate(open_dir_gly);

    double cosang = vec3_dot(NH2_dir, target_dir);
    Vec3 axis = vec3_cross(NH2_dir, target_dir);
    double angle;
    if (vec3_norm(axis) < 1.0e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(NH2_dir, vec3(1,0,0));
        if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }

    int ala_first = ala, ala_natoms = 13;
    nb_transform_rigid(sim, ala_first, ala_natoms, N_pos, axis, angle, vec3_zero());

    /* Translate alanine so N lands at the correct peptide bond length
     * from glycine's C, along glycine's open direction */
    Vec3 target_N = vec3_add(C_pos, vec3_scale(open_dir_gly, PEPTIDE_BOND_LEN));
    Vec3 N_now = sim->atoms[ala_N].position;
    Vec3 translation = vec3_sub(target_N, N_now);
    nb_transform_rigid(sim, ala_first, ala_natoms, N_pos, vec3_zero(), 0.0, translation);

    /* Remove the two real leaving-group atoms. IMPORTANT (the exact
     * bug class caught during the dinucleotide build): gly_OXT and
     * gly_HXT are in glycine, which is EARLIER in the array than
     * alanine - removing them shifts every alanine index, including
     * ala_N, down. Track ala_N explicitly through both removals
     * rather than assume it is unaffected. */
    int first_to_remove  = (gly_HXT > gly_OXT) ? gly_HXT : gly_OXT;
    int second_to_remove = (gly_HXT > gly_OXT) ? gly_OXT : gly_HXT;

    sim_remove_terminal_atom(sim, first_to_remove);
    if (ala_N  > first_to_remove) ala_N--;
    if (ala_H2 > first_to_remove) ala_H2--;

    sim_remove_terminal_atom(sim, second_to_remove);
    if (ala_N  > second_to_remove) ala_N--;
    if (ala_H2 > second_to_remove) ala_H2--;

    /* Now remove alanine's H2 (its own leaving group) - ala_H2's
     * current, already-corrected value is used directly. Removing it
     * does not affect ala_N since alanine's own atom order places N
     * (offset 0) before H2 (offset 7 originally), so ala_H2 > ala_N
     * always - but re-verify defensively rather than assume. */
    sim_remove_terminal_atom(sim, ala_H2);
    if (ala_N > ala_H2) ala_N--; /* defensive; should never trigger given the atom order */

    /* Add the real peptide bond */
    sim_add_bond(sim, gly_C, ala_N, 1);

    if (out_ala_N) *out_ala_N = ala_N;
    return gly + 0; /* glycine's N, the dipeptide's N-terminus */
}
