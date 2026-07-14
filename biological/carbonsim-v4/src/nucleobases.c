#include <math.h>
#include <stdio.h>
#include "../include/nucleobases.h"
#include "../include/sim.h"
#include "../include/forces.h"

/*
 * nucleobases.c
 * See nucleobases.h for full provenance notes on the geometry data.
 *
 * Every coordinate array below is transcribed directly from the
 * "ideal coordinates" columns of the corresponding RCSB PDB Chemical
 * Component Dictionary entry (fetched 2026-06-19):
 *   uracil   -> https://files.rcsb.org/ligands/view/URA.cif
 *   cytosine -> https://files.rcsb.org/ligands/view/CYT.cif
 *   adenine  -> https://files.rcsb.org/ligands/view/ADE.cif
 *   guanine  -> https://files.rcsb.org/ligands/view/GUN.cif
 * Atom order in each array matches the CIF's own pdbx_ordinal order,
 * named in the comment above each row for traceability.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Shared placement helper
 *
 * Centers the molecule (by centroid of all given atoms) at `origin`,
 * adds bonds with r0 set to the EXACT distance in the placed geometry
 * (zero initial strain by construction) while keeping the generic
 * single/double-bond force constant from the BOND_TABLE, then builds
 * angles with theta0 derived from the same placed geometry.
 * ══════════════════════════════════════════════════════════════════════════ */
static int place_molecule(Simulation *sim, Vec3 origin,
                           const double coords[][3], const int Zs[],
                           const double charges[],
                           int n_atoms,
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
        int ia    = first + bonds[i][0];
        int ib    = first + bonds[i][1];
        int order = bonds[i][2];
        int bidx  = sim_add_bond(sim, ia, ib, order);
        if (bidx >= 0) {
            double d = vec3_dist(sim->atoms[ia].position,
                                  sim->atoms[ib].position);
            sim_set_bond_params(sim, bidx, d, sim->bonds[bidx].k);
        }
    }

    /* Generic aromatic-ring bending stiffness, not independently fitted */
    sim_build_angles_geometric(sim, 1.0);

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * URACIL  (RNA pyrimidine base)
 * Atom order: N1 C2 O2 N3 C4 O4 C5 C6 HN1 HN3 H5 H6
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_uracil(Simulation *sim, Vec3 origin) {
    static const double coords[12][3] = {
        { 0.994,  0.000, -1.183},  /* N1  */
        {-0.349, -0.001, -1.135},  /* C2  */
        {-0.986, -0.001, -2.171},  /* O2  */
        {-1.000,  0.003,  0.043},  /* N3  */
        {-0.308, -0.001,  1.200},  /* C4  */
        {-0.896, -0.001,  2.267},  /* O4  */
        { 1.106, -0.000,  1.164},  /* C5  */
        { 1.733,  0.000, -0.031},  /* C6  */
        { 1.445,  0.000, -2.042},  /* HN1 */
        {-1.969,  0.003,  0.059},  /* HN3 */
        { 1.677, -0.000,  2.081},  /* H5  */
        { 2.812,  0.000, -0.078}   /* H6  */
    };
    static const int Zs[12] = {7,6,8,7,6,8,6,6,1,1,1,1};

    /* {atom_i, atom_j, order}  indices into coords[] above */
    static const int bonds[12][3] = {
        {0,1,1}, {0,7,1}, {0,8,1},   /* N1-C2, N1-C6, N1-HN1   */
        {1,2,2}, {1,3,1},            /* C2=O2, C2-N3           */
        {3,4,1}, {3,9,1},            /* N3-C4, N3-HN3          */
        {4,5,2}, {4,6,1},            /* C4=O4, C4-C5           */
        {6,7,2}, {6,10,1},           /* C5=C6, C5-H5           */
        {7,11,1}                    /* C6-H6                  */
    };

    /*
     * Verified RESP partial charges (Aduri et al., J. Chem. Theory
     * Comput. 2007, 3, 1464 - Table 1, "uridine" column), e units.
     * These are nucleoside-context charges; the glycosidic N1 keeps
     * its nucleoside value (the sugar/H distinction barely affects
     * this specific atom's own charge), and the new explicit HN1
     * (replacing the sugar attachment) is set to +0.118186 - derived
     * from charge neutrality, and independently cross-validated: the
     * same 0.118186 residual falls out of cytosine, adenine, and
     * guanine's base-only charge sums too, matching the paper's own
     * stated sugar-charge restraint used during RESP fitting.
     */
    static const double charges[12] = {
         0.1110,  /* N1  */
         0.4539,  /* C2  */
        -0.5407,  /* O2  */
        -0.3681,  /* N3  */
         0.6022,  /* C4  */
        -0.5652,  /* O4  */
        -0.3135,  /* C5  */
        -0.2320,  /* C6  */
         0.118186,/* HN1 (derived, see above) */
         0.3087,  /* HN3 */
         0.1697,  /* H5  */
         0.2557   /* H6  */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 12, bonds, 12);
    /* Polar (N-attached) H's get a small, NONZERO LJ radius - NOT exactly
     * zero like TIP3P water's hydrogens. Exactly-zero LJ was tried first
     * (over-generalizing the TIP3P precedent) and caused a numerical
     * blowup: with eps=0, LJ combining rules give zero repulsion for
     * ANY pair involving that atom, leaving pure 1/r Coulomb attraction
     * with no floor as a donor H approaches an acceptor - in TIP3P this
     * is masked because the donor O's own real LJ provides an indirect
     * floor before collapse; here it wasn't enough. These small values
     * (sigma=0.5 A, eps=0.001 eV) are a stability choice, not a verified
     * literature parameter - generic UFF H radii would be too large
     * (placing real H-bond distances inside the repulsive wall), and
     * exactly zero is too small (no floor at all). Ring C-H's (H5, H6)
     * are not involved in WC pairing and keep generic UFF values. */
    sim_set_atom_lj(sim, first+8, 0.001, 0.5);  /* HN1 */
    sim_set_atom_lj(sim, first+9, 0.001, 0.5);  /* HN3 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * CYTOSINE  (DNA/RNA pyrimidine base)
 *
 * IMPORTANT CORRECTION: the raw PDB CCD "CYT" ligand entry encodes the
 * minor "imino" tautomer of cytosine - independently verified by tracing
 * its canonical SMILES (NC1=CC=NC(=O)N1) around the ring: the exchangeable
 * hydrogen sits on the ring nitrogen DIRECTLY ADJACENT to the amino-
 * bearing carbon. In the dominant, Watson-Crick-relevant amino-oxo
 * tautomer that actually base-pairs in real DNA/RNA, the exchangeable
 * (glycosidic-position) hydrogen sits on N1 - the nitrogen NOT adjacent
 * to the amino carbon (separated from it by the ring on both sides) -
 * while the nitrogen adjacent to the amino carbon (N3) must be a bare
 * lone-pair acceptor for guanine's N1-H to bond to. Using the raw PDB
 * tautomer as-is would put the hydrogen on the wrong nitrogen for any
 * Watson-Crick pairing calculation.
 *
 * Fix: the heavy-atom ring skeleton (positions of N1,N3,C2,C4,C5,C6,O2,N4)
 * is kept exactly as verified from the PDB ideal coordinates (the ring
 * shape itself does not change meaningfully between tautomers). Only the
 * exchangeable hydrogen is relocated: removed from N3, and a new position
 * on N1 is computed via vector math (external angle bisector of C2-N1-C6,
 * at the standard 1.01 A N-H bond length, in-plane by construction since
 * both ring-neighbor directions already lie in the ring plane) - the
 * same construction technique already used for thymine's methyl group.
 *
 * Atom order: N3 C4 N1 C2 O2 N4 C5 C6 HN1(corrected) HN41 HN42 H5 H6
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_cytosine(Simulation *sim, Vec3 origin) {
    /* Heavy atoms + the 4 hydrogens that are NOT affected by the tautomer
     * fix, exactly as verified from the PDB ideal coordinates */
    static const double heavy[12][3] = {
        {-0.356, -2.061, -3.112},  /* N3   */
        { 0.341, -0.883, -3.202},  /* C4   */
        {-0.784, -1.739, -0.794},  /* N1   */
        {-0.931, -2.520, -1.928},  /* C2   */
        {-1.546, -3.587, -1.909},  /* O2   */
        { 0.871, -0.510, -4.410},  /* N4   */
        { 0.498, -0.109, -2.124},  /* C5   */
        {-0.124, -0.631, -0.895},  /* C6   */
        { 1.766, -0.078, -4.418},  /* HN41 */
        { 0.341, -0.687, -5.232},  /* HN42 */
        { 1.030,  0.832, -2.104},  /* H5   */
        { 0.001, -0.001,  0.000}   /* H6   */
    };

    Vec3 N1 = vec3(heavy[2][0], heavy[2][1], heavy[2][2]);
    Vec3 C2 = vec3(heavy[3][0], heavy[3][1], heavy[3][2]);
    Vec3 C6 = vec3(heavy[7][0], heavy[7][1], heavy[7][2]);

    /* External bisector of the C2-N1-C6 ring angle: negative sum of the
     * two unit vectors from N1 toward its ring neighbours. Both inputs
     * lie in the (already-verified-planar) ring plane, so the result
     * does too - no separate planarity enforcement needed. */
    Vec3 to_C2 = vec3_normalize(vec3_sub(C2, N1));
    Vec3 to_C6 = vec3_normalize(vec3_sub(C6, N1));
    Vec3 outward = vec3_normalize(vec3_negate(vec3_add(to_C2, to_C6)));
    Vec3 HN1_pos = vec3_add(N1, vec3_scale(outward, 1.01)); /* N-H bond length */

    /* Assemble final 13-atom array: heavy atoms (with the old, now-unused
     * N3-adjacent H position dropped) plus the new, correctly-placed N1-H */
    double coords[13][3];
    for (int i = 0; i < 12; i++) {
        coords[i][0] = heavy[i][0]; coords[i][1] = heavy[i][1]; coords[i][2] = heavy[i][2];
    }
    coords[12][0] = HN1_pos.x; coords[12][1] = HN1_pos.y; coords[12][2] = HN1_pos.z;

    /* Index map: 0=N3 1=C4 2=N1 3=C2 4=O2 5=N4 6=C5 7=C6
     *            8=HN41 9=HN42 10=H5 11=H6 12=HN1(new) */
    static const int Zs[13] = {7,6,7,6,8,7,6,6,1,1,1,1,1};

    static const int bonds[13][3] = {
        {0,1,1}, {0,3,1},            /* N3-C4, N3-C2 (N3 now bare - WC acceptor) */
        {1,5,1}, {1,6,2},            /* C4-N4, C4=C5           */
        {2,3,1}, {2,7,2}, {2,12,1},  /* N1-C2, N1=C6, N1-HN1(corrected) */
        {3,4,2},                    /* C2=O2                  */
        {5,8,1}, {5,9,1},            /* N4-HN41, N4-HN42       */
        {6,7,1}, {6,10,1},           /* C5-C6, C5-H5           */
        {7,11,1}                    /* C6-H6                  */
    };

    /*
     * Verified RESP partial charges (Aduri et al. 2007, Table 1,
     * "cytidine" column), e units, mapped to the index order used
     * after the tautomer fix above. N1 keeps its nucleoside value;
     * the new HN1 (replacing the sugar) is +0.118186 from charge
     * neutrality - the SAME residual independently derived for
     * uracil, adenine, and guanine (see sim_place_uracil comment).
     */
    static const double charges[13] = {
        -0.8128,   /* N3   */
         0.9020,   /* C4   */
        -0.2152,   /* N1   */
         0.8867,   /* C2   */
        -0.6560,   /* O2   */
        -0.9919,   /* N4   */
        -0.5972,   /* C5   */
         0.1262,   /* C6   */
         0.4251,   /* HN41 */
         0.4251,   /* HN42 */
         0.2023,   /* H5   */
         0.1875,   /* H6   */
         0.118186  /* HN1 (derived, see above) */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 13, bonds, 13);
    sim_set_atom_lj(sim, first+8,  0.001, 0.5);  /* HN41 */
    sim_set_atom_lj(sim, first+9,  0.001, 0.5);  /* HN42 */
    sim_set_atom_lj(sim, first+12, 0.001, 0.5);  /* HN1 (new) */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * ADENINE  (DNA/RNA purine base - fused 5+6 ring)
 * Atom order: N9 C8 N7 C5 C6 N6 N1 C2 N3 C4 HN9 H8 HN61 HN62 H2
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_adenine(Simulation *sim, Vec3 origin) {
    static const double coords[15][3] = {
        {-0.655,  0.000, -2.079},  /* N9   */
        { 0.680,  0.000, -2.352},  /* C8   */
        { 1.360,  0.000, -1.242},  /* N7   */
        { 0.507,  0.005, -0.190},  /* C5   */
        { 0.660,  0.000,  1.205},  /* C6   */
        { 1.919,  0.000,  1.780},  /* N6   */
        {-0.432,  0.000,  1.962},  /* N1   */
        {-1.637,  0.000,  1.423},  /* C2   */
        {-1.829,  0.000,  0.121},  /* N3   */
        {-0.796,  0.000, -0.715},  /* C4   */
        {-1.374,  0.000, -2.731},  /* HN9  */
        { 1.110, -0.001, -3.342},  /* H8   */
        { 2.012, -0.001,  2.745},  /* HN61 */
        { 2.709, -0.004,  1.217},  /* HN62 */
        {-2.498,  0.000,  2.075}   /* H2   */
    };
    static const int Zs[15] = {7,6,7,6,6,7,7,6,7,6,1,1,1,1,1};

    static const int bonds[16][3] = {
        {0,1,1}, {0,9,1}, {0,10,1},  /* N9-C8, N9-C4, N9-HN9   */
        {1,2,2}, {1,11,1},           /* C8=N7, C8-H8           */
        {2,3,1},                    /* N7-C5                  */
        {3,4,1}, {3,9,2},            /* C5-C6, C5=C4           */
        {4,5,1}, {4,6,2},            /* C6-N6, C6=N1           */
        {5,12,1}, {5,13,1},          /* N6-HN61, N6-HN62       */
        {6,7,1},                    /* N1-C2                  */
        {7,8,2}, {7,14,1},           /* C2=N3, C2-H2           */
        {8,9,1}                     /* N3-C4                  */
    };

    /*
     * Verified RESP charges (Aduri et al. 2007, Table 1, "adenosine"
     * column), e units. N9 keeps its nucleoside value; new HN9
     * (replacing the sugar) is +0.118186 from charge neutrality -
     * same cross-validated residual as the other three bases.
     * NOTE: mapped explicitly by atom name to this file's coords
     * order (N9,C8,N7,C5,C6,N6,N1,C2,N3,C4,...) which differs from
     * the paper's own table listing order (N9,C8,N7,C6,N6,C5,C4,
     * N3,C2,N1,...) - a first pass here copied values in the
     * paper's order without remapping, which silently scrambled
     * charges 3 through 9 onto the wrong atoms. Re-verified name by
     * name below.
     */
    static const double charges[15] = {
         0.0172,   /* N9   */
         0.1299,   /* C8   */
        -0.5850,   /* N7   */
         0.0586,   /* C5   */
         0.7111,   /* C6   */
        -0.9386,   /* N6   */
        -0.7536,   /* N1   */
         0.5741,   /* C2   */
        -0.6835,   /* N3   */
         0.3050,   /* C4   */
         0.118186, /* HN9 (derived) */
         0.1749,   /* H8   */
         0.4125,   /* HN61 */
         0.4125,   /* HN62 */
         0.0467    /* H2   */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 15, bonds, 16);
    sim_set_atom_lj(sim, first+10, 0.001, 0.5);  /* HN9  */
    sim_set_atom_lj(sim, first+12, 0.001, 0.5);  /* HN61 */
    sim_set_atom_lj(sim, first+13, 0.001, 0.5);  /* HN62 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * GUANINE  (DNA/RNA purine base - fused 5+6 ring)
 * Atom order: N9 C8 N7 C5 C6 O6 N1 C2 N2 N3 C4 HN9 H8 HN1 HN21 HN22
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_guanine(Simulation *sim, Vec3 origin) {
    static const double coords[16][3] = {
        { 1.510,  0.000, -1.787},  /* N9   */
        { 0.519,  0.000, -2.725},  /* C8   */
        {-0.642,  0.000, -2.139},  /* N7   */
        {-0.466,  0.000, -0.795},  /* C5   */
        {-1.345,  0.001,  0.313},  /* C6   */
        {-2.554,  0.001,  0.152},  /* O6   */
        {-0.812, -0.004,  1.554},  /* N1   */
        { 0.540,  0.000,  1.723},  /* C2   */
        { 1.053,  0.001,  2.996},  /* N2   */
        { 1.367,  0.000,  0.701},  /* N3   */
        { 0.912,  0.000, -0.556},  /* C4   */
        { 2.464, -0.001, -1.961},  /* HN9  */
        { 0.675, -0.000, -3.794},  /* H8   */
        {-1.395, -0.004,  2.330},  /* HN1  */
        { 2.013,  0.000,  3.131},  /* HN21 */
        { 0.455,  0.005,  3.759}   /* HN22 */
    };
    static const int Zs[16] = {7,6,7,6,6,8,7,6,7,7,6,1,1,1,1,1};

    /*
     * Verified RESP charges (Aduri et al. 2007, Table 1, "guanosine"
     * column), e units, mapped name-by-name to this file's coords
     * order. Note: guanine's N1-H is a SEPARATE, always-present
     * feature of the ring (not the glycosidic position - only N9
     * is, for purines), so its charge (0.3546) comes directly from
     * the table, not derived. Only HN9 (replacing the sugar) uses
     * the +0.118186 charge-neutrality residual.
     * Neutrality check: this set sums to -0.000014 (~0, 4-decimal
     * rounding only) - confirms no atom got mismapped.
     */
    static const double charges[16] = {
         0.0268,   /* N9   */
         0.1066,   /* C8   */
        -0.5575,   /* N7   */
         0.1513,   /* C5   */
         0.5316,   /* C6   */
        -0.5483,   /* O6   */
        -0.5287,   /* N1   */
         0.7191,   /* C2   */
        -0.9044,   /* N2   */
        -0.5959,   /* N3   */
         0.1563,   /* C4   */
         0.118186, /* HN9 (derived) */
         0.1767,   /* H8   */
         0.3546,   /* HN1 (table value, not derived) */
         0.3968,   /* HN21 */
         0.3968    /* HN22 */
    };

    static const int bonds[17][3] = {
        {0,1,1}, {0,10,1}, {0,11,1}, /* N9-C8, N9-C4, N9-HN9   */
        {1,2,2}, {1,12,1},           /* C8=N7, C8-H8           */
        {2,3,1},                    /* N7-C5                  */
        {3,4,1}, {3,10,2},           /* C5-C6, C5=C4           */
        {4,5,2}, {4,6,1},            /* C6=O6, C6-N1           */
        {6,7,1}, {6,13,1},           /* N1-C2, N1-HN1          */
        {7,8,1}, {7,9,2},            /* C2-N2, C2=N3           */
        {8,14,1}, {8,15,1},          /* N2-HN21, N2-HN22       */
        {9,10,1}                    /* N3-C4                  */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 16, bonds, 17);
    sim_set_atom_lj(sim, first+11, 0.001, 0.5);  /* HN9  */
    sim_set_atom_lj(sim, first+13, 0.001, 0.5);  /* HN1  */
    sim_set_atom_lj(sim, first+14, 0.001, 0.5);  /* HN21 */
    sim_set_atom_lj(sim, first+15, 0.001, 0.5);  /* HN22 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * THYMINE  (DNA pyrimidine base = 5-methyluracil)
 *
 * Built from uracil's own verified ring geometry directly: the C5-H5
 * bond direction is taken from uracil's real coordinates, extended to
 * the standard aromatic-ring-to-methyl C-C bond length (1.51 A), and a
 * tetrahedral methyl group is placed there using the same construction
 * already validated in sim_place_ch4 (cos/sin of the 109.47 degree
 * tetrahedral angle, three-fold azimuthal symmetry about the bond axis).
 *
 * This is computed via vector math at runtime below, not by hand -
 * the only hand-transcribed numbers are uracil's own already-verified
 * ring coordinates (same ones used in sim_place_uracil).
 *
 * Atom order: N1 C2 O2 N3 C4 O4 C5 C6 HN1 HN3 H6 CM HM1 HM2 HM3
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_thymine(Simulation *sim, Vec3 origin) {
    /* Uracil ring atoms, identical to sim_place_uracil (H5 omitted - it
     * is replaced by the methyl group; its raw position is kept only
     * to determine the correct bond direction for that replacement) */
    static const double ring[10][3] = {
        { 0.994,  0.000, -1.183},  /* N1  */
        {-0.349, -0.001, -1.135},  /* C2  */
        {-0.986, -0.001, -2.171},  /* O2  */
        {-1.000,  0.003,  0.043},  /* N3  */
        {-0.308, -0.001,  1.200},  /* C4  */
        {-0.896, -0.001,  2.267},  /* O4  */
        { 1.106, -0.000,  1.164},  /* C5  */
        { 1.733,  0.000, -0.031},  /* C6  */
        { 1.445,  0.000, -2.042},  /* HN1 */
        {-1.969,  0.003,  0.059}   /* HN3 */
    };
    static const double H6_raw[3] = { 2.812,  0.000, -0.078};
    static const double H5_raw[3] = { 1.677, -0.000,  2.081}; /* direction only */

    Vec3 C5 = vec3(ring[6][0], ring[6][1], ring[6][2]);
    Vec3 H5 = vec3(H5_raw[0], H5_raw[1], H5_raw[2]);

    /* Direction of the original C5-H5 bond, extended to a C-C single
     * bond length for the methyl substituent (toluene-type ring-to-
     * methyl bond, ~1.51 A) */
    Vec3 u  = vec3_normalize(vec3_sub(H5, C5));
    Vec3 CM = vec3_add(C5, vec3_scale(u, 1.51));

    /* Tetrahedral methyl hydrogens about the C5-CM axis u, exactly the
     * same construction used in sim_place_ch4: pick any vector not
     * parallel to u, cross to get an orthonormal frame, place 3 H's at
     * the tetrahedral angle (cos = -1/3) spaced 120 degrees apart. */
    Vec3 v_seed = (fabs(u.y) < 0.9) ? vec3(0,1,0) : vec3(1,0,0);
    Vec3 v = vec3_normalize(vec3_sub(v_seed, vec3_scale(u, vec3_dot(v_seed, u))));
    Vec3 w = vec3_cross(u, v);

    const double r_CH    = 1.09;
    const double cos_tet = -1.0/3.0;
    const double sin_tet = sqrt(1.0 - cos_tet*cos_tet);
    const double TWO_PI_3 = 2.0 * 3.14159265358979323846 / 3.0;

    Vec3 HM[3];
    for (int k = 0; k < 3; k++) {
        double phi = k * TWO_PI_3;
        Vec3 radial = vec3_add(vec3_scale(v, cos(phi)), vec3_scale(w, sin(phi)));
        Vec3 dir    = vec3_add(vec3_scale(u, cos_tet), vec3_scale(radial, sin_tet));
        HM[k] = vec3_add(CM, vec3_scale(dir, r_CH));
    }

    /* Assemble final 15-atom coordinate set (atom order per header comment) */
    double coords[15][3];
    for (int i = 0; i < 10; i++) {
        coords[i][0] = ring[i][0]; coords[i][1] = ring[i][1]; coords[i][2] = ring[i][2];
    }
    coords[10][0] = H6_raw[0]; coords[10][1] = H6_raw[1]; coords[10][2] = H6_raw[2];
    coords[11][0] = CM.x;      coords[11][1] = CM.y;      coords[11][2] = CM.z;
    coords[12][0] = HM[0].x;   coords[12][1] = HM[0].y;   coords[12][2] = HM[0].z;
    coords[13][0] = HM[1].x;   coords[13][1] = HM[1].y;   coords[13][2] = HM[1].z;
    coords[14][0] = HM[2].x;   coords[14][1] = HM[2].y;   coords[14][2] = HM[2].z;

    static const int Zs[15] = {7,6,8,7,6,8,6,6,1,1,1,6,1,1,1};

    /*
     * CHARGES - LOWER CONFIDENCE than the other three bases.
     *
     * Unlike uracil/cytosine/adenine/guanine, this is NOT a verified
     * RESP fit for thymine specifically - no such source was located
     * for the free base. Instead: the 10 ring atoms shared with
     * uracil keep uracil's verified RESP charges unchanged (N1, C2,
     * O2, N3, C4, O4, C6, H6, HN1, HN3), C5 is shifted to
     * -0.143686 (the exact value needed for overall neutrality,
     * given the methyl group below), and the new methyl group (CM +
     * 3 H) uses generic, standard, self-contained-neutral aliphatic
     * values (CM=-0.18, each H=+0.06 - typical AMBER-GAFF-style
     * aromatic-ring-methyl numbers, not independently verified for
     * THIS molecule). Real force fields refit several nearby ring
     * atoms when adding this substituent; this is a charge-
     * conserving approximation, not a verified fit. Since the WC
     * pairing face (N3-H, O4) is far from the methyl and uses
     * unchanged, verified values, this approximation should have
     * minimal impact specifically on pairing energetics.
     */
    static const double charges[15] = {
         0.1110,    /* N1  */
         0.4539,    /* C2  */
        -0.5407,    /* O2  */
        -0.3681,    /* N3  */
         0.6022,    /* C4  */
        -0.5652,    /* O4  */
        -0.143686,  /* C5  (shifted from uracil's -0.3135, see above) */
        -0.2320,    /* C6  */
         0.118186,  /* HN1 (derived) */
         0.3087,    /* HN3 */
         0.2557,    /* H6  */
        -0.18,      /* CM  (generic) */
         0.06,      /* HM1 (generic) */
         0.06,      /* HM2 (generic) */
         0.06       /* HM3 (generic) */
    };

    static const int bonds[15][3] = {
        {0,1,1}, {0,7,1}, {0,8,1},   /* N1-C2, N1-C6, N1-HN1   */
        {1,2,2}, {1,3,1},            /* C2=O2, C2-N3           */
        {3,4,1}, {3,9,1},            /* N3-C4, N3-HN3          */
        {4,5,2}, {4,6,1},            /* C4=O4, C4-C5           */
        {6,7,2},                    /* C5=C6                  */
        {7,10,1},                   /* C6-H6                  */
        {6,11,1},                   /* C5-CM (new methyl bond)*/
        {11,12,1}, {11,13,1}, {11,14,1}  /* CM-HM1/2/3         */
    };

    return place_molecule(sim, origin, coords, Zs, charges, 15, bonds, 15);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Planarity check: reference plane defined by the first 3 given atoms,
 * report max perpendicular deviation of all remaining atoms from it.
 * ══════════════════════════════════════════════════════════════════════════ */
double nb_planarity_deviation(const Simulation *sim,
                               const int *atom_indices, int n) {
    if (n < 3) return 0.0;

    Vec3 p0 = sim->atoms[atom_indices[0]].position;
    Vec3 p1 = sim->atoms[atom_indices[1]].position;
    Vec3 p2 = sim->atoms[atom_indices[2]].position;

    Vec3 normal = vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));

    double max_dev = 0.0;
    for (int i = 0; i < n; i++) {
        Vec3 d = vec3_sub(sim->atoms[atom_indices[i]].position, p0);
        double dev = fabs(vec3_dot(d, normal));
        if (dev > max_dev) max_dev = dev;
    }
    return max_dev;
}

Vec3 nb_ring_normal(const Simulation *sim, const int *atom_indices) {
    Vec3 p0 = sim->atoms[atom_indices[0]].position;
    Vec3 p1 = sim->atoms[atom_indices[1]].position;
    Vec3 p2 = sim->atoms[atom_indices[2]].position;
    return vec3_normalize(vec3_cross(vec3_sub(p1, p0), vec3_sub(p2, p0)));
}

/*
 * Signed angle (radians) to rotate vector v_from onto vector v_to,
 * both measured as rotation AROUND axis `n` (right-hand rule). Both
 * vectors are projected onto the plane perpendicular to n first, so
 * any small out-of-plane component doesn't bias the result. Used to
 * fix the azimuthal (in-plane) rotation left unconstrained by aligning
 * ring normals alone - without this, two coplanar rings can still end
 * up with substituents clashing into each other.
 */
double nb_signed_inplane_angle(Vec3 v_from, Vec3 v_to, Vec3 n) {
    Vec3 a = vec3_sub(v_from, vec3_scale(n, vec3_dot(v_from, n)));
    Vec3 b = vec3_sub(v_to,   vec3_scale(n, vec3_dot(v_to,   n)));
    a = vec3_normalize(a);
    b = vec3_normalize(b);
    double s = vec3_dot(vec3_cross(a, b), n);
    double c = vec3_dot(a, b);
    return atan2(s, c);
}

void nb_transform_rigid(Simulation *sim, int first_atom, int n_atoms,
                         Vec3 pivot, Vec3 axis, double angle,
                         Vec3 translation) {
    for (int i = first_atom; i < first_atom + n_atoms; i++) {
        Vec3 rel = vec3_sub(sim->atoms[i].position, pivot);
        Vec3 rotated = (vec3_norm(axis) > 1.0e-12)
                       ? vec3_rotate_axis_angle(rel, axis, angle)
                       : rel;
        sim->atoms[i].position = vec3_add(vec3_add(pivot, rotated), translation);
    }
}
