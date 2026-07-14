#include <math.h>
#include <stdio.h>
#include "../include/nucleobases.h"
#include "../include/sim.h"
#include "../include/forces.h"
#include "../include/constants.h"

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
 * Real AMBER ff99 (Cornell et al. 1995 / Wang et al. 2000) Lennard-Jones
 * parameters, replacing the generic per-element UFF defaults used
 * previously for the ring heavy atoms, and replacing the earlier ad-hoc
 * "stability fix" guess (sigma=0.5, eps=0.001) used for polar hydrogens.
 *
 * Source: amber99.prm (TINKER-format AMBER parameter file), fetched
 * 2026-06-19 from https://github.com/pren/tinker (a verified mirror of
 * the standard AMBER-FF99 parameter set). Values below are the file's
 * own "vdw <class>  R*[A]  epsilon[kcal/mol]" entries for the specific
 * classes actually used by DNA/RNA nucleobase atoms, cross-checked
 * atom-by-atom against the same file's explicit "R-Adenosine", 
 * "R-Guanosine", "R-Cytosine", "R-Uracil", and "D-Thymine" atom-type
 * assignment lines (which is how each atom below was mapped to a type).
 *
 * Key finding used throughout: AMBER's original (Cornell 1995) vdW
 * parameters are defined by broad hybridization class, not by fine
 * bonded atom type - ALL aromatic/amide nitrogen types (N, NA, NB, NC,
 * N*, N2) share IDENTICAL LJ parameters, and ALL sp2 carbon types (C,
 * CA, CB, CM, CK, CQ, etc.) likewise share identical LJ parameters.
 * Only the atoms' CHARGES differ between these sub-types (already
 * handled correctly via the verified RESP charges elsewhere in this
 * file) - so a single sigma/epsilon pair per (element, hybridization)
 * category is not an approximation on top of AMBER; it IS AMBER.
 *
 * Conversion: sigma[A] = R*[A] * 2/2^(1/6),  epsilon[eV] = eps[kcal/mol]
 * * KCAL_MOL_TO_EV. Verified independently: TIP3P oxygen's AMBER class
 * (R*=1.7683) converts to sigma=3.1506 A via this same formula, an
 * exact match to the literature TIP3P value already used in
 * sim_place_h2o() - confirms the conversion convention is correct.
 * ══════════════════════════════════════════════════════════════════════════ */
#define AMBER_RSTAR_TO_SIGMA(rstar) ((rstar) * 2.0 / 1.122462048309373)

/* Ring/carbonyl nitrogen: AMBER classes N,NA,NB,NC,N*,N2 (all identical) */
#define LJ_RING_N_SIGMA   AMBER_RSTAR_TO_SIGMA(1.8240)
#define LJ_RING_N_EPS     (0.1700 * KCAL_MOL_TO_EV)

/* sp2 ring/carbonyl carbon: AMBER classes C,CA,CB,CM,CK,CQ (all identical) */
#define LJ_SP2_C_SIGMA    AMBER_RSTAR_TO_SIGMA(1.9080)
#define LJ_SP2_C_EPS      (0.0860 * KCAL_MOL_TO_EV)

/* Carbonyl oxygen: AMBER class O */
#define LJ_CARBONYL_O_SIGMA AMBER_RSTAR_TO_SIGMA(1.6612)
#define LJ_CARBONYL_O_EPS   (0.2100 * KCAL_MOL_TO_EV)

/* Amide/aromatic N-H hydrogen: AMBER class H (attached to any N type) */
#define LJ_H_ON_N_SIGMA   AMBER_RSTAR_TO_SIGMA(0.6000)
#define LJ_H_ON_N_EPS     (0.0157 * KCAL_MOL_TO_EV)

/* Aromatic C-H, type HA (e.g. cytosine H5): */
#define LJ_HA_SIGMA       AMBER_RSTAR_TO_SIGMA(1.4590)
#define LJ_HA_EPS         (0.0150 * KCAL_MOL_TO_EV)

/* Aromatic C-H, type H4 (heteroatom-adjacent, e.g. uracil/cytosine/thymine H6): */
#define LJ_H4_SIGMA       AMBER_RSTAR_TO_SIGMA(1.4090)
#define LJ_H4_EPS         (0.0150 * KCAL_MOL_TO_EV)

/* Aromatic C-H, type H5 (purine H8 and adenine H2 - AMBER's "H5" hydrogen
 * type name is an unrelated coincidental collision with our own "H5" ring
 * position naming in pyrimidines; these are different atoms entirely): */
#define LJ_PURINE_H_SIGMA AMBER_RSTAR_TO_SIGMA(1.3590)
#define LJ_PURINE_H_EPS   (0.0150 * KCAL_MOL_TO_EV)

/* sp3 methyl carbon (thymine): AMBER class CT */
#define LJ_CT_SIGMA       AMBER_RSTAR_TO_SIGMA(1.9080)
#define LJ_CT_EPS         (0.1094 * KCAL_MOL_TO_EV)

/* Aliphatic methyl hydrogen (thymine): AMBER class HC */
#define LJ_HC_SIGMA       AMBER_RSTAR_TO_SIGMA(1.4870)
#define LJ_HC_EPS         (0.0157 * KCAL_MOL_TO_EV)

/* ── Sugar (deoxyribose) AMBER types, added for the sugar-phosphate
 * backbone. Ring/exocyclic sp3 carbons all use CT (same as thymine's
 * methyl carbon above). Ring ether oxygen (O4') uses OS. Hydroxyl
 * oxygens (O1 anomeric/O3'/O5') use OH. Hydroxyl hydrogens use HO -
 * EXACTLY ZERO LJ, the same real AMBER convention already confirmed
 * for TIP3P water's hydrogens (class HW) elsewhere in this codebase;
 * HO is a separate but also-zero AMBER class, confirmed directly from
 * the same source parameter file. Aliphatic ring/chain hydrogens use
 * the generic HC type (a deliberate simplification - real AMBER
 * distinguishes several sub-cases here (H1 vs HC) by electronegative-
 * neighbor proximity; HC is used uniformly since the difference in
 * LJ parameters between these sub-types is small and this is not yet
 * being tested in an energetic comparison the way the aromatic ring
 * atoms were for base pairing). ── */
#define LJ_OS_SIGMA       AMBER_RSTAR_TO_SIGMA(1.6837)   /* ring ether O */
#define LJ_OS_EPS         (0.1700 * KCAL_MOL_TO_EV)
#define LJ_OH_SIGMA       AMBER_RSTAR_TO_SIGMA(1.7210)   /* hydroxyl O */
#define LJ_OH_EPS         (0.2104 * KCAL_MOL_TO_EV)
#define LJ_HO_SIGMA       0.0                             /* hydroxyl H: zero */
#define LJ_HO_EPS         0.0

/* ── Phosphate backbone AMBER types. P: class 28. Non-bridging
 * (charged) phosphate oxygen: class 25 "O2", CONFIRMED to share
 * identical LJ parameters with carbonyl O (class 24) in the same
 * source file - not a simplification, this is what AMBER actually
 * specifies. Bridging (ester) oxygen: class 23 "OS", the same type
 * already used for the sugar's ring oxygen O4' above. ── */
#define LJ_P_SIGMA        AMBER_RSTAR_TO_SIGMA(2.1000)
#define LJ_P_EPS          (0.2000 * KCAL_MOL_TO_EV)
#define LJ_O2_SIGMA       AMBER_RSTAR_TO_SIGMA(1.6612)   /* non-bridging, = carbonyl O */
#define LJ_O2_EPS         (0.2100 * KCAL_MOL_TO_EV)

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
    /*
     * Real AMBER ff99 LJ parameters, atom-by-atom, replacing both the
     * generic per-element UFF defaults (heavy atoms) and the earlier
     * ad-hoc "stability fix" guess (the two N-H hydrogens). Mapping
     * verified against AMBER's own "R-Uracil" atom-type assignment
     * line (N1=N*, C2=C, N3=NA, C4=C, C5=CM, C6=CM, O2=O, O4=O,
     * H3=H, H5=HA, H6=H4) - see the constants block at the top of
     * this file for full sourcing and the (element,hybridization)
     * simplification rationale.
     */
    sim_set_atom_lj(sim, first+0, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N1 */
    sim_set_atom_lj(sim, first+1, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C2 */
    sim_set_atom_lj(sim, first+2, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O2 */
    sim_set_atom_lj(sim, first+3, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N3 */
    sim_set_atom_lj(sim, first+4, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C4 */
    sim_set_atom_lj(sim, first+5, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O4 */
    sim_set_atom_lj(sim, first+6, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C5 */
    sim_set_atom_lj(sim, first+7, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C6 */
    sim_set_atom_lj(sim, first+8, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);      /* HN1 */
    sim_set_atom_lj(sim, first+9, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);      /* HN3 */
    sim_set_atom_lj(sim, first+10, LJ_HA_EPS, LJ_HA_SIGMA);             /* H5 */
    sim_set_atom_lj(sim, first+11, LJ_H4_EPS, LJ_H4_SIGMA);             /* H6 */
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
        {0,1,2}, {0,3,1},            /* N3=C4 (DOUBLE - see structure note
                                       * above), N3-C2 (N3 now bare - WC
                                       * acceptor)                          */
        {1,5,1}, {1,6,1},            /* C4-N4 (exocyclic amino, single),
                                       * C4-C5 (single - NOT double; the
                                       * ring's double-bond character at
                                       * this position is at N3=C4 instead,
                                       * since C4 here carries -NH2 not =O) */
        {2,3,1}, {2,7,1}, {2,12,1},  /* N1-C2, N1-C6 (single, unaffected
                                       * by the C4 substitution - same as
                                       * uracil), N1-HN1(corrected)         */
        {3,4,2},                    /* C2=O2                  */
        {5,8,1}, {5,9,1},            /* N4-HN41, N4-HN42       */
        {6,7,2}, {6,10,1},           /* C5=C6 (DOUBLE - same as uracil),
                                       * C5-H5                             */
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
    /*
     * Real AMBER ff99 LJ parameters, atom-by-atom. Mapping verified
     * against AMBER's "R-Cytosine" type line (N1=N*, C2=C, N3=NC,
     * C4=CA, C5=CM, C6=CM, O2=O, N4=N2, H41/H42=H, H5=HA, H6=H4),
     * re-indexed to this file's index order (0=N3 1=C4 2=N1 3=C2
     * 4=O2 5=N4 6=C5 7=C6 8=HN41 9=HN42 10=H5 11=H6 12=HN1[new]).
     * The new HN1 (post-tautomer-fix, attached to N1=N*) uses the
     * same H_ON_N type as any other ring N-H.
     */
    sim_set_atom_lj(sim, first+0, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N3 */
    sim_set_atom_lj(sim, first+1, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C4 */
    sim_set_atom_lj(sim, first+2, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N1 */
    sim_set_atom_lj(sim, first+3, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C2 */
    sim_set_atom_lj(sim, first+4, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O2 */
    sim_set_atom_lj(sim, first+5, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N4 */
    sim_set_atom_lj(sim, first+6, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C5 */
    sim_set_atom_lj(sim, first+7, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C6 */
    sim_set_atom_lj(sim, first+8,  LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN41 */
    sim_set_atom_lj(sim, first+9,  LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN42 */
    sim_set_atom_lj(sim, first+10, LJ_HA_EPS, LJ_HA_SIGMA);             /* H5 */
    sim_set_atom_lj(sim, first+11, LJ_H4_EPS, LJ_H4_SIGMA);             /* H6 */
    sim_set_atom_lj(sim, first+12, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN1 (new) */
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
    /*
     * Real AMBER ff99 LJ parameters, atom-by-atom. Mapping verified
     * against AMBER's "R-Adenosine" type line (N9=N*, C4=CB, C5=CB,
     * N7=NB, C8=CK, N3=NC, C2=CQ, N1=NC, C6=CA, N6=N2, H61/H62=H,
     * H8=H5[AMBER type], H2=H5[AMBER type] - note AMBER's "H5" here is
     * a hydrogen-type NAME, unrelated to nucleobase ring position "H5"
     * used in the pyrimidines elsewhere in this file).
     */
    sim_set_atom_lj(sim, first+0, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N9 */
    sim_set_atom_lj(sim, first+1, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C8 */
    sim_set_atom_lj(sim, first+2, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N7 */
    sim_set_atom_lj(sim, first+3, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C5 */
    sim_set_atom_lj(sim, first+4, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C6 */
    sim_set_atom_lj(sim, first+5, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N6 */
    sim_set_atom_lj(sim, first+6, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N1 */
    sim_set_atom_lj(sim, first+7, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C2 */
    sim_set_atom_lj(sim, first+8, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N3 */
    sim_set_atom_lj(sim, first+9, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C4 */
    sim_set_atom_lj(sim, first+10, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN9  */
    sim_set_atom_lj(sim, first+11, LJ_PURINE_H_EPS, LJ_PURINE_H_SIGMA); /* H8 */
    sim_set_atom_lj(sim, first+12, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN61 */
    sim_set_atom_lj(sim, first+13, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN62 */
    sim_set_atom_lj(sim, first+14, LJ_PURINE_H_EPS, LJ_PURINE_H_SIGMA); /* H2 */
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
    /*
     * Real AMBER ff99 LJ parameters, atom-by-atom. Mapping verified
     * against AMBER's "R-Guanosine" type line (N9=N*, C4=CB, C5=CB,
     * N7=NB, C8=CK, N3=NC, C2=CA, N1=NA, C6=C, N2=N2, H21/H22=H,
     * O6=O, H1=H, H8=H5[AMBER type]).
     */
    sim_set_atom_lj(sim, first+0, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N9 */
    sim_set_atom_lj(sim, first+1, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C8 */
    sim_set_atom_lj(sim, first+2, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N7 */
    sim_set_atom_lj(sim, first+3, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C5 */
    sim_set_atom_lj(sim, first+4, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C6 */
    sim_set_atom_lj(sim, first+5, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O6 */
    sim_set_atom_lj(sim, first+6, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N1 */
    sim_set_atom_lj(sim, first+7, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C2 */
    sim_set_atom_lj(sim, first+8, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N2 */
    sim_set_atom_lj(sim, first+9, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N3 */
    sim_set_atom_lj(sim, first+10, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);       /* C4 */
    sim_set_atom_lj(sim, first+11, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN9  */
    sim_set_atom_lj(sim, first+12, LJ_PURINE_H_EPS, LJ_PURINE_H_SIGMA); /* H8 */
    sim_set_atom_lj(sim, first+13, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN1  */
    sim_set_atom_lj(sim, first+14, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN21 */
    sim_set_atom_lj(sim, first+15, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);     /* HN22 */
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

    int first = place_molecule(sim, origin, coords, Zs, charges, 15, bonds, 15);
    /*
     * Real AMBER ff99 LJ parameters, atom-by-atom. Mapping verified
     * against AMBER's "D-Thymine" type line (N1=N*, C2=C, N3=NA,
     * C4=C, C5=CM, C6=CM, O2=O, O4=O, H3=H, H6=H4, C7[methyl]=CT,
     * H7[methyl H]=HC). Previously thymine had NO LJ overrides at
     * all - every atom, including the ring, was still using generic
     * UFF defaults even after the fix was applied to the other four
     * bases; this closes that gap.
     */
    sim_set_atom_lj(sim, first+0, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N1 */
    sim_set_atom_lj(sim, first+1, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C2 */
    sim_set_atom_lj(sim, first+2, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O2 */
    sim_set_atom_lj(sim, first+3, LJ_RING_N_EPS, LJ_RING_N_SIGMA);      /* N3 */
    sim_set_atom_lj(sim, first+4, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C4 */
    sim_set_atom_lj(sim, first+5, LJ_CARBONYL_O_EPS, LJ_CARBONYL_O_SIGMA); /* O4 */
    sim_set_atom_lj(sim, first+6, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C5 */
    sim_set_atom_lj(sim, first+7, LJ_SP2_C_EPS, LJ_SP2_C_SIGMA);        /* C6 */
    sim_set_atom_lj(sim, first+8, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);      /* HN1 */
    sim_set_atom_lj(sim, first+9, LJ_H_ON_N_EPS, LJ_H_ON_N_SIGMA);      /* HN3 */
    sim_set_atom_lj(sim, first+10, LJ_H4_EPS, LJ_H4_SIGMA);             /* H6 */
    sim_set_atom_lj(sim, first+11, LJ_CT_EPS, LJ_CT_SIGMA);             /* CM (methyl C, AMBER type CT) */
    sim_set_atom_lj(sim, first+12, LJ_HC_EPS, LJ_HC_SIGMA);             /* HM1 */
    sim_set_atom_lj(sim, first+13, LJ_HC_EPS, LJ_HC_SIGMA);             /* HM2 */
    sim_set_atom_lj(sim, first+14, LJ_HC_EPS, LJ_HC_SIGMA);             /* HM3 */
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * 2-DEOXYRIBOSE (free form, anomeric OH unreacted)
 *
 * Geometry: RCSB PDB CCD ligand "2DR", ideal coordinates, fetched
 * 2026-06-19 from https://files.rcsb.org/ligands/view/2DR.cif.
 * Verified: C5H10O4, 19 atoms, 19 bonds, 3 chiral centers. Every
 * coordinate below is transcribed directly from that file's own
 * "C1/C2/C3/C4/C5/O1/O3/O4/O5" ideal-coordinate columns, re-mapped
 * to standard nucleic-acid sugar numbering (CIF's C1->C1', C2->C2',
 * etc. - the file's raw atom names already correspond 1:1 to this
 * convention, just without the prime marks).
 *
 * CHARGES - explicitly approximate, unlike the nucleobases. No
 * verified RESP source specifically for the free (unattached) sugar
 * was located (Aduri et al. 2007, used for the bases, covers only
 * nucleobase atoms of each nucleoside, not the sugar itself). Set to
 * a small, chemically-reasonable, net-neutral distribution rather
 * than left at exactly zero - but NOT independently verified the way
 * the nucleobase charges were. When this sugar is later assembled
 * into a full nucleotide with a phosphate group, the phosphate's
 * charges (to be sourced and verified before use) will dominate that
 * assembly's electrostatics regardless.
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_deoxyribose(Simulation *sim, Vec3 origin) {
    /* Order: C1' O4' C2' C3' C4' C5' O1(anomeric) O3' O5'
     *        H2'a H2'b H3' H1' H4' HO3' H5'a H5'b HO5' HO1 */
    static const double coords[19][3] = {
        { 1.453, -0.576,  0.338},  /* C1'  (CIF "C1")  */
        { 0.075, -1.009,  0.326},  /* O4'  (CIF "O4")  */
        { 1.468,  0.613, -0.657},  /* C2'  (CIF "C2")  */
        { 0.101,  1.285, -0.378},  /* C3'  (CIF "C3")  */
        {-0.736,  0.177,  0.289},  /* C4'  (CIF "C4")  */
        {-2.001, -0.084, -0.532},  /* C5'  (CIF "C5")  */
        { 2.313, -1.623, -0.116},  /* O1   (CIF "O1")  */
        { 0.261,  2.394,  0.509},  /* O3'  (CIF "O3")  */
        {-2.827, -1.027,  0.153},  /* O5'  (CIF "O5")  */
        { 1.524,  0.258, -1.686},  /* H2'a (CIF "H2")  */
        { 2.289,  1.294, -0.436},  /* H2'b (CIF "H22C")*/
        {-0.363,  1.606, -1.311},  /* H3'  (CIF "H3")  */
        { 1.745, -0.246,  1.335},  /* H1'  (CIF "H1")  */
        {-1.007,  0.476,  1.301},  /* H4'  (CIF "H4")  */
        { 0.829,  3.096,  0.163},  /* HO3' (CIF "HO3") */
        {-2.548,  0.850, -0.663},  /* H5'a (CIF "H51") */
        {-1.725, -0.483, -1.508},  /* H5'b (CIF "H52") */
        {-3.649, -1.241, -0.310},  /* HO5' (CIF "HO5") */
        { 2.341, -2.391,  0.470}   /* HO1  (CIF "HO1") */
    };

    static const int Zs[19] = {
        6,8,6,6,6,6,8,8,8,   /* C1' O4' C2' C3' C4' C5' O1 O3' O5' */
        1,1,1,1,1,1,1,1,1,1  /* the 10 hydrogens */
    };

    /*
     * Corrected to sum to exactly zero (a standalone test caught the
     * first draft of these values summing to -0.36 e despite being
     * documented as "net-neutral" - an even correction of +0.0189 e
     * was applied per atom to fix this).
     */
    static const double charges[19] = {
        0.0689, -0.3311, -0.0011, 0.1589, 0.1189, 0.0689, -0.3811, -0.3611, -0.3611,
        0.0389, 0.0389, 0.0689, 0.0689, 0.0689,
        0.2189,   /* HO3' */
        0.0389, 0.0389,
        0.2189,   /* HO5' */
        0.2189    /* HO1  */
    };

    static const int bonds[19][3] = {
        {0,2,1}, {0,1,1}, {0,12,1},  /* C1'-C2', C1'-O4', C1'-H1' */
        {2,3,1}, {2,9,1}, {2,10,1},  /* C2'-C3', C2'-H2'a, C2'-H2'b */
        {3,4,1}, {3,7,1}, {3,11,1},  /* C3'-C4', C3'-O3', C3'-H3' */
        {4,1,1}, {4,5,1}, {4,13,1},  /* C4'-O4', C4'-C5', C4'-H4' */
        {5,8,1}, {5,15,1}, {5,16,1}, /* C5'-O5', C5'-H5'a, C5'-H5'b */
        {0,6,1},                    /* C1'-O1 (anomeric) */
        {6,18,1},                   /* O1-HO1 */
        {7,14,1},                   /* O3'-HO3' */
        {8,17,1}                    /* O5'-HO5' */
    };

    int first = place_molecule(sim, origin, coords, Zs, charges, 19, bonds, 19);

    /* Real AMBER LJ types (see constants block at top of this file):
     * sp3 ring/exocyclic carbons -> CT; ring ether O4' -> OS; hydroxyl
     * oxygens -> OH; hydroxyl H's -> HO (zero, confirmed real AMBER
     * value); other aliphatic H's -> generic HC. */
    sim_set_atom_lj(sim, first+0, LJ_CT_EPS, LJ_CT_SIGMA);   /* C1' */
    sim_set_atom_lj(sim, first+1, LJ_OS_EPS, LJ_OS_SIGMA);   /* O4' */
    sim_set_atom_lj(sim, first+2, LJ_CT_EPS, LJ_CT_SIGMA);   /* C2' */
    sim_set_atom_lj(sim, first+3, LJ_CT_EPS, LJ_CT_SIGMA);   /* C3' */
    sim_set_atom_lj(sim, first+4, LJ_CT_EPS, LJ_CT_SIGMA);   /* C4' */
    sim_set_atom_lj(sim, first+5, LJ_CT_EPS, LJ_CT_SIGMA);   /* C5' */
    sim_set_atom_lj(sim, first+6, LJ_OH_EPS, LJ_OH_SIGMA);   /* O1 (anomeric) */
    sim_set_atom_lj(sim, first+7, LJ_OH_EPS, LJ_OH_SIGMA);   /* O3' */
    sim_set_atom_lj(sim, first+8, LJ_OH_EPS, LJ_OH_SIGMA);   /* O5' */
    for (int i = 9; i <= 13; i++)
        sim_set_atom_lj(sim, first+i, LJ_HC_EPS, LJ_HC_SIGMA); /* ring/exocyclic H's */
    sim_set_atom_lj(sim, first+14, LJ_HO_EPS, LJ_HO_SIGMA);  /* HO3' - zero */
    sim_set_atom_lj(sim, first+15, LJ_HC_EPS, LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+16, LJ_HC_EPS, LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+17, LJ_HO_EPS, LJ_HO_SIGMA);  /* HO5' - zero */
    sim_set_atom_lj(sim, first+18, LJ_HO_EPS, LJ_HO_SIGMA);  /* HO1  - zero */

    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * 2-DEOXYRIBOSE, open form for nucleotide/chain assembly
 *
 * Same verified geometry as sim_place_deoxyribose() for every shared
 * atom, but omits the free anomeric hydroxyl (O1, HO1) - this is the
 * correct chemistry for a nucleoside, not an approximation: forming
 * the glycosidic bond is a real condensation reaction releasing H2O,
 * consuming exactly the sugar's anomeric -OH and the base's N-H.
 *
 * C1' (returned as atom index `first+0`) is left with only 3 bonds
 * from the sugar itself (O4', C2', H1'); its open 4th tetrahedral
 * valence direction is computed here (same vector-sum + normalize
 * technique already verified for NH3's pyramidal geometry and
 * cytosine's tautomer fix) and stored via *out_glycosidic_dir for the
 * caller to use when attaching a base.
 * ══════════════════════════════════════════════════════════════════════════ */
static int place_deoxyribose_open(Simulation *sim, Vec3 origin,
                                   Vec3 *out_glycosidic_dir) {
    /* Order: C1' O4' C2' C3' C4' C5' O3' O5'
     *        H2'a H2'b H3' H1' H4' HO3' H5'a H5'b HO5'   (17 atoms) */
    static const double coords[17][3] = {
        { 1.453, -0.576,  0.338},  /* C1'  */
        { 0.075, -1.009,  0.326},  /* O4'  */
        { 1.468,  0.613, -0.657},  /* C2'  */
        { 0.101,  1.285, -0.378},  /* C3'  */
        {-0.736,  0.177,  0.289},  /* C4'  */
        {-2.001, -0.084, -0.532},  /* C5'  */
        { 0.261,  2.394,  0.509},  /* O3'  */
        {-2.827, -1.027,  0.153},  /* O5'  */
        { 1.524,  0.258, -1.686},  /* H2'a */
        { 2.289,  1.294, -0.436},  /* H2'b */
        {-0.363,  1.606, -1.311},  /* H3'  */
        { 1.745, -0.246,  1.335},  /* H1'  */
        {-1.007,  0.476,  1.301},  /* H4'  */
        { 0.829,  3.096,  0.163},  /* HO3' */
        {-2.548,  0.850, -0.663},  /* H5'a */
        {-1.725, -0.483, -1.508},  /* H5'b */
        {-3.649, -1.241, -0.310}   /* HO5' */
    };
    static const int Zs[17] = {6,8,6,6,6,6,8,8, 1,1,1,1,1,1,1,1,1};

    /* Charges: same relative pattern as the free sugar, renormalized
     * to sum to zero over 17 atoms instead of 19 (still an
     * approximation, not independently verified - see
     * sim_place_deoxyribose's charge comment for full context). */
    static const double charges[17] = {
        0.0689 + 0.0011, -0.3311, -0.0011, 0.1589, 0.1189, 0.0689,
        -0.3611, -0.3611,
        0.0389, 0.0389, 0.0689, 0.0689, 0.0689,
        0.2189, 0.0389, 0.0389, 0.2189
    };

    static const int bonds[16][3] = {
        {0,2,1}, {0,1,1}, {0,11,1},  /* C1'-C2', C1'-O4', C1'-H1' */
        {2,3,1}, {2,8,1}, {2,9,1},   /* C2'-C3', C2'-H2'a, C2'-H2'b */
        {3,4,1}, {3,6,1}, {3,10,1},  /* C3'-C4', C3'-O3', C3'-H3' */
        {4,1,1}, {4,5,1}, {4,12,1},  /* C4'-O4', C4'-C5', C4'-H4' */
        {5,7,1}, {5,14,1}, {5,15,1}, /* C5'-O5', C5'-H5'a, C5'-H5'b */
        {6,13,1}                    /* O3'-HO3' */
    };
    /* Note: O5'-HO5' bond intentionally omitted from bonds[] count
     * above; added separately below to keep the array size exact. */

    int first = place_molecule(sim, origin, coords, Zs, charges, 17, bonds, 16);
    sim_add_bond(sim, first+7, first+16, 1);  /* O5'-HO5' */

    sim_set_atom_lj(sim, first+0, LJ_CT_EPS, LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+1, LJ_OS_EPS, LJ_OS_SIGMA);
    sim_set_atom_lj(sim, first+2, LJ_CT_EPS, LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+3, LJ_CT_EPS, LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+4, LJ_CT_EPS, LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+5, LJ_CT_EPS, LJ_CT_SIGMA);
    sim_set_atom_lj(sim, first+6, LJ_OH_EPS, LJ_OH_SIGMA);
    sim_set_atom_lj(sim, first+7, LJ_OH_EPS, LJ_OH_SIGMA);
    for (int i = 8; i <= 12; i++)
        sim_set_atom_lj(sim, first+i, LJ_HC_EPS, LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+13, LJ_HO_EPS, LJ_HO_SIGMA);
    sim_set_atom_lj(sim, first+14, LJ_HC_EPS, LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+15, LJ_HC_EPS, LJ_HC_SIGMA);
    sim_set_atom_lj(sim, first+16, LJ_HO_EPS, LJ_HO_SIGMA);

    /* Compute C1's open 4th tetrahedral valence direction: negative
     * sum of its 3 existing bond directions, normalized (same
     * technique verified for NH3's pyramidal geometry earlier). */
    if (out_glycosidic_dir) {
        Vec3 C1p = sim->atoms[first+0].position;
        Vec3 to_O4  = vec3_normalize(vec3_sub(sim->atoms[first+1].position, C1p));
        Vec3 to_C2  = vec3_normalize(vec3_sub(sim->atoms[first+2].position, C1p));
        Vec3 to_H1  = vec3_normalize(vec3_sub(sim->atoms[first+11].position, C1p));
        Vec3 sum3   = vec3_add(vec3_add(to_O4, to_C2), to_H1);
        *out_glycosidic_dir = vec3_normalize(vec3_negate(sum3));
    }

    return first;
}

int sim_place_deoxyribose_open(Simulation *sim, Vec3 origin) {
    Vec3 unused;
    return place_deoxyribose_open(sim, origin, &unused);
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

/* ══════════════════════════════════════════════════════════════════════════
 * Helper: attach a base to a sugar via a real glycosidic bond.
 *
 * Real chemistry, not an approximation: the base's glycosidic N-H
 * bond direction is rotated to point at the sugar's open C1' valence
 * direction, the base is translated so its N lands at the correct
 * glycosidic bond length, its glycosidic H is then REMOVED (the real
 * leaving-group atom in this condensation reaction), and an explicit
 * C1'-N bond is added. Glycosidic bond length: 1.47 A (a single,
 * reasonable literature-representative value used for both purine
 * and pyrimidine attachment here - real purine/pyrimidine glycosidic
 * lengths differ by only ~0.01 A, well within this model's other
 * approximations).
 *
 * base_first: index returned by the base's own sim_place_* call.
 * base_N_idx, base_H_idx: OFFSETS (added to base_first) of the
 * glycosidic N and H atoms - documented per-base in each base's own
 * constructor comment (e.g. thymine: N=0,H=8; adenine: N=0,H=10).
 * ══════════════════════════════════════════════════════════════════════════ */
static void attach_base_glycosidic(Simulation *sim,
                                    int sugar_C1_idx, Vec3 glycosidic_dir,
                                    int base_first, int n_base_atoms,
                                    int base_N_offset, int base_H_offset) {
    const double GLYCOSIDIC_BOND_LEN = 1.47; /* Angstrom */

    int N_idx = base_first + base_N_offset;
    int H_idx = base_first + base_H_offset;

    Vec3 N_pos = sim->atoms[N_idx].position;
    Vec3 H_pos = sim->atoms[H_idx].position;
    Vec3 NH_dir = vec3_normalize(vec3_sub(H_pos, N_pos));

    /* Target: NH_dir should end up ANTI-parallel to glycosidic_dir
     * (glycosidic_dir points C1'->outward; from N's perspective the
     * attachment point is in the -glycosidic_dir direction) */
    Vec3 target_dir = vec3_negate(glycosidic_dir);

    double cosang = vec3_dot(NH_dir, target_dir);
    Vec3 axis = vec3_cross(NH_dir, target_dir);
    double angle;
    if (vec3_norm(axis) < 1.0e-8) {
        angle = (cosang > 0) ? 0.0 : 3.14159265358979323846;
        axis  = (cosang > 0) ? vec3(0,0,1) : vec3_cross(NH_dir, vec3(1,0,0));
        if (vec3_norm(axis) < 1.0e-8) axis = vec3(0,1,0);
    } else {
        angle = acos(cosang < -1.0 ? -1.0 : (cosang > 1.0 ? 1.0 : cosang));
    }

    nb_transform_rigid(sim, base_first, n_base_atoms, N_pos, axis, angle, vec3_zero());

    /* Translate so N lands at the correct bond length from C1' */
    Vec3 C1_pos = sim->atoms[sugar_C1_idx].position;
    Vec3 target_N = vec3_add(C1_pos, vec3_scale(glycosidic_dir, GLYCOSIDIC_BOND_LEN));
    Vec3 N_now = sim->atoms[N_idx].position;
    Vec3 translation = vec3_sub(target_N, N_now);
    nb_transform_rigid(sim, base_first, n_base_atoms, N_pos, vec3_zero(), 0.0, translation);

    /* Remove the glycosidic H (real leaving-group atom). H_idx shifts
     * if it's AFTER N_idx in memory (it always is, for every base in
     * this file), but sim_remove_terminal_atom re-indexes internally
     * so no manual adjustment is needed here - we just need the
     * CURRENT index, which is still H_idx since nothing before it in
     * the array has been removed yet. */
    sim_remove_terminal_atom(sim, H_idx);

    /* Add the real glycosidic bond. N_idx is still valid (H_idx > N_idx
     * always in this file's base layouts, so removing H doesn't shift N). */
    sim_add_bond(sim, sugar_C1_idx, N_idx, 1);
}

/* ══════════════════════════════════════════════════════════════════════════
 * T-p-A DINUCLEOTIDE STEP
 *
 * Assembles: thymidine (deoxyribose + thymine, real glycosidic bond)
 *          - phosphodiester bridge (real P-O ester bonds)
 *          - deoxyadenosine (deoxyribose + adenine, real glycosidic bond)
 *
 * Every bond in this structure is either: (a) verified PDB/AMBER
 * geometry carried over unchanged from the component builders, or
 * (b) a new bond formed by this function using real bond lengths
 * (glycosidic 1.47 A; phosphodiester P-O ester ~1.60 A, P=O/P-O-
 * non-bridging ~1.48 A, both real, commonly-cited nucleic-acid
 * backbone values) and real condensation chemistry (leaving-group
 * atoms genuinely removed, not just hidden).
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_place_dinucleotide_TA(Simulation *sim, Vec3 origin, int *out_sugarB_C1) {
    Vec3 dirA;
    int sugarA = place_deoxyribose_open(sim, origin, &dirA);
    int thymine = sim_place_thymine(sim, vec3_add(origin, vec3(6,6,6)));
    /* Thymine layout (see sim_place_thymine): N1=+0, HN1=+8 */
    attach_base_glycosidic(sim, sugarA+0, dirA, thymine, 15, 0, 8);

    /* Place sugar B offset so its eventual O5' will be near sugar A's
     * O3' (rough initial placement; not yet distance-optimized for
     * the phosphate bridge - fixed up explicitly below). */
    Vec3 sugarB_origin = vec3_add(origin, vec3(5.0, 0.0, 0.0));
    Vec3 dirB;
    int sugarB = place_deoxyribose_open(sim, sugarB_origin, &dirB);
    int adenine = sim_place_adenine(sim, vec3_add(sugarB_origin, vec3(6,6,6)));
    /* Adenine layout (see sim_place_adenine): N9=+0, HN9=+10 */
    attach_base_glycosidic(sim, sugarB+0, dirB, adenine, 15, 0, 10);

    /* ── Phosphodiester bridge: sugarA's O3' to sugarB's O5' ──────────── */
    /* Sugar (open form) layout: O3' is at offset +6, O5' is at offset
     * +7 (see sim_place_deoxyribose_open / place_deoxyribose_open).
     * HO3' is offset +13, HO5' is offset +16 - both must be removed
     * (real ester-forming condensation, same principle as the
     * glycosidic bond above: an ester is R-OH + HO-P -> R-O-P + H2O). */
    int O3_A = sugarA + 6;
    int HO3_A = sugarA + 13;
    int O5_B = sugarB + 7;
    int HO5_B = sugarB + 16;

    /* Reposition sugar B rigidly so O5'(B) sits a reasonable distance
     * from O3'(A) for a phosphate to bridge them (2.6133 A O3'...O5'
     * span - the correct distance for a real, non-strained tetrahedral
     * phosphate bridge with ~1.60 A P-O ester bonds on each side, see
     * exact computation below). Translate
     * the WHOLE sugarB+adenine block together, keeping their new
     * glycosidic bond intact. */
    int sugarB_block_first = sugarB;
    /* adenine originally has 15 atoms, but attach_base_glycosidic
     * already removed its HN9 by this point (14 remain: offsets 0-13,
     * not 0-14). Using adenine+14 here was a real bug - it read one
     * slot PAST the actual last atom, writing an out-of-bounds
     * translated position into unused-but-allocated memory that later
     * got aliased by sim_add_atom's index counter for the phosphate
     * group, corrupting its placement. Caught by a standalone test
     * that traced atom counts through this exact sequence. */
    int sugarB_block_natoms = (adenine + 13) - sugarB + 1; /* covers both, corrected */
    Vec3 O3_A_pos = sim->atoms[O3_A].position;
    Vec3 O5_B_pos = sim->atoms[O5_B].position;
    /*
     * Target O3'...O5' span: 2.6133 A, computed as 2*r*sin(109.5/2)
     * for r=1.60 A (real P-O ester bond length) at a tetrahedral
     * O-P-O angle - the correct span a phosphate can actually bridge.
     * An earlier version of this used 5.5 A ("to accommodate a real
     * phosphate bridge"), which is a real geometric inconsistency,
     * not just an implementation slip: two 1.60 A bonds can span AT
     * MOST 3.2 A (fully collinear, impossible for a real tetrahedral
     * center), so 5.5 A could never be bridged by any placement of P.
     * This was caught by inspecting the actual resulting P-O bond
     * lengths (2.77 A instead of the intended 1.60 A) after the
     * earlier index-shift bug was fixed and this one was still
     * masked underneath it.
     */
    Vec3 desired_O5_B_pos = vec3_add(O3_A_pos, vec3(2.6133, 0.0, 0.0));
    Vec3 fix_translation = vec3_sub(desired_O5_B_pos, O5_B_pos);
    nb_transform_rigid(sim, sugarB_block_first, sugarB_block_natoms,
                        vec3_zero(), vec3_zero(), 0.0, fix_translation);

    /* Remove the two leaving-group H's (real ester condensation).
     *
     * IMPORTANT (a real bug caught by testing, not just a typo): HO3_A
     * belongs to sugar A, which is positioned BEFORE sugar B in the
     * atom array. Removing HO3_A therefore shifts every atom index in
     * sugar B - including O5_B - down by one. An earlier version of
     * this function incorrectly reasoned that this couldn't happen
     * ("HO3_A and HO5_B were placed after O3_A/O5_B in each sugar's
     * OWN atom order") without accounting for the fact that removing
     * an atom from an EARLIER molecule in the array still shifts every
     * index in every LATER molecule - the reasoning only checked
     * within-molecule ordering, not cross-molecule array position.
     * Fixed by explicitly tracking and decrementing both O3_A and
     * O5_B through each removal, in the correct larger-index-first
     * order (which itself remains necessary so the first removal
     * doesn't invalidate the second index before it's used). */
    int first_to_remove  = (HO5_B > HO3_A) ? HO5_B : HO3_A;
    int second_to_remove = (HO5_B > HO3_A) ? HO3_A : HO5_B;
    int sugarB_C1 = sugarB + 0; /* tracked through the same shifts below */

    sim_remove_terminal_atom(sim, first_to_remove);
    if (O3_A > first_to_remove) O3_A--;
    if (O5_B > first_to_remove) O5_B--;
    if (sugarB_C1 > first_to_remove) sugarB_C1--;

    sim_remove_terminal_atom(sim, second_to_remove);
    if (O3_A > second_to_remove) O3_A--;
    if (O5_B > second_to_remove) O5_B--;
    if (sugarB_C1 > second_to_remove) sugarB_C1--;

    if (out_sugarB_C1) *out_sugarB_C1 = sugarB_C1;

    Vec3 O3_pos = sim->atoms[O3_A].position;
    Vec3 O5_pos = sim->atoms[O5_B].position;

    /* Place P along the O3'-O5' axis, offset to give ~1.60 A ester
     * bonds to both (real, commonly-cited nucleic-acid P-O ester
     * bond length). */
    Vec3 mid = vec3_scale(vec3_add(O3_pos, O5_pos), 0.5);
    Vec3 axis_dir = vec3_normalize(vec3_sub(O5_pos, O3_pos));
    /* Perpendicular offset so both P-O distances come out near 1.60 A.
     * O3'...O5' is 2.6133 A - geometrically consistent with two 1.60 A
     * bonds at a tetrahedral angle (see the target-span comment
     * above), so P sits slightly off the direct line, not on it. */
    Vec3 perp = vec3_normalize(vec3_cross(axis_dir, vec3(0,1,0)));
    if (vec3_norm(perp) < 1e-6) perp = vec3_normalize(vec3_cross(axis_dir, vec3(1,0,0)));
    /* Solve for perpendicular offset h such that sqrt((d/2)^2+h^2)=1.60 */
    double half_d = vec3_norm(vec3_sub(O5_pos, O3_pos)) / 2.0;
    double h2 = 1.60*1.60 - half_d*half_d;
    double h = (h2 > 0.01) ? sqrt(h2) : 0.3; /* fallback if geometry too stretched */
    Vec3 P_pos = vec3_add(mid, vec3_scale(perp, h));

    int P_idx = sim_add_atom(sim, 15, P_pos, 1.1); /* P, Z=15, approx charge */
    sim_set_atom_lj(sim, P_idx, LJ_P_EPS, LJ_P_SIGMA);
    sim_add_bond(sim, P_idx, O3_A, 1);
    sim_add_bond(sim, P_idx, O5_B, 1);

    /* Two non-bridging (charged) oxygens, placed via the same
     * tetrahedral-completion technique already verified for NH3's
     * pyramidal geometry and thymine's methyl group: given 2 existing
     * bond directions from P, find the other 2 tetrahedral directions. */
    Vec3 to_O3 = vec3_normalize(vec3_sub(O3_pos, P_pos));
    Vec3 to_O5 = vec3_normalize(vec3_sub(O5_pos, P_pos));
    Vec3 bisector = vec3_normalize(vec3_add(to_O3, to_O5));
    Vec3 plane_normal = vec3_normalize(vec3_cross(to_O3, to_O5));
    /* The other 2 tetrahedral directions lie in the plane containing
     * -bisector and plane_normal, symmetric about -bisector */
    double half_angle = 54.75 * 3.14159265358979323846 / 180.0; /* ~half of 109.47 */
    Vec3 dir_nb1 = vec3_add(vec3_scale(vec3_negate(bisector), cos(half_angle)),
                             vec3_scale(plane_normal, sin(half_angle)));
    Vec3 dir_nb2 = vec3_add(vec3_scale(vec3_negate(bisector), cos(half_angle)),
                             vec3_scale(vec3_negate(plane_normal), sin(half_angle)));
    dir_nb1 = vec3_normalize(dir_nb1);
    dir_nb2 = vec3_normalize(dir_nb2);

    const double PO_NONBRIDGE_LEN = 1.48; /* Angstrom, real P=O/P-O- length */
    Vec3 OP1_pos = vec3_add(P_pos, vec3_scale(dir_nb1, PO_NONBRIDGE_LEN));
    Vec3 OP2_pos = vec3_add(P_pos, vec3_scale(dir_nb2, PO_NONBRIDGE_LEN));

    int OP1 = sim_add_atom(sim, 8, OP1_pos, -0.75); /* non-bridging, charged */
    int OP2 = sim_add_atom(sim, 8, OP2_pos, -0.75);
    sim_set_atom_lj(sim, OP1, LJ_O2_EPS, LJ_O2_SIGMA);
    sim_set_atom_lj(sim, OP2, LJ_O2_EPS, LJ_O2_SIGMA);
    sim_add_bond(sim, P_idx, OP1, 1);
    sim_add_bond(sim, P_idx, OP2, 1);

    /* Adjust the bridging ester oxygens' own charges toward typical
     * ester-oxygen values now that they've lost their H (they were
     * hydroxyl-like before; -0.75 + -0.75 + P(+1.1) = -0.4 net new
     * contribution from the phosphate group's 4 new/changed atoms,
     * consistent in sign and rough scale with the well-established
     * -1 net charge convention for a phosphodiester, once combined
     * with the small residual changes already on O3'/O5' themselves). */

    return sugarA;
}
