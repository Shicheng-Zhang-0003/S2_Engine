#ifndef NUCLEOBASES_H
#define NUCLEOBASES_H

#include "types.h"

/*
 * nucleobases.h
 *
 * Constructors for the five DNA/RNA nucleobases: uracil, cytosine,
 * thymine, adenine, guanine.
 *
 * PROVENANCE OF GEOMETRY:
 * Seed atomic positions for uracil, cytosine, adenine, and guanine are
 * taken from the RCSB PDB Chemical Component Dictionary's "ideal
 * coordinates" (ligand codes URA, CYT, ADE, GUN respectively), fetched
 * directly from https://files.rcsb.org/ligands/view/<CODE>.cif on
 * 2026-06-19. These are algorithmically-generated (OpenEye OMEGA /
 * Molecular Networks Corina), chemically-sensible 3D embeddings of the
 * verified 2D connectivity - they are NOT raw spectroscopic or
 * crystallographic measurements, but they do encode correct bond
 * connectivity, correct ring topology, and reasonable bond lengths/
 * angles for every atom, independently of any derivation in this file.
 *
 * Thymine has no standalone PDB ligand entry under the code "THY" (that
 * code is taken by an unrelated thiamine diphosphate derivative - a
 * PDB code collision). Thymine is instead constructed by taking
 * uracil's own verified C5-H5 bond vector and replacing that hydrogen
 * with a tetrahedral methyl group at the standard aromatic-ring-to-
 * methyl C-C bond length (1.51 Angstrom) - this is the correct
 * chemical relationship (thymine = 5-methyluracil) and avoids
 * depending on a second, harder-to-verify external coordinate source.
 *
 * EQUILIBRIUM GEOMETRY CONSTRUCTION:
 * Every bond's r0 and every angle's theta0 are set to EXACTLY match the
 * placed Cartesian seed geometry (computed at construction time via
 * vec3_dist / vec3_angle, not hardcoded) - so these molecules start at
 * zero strain energy, and a static geometry report computed from the
 * positions will exactly reproduce the seed bond lengths/angles.
 * Force constants (k) are reasonable generic single/double-bond and
 * aromatic-ring-bending values (same order of magnitude as the rest of
 * this codebase's AMBER-derived table), not independently fitted to
 * spectroscopic data for these specific rings. Partial charges are left
 * at 0.0 - accurate RESP/AMBER nucleobase charges (needed for
 * quantitative base-pairing energetics) are deferred to a follow-up;
 * fabricating precise-looking but unverified per-atom charges here
 * would be worse than leaving them at a clearly-flagged placeholder.
 *
 * All constructors return the index of the first atom added, atoms
 * follow in PDB atom-naming order as commented in each function.
 */

int sim_place_uracil   (Simulation *sim, Vec3 origin);
int sim_place_cytosine (Simulation *sim, Vec3 origin);
int sim_place_thymine  (Simulation *sim, Vec3 origin);
int sim_place_adenine  (Simulation *sim, Vec3 origin);
int sim_place_guanine  (Simulation *sim, Vec3 origin);

/*
 * 2-deoxyribose (the DNA sugar), free form (anomeric OH unreacted).
 *
 * Geometry source: RCSB PDB Chemical Component Dictionary, ligand code
 * "2DR" (2-deoxy-beta-D-erythro-pentofuranose), ideal coordinates,
 * fetched 2026-06-19 from https://files.rcsb.org/ligands/view/2DR.cif.
 * Verified: C5H10O4, 19 atoms, 19 bonds, 3 chiral centers (correct -
 * deoxyribose has stereocenters at C1', C3', C4').
 *
 * Atom order matches standard nucleic-acid sugar numbering directly
 * (the PDB file's own atom names already use this convention):
 *   C1' O4' C2' C3' C4' C5' O1(anomeric OH) O3' O5'
 *   H2'(x2) H3' H1' H4' HO3' H5'(x2) HO5' HO1(anomeric H)
 *
 * In a real nucleotide: O1 (anomeric OH on C1') is replaced by the
 * glycosidic bond to a base's N9 (purine) or N1 (pyrimidine); O3'-H
 * and O5'-H are replaced by phosphodiester linkages when part of a
 * chain (O5' typically carries the 5'-phosphate, O3' links to the
 * next nucleotide's 5'-phosphate). C2' bears only 2 H's and no
 * oxygen substituent - this absence (vs. ribose's O2'-H) is exactly
 * what "2-deoxy" means.
 */
int sim_place_deoxyribose(Simulation *sim, Vec3 origin);

/*
 * Sugar variant for nucleotide/chain assembly: identical geometry to
 * sim_place_deoxyribose() for all shared atoms, but OMITS the free
 * anomeric hydroxyl (O1, HO1) entirely, leaving C1' with an open
 * valence. This is the correct chemistry for a nucleoside/nucleotide
 * (not an approximation): forming the glycosidic bond is a real
 * condensation reaction that releases H2O, consuming exactly the
 * sugar's anomeric -OH and the base's N-H. 17 atoms (19 minus O1,
 * HO1), same atom order as sim_place_deoxyribose with those two
 * dropped.  C1' is atom index 0 of the returned block, with its 4th
 * bond left for the caller to complete via an explicit glycosidic
 * bond to a base's N9 (purine) or N1 (pyrimidine).
 */
int sim_place_deoxyribose_open(Simulation *sim, Vec3 origin);

/*
 * Builds a real, verified two-nucleotide DNA fragment (a "T-p-A step"):
 * two complete nucleosides (thymidine and deoxyadenosine, each formed
 * via a real glycosidic condensation bond - sugar loses its anomeric
 * -OH, base loses its glycosidic -H, exactly as the real reaction
 * does) linked by one real phosphodiester bridge between the first
 * sugar's O3' and the second sugar's O5' - the actual chain-forming
 * bond of the DNA backbone.
 *
 * This is a concrete, fully worked example rather than a generic
 * "any base combination" builder: each base's internal atom layout
 * places its glycosidic N and H at different indices (documented
 * per-base in nucleobases.c), and generalizing this into a single
 * polymorphic function would need an additional per-base metadata
 * table not yet built - a reasonable next extension, not implemented
 * here.
 *
 * Returns the atom index of sugar A's C1'. Sugar B's C1' index is
 * written to *out_sugarB_C1 if that pointer is non-NULL - this is
 * NOT simply "return value + a fixed offset": internal condensation
 * chemistry (real atoms genuinely removed, not hidden) shifts every
 * index after each removal, and sugar B's own atoms shift when an
 * earlier atom (sugar A's HO3', during phosphodiester formation) is
 * removed. Callers needing to inspect sugar B's structure MUST use
 * *out_sugarB_C1 rather than assume any fixed relationship to the
 * return value.
 */
int sim_place_dinucleotide_TA(Simulation *sim, Vec3 origin, int *out_sugarB_C1);

/*
 * Fit a least-squares plane through the given atom indices and return
 * the maximum perpendicular distance (Angstrom) of any atom from that
 * plane. Aromatic rings should be very close to planar (a few
 * hundredths of an Angstrom at most) - this is an independent physical
 * sanity check on the placed geometry that doesn't depend on any
 * external reference numbers.
 */
double nb_planarity_deviation(const Simulation *sim,
                               const int *atom_indices, int n);

/*
 * Compute the unit normal vector of the plane defined by the first 3
 * given atom indices (cross product of the two edge vectors from the
 * first atom). For a near-planar ring, this is a good approximation
 * of the ring's plane normal, useful for rigid-body alignment of one
 * planar molecule against another (e.g. orienting a base pair).
 */
Vec3 nb_ring_normal(const Simulation *sim, const int *atom_indices);

/*
 * Signed angle to rotate v_from onto v_to, as a rotation around axis n
 * (right-hand rule), after projecting both onto the plane perpendicular
 * to n. Needed to fix the azimuthal orientation left unconstrained by
 * plane-normal alignment alone.
 */
double nb_signed_inplane_angle(Vec3 v_from, Vec3 v_to, Vec3 n);

/*
 * Rigidly rotate (around `pivot`, by `angle` radians about `axis`) and
 * then translate (by `translation`) every atom in [first_atom,
 * first_atom + n_atoms). Since rotation+translation is an isometry,
 * this does not disturb any bond length or angle already set up for
 * this group of atoms - safe to call on an already-bonded molecule.
 */
void nb_transform_rigid(Simulation *sim, int first_atom, int n_atoms,
                         Vec3 pivot, Vec3 axis, double angle,
                         Vec3 translation);

#endif /* NUCLEOBASES_H */
