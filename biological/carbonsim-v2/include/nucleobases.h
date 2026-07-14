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
 * Fit a least-squares plane through the given atom indices and return
 * the maximum perpendicular distance (Angstrom) of any atom from that
 * plane. Aromatic rings should be very close to planar (a few
 * hundredths of an Angstrom at most) - this is an independent physical
 * sanity check on the placed geometry that doesn't depend on any
 * external reference numbers.
 */
double nb_planarity_deviation(const Simulation *sim,
                               const int *atom_indices, int n);

#endif /* NUCLEOBASES_H */
