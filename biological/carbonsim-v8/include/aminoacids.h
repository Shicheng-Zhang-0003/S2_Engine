#ifndef AMINOACIDS_H
#define AMINOACIDS_H

#include "types.h"

/*
 * aminoacids.h
 *
 * Constructors for free amino acids and real peptide-bond assembly,
 * following the same verification discipline as nucleobases.h/.c.
 *
 * GEOMETRY PROVENANCE: seed atomic positions for glycine and alanine
 * are taken from the RCSB PDB Chemical Component Dictionary's "ideal
 * coordinates" (ligand codes GLY, ALA), fetched directly from
 * https://files.rcsb.org/ligands/view/<CODE>.cif. Both are the
 * standard, non-terminal-modified, neutral (non-zwitterionic) free
 * amino acid forms - confirmed by each CIF's own SMILES string
 * (glycine: "NCC(O)=O"; alanine: "C[C@H](N)C(O)=O", explicitly showing
 * unionized -NH2 and -COOH groups, not the -NH3+/-COO- zwitterion).
 *
 * PEPTIDE BOND CHEMISTRY: real condensation reaction, not an
 * approximation. The PDB CCD itself flags which specific atoms are
 * "leaving atoms" (pdbx_leaving_atom_flag) when a residue is
 * incorporated into a polymer chain as a non-terminal residue: one
 * amine hydrogen (named "H2" in both GLY and ALA's CCD entries) and
 * the entire carboxyl hydroxyl (OXT + HXT). This is authoritative,
 * PDB-defined guidance, not an inferred convention - the N-terminal
 * residue in a peptide keeps both its amine hydrogens and loses only
 * OXT/HXT (bonding via its carbonyl C to the next residue's N); the
 * C-terminal residue keeps OXT/HXT and loses only H2 (bonding via its
 * N to the previous residue's carbonyl C).
 *
 * CHARGES: explicitly approximate, same status as the sugar/phosphate
 * charges in nucleobases.c. No sufficiently authoritative primary
 * source (an actual AMBER parameter file or a peer-reviewed RESP
 * table) was located for standard protein backbone charges during
 * this session - only a third-party software package's claimed
 * values, which was not trusted as a citation-worthy source given
 * this codebase's standard of verifying against primary data. Values
 * here are a physically-reasonable, charge-balanced approximation,
 * not independently verified the way the nucleobase RESP charges are.
 */

/* Free glycine (NH2-CH2-COOH), neutral form, 10 atoms.
 * Atom order: N CA C O OXT H H2 HA2 HA3 HXT */
int sim_place_glycine(Simulation *sim, Vec3 origin);

/* Free L-alanine (NH2-CH(CH3)-COOH), neutral form, 13 atoms.
 * Atom order: N CA C O CB OXT H H2 HA HB1 HB2 HB3 HXT */
int sim_place_alanine(Simulation *sim, Vec3 origin);

/*
 * Assembles a real Gly-Ala dipeptide: glycine as the N-terminal
 * residue (keeps both amine H's, loses its carboxyl OH), alanine as
 * the C-terminal residue (keeps its carboxyl OH, loses one amine H),
 * linked by a genuine peptide (amide) bond formed via the same real
 * condensation chemistry already used for the DNA backbone's
 * phosphodiester bond (sim_remove_terminal_atom removes the real
 * leaving-group atoms; nothing is hidden).
 *
 * Returns the atom index of glycine's N (the dipeptide's N-terminus).
 * Alanine's N index is written to *out_ala_N if non-NULL - like the
 * dinucleotide function, internal atom removal shifts indices, so
 * callers must use this rather than assume a fixed offset.
 */
int sim_place_dipeptide_GlyAla(Simulation *sim, Vec3 origin, int *out_ala_N);

/*
 * Tracks one residue's key backbone/leaving-group atom indices as a
 * peptide chain is assembled. Used with aa_adjust_residues_after_removal
 * below - a centralized, defensive alternative to manually reasoning
 * about which specific indices are affected by each atom removal
 * (exactly the reasoning that produced real bugs earlier in this
 * project). Every tracked index across every residue is unconditionally
 * checked and decremented if it was greater than whatever was just
 * removed - correct by construction, regardless of how many residues
 * or removals are involved.
 */
typedef struct {
    int N, CA, C, O, OXT, HXT, H2;
    int H; /* amide backbone N-H (the hydrogen that STAYS after H2 is
            * removed) - tracked explicitly rather than assumed at a
            * fixed offset from N, since OXT sits BETWEEN N and H in
            * alanine's own atom order and is removed for every
            * non-C-terminal residue, shifting H's relative offset. */
} AAResidue;

/* Call after every sim_remove_terminal_atom(sim, removed_idx) while
 * assembling a chain: decrements every tracked field across all
 * `count` residues that was greater than removed_idx. */
void aa_adjust_residues_after_removal(AAResidue *residues, int count,
                                       int removed_idx);

#define MAX_CHAIN_RESIDUES 16

/*
 * Assembles a chain of n_residues alanine residues, linked by
 * n_residues-1 real peptide bonds (same condensation chemistry as
 * sim_place_dipeptide_GlyAla, generalized to a sequential chain, with
 * the angle terms spanning each peptide bond junction included).
 * residues_out (caller-provided array, size >= n_residues) receives
 * each residue's final, post-assembly backbone atom indices for use
 * in adding dihedral restraints or checking hydrogen bonding
 * afterward. n_residues is capped at MAX_CHAIN_RESIDUES.
 * Returns the first residue's N index (the chain's N-terminus).
 */
int sim_place_polyalanine(Simulation *sim, Vec3 origin, int n_residues,
                           AAResidue *residues_out);

#endif /* AMINOACIDS_H */
