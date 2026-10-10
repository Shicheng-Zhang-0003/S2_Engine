#ifndef AMBER_LJ_H
#define AMBER_LJ_H

/*
 * amber_lj.h — AMBER ff99 Lennard-Jones types for biomolecular atoms.
 *
 * WHY THIS FILE EXISTS (streamlining pass)
 *
 * Before this, three modules each carried their own private copy of the
 * same AMBER ff99 numbers, and main.c carried a fourth spelling of the
 * same carbonyl-oxygen sigma inline in seven separate functions:
 *
 *     src/nucleobases.c   AMBER_RSTAR_TO_SIGMA + LJ_RING_N_*, LJ_SP2_C_*, ...
 *     src/aminoacids.c    AA_RSTAR_TO_SIGMA   + AA_LJ_N_*,  AA_LJ_C_*,  ...
 *     src/kcsa_filter.c   its own copies
 *     src/main.c          1.6612 * 2.0 / 1.122462048309373
 *                         (14 occurrences across 7 functions)
 *
 * The literal 2^(1/6) appeared NINETEEN times across four files before
 * this consolidation (verified by grep against the pre-refactor tree:
 * sixteen in main.c - fourteen carbonyl-sigma expressions plus the two
 * JC ion-sigma macros - and one each in nucleobases.c, aminoacids.c and
 * kcsa_filter.c). aminoacids.c
 * admitted the duplication in its own comment — "duplicated here rather than
 * shared via a header refactor to avoid touching already-validated, working
 * code under time pressure" — which is an honest note about a debt that was
 * never paid. This file pays it.
 *
 * THE CONVERSION, ONCE
 *
 * AMBER/TINKER prm files tabulate a parameter R* (Rstar = Rmin/2: HALF
 * the distance of the 12-6 well minimum) and write the potential as
 *
 *     E = eps * [ (Rmin/r)^12 - 2 (Rmin/r)^6 ],   Rmin = 2 Rstar
 *
 * whose MINIMUM sits at r = Rmin. The standard 12-6 form used by
 * pair_nonbonded_core() is
 *
 *     E = 4 eps [ (sigma/r)^12 - (sigma/r)^6 ]
 *
 * whose minimum sits at r = 2^(1/6) sigma. Matching the two minima,
 * 2 Rstar = 2^(1/6) sigma, gives
 *
 *     sigma = 2 Rstar / 2^(1/6) = Rstar * 2^(5/6).
 *
 * The factor of 2 is required because AMBER's Rstar is a HALF-distance:
 * dropping it (sigma = Rstar / 2^(1/6)) would halve every sigma. This
 * macro computes the full expression and reproduces AMBER's published
 * sigmas exactly:
 *
 *     N  amide   Rstar 1.8240 -> 3.25000  (AMBER 3.2500)
 *     C  sp2     Rstar 1.9080 -> 3.39967  (AMBER 3.3997)
 *     O  carbonyl Rstar 1.6612 -> 2.95992  (AMBER 2.9600)
 *     CT sp3     Rstar 1.9080 -> 3.39967  (AMBER 3.3997)
 *
 * FULL-AUDIT Q5: the worked examples above used to read 3.24979 and
 * 2.96000. Both were wrong in the fifth decimal: 2*1.8240/2^(1/6) =
 * 3.2499985 and 2*1.6612/2^(1/6) = 2.9599219. The macro was always
 * right — the arithmetic printed beside it was not, which is exactly the
 * kind of number a reader copies into their own notes.
 *
 * The periodic table converts with a DIFFERENT expression because UFF
 * tabulates a different quantity: UFF's x1 is the FULL distance of the
 * potential minimum, not half of it, so its correct conversion is
 * sigma = x1 / 2^(1/6). The v9R4 audit found periodic_table.c storing
 * UFF's x1 directly in a field declared as sigma, putting every
 * effective LJ well 12.2% too far out; it now divides by 2^(1/6)
 * itself. The two expressions differ by the factor of 2 precisely
 * because the tabulated inputs differ - neither may be swapped for the
 * other.
 */

#include "constants.h"

/* 2^(1/6) to full double precision. ONE definition, everywhere. */
#define TWOPOW_SIXTH   1.122462048309373

/* AMBER Rstar (a half-distance) -> standard 12-6 collision diameter. */
#define AMBER_RSTAR_TO_SIGMA(rstar) ((rstar) * 2.0 / TWOPOW_SIXTH)

/* ── Nitrogen ─────────────────────────────────────────────────────────── */
/* AMBER classes N, NA, NB, NC, N*, N2 — all identical in ff99. */
#define LJ_AMBER_N_SIGMA      AMBER_RSTAR_TO_SIGMA(1.8240)
#define LJ_AMBER_N_EPS        (0.1700 * KCAL_MOL_TO_EV)

/* ── sp2 carbon ───────────────────────────────────────────────────────── */
/* AMBER classes C, CA, CB, CM, CK, CQ. */
#define LJ_AMBER_C2_SIGMA     AMBER_RSTAR_TO_SIGMA(1.9080)
#define LJ_AMBER_C2_EPS       (0.0860 * KCAL_MOL_TO_EV)

/* ── sp3 carbon (methyl / aliphatic) ──────────────────────────────────── */
/* AMBER class CT. */
#define LJ_AMBER_CT_SIGMA     AMBER_RSTAR_TO_SIGMA(1.9080)
#define LJ_AMBER_CT_EPS       (0.1094 * KCAL_MOL_TO_EV)

/* ── Oxygen ───────────────────────────────────────────────────────────── */
/* Carbonyl / backbone O. */
#define LJ_AMBER_O_SIGMA      AMBER_RSTAR_TO_SIGMA(1.6612)
#define LJ_AMBER_O_EPS        (0.2100 * KCAL_MOL_TO_EV)

/* Hydroxyl O (serine, threonine, tyrosine, sugar hydroxyls). */
#define LJ_AMBER_OH_SIGMA     AMBER_RSTAR_TO_SIGMA(1.7210)
#define LJ_AMBER_OH_EPS       (0.2104 * KCAL_MOL_TO_EV)

/* Ring ether O: the furanose/deoxyribose ring oxygen. Distinct from the
 * hydroxyl above — same element, different bonding, different AMBER class. */
#define LJ_AMBER_OS_SIGMA     AMBER_RSTAR_TO_SIGMA(1.6837)
#define LJ_AMBER_OS_EPS       (0.1700 * KCAL_MOL_TO_EV)

/* Phosphorus: the phosphate backbone of nucleic acids. */
#define LJ_AMBER_P_SIGMA      AMBER_RSTAR_TO_SIGMA(2.1000)
#define LJ_AMBER_P_EPS        (0.2000 * KCAL_MOL_TO_EV)

/* ── Hydrogens ────────────────────────────────────────────────────────── */
/* H on nitrogen: AMBER class H. */
#define LJ_AMBER_HN_SIGMA     AMBER_RSTAR_TO_SIGMA(0.6000)
#define LJ_AMBER_HN_EPS       (0.0157 * KCAL_MOL_TO_EV)

/* Aliphatic C-H: AMBER class HC. */
#define LJ_AMBER_HC_SIGMA     AMBER_RSTAR_TO_SIGMA(1.4870)
#define LJ_AMBER_HC_EPS       (0.0157 * KCAL_MOL_TO_EV)

/* H on sp3 carbon bearing one electronegative neighbour (CA-H): H1. */
#define LJ_AMBER_H1_SIGMA     AMBER_RSTAR_TO_SIGMA(1.3870)
#define LJ_AMBER_H1_EPS       (0.0157 * KCAL_MOL_TO_EV)

/* Aromatic C-H, heteroatom-adjacent: AMBER class H4 (H6 of pyrimidines). */
#define LJ_AMBER_H4_SIGMA     AMBER_RSTAR_TO_SIGMA(1.4090)
#define LJ_AMBER_H4_EPS       (0.0150 * KCAL_MOL_TO_EV)

/* Aromatic C-H: AMBER class HA (H5 of pyrimidines). */
#define LJ_AMBER_HA_SIGMA     AMBER_RSTAR_TO_SIGMA(1.4590)
#define LJ_AMBER_HA_EPS       (0.0150 * KCAL_MOL_TO_EV)

/* Purine C-H: AMBER class H5. The name collides with this codebase's own
 * "H5" ring-position label in pyrimidines; they are unrelated. */
#define LJ_AMBER_H5_SIGMA     AMBER_RSTAR_TO_SIGMA(1.3590)
#define LJ_AMBER_H5_EPS       (0.0150 * KCAL_MOL_TO_EV)

/*
 * Hydroxyl hydrogen (AMBER class HO).
 *
 * ff99 gives this EXACTLY zero, and it is deliberately NOT used here.
 * Pair_nonbonded_core() combines epsilons by geometric mean, so a single
 * zero-epsilon atom drives a whole pairwise interaction to exactly zero
 * regardless of the partner's real LJ — leaving nothing but Coulomb
 * between a hydroxyl H and any oppositely charged atom it approaches.
 * That was confirmed dangerous by direct testing before these small
 * nonzero values were adopted. This is a deviation from ff99 and is
 * labelled as one wherever it is used.
 */
/* These are NOT derived from a table: both nucleobases.c and aminoacids.c
 * independently settled on 0.5 A / 0.001 kcal/mol after the failure
 * described above, and the consolidation confirmed the two copies agreed
 * rather than assuming they did.
 *
 * FULL-AUDIT M41: epsilon is 0.001 KCAL/MOL, not 0.001 eV. The bare
 * literal overstated the well depth 23.06x (0.001 eV = 0.02306 kcal/mol)
 * and propagated sqrt(23)~4.8x into every HO-X pair via the geometric
 * mean. The conversion is now explicit like every sibling macro. */
#define LJ_AMBER_HO_SIGMA     0.5
#define LJ_AMBER_HO_EPS       (0.001 * KCAL_MOL_TO_EV)

#endif /* AMBER_LJ_H */
