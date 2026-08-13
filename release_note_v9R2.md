# carbonsim v9A1 — dihedral forces, emergent helix geometry, and an honest KcsA negative result

The first version with genuine 3D structure capability: torsion forces, energy minimization, and multi-residue chain assembly — plus the first demo that tests the engine against real ion-channel physics and reports a genuine negative result.

## Highlights

- **Demo 11 — alpha helix:** given only correct local backbone torsion geometry (textbook φ/ψ/ω applied as restraints), the defining i,i+4 backbone hydrogen bond forms from the same validated Coulomb+LJ physics proven on the water trimer and G-C/A-U pairing (H···O = 2.16 Å, N···O = 3.13 Å). Not a spontaneous-folding test, and the output says so.
- **Demo 12 — KcsA selectivity filter:** K⁺ vs Na⁺ across fixed-radius site tests, a radius scan, and a two-ring antiprism construction, using real backbone-carbonyl charges (−0.55 e) and amino-acid-specific carbonyl LJ typing. **Result: Na⁺ favored at every test — the wrong direction.** Systematic elimination rules out generic oxygen typing and geometry-fit; the diagnosis points at missing electronic polarizability. Reported as a finding, not hidden.
- **Demo 12's antiprism block is a placeholder, not a result** — the ring z-separation was never sourced, and its numbers are dominated by an artificial O–O overlap. It says so in the output and must not be cited as selectivity energetics.

## New in this version

- Dihedral (torsion) forces in AMBER/CHARMM form, gradients by central finite difference — correctness-by-construction over a hand-derived analytical chain rule
- Steepest-descent energy minimizer with adaptive step, per-atom displacement cap, divergence guard, and a frozen-atom variant
- Poly-alanine chain builder (up to 16 residues) with centralized, defensive index tracking across every atom removal
- v9 masterplan Step 0 resolved in-file: potassium's `periodic_table.c` row verified against tabulated data, stale comment removed; LJ parameters remain flagged UFF-sourced but not independently re-verified

## Honest scope

- No explicit solvent anywhere — vacuum + dielectric divisor; this is the root of Demo 7's quantitative overshoot and the reason Demo 12 cannot express the dehydration term real selectivity depends on
- Base-pairing magnitudes are qualitative only (ordering validated, numbers are not)
- Sugar/phosphate/amino-acid charges are charge-balanced approximations, not verified RESP fits (nucleobase charges are — see the confidence tiers in the readme)
- O(N²) pair loop, no 1-4 scaling, steepest descent only

## Build

```bash
cd biological/v9A1
make
./carbonsim > output.txt
```

Requires a C11 compiler and `make`. ASan-clean; fast-math and no-fast-math outputs verified byte-identical (non-fast-math).

**Full record** — all twelve demo outputs verbatim, complete parameter provenance, verification discipline, and known limitations: `readme.md`.
