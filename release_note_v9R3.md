# carbonsim v9R3 — the audit & correctness release

Not a new-physics release — a systematic, sector-by-sector correctness
audit of the whole engine, held to the project's own standard: every
number sourced, every approximation flagged, every silent behavior made
explicit. Every finding was fixed in code or documented as a deliberate
limitation; the full 12-demo sequence was regenerated and re-asserted
after every behavior-affecting change.

## Highlights

- **The record is clean and reproducible.** Every demo re-ran, every
  assertion re-passed, and the captured `output.txt` is the verbatim
  post-audit record.
- **Demo 11's emergent helix H-bond survived the audit** (H···O ≈ 2.16 Å,
  N···O ≈ 3.13 Å) — the clean hard-cutoff baseline, restored once a
  silently-active cutoff switch was gated off.
- **Demo 12's KcsA negative result stands** (Na⁺ favored, wrong direction,
  in every test) — the audit changed no physics that was already honest.

## What the audit fixed

- **Constants derived in-line (C1–C4)**: `KCAL_MOL_TO_EV`, `COULOMB_MD`,
  `J_TO_EV`, `EV_TO_HARTREE`, `EV_TO_KCAL_MOL` are computed from 2019
  CODATA primaries at the point of definition; every reciprocal pair now
  multiplies to exactly 1.0 by construction.
- **BOND_TABLE provenance corrected (F2)**: the false "AMBER ff14SB /
  CHARMM36" claim replaced with the honest tier-2/3 description (generic
  spectroscopic-order stiffness constants) plus spot-check numbers.
- **Fallback warnings (F3)**: bond/angle parameter fallbacks now print a
  one-line warning when taken.
- **Topology & integrator hygiene (S1–S3, I1–I2)**: angle-rebuild
  functions renamed to make their reset semantics explicit; dead fields
  removed; the silent no-op Nosé-Hoover enum removed; the 3N−3
  temperature comment fixed; the global RNG replaced with a
  per-simulation seedable state.
- **Potassium verified (T1)** against Rappé 1992 UFF, fully closing the
  masterplan Step 0 flag.
- **Makefile ASan path (B1)**: the link rule now passes `$(CFLAGS)`, so
  the sanitizer build works through make.
- **k_restraint macro (M1)**: the one hardcoded kcal→eV literal now uses
  `KCAL_MOL_TO_EV` — one conversion value in the whole tree.

## Investigated and deliberately rejected

- **1-4 non-bonded scaling (F1)**: implemented, tested, and found to
  break the validated helix i,i+4 H-bond (Demo 11) — AMBER's scaling
  presupposes co-fitted torsions this force field doesn't have. Reverted;
  1-4 pairs stay at full strength; reasoning kept in-source.

## Added but gated off

- **Cutoff switching (F5)**: a CHARMM-style smooth switch exists for
  future condensed-phase use but is opt-in (default off). The audit found
  it had been silently active in the gas-phase demos (attenuating the
  9.6–12 Å band), perturbing Demos 7 and 11; gating it off restored the
  clean baseline.

## Behavior changes vs v9R2

- Demo 7 and Demo 11 numbers shift slightly (F5 switching removed + C1/C4
  constant refinement). All qualitative results are unchanged: G-C > A-U,
  water trimer bound, helix i,i+4 H-bond forms, HH worked-example values,
  KcsA still honestly Na⁺-favored.

## Deferred (P2, flagged not failed)

- Analytic dihedral gradients (F4) — finite differences remain the oracle,
  correct by construction.

## Record (SHAs)

- Normal build output (`output.txt`): `875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935`
- ASan build output (`output.asan.txt`): `875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935`
- ASan stderr: empty (no memory errors). ASan output is byte-for-byte identical to the normal-build output (same SHA-256) - ASan instrumentation changes no physics.

## Build

```bash
cd biological/v9R3
make
./carbonsim > output.txt
```

Requires a C11 compiler and `make`. ASan-clean (empty stderr on the full
demo run, buildable through make).

**Full record** — all twelve demo outputs verbatim, complete parameter
provenance, verification discipline, and known limitations: `readme.md`.
