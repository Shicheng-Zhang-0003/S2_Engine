# carbonsim v9R4 — release notes

v9R4 is a physics release: the quantum foundation (v4), the KcsA
selectivity program, the AGTC duplex system, and a reinforced classical
core, all verified end-to-end (record SHA below, s01 green, ASan+UBSan
clean, three selftests). It supersedes the v9R3 audit posture of
behavior-neutral hygiene: the record moved, on purpose, with every move
derived and documented.

## Highlights

- **K⁺ selectivity computed, not assumed.** Single-ion vacuum filter
  statics still favor Na⁺ (+0.2…+0.4 eV, six legs — kept as the honest
  baseline); K⁺ wins the deterministic rigid-filter exchange (−0.578 eV
  from our own explicit-water clusters), the SCF-polar transit barrier
  (−0.36 eV), and the valid-subset knock-on landscape (−1.62 eV).
  Sampled polar kinetics favors Na⁺ (+0.33 ± 0.06) — reported as an
  open bracket against the single-point barrier, not averaged away.
- **Quantum foundation v4.** Dual-sourced Clementi–Raimondi exponents
  for spatial ranges; exact H-like expectations; real f harmonics;
  lobe-max σ + sum-rule π overlap; charge-responsive screening live in
  overlap/dispersion/α; self-consistent dipoles with Hellmann–Feynman
  forces; Slater–Kirkwood dispersion; 1-pin and 2-pin SCF charge loops.
- **Systems built.** Knock-on ion pairs (relaxed + conductive +
  landscape), 6-vs-8 coordination probe, forensic table printing every
  leg side by side; duplex Phase 2 with four real glycosidic sugar
  tethers (4× C1′–N at exactly 1.47 Å).
- **Classical core reinforced.** Analytic dihedral gradients (audit P2
  closed, 22-check FD-oracle selftest); Andersen canonical thermostat
  (WHAM legs); FIRE minimizer (KcsA relax legs; steepest kept for
  clash relief after head-to-head evidence); r2 micro-opts.
- **Display split.** Deterministic stdout (the record) vs TTY-gated
  stderr (timings, heartbeats, leg progress); version banner; 13-line
  result recap; `verify_scripts.sh` gates record, both selftests plus
  fire test, datastream seal, and key/unit compliance.

## What the v9R3 audit had fixed (history, kept)

- **Q1 — quantum normalization.** `(n+l)!` cubed → first power in the
  radial normalization (2× for 2s, 6× for 2p; H 1s exact either way).
  Behavior-neutral by construction (shape-normalized plots, argmax).
- **B1–B3 — initializer braces.** Fully explicit `{{{0}}, 0, 0}`
  initializers; dead `-Wno-missing-braces` suppression removed (M1).
- **T1 — potassium row verified** against Rappé et al. 1992 UFF.
- **F1/F5 — 1-4 scaling rejected** (broke Demo 11), cutoff switching
  gated opt-in/off-default. **F4/P2** — analytic dihedrals, deferred
  then, closed now (see above).
- **N1 — naming hygiene**, v9R4 operational references.

## Gate discipline — misses the gates caught (this arc)

- Finite-difference cross-checks rejected the first analytic-dihedral
  draft (real sign bug in the m1·dn2 terms, localized per-atom with
  sympy ground truth) before it could touch dynamics.
- Undamped induction collapsed hydration minimization (−28000 eV);
  Thole damping fixed it. SCF charge feedback collapsed Na⁺ (−51 eV);
  gain/mixing/freezing fixed it. Soft-scaffold vacuum cages
  dissociated (<d>→5.5 Å); protein-like stiffness + free-ion framing
  fixed the question being asked.
- A scratch trajectory edit overflowed a fixed-size sample array; the
  stack protector aborted it, and sampling is now stride-based with
  bounds guards. A capacity audit caught the duplex Phase-2 atom count
  (127 > 96) as a segfault before release; sim sized to 160/128.
- WHAM barrier estimates flipped across builds (0.21 vs 0.61 eV) from
  tail-bin noise; count-thresholded barriers over 3 seed repeats fixed
  the estimator, and the fixed-charge reference is kept alongside.
- Every behavior-affecting change re-ran the full record through s01
  (normal + ASan builds, output SHAs, both stderrs).

## Open items (carried, precisely)

1. Bulk solvent (PBC box, Ewald/PME); 6-water clusters feed exchange
   but are not bulk. 2. Mixed-occupancy knock-on (KNa/NaK).
3. Charge-transfer covalency, exchange/correlation beyond Slater–
   Kirkwood/Pauli estimates. 4. Protein↔electrophysiology loop
   (gating from structure). 5. Full phosphodiester duplex polymer.

## Record (SHAs)

- Normal build output (`output.txt`):
  `58bb9f69f15779373cbc429998171f90deb61abc6bbd08c19c8e4c431a94134e`
  (v9R4-release: banner + all legs above.)
- ASan build output (`output.asan.txt`): byte-for-byte identical to the
  normal build (verified end-to-end: full ASan+UBSan run, empty stderr,
  SHA match). The readme's embedded output block is regenerated and
  matches this record.
- Datastream (`kcsa.cvmds`, schema 1): 100+ claims, seal verified.

## Build

```bash
cd biological/v9R4
make            # warning-clean -Wall -Wextra, no suppressions
./carbonsim > output.txt 2> stderr.txt   # + kcsa.cvmds side effect
make selftest selftest-forces selftest-fire
bash verify_scripts.sh
bash ../s01_verify_record.sh
```

Requires a C11 compiler and `make`.

---

## v9R4 Audit & Hardening Addendum (2026-09-28)

This addendum documents the comprehensive audit and hardening performed post-release:

### Mathematical Verification
- All Slater screening values verified against Slater's 1930 worked examples
- Radial wavefunction normalization verified to machine precision for n=1..3
- Most probable radii match analytical values for all hydrogen-like orbitals
- UFF LJ parameters (H–Ca) match Rappé et al. 1992 reference table exactly
- Clementi-Raimondi exponents for K dual-sourced (Wikipedia + WebElements)

### Code Quality Improvements
- Fallback warnings in `forces.c` now print once per unique missing combination
- Compile-time static assertions for all struct layouts, enums, and constants
- PCG64 RNG available as opt-in (LCG retained for record reproducibility)
- Automatic dependency generation (`-MMD -MP`) in makefile
- Static analysis clean: ASan+UBSan zero findings

### Testing
- All selftests pass (17 datastream + 22 forces + 7 fire checks)
- Record byte-identical across normal and ASan builds
- Zero memory errors, zero UB findings under sanitizers
