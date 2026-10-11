# v9R4 / V0.9RC4 — Deep Audit Release

**Internal:** `v9R4`  **External tag:** `V0.9RC4`
**Repository:** https://github.com/Shicheng-Zhang-0003/S2_Engine

---

## The one-line version

The tree did not compile, the KcsA model was not KcsA, and the previous
release notes documented a compile-time assertion that failed to compile, a
PCG64 random-number generator that did not exist, and a set of Lennard-Jones
parameters that were verified only in their epsilon column. All of that is
repaired, and the KcsA filter is now built from deposited coordinates rather
than assembled by hand.

---

## 1. The build was broken, and the tracked artifacts hid it

Two defects in sequence:

1. The first audit pass added static assertions to `include/types.h` that were
   wrong four ways — four asserted struct sizes computed without padding, and
   four that referenced the `SimError` enum about twenty lines before it was
   declared. `make` failed on every source file.
2. The repair removed that block — and dropped a closing brace in
   `src/forces.c`, leaving `main` unbuildable from a clean tree, with every
   function after line 170 swallowed into an unclosed conditional.

Neither was visible because the repository **tracks `build/*.o` and the
`carbonsim` binary**. An in-place `make` sees the committed objects as
up-to-date, rebuilds nothing, and exits 0. `verify_scripts.sh` counted build
warnings with exactly that bare `make`, so a source tree that could not
compile reported "0 warnings". Both now run `make clean` first, and a failed
clean build is a hard gate.

This is worth stating plainly: **the previously published record digest was
produced by a stale binary, not by the source in the repository.** Every
digest before this release should be read with that in mind.

## 2. Physics corrections

Each was found by measurement — analytic forces against finite differences,
constants against their primary source, tabulated data against published
values — and each is pinned by a regression test that fails on the old code.

### Critical

* **The dispersion force had the wrong sign.** The energy `−C6 f6(r)/r6` is
  attractive and every printed `E_disp` was correct, but the force was
  `−(dE/dr)(d/r)` where the convention requires `+`. Finite differences found
  the analytic force exactly anti-parallel to the true gradient — relative L2
  error **2.000**, the signature of `F = −F_FD`, and nothing else. Correcting
  the sign takes it to 5.9e-11. The v3 headline term was an effective van der
  Waals *repulsion* in the dynamics while the energy printout looked perfect.
* **The self-consistent polarization energy was a discontinuous function of
  geometry.** The Hellmann-Feynman shortcut does not apply once the dipoles
  are self-consistent (`dU/dE0 = −mu (I − T alpha)^-1`, not `−mu/2`), and on
  solver failure the code silently substituted a *different* energy
  functional. Measured: a **5.92 eV jump over a 0.01 Å displacement**, so
  nothing conserved energy and the force was the gradient of no potential.
* **The "SCF" was a fixed-point iteration for a linear system.** It consumed
  96–97 of its 100 iterations for a *four-atom* system; a
  tuned-acceleration replacement still failed on 121 of 401 geometries.
  Replaced by an exact 128-bit Gaussian elimination on
  `(I − A) mu = alpha (*) E0` — no iteration count, no convergence test, no
  failure mode.
* **Induced-dipole catastrophe cancellation on bonded pairs.** With no dipole
  hard core, a 1.3 Å C=O bond with alpha_C = 1.76 gives an off-diagonal
  coupling near 2.8, pushing an eigenvalue of `(I − A)` past 1 so that no
  solution exists. The direct solve correctly reported "singular" on all 401
  geometries. 1-2 and 1-3 pairs now carry the same topological exclusion the
  QEq solve already used.

### Major

* **LJ sigma stored UFF's Rmin, not sigma**, putting every effective
  Lennard-Jones well **12.2% too far out** — carbon's minimum at 4.32 Å
  instead of 3.85 Å. Caught by recomputing the well from the engine's own
  shipped parameters, not by reading the table.
* **Carbonyl oxygen and amide nitrogen were labelled sp3.** The rule
  `steric >= 4 || (s_count > 0 && p_count >= 3)` fires for every nitrogen and
  oxygen, overriding the steric logic — on precisely the two functional groups
  the whole biomolecular model is built from. Hybridization now decides by
  steric number plus a conjugation test, which is what separates an amine
  (sp3) from an amide (sp2) and which no bare `Atom` can see, so a
  topology-aware entry point was added.
* **`r_mp^3` was being used as a polarizability.** It is a screening-radius
  volume, and because the branch was selected by atomic number the induction
  energy jumped **13.6× between adjacent elements** (H 0.420 against He
  0.031) — and the same number also fed Pauli screening and the
  Slater-Kirkwood C6, so the error reached three terms. Replaced with a table
  of measured atomic polarizabilities for H–Kr, stored in atomic units with
  the conversion in code. A Lorentz-Lorenz derivation was tried and
  **rejected**: it gives F = 87 and Ne = 143 Å³ against measured 0.557 and
  0.395.
* **QEq had neither a hard core nor a clamp.** Standard QEq zeroes `A_ij` for
  1-2 and 1-3 pairs; with full 1/r on bonded neighbours the electrostatics
  beat the hardness and the solve ran away to **q(C) = +4.82 e, q(O) = −5.63 e**
  on a bonded C-O-O fragment at real bond lengths.
* **The minimiser divergence floor did not scale with the system.** A hard
  −50000 eV threshold silently disables minimization for anything below it.
  No shipped demo comes near it, but a protein at a routine −15 eV/atom
  crosses it near 3300 atoms, `s10_1K4C.pdb` is in this repository, and the
  result would have been an unminimized structure returned as a success. Now
  −60 eV per atom.
* **Free-energy sampling had no statistical treatment.** 300 nominal samples
  per umbrella window, autocorrelation never measured, barriers printed to
  two decimals. Trajectories are 8× longer, and the integrated
  autocorrelation time and effective sample size are now **measured and
  printed**: tau = 50–96 samples, **N_eff ≈ 15–30 independent samples per
  window** out of 2400 nominal. Bins dropped by the `min_count` filter are
  counted and reported.
* **PCG64 was documented but did not exist.** The previous notes stated "PCG64
  RNG added as opt-in (LCG retained for record reproducibility)"; there was no
  PCG64 anywhere in the tree. It is implemented now — 128-bit state, xsl-rr
  output permutation, 2⁶⁴ period — tested, and selectable via
  `sim->rng_kind`. The LCG remains the default so the record stays
  byte-identical. The uniform conversion also moved from 31 to 53 bits.

### Memory safety

Two out-of-bounds **writes**, both confirmed under UBSan:

```
forces.c:147  index 4 out of bounds for type 'int [4]'
forces.c:199  index 118 out of bounds for type 'int [118][118][118]'
```

The de-duplication arrays used caller-supplied bond orders and atomic numbers
as indices, and consumed **6.9 MB of BSS** to implement "print once". Replaced
with a bounded key set, so no caller-supplied value is ever an index. Also
fixed: five unguarded `->element->` dereferences (the periodic table stops at
Kr while `MAX_ELEMENTS` advertises 118, so `sim_add_atom(s, 37, ...)` yields
`element == NULL` and then segfaults), `vec3_pbc` producing NaN on a
degenerate box (`0 * round(dr/0)`), `nb_planarity_deviation` reading past the
atom array, `qm_alpha` dereferencing before its NULL check, a fixed 32-element
stack array indexed by an atom index, and silent truncation of the orbital
table.

## 3. KcsA: the filter is now the real thing

`src/kcsa_filter.c` builds the **actual KcsA selectivity filter from
deposited PDB 1K4C coordinates** — chain C, THR75–VAL76–GLY77–TYR78–GLY79,
the TVGYG signature motif — with the other three tetramer subunits generated
by the C4 rotation about the pore axis. 1K4C was already in the repository as
`s10_1K4C.pdb`, so nothing is reconstructed from a secondary description.

**The symmetry was derived, not assumed.** Only the correct rotation axis
reproduces a coordination number of 8 at every deposited K⁺; a wrong axis or
centre gives 2. That is what makes the identification unambiguous, and it now
runs as a regression check — the geometry is a checksum on the symmetry.

```
164 atoms (41 per subunit x 4), net charge +1.8e-15 e
site 1  CN=8  <ion-O>=2.932 A      site 3  CN=8  <ion-O>=2.773 A
site 2  CN=8  <ion-O>=2.777 A      site 4  CN=8  <ion-O>=2.912 A
C=O 1.219-1.242 A   C(i)-N(i+1) 1.327-1.335 A   CB-OG1 1.436 A
```

Two structural facts emerged from the coordinates that a hand-built cage
could never have produced, and both are load-bearing:

* **The innermost gate is the Thr75 side-chain hydroxyl (OG1)**, not a
  backbone carbonyl. The poly-alanine model replaced it with a methyl, so the
  single most important ligand in the real filter was missing entirely.
* **GLY79's backbone oxygen sits 4.82 Å off the pore axis and coordinates no
  ion.** GLY79 is part of the motif but is not a filter ligand. Tyr78's phenol
  OH is 9.2 Å from the nearest ion — it points into the pore *wall*. Tyr79
  (canonical numbering) is a gating residue, not a selectivity residue.

### What the K⁺/Na⁺ comparison shows, and what it does not

The site is K⁺-sized: the deposited contacts (2.77–2.93 Å) match K⁺
(1.51 VIII + 1.40 = 2.91 Å; full-audit M3 corrects previous 1.38/2.78 VI)
and miss Na⁺ (1.18 VIII + 1.40 = 2.58 Å; was 1.02/2.42 VI) by 0.19–0.35 Å on
all eight ligands. The dominant energetic term is desolvation, a measured bulk
quantity rather than a force-field term:

```
K+   -295.3 kJ/mol  ->  3.061 eV to dehydrate
Na+  -365.3 kJ/mol  ->  3.786 eV
K+ enters 0.726 eV cheaper
```

*(Second-pass audit N1: this block previously carried -322/-454 kJ/mol and
1.368 eV. Those are the old absolute-scale hydration ENTHALPIES, not the
TATB FREE energies; the TATB free energies are Marcus 1991's -295.3/-365.3.
The enthalpy set overstated the K⁺ advantage by 62 kJ/mol and moved the
two-leg sum away from the experimental -0.1786 eV reference; with ΔG the
deviation is -0.547 eV instead of -1.190 eV. The code comment, the record,
and every document that quoted 1.368 eV now carry the free-energy value.)*

But the honest result is that the **rigid cage separates the two cations by
only 0.0098 eV** (was 0.0130 eV with VI radii; VIII gives 0.0393/4),
because a symmetric eight-oxygen cage pulls both onto the
axis regardless of radius. The K⁺ preference at the *site* is therefore
carried essentially entirely by the measured hydration free energy.

*(Audit fix V3: this figure previously read 0.05 eV, which matched neither the
engine's per-site value nor any per-atom figure — it overstated the cage
separation by a factor of about four. The engine printed
0.0130 eV with VI radii; full-audit M3 prints 0.0098 eV with VIII.)*

**No quantitative selectivity free energy is claimed.** Turning that 0.19–0.35 Å
geometric mismatch into a binding free energy requires filter flexibility and
a solvation model, neither of which is here. A flexible-filter calculation was
implemented, found to collapse — a restraint weak enough to let the filter
respond cannot hold a +1 ion's Coulomb field — and is documented rather than
tuned away. Demo 12's constructed cage is retained and explicitly labelled as
the legacy artefact it is, so the before and after can be compared.

*(Second-pass audit: the sentence this paragraph defended — "the protein is
poly-alanine, so the real TVGYG selectivity chemistry is not present" — is
no longer true and has been removed from the limitations. Demo 12's legacy
constructed cage never contained a poly-alanine chain in its final form, and
Demo 12b is the real deposited TVGYG filter, held rigid. The limitation that
survives is the one that matters: no quantitative selectivity free energy is
claimed.)*

## 4. Streamlining

* **`include/amber_lj.h`** now holds the AMBER ff99 Lennard-Jones types for
  every biomolecular atom, and the `2^(1/6)` conversion that turns an AMBER
  Rstar into a standard collision diameter. Before this pass those numbers
  existed in four private copies, and the literal `1.6612 * 2.0 /
  1.122462048309373` was written out in **sixteen** places in `main.c`
  (fourteen carbonyl-sigma expressions plus two JC macros; the second-pass
  audit measured nineteen literal occurrences across the four files, where
  the header's own note had said thirteen). `2^(1/6)` is now written exactly
  once in executable code (`TWOPOW_SIXTH`); the second pass also removed the
  one straggler copy in `tui.c`. `aminoacids.c` had
  admitted the duplication in its own comment — "duplicated here rather than
  shared via a header refactor to avoid touching already-validated, working
  code under time pressure" — a debt never paid.
* **The consolidation immediately paid for itself.** The KcsA filter module
  had hardcoded the carbonyl-oxygen sigma as **3.06615** when the AMBER ff99
  value is **2.95992**, and the regression test had copied the same wrong
  literal — so the module and its test agreed with each other and both were
  wrong. The shared constant fixed both, moving the filter energies by
  ~0.0014 eV.
* Release identity is defined once, as `S2_VERSION_INTERNAL` and
  `S2_VERSION_EXTERNAL`, and consumed by the program banner and the datastream
  header. The previous release carried a hand-typed source string in
  `main.c` that had been stale for several releases.
* `verify_scripts.sh` now cleans before building, and its regression count no
  longer counts the suite's own summary line as one of its own checks.

## 5. Confirmed correct

Equally worth recording, because it bounds where the doubt lies.

* **The hydrogenic quantum mechanics is exact** — `∫r²|R_nl|²dr`, `<r>`, `<r²>`
  for n = 1…3, l = 0…n−1 agree with Simpson quadrature to **1e-13**.
* **Slater's rules are correct**, including the single 0.85 tier spanning
  `(n−1)s`, `(n−1)p` and `(n−1)d`; reproduces Slater's own Fe 4s and Fe 3d
  worked examples exactly.
* **The Clementi–Raimondi table is sound.** O 1s was suspected of a digit
  transposition (7.6579 vs 7.7579) and is correct: the 1s series rises a
  uniform 0.9925 per element from Li to Kr, and 7.6579 sits on that line. The
  suite now asserts that smoothness, so a real transposition would be caught.
* **All 49 real spherical harmonics are orthonormal** for l ≤ 3, including the
  hand-normalized f₂ coefficient.
* **SHA-256 matches all three NIST vectors**; the datastream verifier resists
  parser confusion and detects tampering.
* **Velocity Verlet conserves energy** below 1e-7 eV/step over 20 000 steps,
  flat in `dt`.
* **CHARMM cutoff switching, the angle gradient, the first-order induction
  gradient, and Tang-Toennies damping** (including `df6/dr`) are all correct.
* **MBAR is implemented correctly** — the problem was the sampling, not the
  estimator. *(Full-audit M9: solver is WHAM iteration, not MBAR; fixed in prose.)*
* **`neuron.c` is a clean Hodgkin–Huxley implementation**, with the removable
  singularities in alpha_m and alpha_n handled by a correct Taylor expansion.

## 6. Verification

```
clean build (-Werror -ffp-contract=off)  0 errors, 0 warnings
selftest / forces / fire             17 / 22 / 7 checks green
selftest-regression                  165 checks green
selftest-external                    51 checks green (49 + HO-eps unit pins)
selftest-loop                        60 checks green
ASan + UBSan                         0 memory errors, 0 UB, empty stderr
stdout byte-identical across runs    yes
stdout byte-identical under ASan     yes
record SHA-256                       see CURRENT_BASELINE_SHA.txt (do not duplicate digests in prose)
```

`kcsa.cvmds` is regenerated with each run and is deliberately **not**
byte-stable: `ds_open` writes a live wallclock unless `SOURCE_DATE_EPOCH` is
set, which its own comment admits. The baseline digest is the *stdout* digest,
which `verify_scripts.sh` checks — it is not the `.cvmds` digest, and the two
sit next to each other in a way that invites the wrong conclusion.

## 7. Open

* No bulk solvent, no membrane potential, no ion concentrations. A
  relative-permittivity divisor stands in for condensed-phase screening.
* 1-4 non-bonded scaling is deliberately absent — a non-standard choice,
  documented in place.
* Harmonic bonds cannot break.
* The three tracks (protein / nucleic acid / electrophysiology) remain
  unconnected.
* Base-pairing energetics are qualitative: the G–C > A–U ordering is
  validated, the absolute magnitudes overshoot gas-phase references ~4×.
* `qm_overlap` is a heuristic, not an overlap integral — the Slater polynomial
  prefactor is dropped and `m` is maximised independently per atom.
* The umbrella barriers are reported to four decimals with a three-seed
  spread, and the measured N_eff of 15–30 per window means they are not
  determined to that precision. The program says so in its own output.

---

## 8. Full-audit follow-through (post-v9R4 operation)

A second, deeper pass (118 findings: 47 math/physics, 34 programming, 37
operational) verified every equation against its canonical external source
and repaired what it found. Truth changes, all pinned by tests:

* **HO epsilon unit fix (CRITICAL).** Hydroxyl-H LJ was `0.001` bare-eV
  against a documented `0.001 kcal/mol` (23.06× too deep, ~4.8× into every
  HO–X pair via the geometric mean). Now `(0.001 * KCAL_MOL_TO_EV)` with a
  `test_external` pin. Record moves honestly: helix PE 9.689 → 9.596 eV,
  H-bond 2.1634 → 2.1204 Å, polar-WHAM gap −0.0422 ± 0.2332 → +0.1155 ±
  0.0886 eV (sign flips inside spread — still UNDECIDED, still no
  selectivity claim; verdicts HELD/BREATHING unchanged).
* **PCG64 claim narrowed.** The stream is O'Neill-inspired, not reference
  pcg64 (64-bit increment, non-reference output permutation). Renamed
  S2-PCG64-like in code and docs; LCG stays the record default.
* **QEq/Thole relabelled.** QEq-lite (topological hard core, not Rappe
  shielding) and absolute-width Thole damping are now named as what they are.
* **Force-layer trust boundaries.** Bond/angle/dihedral/pair entry points now
  take `num_atoms` and reject upper-bound OOB (wild read + wild force write
  on corrupt topology); LJ-energy diagnostic mirrors the P9 guard; angle
  forces gate finiteness; `kick_clamped` zeroes non-finite kicks instead of
  integrating inf; per-sim temperature-cap latch; constrained-DOF validation.
* **Fail loud.** QM caps (SCF 128, polar 256, coupled-SCF 64) warn-once on
  stderr instead of substituting U=0 at exit 0; `main` returns 1 on demo
  failure; temperature-cap diagnostics are TTY-gated so the record stderr
  stays empty.
* **TUI hardening.** Create-then-swap `new`, full-consumption numeric parsing
  with INT range checks, unsigned `$(( ))` arithmetic (no UB), lossy-save
  WARNING (S2SAVE1 drops dihedrals/angles/box/flags/RNG), velocity finiteness,
  restraint-cap-bounded loader, `S2_NO_HOST` restricted mode with documented
  trust boundary, TMPDIR validation, overlap-safe `shift`, man stubs for
  `test/help/tput/fc/ps-keys/vi-keys`.
* **Build.** `-ffp-contract=off` (FMA fusion pinned, not just `-march` banned)
  and `-Werror`; datastream `record-tree-<git-sha>` source-hash, `\r`/bracket
  hygiene, `fseeko`/`ftello`, NULL-safe seal search; `verify_cvmds` C tool
  calls `ds_verify_file` as primary with python as second opinion.
* **Operations.** `./run` stages `runs/<stamp>.{txt,cvmds}` and promotes
  atomically, defaults `SOURCE_DATE_EPOCH` to the tree commit time, and
  refuses `--accept` unless `make test` is green; `make run` execs `./run`;
  `verify_scripts.sh` asserts exact gate counts (17/22/7/165/51/60), header
  presence, provenance/unit hygiene, no-digest-in-prose, and a TUI smoke;
  spec schema-1 frozen with `computed-<leg>` grammar and the sign convention
  carried once in the Demo 12 header; `run_audit.sh` builds hermetically and
  requires ≥4 red checks on revert; history scripts carry HISTORICAL headers.

Verification at promotion: `make test` 6/6 green, `./verify_scripts.sh`
VERIFY PASSED (incl. TUI smoke + seal), record digest `see
CURRENT_BASELINE_SHA.txt`, stderr empty, `s01` ASan/UBSan byte-identity per
`s01_verify_record.sh`.

---

## 9. EHT semiempirical SCF (post-full-audit operation)

The quantum layer outgrows its heuristic overlap (`qm_overlap`,
abstraction #7) with a real electronic-structure method: Hoffmann-1963
Extended Hückel + self-consistent charge (`qm_eht.h`, Demo 18, `test_eht`
31/31 green as the seventh gate).

* **Exact overlaps.** Prolate-spheroidal STO machinery (Rosen/Roothaan):
  closed A-sums + GL quadrature, heteronuclear-exact, p-pi by documented
  2D quadrature. Every channel validated against brute-force 3D Cartesian
  grid integration of independently-coded STOs; the grid caught two real
  bugs pre-record (an `a^{2k-n-1}` exponent slip at 2.5x overshoot, a
  double-counted radial power on p-channels).
* **Honest SCF.** VSIP Hamiltonian, WH K=1.75, Löwdin orthogonalization
  (Jacobi eigensolver, closed-2x2 oracle), Aufbau+Fermi fill, Mulliken
  charges, Mayer bond orders, gaps. H2: 2x2-exact MOs, BO ~1, neutrality
  to 1e-9, band dissociates to 2xH. H2O: O −0.54/H +0.27, gap 14.8 eV,
  rotation-invariant to 1e-6.
* **Three measured failures kept as findings, not tuned away.** E10: the
  SCC sum initially excluded on-site Hubbard terms and ran away to full
  ±2.0 e transfer at ANY separation (fixed: U on-site restores hardness;
  eq water went −1.86 → −0.54, SCC 31 → 6 iters). Fermi annealing for
  stretched-heteronuclear flips was implemented and REVERTED same-day
  (smoothed roots drift with kT; lowering re-triggers flips, once AT
  resid 1.4e-07). E11: stretched heteronuclear level-crossing is rc=-3
  fail-loud, gated as such. E5/E7: band-only energies (bare 1/R sums
  unclaimed — water reads +50 eV otherwise), EHT-large gaps, no interior
  H2 minimum on the grid (argmin at R=4.0 is the asymptote, not a well).
* **Reactive milestone.** Overlap-gated Mayer switch (E8/E12): H-H factor
  1.0 → 0.007 from equilibrium to 4 A; at 2.0 A the BO-scaled bond
  carries 7.6 eV against the harmonic 28.5 eV wall. `use_eht_bo` is
  opt-in, default off; caps fail loud (C9 pattern).
* **Build hygiene (E9).** Test `.d` files were written but never
  `-include`d, so a `types.h` field left `test_regression.o` stale and
  silently skewed struct layouts (PCG64 reads degraded to LCG, one red
  gate, zero source errors — caught only because a test compares two
  streams). Test deps are included now; clean-first remains the rule.

---

## 10. Third audit pass (2026-10-10) — independent oracles, four defects fixed

A full mathematical / programming / operational re-audit of the whole tree.
The method that matters: it did **not** reuse the project's own oracles. A
separate verification harness was built outside the repository (in `/tmp`,
linked against the shipped engine objects, 78 checks) using independent
derivations — central finite differences, closed-form hydrogenic results,
statistical-ensemble identities, and machine-readable transcriptions of the
primary parameter tables. A shared mistaken assumption could therefore not
pass twice. Where the harness and the tree disagreed, the harness was
checked against the primary source before the tree was judged.

**Confirmed clean (the majority of the tree).**

* Every SI literal in `include/constants.h` matches CODATA 2019 exactly,
  and every derived conversion re-derives in Python from its primaries.
* UFF Lennard-Jones: **36/36** σ and ε verified cell-by-cell against an
  independent transcription of Rappé Table II. Two values that secondary
  sources commonly mis-report are **correct here**: Si σ = 3.82641 Å
  (x₁ = 4.295, not 4.195) and Ar σ = 3.44600 Å (x₁ = 3.868, not 3.912).
* AMBER ff99: **13/13** atom types match `parm99.dat` `MOD4 RE` exactly.
* Analytic forces vs central differences on every term — including CHARMM
  cutoff switching, flat-bottom restraints, PBC minimum-image pairs, the
  coupled 3N×3N dipole solve, Pauli and dispersion, none of which the
  shipped suites covered. Worst deviation **5.2e-10 eV/Å**.
* Integrator is symplectic: harmonic oscillator conserves E exactly over
  20 000 steps; dimer energy error is bounded (late/early ratio 0.98), not
  secular; Maxwell–Boltzmann, Andersen and Langevin all reproduce
  ⟨v²⟩ = 3kT/m; Berendsen relaxes 1000 K → 300 K.
* Hodgkin–Huxley: resting gating matches the published 1952 worked example
  to four decimals, both removable singularities take their L'Hôpital limits
  exactly, the action potential peaks at +40.3 mV.
* KcsA: Marcus 1991 TATB free energies, Shannon 1976 VIII radii, and
  CN = 8 at all four deposited 1K4C sites verify.
* SHA-256: FIPS 180-4 known-answer vectors pass; a tampered payload byte is
  rejected.

**Four real defects found and fixed.**

### Q1 — QEq could invert electronegativity silently *(physics)*

`qm_qeq` / `qm_qeq_pinned` build the augmented system with a **bare 1/r**
off-diagonal (the Rappé–Goddard screened Coulomb integral is not
implemented). With this engine's own parameters the A-matrix goes indefinite
once a non-excluded pair is closer than ≈ 2.2 Å — for an H/O pair the
smallest eigenvalue is −8.15 eV at 1.0 Å, −3.35 eV at 1.5 Å, −0.95 eV at
2.0 Å, +0.49 eV at 2.5 Å. Indefinite is not singular, so the Gaussian
elimination returned a perfectly well-conditioned "solution" that satisfied
`sum(q) = total_q`, respected the ±2 e bound, and placed **positive charge
on oxygen**:

```
O...H at 2.0 A:  q_O = +0.191402 e   q_H = -0.191402 e
```

at exit status 0, flowing straight into every downstream Coulomb term.

A positive-definiteness requirement was implemented and measured against the
record's own ion cage, where the free-oxygen block has a Cholesky pivot of
**−4.39 eV** — genuinely indefinite, and yet it produces the correct
`q_O = −0.5462 e`, because a many-oxygen shell is dominated by the diagonal
hardness even when the block is indefinite. Demanding PD would have deleted
real, shipped physics. The gate that ships instead tests the quantity that
actually broke: **electronegativity ordering**, as a group-mean
monotonicity check. A linear O–C–O–C–O chain drove its central oxygen to
the +2 e bound while the terminal oxygens sat at −1.38 e, which a *pairwise*
rule misreads as an inversion; the group means (O −0.25 e, C +0.38 e) are
correctly ordered, and the group-mean form passes it. Both the refusal and
the non-refusal are pinned. `qm_scf_run` is deliberately **not** guarded: it
is a documented semi-empirical under-relaxed heuristic, and gating it
refuses the record's own cage.

### Q2 — the quantum layer ignored periodic boundaries *(physics)*

`forces.c` applied the minimum-image convention in its pair loop; `qm.c`
contained **zero** calls to `vec3_pbc_box`. QEq, the Thole-damped field,
the coupled-dipole tensor, Pauli and dispersion therefore all used raw
Cartesian separations. In vacuum the two are identical, which is why it
survived; in a periodic box (the TUI petri dish is a 48 Å periodic world)
every quantum-enhancement term was wrong by a geometry-dependent factor, and
— worse — the two halves of one potential energy disagreed about a pair's
distance, so the reported PE could not be differentiated into the reported
forces. All pair terms now go through one `qm_sep()` helper. Two pinned
checks place unit charges on opposite faces of a 6 Å box and require the
QM field to see the 1 Å separation. **No recorded number changed.**

### Q3 — the datastream seal did not reach the end of the file *(integrity)*

`ds_verify_file`'s comment claimed the hash line was required to be final
"no trailing claims", but the code only checked that 64 hex digits were
followed by a newline or EOF. Content appended after that newline is
outside the sealed payload *and* outside every structural check, so it
verified clean — a claim that no digest covers. Measured before the fix:
append `"injected.claim 999 eV computed"` to a sealed file →
`ds_verify_file() = 0` (accepted). The seal is now required to be the last
content, with the rejection pinned by a new datastream check.

### Q4 — scratch paths leaked a tool-specific directory *(operational)*

`run`, `verify_scripts.sh` and `audit/external/hydro.py` hard-coded
`/tmp/opencode` — a scratch directory belonging to one tool, not a general
temporary location. On any machine without it, `./run` and the verifier
silently created a stray `opencode` tree under `/tmp`, and `hydro.py` died
with `FileNotFoundError` on its very first `open()` because it never created
the directory. Now `${TMPDIR:-/tmp}` with an explicit `mkdir -p`, and a
`tempfile.mkdtemp` the script owns and removes.

**Gates.** `test_datastream` 17 → 19, `test_regression` 165 → 170, for a
total of **360** pinned checks across seven suites, all green. The record
reproduces **byte-identically** (`a20cb3a7…`) with empty stderr, because
none of the four fixes changes a number the record contains: the QEq guard
only refuses geometries the record never visits, minimum image is a no-op in
vacuum, the seal change only rejects files that were never valid, and the
path change is invisible to the build.

---

## Part 2 — the KcsA ionic layer: Born self-energy, the contact mismatch, and the entropic term

This part answers a standing request — search everything relevant outside the
tree, then build the quantum parts far enough to settle KcsA. The short
version is that the settling came as a **negative result with a positive
mechanism attached**, which is the most this engine can honestly deliver, and
here is what it is.

### The literature, and what it changed

The decisive source is Dixit & Asthagiri, *J. Chem. Theory Comput.* **5**,
1037 (2009), with Asthagiri, Roux & Cross, *Biophys. J.* **95**, 1269 (2008).
Their result is that the previous KcsA story was the wrong kind of story. The
**mean** binding energy slightly *favours Na⁺* — Na⁺ is simply the better fit
for an eight-carbonyl cage — and selectivity comes from the **fluctuation**
term `Δ(βσ²/2) ≈ 5–6 kcal/mol`, whose physical source is carbonyl–carbonyl
repulsion strain of roughly 7.3 kcal/mol. Na⁺ reaches its favourable
contacts only by squeezing the site, and pays for that entropically.

So a model that reports a binding energy and a hydration free energy and
subtracts them has been describing the part of the mechanism that is
approximately *zero*. That is the single most important thing found in this
pass, and it is why the new module exists.

Supporting external facts now used and cited: Doyle et al., *Science* **280**,
69 (1998) (K⁺ at least 10⁴ more permeant than Na⁺ in the generic channel);
LeMasurier et al., *J. Gen. Physiol.* **118** (2001) (KcsA `P_Na/P_K ≤ 0.006`);
B(h)yskov, Asthagiri et al., *JCTC* (2010) QM/MM coordination numbers — K⁺ 6.6,
Na⁺ 5.2, first-peak ion–O 2.93 vs 2.55 Å; Latimer, Pitzer & Slansky (1936)
Born radii, corroborated for Na⁺ by Hummer, Pratt & García,
arXiv:physics/9709042 (1.80 Å there, 1.88 Å in Marcus's later compilation);
Pauling's free-ion polarizabilities; Gong, Hocky & Freed, *PNAS* **105**
(2008) for the protein-interior ε ≈ 4.

### N1 — the Born ionic layer *(new: `src/qm_born.c`)*

Until now the desolvation leg was a transcribed table with nothing in the code
able to check it, vary it, or explain it. There is now a Born self-energy to
compare it against, derived from the CODATA primaries the engine's own
electrostatics uses.

* Differential ΔG(Na⁺−K⁺): Born −51.3 kJ/mol against the measured TATB
  −70.0, a deficit of **+18.7 kJ/mol**. That is the known Born limitation —
  a linear-response continuum has no room for cavity and dispersion — and it
  is reported as a model property, not tuned away.
* The Born radius is not unique. The fitted alternative set (Valisko, reported
  by Antonijevic, *Electrochem. Sci.* **15**, 1511 (2010)) shifts the
  differential by a further −20.4 kJ/mol. That is now a returned number,
  not a caveat.
* Induction had to be **screened**, and the screening is load-bearing. At the
  ~5.2 e/Å a cation sees in the cage, unscreened induction is −42 eV for K⁺
  and −9 eV for Na⁺ — four orders of magnitude too large to be an ion–ligand
  binding energy, and it would have silently dominated any total it entered.
  With ε = 4 it is −0.70 and −0.15 eV: real, correctly signed, still too
  small to decide anything alone.

### N2 — the contact mismatch *(new, and the sharpest KcsA result here)*

A geometric question with a geometric answer. On the pore axis the mean
ion–oxygen distance over the eight coordinating oxygens is **2.777 Å**, and
that is the *shortest the deposited geometry permits*. Against each ion's own
preferred contact:

| ion | prefers | achieves | mismatch |
|---|---|---|---|
| K⁺ | 2.910 Å | 2.777 Å | −0.133 Å (squeezed inside) |
| Na⁺ | 2.580 Å | 2.777 Å | +0.197 Å (held outside) |

The cage is built for K⁺. The larger ion is squeezed 0.13 Å inside its
preference; the smaller is held 0.20 Å outside its preference and cannot get
closer anywhere on the axis. **That is the exclusion mechanism, and it needs
no force field at all.**

Three placements were tried on the way to this and each is instructive:

1. *Relax the ion's radius to the potential minimum* — returns **on-axis for
   both ions**, same 2.777 Å, because the ion–O interaction here is a
   0.05 kcal/mol Lennard-Jones wall with no minimum and the Coulomb term
   simply parks the cation on the axis. Strains then agree to 1%: one
   measurement taken twice. This also revealed that **the 0.36 Å K/Na contact
   difference `kcsa_set_ion_radius` exists to express is never realised
   anywhere in this model** — relax the ion freely and it leaves at its own
   contact distance immediately.
2. *Impose each ion's Shannon contact by bisecting* — impossible on the axis,
   and the bisection silently pinned both ions back at r = 0, losing the
   mismatch entirely. It also produced a two-ligand complex reporting a
   perfect eight-ligand contact distance, because the mean was taken over
   "whatever is in range now" while the ion was moving. The ligand set has to
   be frozen at the deposited geometry.
3. *Allow the ions off-axis* — they both stay on-axis, and strains still agree
   to 1.2%, a ratio of **1.0122 invariant across a 25× stiffness scan**
   because strain goes as 1/k and the ratio does not.

Enabling the exchange (overlap-Pauli) and dispersion terms moves the
difference by under 0.1%. So the sharp negative result: **in this model the
entire K⁺/Na⁺ discrimination is the tabulated hydration free energy.** The
geometric mismatch is real, measurable, and does not convert into an energy.

### N3 — the entropic machinery, and its honest verdict *(new: `src/qm_fluctuation.c`)*

`qm_fluctuation` runs the deposited filter as a breathing peptide with the ion
inside and measures the variance of the ion's binding energy over an explicit
Langevin NVT trajectory — the `Δ(βσ²/2)` term that the whole selectivity
argument turns on and that a rigid cage cannot produce.

**It does not converge.** The filter holds CN = 8 through minimisation
(2.777 Å) and loses it within a few hundred steps of Langevin dynamics at
every temperature from 50 K to 300 K and every Cα restraint stiffness from
0.5 to 16 eV/Å²: CN falls to 3–5 and the O–O diagonals close from the
crystallographic ~4.5–5.5 Å to ~3.5 Å. The cause is structural, not
numerical — **the engine has no side-chain torsion terms, and four of the
eight filter oxygens belong to Thr75's side chain.** Without torsions the
hydroxyl leaves the shell. No thermostat and no restraint recovers it.

The leg is therefore reported **UNRESOLVED**, with the evidence attached and a
`site_maintained` gate that fails closed rather than printing a number from a
filter that came apart. Closing it needs a real protein force field and
explicit solvent, not a smaller timestep.

Three implementation defects were found and fixed on the way, each of which
would have produced a confidently wrong number:

* **The variance was silently zero** whenever the caller did not request the
  series buffers — the moments were computed in a second pass over an
  optional array and fell back to the running mean, turning the module's
  entire reason for existing into a constant. Now Welford streaming, which
  is also the numerically better form here given binding energies of tens of
  eV fluctuating by tenths.
* **The coordination diagnostic counted hydrogens.** With the side-chain H
  completed, an unrestricted count inside the same cutoff returns CN = 10.5
  (the Thr75 OG1–H and Tyr78 OH–H are not ligands), so it now counts
  oxygens, matching `kcsa_coord_stats`, without which the two numbers cannot
  be compared.
* **`kcsa_add_hydrogens` needed an undocumented bond capacity.** Asked for the
  filter's own `KCSA_FILTER_BONDS*4 + 16`, it fills the 164 heavy-atom bonds
  and then fails on the 18th hydrogen, having already appended 17 — a
  Thr/Val/Gly fully hydrogenated and a TYR missing one aromatic hydrogen. It
  fails closed, but the requirement is now stated in the header.

### N4 — `sim_zero_strain_bonds` *(engine fix, `src/sim.c`)*

The bond and angle tables are keyed on element and bond order and are
parameterised for the small gas-phase fragments this engine was originally
built for. Turning on `use_bonds` against a *deposited* structure charges it
bond-length strain it never had, in a direction set by gas-phase constants
rather than by anything about the ion. On the KcsA filter that strain was
worth **11.7 eV** of potential energy and crushed the 8-fold site to CN 4
with the shell diagonals closing from ~4.5–5.5 Å to ~3.7 Å — every one of
those numbers the deposit being deformed by the force field rather than by
the ion.

`kcsa_site_binding` solves this by disabling the bonded terms entirely, which
is right for a rigid cage and useless where the strain mechanism is the
subject: the filter has to be able to be strained. `sim_zero_strain_bonds`
is the third option and the one this work needed — set every bond's `r0` to
the distance its atoms actually have, leaving `k` alone. The bonded terms then
are live, and zero at the reference geometry, so what they report afterwards
is deviation from experiment.

**Gates.** `test_external` 60 → 105, `test_regression` 170 → 202, for a
total of **392** pinned checks across seven suites, all green. Every new
external check names its source: the Born coefficient against the closed form
built from CODATA, the Born radii against Latimer–Pitzer–Slansky with the
Na⁺ value corroborated by Hummer/Pratt/García, the Pauling polarizabilities
cross-tabulated, and the 2.777 Å site mean re-derived from the deposition
rather than read back from the engine's own bookkeeping. The regression suite
additionally pins the defect class this whole part is about — a variance that
is *identically* zero because a code path skipped the sampling would pass
every check that only looked at a mean.

---

## Part 3 — two defects found by reading the record rather than trusting it

Part 2 shipped a working demo and a green tree. Both of the defects below
were present in it, and neither was caught by any gate that existed, which
is the part worth recording.

### V1 — the verdict slot collided, and the recap silently lied *(main.c)*

Demo 12c was given verdict slot 13, because 12 was the last slot anyone had
looked at. **Demo 17 already owned slot 13.** Demo 17 runs after 12c, so
12c's summary was overwritten before the recap ever read it, and the row
labelled `17 duplex` was the only thing visible.

This is precisely the failure `AUDIT FIX V2` in `src/main.c` was written to
prevent. Its comment says the next demo added cannot quietly disappear from
the recap — and the mechanism does not deliver that. The bounds check stops
an out-of-range **slot**; two demos writing the same **in-range** slot is
perfectly in range, and no amount of array sizing catches it. That is what
"necessary but not sufficient" means, and demo 12c demonstrated it inside
the same commit series that added the protection.

The symptom was invisible by construction. The recap table printed
perfectly, every row was populated, and one row simply showed the wrong
demo's text. Nothing looked broken.

Fixed: `DEMO_VERDICT_COUNT` 15 → 16, `names[]` gains `12c ionic layer`, and
duplex and EHT move to slots 14 and 15. Every demo now owns a distinct slot
and the label sits in the same position as the slot.

**And gated**, because the fix that was supposed to prevent this
demonstrably did not. `verify_scripts.sh` now parses `src/main.c`, resolves
each `set_verdict` call to its enclosing demo function, and fails if any
slot has more than one owner, if `names[]` is not exactly as long as
`DEMO_VERDICT_COUNT`, or if the distinct slot count differs from the
declared count. The three together mean every declared slot has exactly
one writer.

The gate was checked against the bug, not only against the fix:
reintroducing the 12c-into-13 collision in a scratch copy of `main.c` turns
it from 0 to 1. A gate that has only ever seen a passing tree is not known
to work.

### V2 — demo 12c contributed no claims to `kcsa.cvmds` *(main.c)*

The sidecar is the sealed provenance ledger for the KcsA program and held
**zero** of demo 12c's quantities. Nothing detected it: `verify_scripts.sh`
checks the seal, the schema, the key case, the units and the provenance
vocabulary, but there was no claim-*coverage* check, because the writer and
the demo it was supposed to be recording for were never compared.

**30 claims added, 117 → 147:** the Born radii and self-energies, the Born
differential with its deficit against the measured table, the radius-choice
spread, both transfer energies, the screened induction pair, the full shell
set (CN, preferred contacts, achieved distances, mismatches, forces,
strains, and the stiffness-independent ratio), and the Dixit entropic
targets. Provenance tags follow the file's own convention — published
transcriptions are tagged as transcriptions (`Latimer-Pitzer-Slansky1936`,
`Pauling-free-ion`, `Shannon1976-VIII+1.40A`) rather than as `computed`,
because a reader checking this file against outside references needs to
know which numbers to check where.

A claim that **cannot** be computed is sealed explicitly rather than
omitted: `kcsa.entropic.computed` reads `0 kcal/mol` tagged
`UNRESOLVED-no-sidechain-torsions`, beside its published target
`kcsa.entropic.target` of 5.5 kcal/mol. An absent claim reads as an
oversight; a zeroed claim with a reason reads as a measurement that failed,
which is what happened.

### The structural problem underneath V2

`ds_open` buffers the entire payload and truncates at `ds_close` — there is
no append path and no reopen — and Demo 12 runs *before* Demos 12b/12c. The
writer therefore cannot reach a value demo 12c has not computed yet. The
obvious fix, hoisting demo 12's hundred-odd locals to file scope so the
whole block can move later, is a large refactor of working code for no
benefit.

The cheap correct answer is that demo 12c's measurement is a **pure
deterministic function of the deposited coordinates**. Factored into
`kcsa_ionic_measure` and called from both the printer and the writer, the
two call sites cannot disagree even in principle, at a cost of ~40 ms
against a record run measured in minutes. The printed contact table now
reads from the same struct the claims are written from, rather than
recomputing beside it — which is how the two could have drifted in the
first place.

That determinism is load-bearing, so it is pinned: 22 checks in
`test_external` assert every value both call sites share, against the same
outside references the claims carry provenance for. If the measurement ever
stops being deterministic, it breaks a gate instead of producing two
plausible, disagreeing numbers that nobody compares.

One of those 22 is deliberately a tripwire rather than a regression:
`kcsa.shell.strain_ratio` must stay within 5% of unity. The ratio being
near 1 **is** the negative result of this release, so a future change that
quietly starts discriminating has to update that line on purpose and say so
in the commit message.

### The `-dirty` stamp, and what it was telling us

The previous `kcsa.cvmds` read `source-hash: record-tree-07af540…-dirty`.
That suffix was honest — it says the file was produced from a working tree
with uncommitted changes, which it was — but it was avoidable. The code
fixes were committed first, the tree was brought fully clean, and only then
was the record regenerated, so the new sidecar reads
`record-tree-1702c9b…` with no suffix. The inherent one-commit lag remains
and is a property of the design rather than a defect: a stamp cannot name
the commit that contains its own output.

**Gates.** `test_external` 105 → 127 for **468** pinned checks across seven
suites, all green, with three new structural checks in `verify_scripts.sh`
covering verdict-slot ownership. The record's stdout diff against part 2 is
a **single added line** — the recap row demo 12c was owed — so no physics in
the record moved; only its self-description was incomplete.

---

## Part 4 — auditing what IS there, not what is missing

Parts 1-3 went after absent physics and against the record's own
self-description. Part 4 asked the opposite question: of the KcsA
machinery already in the tree, what is wrong, or wrong in some parameter
regime. Everything below was found by probing the built engine rather than
by reading it, and two of the initial suspicions turned out to be wrong and
are recorded as such.

### N2 — `kcsa_cation_radius` failed OPEN *(src/kcsa_filter.c)*

`default: return 1.51` handed **potassium's** radius to any ion not on the
list. Measured: Z = 1, 8, 26, 27, 54 and 86 all returned 1.51 Å.

The failure was silent *and* plausible, and it was compound. The two
functions that together define what this program thinks an ion is —
`kcsa_cation_radius` and `kcsa_hydration_free_energy_kJmol` — failed in
**opposite directions**. The radius fell back to potassium; hydration
correctly returned 0.0. An unlisted cation therefore got a full-size
potassium Lennard-Jones site *and* a zero desolvation cost, which is a
complete, coherent-looking, wrong answer rather than an error.

The pair now fails the same way: 0.0 for an unlisted ion from both, with
every caller required to refuse rather than propagate. `kcsa_set_ion_radius`
needed a guard of its own for a subtler reason — with a zero radius the
contact computes to 1.40 Å, which drives `sigma_ion` negative and was
silently clamped to a 0.05 Å floor. A value that looks like a deliberate
number and is not one. `qm_kcsa_elastic_strain` got the same guard.

### N3 — the Joung–Cheatham parameters had no oracle *(kcsa_filter.h)*

Four numbers drove a full row of the KcsA forensic table, carried the
external provenance tag `jc2008-params` in the sealed `kcsa.cvmds`, and
were validated by **no gate anywhere in the tree** — confirmed by grepping
all seven suites and the audit harness. Every other outside number in this
program has an independent check: AMBER against the vendored `parm99.dat`,
Marcus hydration against its primary table, Shannon radii against the
deposited geometry, NIST against CODATA.

They are now `KCSA_JC_*` macros in `include/kcsa_filter.h`, beside the
Shannon radii and the hydration table they are compared against, and 23
checks in `test_external` pin all four against Joung & Cheatham, *J. Phys.
Chem. B* **112**, 9020 (2008), doi:10.1021/jp8001614. The checks include the
orderings a transcription slip would break — K⁺ well deeper than Na⁺, K⁺
R\* larger — and the K/Na gap pinned independently of the absolutes.

**A caveat that was implicit and is now written down.** These put the
ion–oxygen minimum at 3.366 Å (K⁺) and 3.030 Å (Na⁺), while the
Shannon-based contacts are 2.910 and 2.580 Å — a ~0.45 Å disagreement in
absolute distance. They agree on what decides selectivity (gap 0.336 Å vs
0.330 Å, within 2%), so it is a constant offset rather than a disagreement
about the ions. But the JC leg and the real-filter leg of the table are
measured at **different geometries**, and both that fact and the 2% are now
pinned.

### N4 — the 1000:1 reference scale was asserted, not sourced *(main.c)*

`KCSA_EXPT_RATIO` was commented "Experimental K+/Na+ selectivity ~1000:1 for
KcsA" with no citation, and the sealed file tagged the result
`expt-1000:1@300K`, which reads as a measurement. It is not one. The
published numbers disagree with each other and neither is 1000:

* LeMasurier, Hille & Horn, *J. Gen. Physiol.* **118** (2001) 341 — KcsA
  `P_Na/P_K ≤ 0.006`, i.e. **≥ 167-fold**. Channel-specific, conservative.
* Doyle et al., *Science* **280** (1998) 69 — K⁺ "at least **10⁴** times
  more permeant than Na⁺". A potassium channel generally, and a lower bound.

The defensible statement is a bracket of ~1.7×10² to 10⁴. **The constant
keeps its value** — moving every reference-scaled number in the record for
no evidential gain, when a bracket is what the literature supports and a
point is not. What changed is the label (now a chosen representative), the
sealed tag (`representative-in-bracket-167..10000@300K`), and a printed
span showing the bracket runs from -0.238 to -0.132 eV and is 59% wide
relative to the representative point in ln(ratio), so no leg is read
against that scale more tightly than the evidence allows.

### N5 — two different oxygen polarizabilities, one unstated *(main.c)*

`KCSA_POL_ALPHA_O = 0.84` is a **carbonyl**-oxygen value (Applequist),
deliberately not the engine's free-oxygen 0.7854. Different physical
quantities, 7% apart — not an error in either. But nothing said so, which
left a reader comparing 0.84 against the engine's own table with no way to
tell whether one was stale.

It stays 0.84: this is an induction *proxy* for a carbonyl oxygen, and
substituting the free-atom value would make it consistently wrong instead
of correctly-sourced. A gate now pins both numbers **and** pins that they
differ by less than 10%, so a future edit to either surfaces the question
rather than letting the two drift apart silently.

### N6 — the coordination cutoff, stated as a threshold *(kcsa_filter.h)*

Callers pass 3.4 Å, and that single number is what "the filter is
8-coordinate" means throughout the program. It is a threshold, not a
measurement, and it reclassifies silently as the shell deforms — which is
exactly what happens in the flexible legs, where CN falls from 8 to 3–5.
The header now says so, and `test_external` pins both the CN = 8 at the
deposited sites and the fact that the on-axis distance is a **minimum**,
which is the part the exclusion argument actually rests on.

### Checked and CLEARED — recorded so they are not re-audited

* **`kcsa_set_ion_radius` places the LJ minimum exactly on target** —
  2.3200 / 2.5800 / 2.9100 Å, error 0.0000. My first probe reported it off
  by 0.23 Å; I had written `σ^1.2` where the LJ minimum is
  `2^(1/6)·σ`. The probe was wrong, not the code.
* **`qm_alpha` is wrong for cations and it does not matter.** Na⁺ 6.5× too
  small, Cs⁺ 4.3× too *large*, K⁺/Na⁺ ratio 18.6 against Pauling's 4.61 —
  but it is diagnostic-only. Every energy path uses `qm_polarizability`,
  which already handles cations correctly (K⁺ 0.83, Na⁺ 0.18).
* **Cα indices {1,8,15,19,31}** are all carbons with ≥2 bonds.
* **The filter is exactly neutral** per segment and as a tetramer
  (+0.000000 e).
* **The golden section returning r = 0** is the physical on-axis answer,
  not an ion escaping a bracket edge.
* **The vendored `parm99.dat` cannot oracle the JC values** — it carries
  ion hydration radii, not ion Lennard-Jones. The oracle has to be the
  paper, which is why the gate transcribes it.

### Negative testing

Every gate above was checked against a mutation, not only against the fix.
Reintroducing the exact N2 fail-open in `kcsa_cation_radius` — and
verifying the build still succeeded, which is the step that made the first
attempt meaningless — turns **five** checks red, including the
`set_ion_radius` guard firing on the observed σ change from 2.5587 to the
2.2251 clamp. A gate that has only ever seen a passing tree is not known to
work.

**Gates.** `test_external` 127 → 150, for **491** pinned checks across
seven suites, all green; `verify_scripts.sh` passes. The record's stdout
diff against part 3 is the new bracket line printed beside the 1000:1
reference and nothing else — no physics in the record moved.
