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
(1.38 + 1.40 = 2.78 Å) and miss Na⁺ (1.02 + 1.40 = 2.42 Å) by 0.35 Å on all
eight ligands. The dominant energetic term is desolvation, a measured bulk
quantity rather than a force-field term:

```
K+   -322 kJ/mol  ->  3.337 eV to dehydrate
Na+  -454 kJ/mol  ->  4.705 eV
K+ enters 1.368 eV cheaper
```

But the honest result is that the **rigid cage separates the two cations by
only 0.0130 eV**, because a symmetric eight-oxygen cage pulls both onto the
axis regardless of radius. The K⁺ preference at the *site* is therefore
carried essentially entirely by the measured hydration free energy.

*(Audit fix V3: this figure previously read 0.05 eV, which matched neither the
engine's per-site value nor any per-atom figure — it overstated the cage
separation by a factor of about four. The engine has always printed
0.0130 eV here.)*

**No quantitative selectivity free energy is claimed.** Turning that 0.35 Å
geometric mismatch into a binding free energy requires filter flexibility and
a solvation model, neither of which is here. A flexible-filter calculation was
implemented, found to collapse — a restraint weak enough to let the filter
respond cannot hold a +1 ion's Coulomb field — and is documented rather than
tuned away. Demo 12's constructed cage is retained and explicitly labelled as
the legacy artefact it is, so the before and after can be compared.

The poly-alanine chain in Demo 12 is still not a protein. That sentence stays
in the limitations section, because it is still true.

## 4. Streamlining

* **`include/amber_lj.h`** now holds the AMBER ff99 Lennard-Jones types for
  every biomolecular atom, and the `2^(1/6)` conversion that turns an AMBER
  Rstar into a standard collision diameter. Before this pass those numbers
  existed in four private copies, and the literal `1.6612 * 2.0 /
  1.122462048309373` was written out in **fourteen** places in `main.c`.
  `2^(1/6)` is now written exactly once in the tree. `aminoacids.c` had
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
  estimator.
* **`neuron.c` is a clean Hodgkin–Huxley implementation**, with the removable
  singularities in alpha_m and alpha_n handled by a correct Taylor expansion.

## 6. Verification

```
clean build                          0 errors, 0 warnings
selftest / forces / fire             17 / 22 / 7 checks green
selftest-regression                  161 checks green
ASan + UBSan                         0 memory errors, 0 UB, empty stderr
stdout byte-identical across runs    yes
stdout byte-identical under ASan     yes
record SHA-256                       9eb32e5489e421f89e80f85c439f4e7b5880f83c7b1c3db40c12ce31f9d53ec4
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
