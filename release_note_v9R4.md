# v9R4 — Deep Audit Release

## What this release is

A second, adversarial audit of the whole tree, in which nothing the code
claimed about itself was taken at face value. Every analytic force was checked
against an independent finite-difference oracle, every physical constant
against its primary source, every tabulated dataset against published reference
values, and every documentation claim against the code.

The headline: **the tree did not compile**, and two independent defects in
sequence had hidden that fact. The physics underneath also had a
sign error and a functional discontinuity that were invisible in every printed
number.

## The build was broken, and the artifacts hid it

1. The first v9R4 pass added static assertions to `include/types.h` that were
   wrong: four asserted struct sizes that ignore padding, and four that
   referenced the `SimError` enum about twenty lines before it was declared.
   `make` failed on every source file.
2. The repair removed that block — and dropped a closing brace in
   `src/forces.c`, leaving `origin/main` still unable to build from a clean
   tree. Everything after line 170 was swallowed into an unclosed `if`.

Neither was visible because the repository **tracks `build/*.o` and the
`carbonsim` binary**. An in-place `make` sees the committed objects as
up-to-date, rebuilds nothing, and succeeds. `verify_scripts.sh` counted build
warnings with a bare `make`, so it reported "0 warnings" on a source tree that
could not compile. The verify script now runs `make clean` first and treats a
failed clean build as a hard failure.

This is worth stating plainly because it means the previously published record
SHA was produced by a stale binary, not by the source in the repository.

## Physics corrections

Every item was found by measurement, and each is pinned by a regression test
that fails on the old code.

**Critical**

* **The dispersion force had the wrong sign.** The energy `-C6 f6(r)/r6` is
  attractive and every printed `E_disp` was correct, but the force was applied
  as `-(dE/dr)(d/r)` instead of `+(dE/dr)(d/r)` — so the van der Waals term
  acted as a repulsion in the dynamics. Finite differences showed the analytic
  force exactly anti-parallel to the true gradient (relL2 = 2.000, the
  signature of `F = -F_FD`); correcting the sign drops the error to 5.9e-11.
* **The self-consistent polarization energy was a discontinuous function of
  geometry.** The code used a Hellmann-Feynman weight that does not apply once
  the dipoles are self-consistent (`dU/dE0 = -mu (I - T alpha)^-1`, not
  `-mu/2`), and on solver non-convergence it silently substituted a *different*
  energy functional. Measured result: a **5.92 eV jump over a 0.01 Å
  displacement**. Nothing conserved energy on that path.
* **The "SCF" was a fixed-point iteration for a linear system.** It consumed
  96–97 of its 100 iterations for a *four-atom* system; a tuned-acceleration
  replacement still failed on 121 of 401 sampled geometries. Replaced with an
  exact 128-bit linear solve.
* **Induced-dipole catastrophe cancellation on bonded pairs.** With no dipole
  hard core, a 1.3 Å C=O bond with α_C = 1.76 gives an off-diagonal coupling
  of about 2.8, pushing an eigenvalue of `(I - A)` past 1 so that no solution
  exists. Standard 1-2/1-3 exclusion applied, matching what the QEq solve
  already did.

**Major**

* **LJ σ stored UFF's Rmin rather than σ**, putting every effective Lennard-Jones
  well 12.2% too far out (carbon minimum at 4.32 Å instead of 3.85 Å). All 36
  elements converted, and now consistent with the AMBER conversions already
  used for amino acids and nucleobases.
* **Carbonyl oxygen and amide nitrogen were labelled sp3.** The classifier's
  `p_count >= 3` clause fired for every nitrogen and oxygen, overriding the
  steric logic — on exactly the functional groups the whole biomolecular model
  is made of. Hybridization now decides by steric number with a proper
  conjugation test, which is what separates an amine (sp3) from an amide
  (sp2); that distinction is invisible to a bare `Atom`, so a topology-aware
  entry point was added.
* **`r_mp³` was being used as a polarizability.** It is a screening-radius
  volume, and because the branch was selected by atomic number the induction
  energy jumped 13.6× between two adjacent elements (H 0.420 vs He 0.031). It
  also fed Pauli screening and the Slater-Kirkwood C6. Replaced with a
  tabulated set of measured atomic polarizabilities for H–Kr, stored in atomic
  units and converted in code.
* **QEq had neither a hard core nor a clamp.** Standard QEq zeroes `A_ij` for
  1-2 and 1-3 pairs; including the full 1/r term let the solve run away to
  q(C) = +4.82 e and q(O) = −5.63 e on a bonded C-O-O fragment at real bond
  lengths.
* **The minimizer divergence floor did not scale with the system.** A hard
  −50000 eV threshold silently disables minimization for anything below it — a
  protein at a routine −15 eV/atom crosses it near 3300 atoms, and
  `s10_1K4C.pdb` is in this repository. The result would have been an
  unminimized structure returned as a success. Now scaled per atom.
* **Free-energy sampling had no statistical treatment.** 300 nominal samples
  per umbrella window, autocorrelation ignored, barriers quoted to two decimals.
  Trajectories lengthened 8×, and the integrated autocorrelation time and
  effective sample size are now **measured and printed**: τ = 50–96 samples,
  N_eff ≈ 15–30 independent samples per window out of 2400 nominal. Bins
  dropped by the `min_count` filter are counted and reported rather than
  silently discarded.
* **PCG64 was documented but did not exist.** The previous notes stated "PCG64
  RNG added as opt-in (LCG retained for record reproducibility)". There was no
  PCG64 anywhere in the tree. It is implemented now, tested, and selectable via
  `sim->rng_kind`; the LCG remains the default so the record stays
  byte-identical. The uniform conversion also moved from 31 to 53 bits.

**Moderate and minor**

* Silent orbital truncation now reported through an out-parameter instead of
  quietly dropping electrons for heavy elements.
* Unclamped QEq charges on the `forces_calculate` path.
* `qm_alpha` dereferenced before its NULL check; `integrator_remove_com_velocity`
  had no guard at all.
* Five unguarded `->element->` dereferences. The periodic table stops at Kr
  (Z = 36) while `MAX_ELEMENTS` advertises 118, so `sim_add_atom(s, 37, ...)`
  yields `element == NULL` and then segfaults in five places.
* `vec3_pbc` produced NaN on a degenerate box dimension; `nb_planarity_deviation`
  read past the atom array on a bad index.
* A dead per-atom displacement clamp in both steepest-descent minimizers.
* The induction constant was hand-typed as a 2.2e-6 truncation of
  `1/COULOMB_MD`, against the codebase's own derive-in-line rule.
* `frozen[32]` stack array indexed by an atom index and read for every atom.

## Memory safety

Two out-of-bounds **writes**, both confirmed under UBSan:

```
forces.c:147  index 4 out of bounds for type 'int [4]'
forces.c:199  index 118 out of bounds for type 'int [118][118][118]'
```

The de-duplication arrays used caller-supplied bond orders and atomic numbers as
indices, and consumed 6.9 MB of BSS to implement "print once". Replaced with a
bounded key set, so no caller-supplied value is ever an index.

## Confirmed correct

Equally worth recording, because it bounds where the doubt lies.

* **The hydrogenic quantum mechanics is exact** — `∫r²|R_nl|²dr`, `<r>`, `<r²>`
  for n = 1…3, l = 0…n−1 agree with quadrature to **1e-13**.
* **Slater's rules are correct**, including the single 0.85 tier for
  `(n−1)s/(n−1)p/(n−1)d`; reproduces Slater's own Fe 4s and Fe 3d worked
  examples exactly.
* **The Clementi–Raimondi table is sound.** O 1s was suspected of a digit
  transposition (7.6579 vs 7.7579) and is correct: the 1s series rises a
  uniform 0.9925 per element from Li to Kr, and 7.6579 sits on that line. The
  suite now asserts that smoothness so a real transposition would be caught.
* **All 49 real spherical harmonics are orthonormal** for l ≤ 3, including the
  hand-normalized f₂ coefficient.
* **SHA-256 matches all three NIST vectors**; the datastream verifier correctly
  resists parser confusion and detects tampering.
* **Velocity Verlet conserves energy** below 1e-7 eV/step over 20 000 steps,
  flat in `dt`.
* **CHARMM cutoff switching, the angle gradient, the first-order induction
  gradient, and Tang-Toennies damping (including `df6/dr`) are all correct.**
* **MBAR is implemented correctly** — the problem was sampling, not the
  estimator.
* **`neuron.c` is a clean Hodgkin–Huxley implementation**, with the removable
  singularities in α_m and α_n handled by a correct Taylor expansion.

## Verification

New `make selftest-regression`: **138 checks**, one per defect, on independent
oracles so it can fail — including the KcsA filter's coordination number,
deposited distances, exact neutrality and the derived C4 symmetry assignment.
`make test` runs every gate.

```
clean build                         0 warnings, 0 errors
selftest / forces / fire            17 / 22 / 7 checks green
selftest-regression                 138 checks green
ASan + UBSan                         0 memory errors, 0 UB, empty stderr
stdout byte-identical across runs    yes
stdout byte-identical under ASan     yes
record SHA-256                       462a7be59ac55d87d0199b75260ee9761d370e9614f66ab3f06014a3d821eafa
```

The record digest changed. That is the honest consequence of fixing the
physics: a force that pointed the wrong way, a corrected LJ well, a tabulated
polarizability, a hard-core QEq, corrected hybridization labels, a direct
dipole solve and a wider uniform conversion all change trajectories. The
previous digest was produced by a binary that the committed source could not
rebuild.

## KcsA: the filter is now the real thing

This release also fixes the biological claim, which the previous note
conceded but did not repair. `src/kcsa_filter.c` builds the **actual KcsA
selectivity filter from deposited PDB 1K4C coordinates** — chain C,
THR75–VAL76–GLY77–TYR78–GLY79, the TVGYG signature motif — with the other
three tetramer subunits generated by the C4 rotation about the pore axis.

**The symmetry was derived, not assumed.** Only the correct rotation axis
reproduces a coordination number of 8 at every deposited K⁺; a wrong axis gives
2. That is the check, and it is asserted in the regression suite.

The build reproduces the deposited geometry exactly:

```
164 atoms (41 per subunit x 4), net charge +1.8e-15 e
site 1  CN = 8   <ion-O> = 2.932 A
site 2  CN = 8   <ion-O> = 2.777 A
site 3  CN = 8   <ion-O> = 2.773 A
site 4  CN = 8   <ion-O> = 2.912 A
C=O 1.219-1.242 A    C(i)-N(i+1) 1.327-1.335 A    CB-OG1 1.436 A
```

Two structural facts emerged from the coordinates that a hand-built cage
could never have produced, and both are load-bearing:

* **The innermost gate is the Thr75 side-chain hydroxyl (OG1)**, not a
  backbone carbonyl. The poly-alanine model had no such atom, so the most
  important single ligand in the real filter was absent entirely.
* **GLY79's backbone oxygen sits 4.82 Å off the pore axis and coordinates no
  ion.** GLY79 is in the motif but is not a filter ligand. Tyr78's phenol OH
  is 9.2 Å from the nearest ion — it points into the pore *wall*. Tyr79
  (canonical numbering) is a gating residue, not a selectivity residue.

**What the K⁺/Na⁺ comparison now shows.** The site is K⁺-sized: the deposited
contacts (2.77–2.93 Å) match K⁺ (1.38 Å + 1.40 Å = 2.78 Å) and miss Na⁺
(1.02 + 1.40 = 2.42 Å) by 0.35 Å on all eight ligands. The dominant energetic
term is the cost of desolvating the ion, which is a measured bulk quantity
rather than a force-field term:

```
K+   -322 kJ/mol  ->  3.337 eV to dehydrate
Na+  -454 kJ/mol  ->  4.705 eV
K+ therefore enters 1.368 eV cheaper
```

But the honest result is that in this model the rigid cage separates the two
cations by only **0.05 eV** — a symmetric eight-oxygen cage pulls both to the
same axis position regardless of radius. The K⁺ preference in Demo 12b is
therefore carried essentially entirely by the measured hydration free energy.

**What is still not claimed.** No quantitative selectivity free energy. The
real filter's 0.35 Å geometric mismatch is genuine and now measured, but
turning it into a binding-energy difference requires filter flexibility and a
solvation model, neither of which is here. A flexible-filter calculation was
implemented, found to collapse — a restraint weak enough to let the filter
respond cannot hold a +1 ion's Coulomb field — and is documented rather than
tuned away. The poly-alanine chain is still not a protein. So: the geometry is
genuinely KcsA, the mechanism is correctly identified, and the number is
reported with its provenance and its limits rather than presented as a
prediction.

Also unchanged and still open: no bulk solvent, no membrane potential, no ion
concentrations, harmonic bonds that cannot break, the three tracks
(protein / nucleic acid / electrophysiology) still unconnected, qualitative
base-pairing energetics that overshoot gas-phase references by about 4x, and
`qm_overlap` which is a heuristic rather than an overlap integral.

The `±` on the umbrella barriers is a sample standard deviation over three
seed repeats — not a standard error, not a confidence interval — and the
measured `N_eff` of 15–30 means the barriers are not determined to the digits
at which they are printed. The program now says so in its own output.
