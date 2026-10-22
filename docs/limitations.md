# Known limitations

*S2 Engine v9R4 · stated plainly, because the numbers above cannot
be read correctly without them.*

<a id="known-limitations"></a>
## Known limitations

Stated plainly, because the numbers above cannot be read correctly without them.

* **No bulk solvent.** Everything runs in vacuum or with a relative-permittivity
  divisor; the ion program adds explicit octahedral 6-water clusters but no
  periodic water box, Ewald/PME, or membrane potential.
* **No quantitative KcsA selectivity free energy.** The real filter is rigid,
  the legs disagree, and the umbrella `±` is a three-seed spread, not a
  standard error. Read every K⁺/Na⁺ number as "in this reduced model".
* **The KcsA entropic term is implemented but does not converge, and the
  site does not survive long enough to sample.** The machinery to measure
  `Δ(βσ²/2)` exists, with side-chain + backbone dihedral torsions
  (`kcsa_add_dihedrals`) and an ion axis constraint. Measured: the 8-fold
  site holds for the first ~50 fs (CN = 8.000) and then the ion's Coulomb
  pull closes the shell — CN 6.2 at 100 fs, CN 4.0 with O–O at 3.28 Å at
  200 fs — at every temperature (50–300 K), equilibration length
  (100–2000 steps) and tether stiffness (0.5–8 eV/Å²) tested. Stiffening
  the angle term makes it worse, not better. The sampler would therefore be
  measuring the variance of a **collapsing** cage, so the site-integrity gate
  now **refuses** to return a variance for a site that was not maintained
  (`qm_kcsa_binding_series()` returns −1 and clears `ok`). Demo 12c reports
  this as UNRESOLVED with the evidence. It is *not* filled with a plausible
  number, and it cannot be filled with one by accident either.
* **Shell strain does not discriminate in this model.** The geometric contact
  mismatch is real (K⁺ −0.133 Å, Na⁺ +0.197 Å against their preferred
  distances) but the ion–oxygen Lennard-Jones wall is ~4 orders of magnitude
  too weak to turn it into an energy difference: the two ions' strains agree
  to 1.2% and that ratio is exactly invariant across a 25× stiffness scan.
* **The Born radius is not unique.** It shifts the hydration differential by
  −20.4 kJ/mol between the Latimer–Pitzer–Slansky tabulation and a fitted
  set. The Born self-energy is a cross-check on the desolvation leg, not a
  replacement for the measured TATB values, and its +18.7 kJ/mol deficit is a
  known property of the model rather than a bug to be tuned away.
* **The three tracks are unconnected.** Nucleic acids, proteins and
  electrophysiology each run standalone; gating is not derived from protein
  structure.
* **Base-pairing energetics are qualitative.** The G–C > A–T ordering is
  correct; monomer-subtracted interactions underbind gas-phase ab initio
  references (~3.5x GC, ~2.4x AT; ratio 1.5 vs 2.2 real), and no single
  dielectric fixes both pairs at once. Prior 4x-overshoot read compared
  dimer totals to interaction refs; Demo 7 now prints both.
* **Harmonic bonds cannot break — unless the EHT switch is on.** Plain
  bonds cost infinite energy at separation. With `use_eht_bo`, each bond
  scales by its overlap-gated Mayer index (1 at equilibrium, →0
  dissolved), so stretched bonds let go. Default off; Demo 18 prints both.
* **Non-nucleobase charges are approximations.** Sugar, phosphate and
  amino-acid partial charges are charge-balanced but not verified RESP fits;
  the glycine table is now exactly neutral.
* **The Pauli term still uses the heuristic overlap.** `qm_overlap()`
  dispatches to the exact kernel, but the charge-responsive overlap
  `qm_overlap_q()` — the `exp(−ζR)·ang` estimate — is what the Pauli
  repulsion is built from, and `qm_bond_order` inherits it. The EHT layer
  does NOT use it: it integrates exact STO overlaps (prolate machinery +
  documented p-pi quadrature), validated against brute-force 3D grids by
  `test_eht`.
* **QEq-lite is unscreened, and that has a hard edge.** The charge
  equilibration here zeroes the Coulomb coupling for bonded 1-2/1-3 pairs
  instead of using the Rappé–Goddard screened Coulomb integral
  `J_AB(R) = ∫ φ_A² φ_B² / |r−r′|`. The bare `1/r` that remains on every
  other pair dominates the atomic hardness `J ≈ 6 eV` once a pair is closer
  than ≈ 2.2 Å, the A-matrix goes indefinite, and the "minimum" the solver
  then returns is a saddle that puts **positive** charge on the more
  electronegative atom. Measured: an O···H pair at 2.0 Å returned
  `q_O = +0.1914 e`. The engine now **refuses** such a solution
  (`qm_qeq`/`qm_qeq_pinned` return −1, charges left untouched, one warning)
  rather than using it — see §8 and the `qm.c` Q1 note. The record is
  unaffected because every shipped QEq leg solves over a single element
  (the eight coordinating oxygens) with a pinned cation, which has no
  electronegativity ordering to violate. The guard is a *detector*, not a
  cure: the honest fix is the screened Coulomb, which needs per-element
  screening radii from the QEq parameter table that this tree does not
  carry. Inventing them would be an unsourced parameter.
* **`qm_scf_run` is deliberately not ordering-guarded.** The SCF charge
  driver uses bare `1/r` with no topological exclusion at all and a 30/70
  under-relaxation, and is a documented semi-empirical heuristic rather
  than a strict minimisation. Per-iterate and at-convergence ordering
  checks were both measured and both refuse the record's own KcsA cage
  (a Cholesky pivot of −4.39 eV on the free-oxygen block, yet
  `q_O = −0.5462 e`, which is correct). Refusing there would delete real,
  shipped physics; the ordering guard therefore applies only to the strict
  `qm_qeq`/`qm_qeq_pinned` entry points.
* **1-4 non-bonded scaling is deliberately absent** — AMBER's scaling
  presupposes its co-fitted torsions, which this force field does not use; the
  choice is documented in `forces.c` and was tested against the helix.
* **A diagnostic velocity cap is ON by default and is not NVE.**
  `max_temperature` defaults to 2000 K and rescales all velocities whenever
  the instantaneous temperature exceeds it — correct for its purpose
  (telling a heat artefact from real dynamics) but easy to trip over when
  building a test: a hand-seeded kick of 0.05 Å/fs on Mg is ≈ 24 000 K, so
  the cap fires and removes ~90 % of the kinetic energy in one step, and the
  run silently stops being NVE. Verified in the third-pass harness: two free
  particles at 0.1 Å/fs (a 8083 K state) lost 75 % of their KE in one step
  until the cap was disabled with `maxtemp 0`. Set `max_temperature <= 0` to
  turn it off; any conservation test must.
* **Cutoff switching is implemented but off by default.** Every current demo
  is gas-phase and wants the plain hard cutoff; enable the switch before any
  condensed-phase system.
* **O(N²) non-bonded loop with cell-list opt-in.** Fine at the current
  system sizes (≤ hundreds of atoms); pair exclusion is O(coordination).
  `forces_calculate` warns once above 2000 atoms. **`make perf` builds with
  cell lists** (spatial grid, cell size ≥ cutoff, O(N) pair enumeration) and
  OpenMP parallelisation over cell pairs, induction atoms, Pauli pairs, and
  dispersion pairs. Cell lists are used automatically for vacuum systems
  with N > 512 atoms when `sim->use_cell_list` is set; PBC boxes use the
  serial path (cell-list PBC minimum-image is future work).
  The PCG64 is the **reference** pcg64 of pcg-random.org (O'Neill 2014,
  XSL-RR 128/64 output) — it used to be an "S2-PCG64-like" invention with a
  64-bit increment and a made-up output permutation, and it is now pinned
  against the upstream `pcg-c` expected-output file. The LCG remains the
  record default, so the record is unaffected by the choice. QEq is QEq-lite (topological 1-2/1-3
  hard core, not Rappe shielded Coulomb) with a charge-conserving ±2e clamp.
  Thole damping uses an absolute 2.0 Å width (not AMOEBA dimensionless a).
  Hydroxyl-H LJ (0.5 Å, 0.001 kcal/mol) is a documented ff99 deviation whose
  epsilon was a 23× unit bug before the full audit (bare eV); now pinned.
* **Filter flexibility was tried and collapses** without a real protein force
  field: a restraint weak enough to let the filter respond cannot hold a +1
  ion's Coulomb field. The failure is documented, not tuned away.
* **The Roothaan–Hartree–Fock layer is qualitative.** Its overlap matrix
  is exact and its forces are a clean gradient, but the core Hamiltonian
  is the EHT VSIP diagonal, the off-diagonal coupling is a
  Wolfsberg–Helmholtz resonance integral calibrated for a different
  energy formulation, and the two-electron term is a ZDO approximation.
  The H₂ curve therefore has its minimum at 1.300 Å against the
  experimental 0.7414 Å, and **water does not produce a bound solution at
  all** in this minimal basis. The demos report both facts; the
  not-bound gate (−4) refuses to report an unbound solution. The
  *validated* electronic-structure module is the EHT one (Demo 18).
* **The screened-Coulomb kernel is not Mataga–Nishimoto.** It is a
  Hubbard-screened Coulomb that agrees with a grid oracle for the
  two-centre integral to 1.0 %. The Mataga–Nishimoto and DFTB
  interpolations would be 19 % and 26 % further from that oracle
  respectively, so the textbook substitution is *worse*, not better. The
  label is what was wrong, and it has been corrected.
* **The free-ion polarizability table is a compilation, not Pauling's own
  numbers.** The values in the tree (Li⁺ 0.03, Na⁺ 0.18, K⁺ 0.83, Rb⁺
  1.40, Cs⁺ 2.40 Å³) give the K⁺/Na⁺ ratio 4.61 that the induction leg
  uses. Pauling's own third-edition table gives 0.028 / 0.147 / 0.810 /
  1.35 / 1.94 Å³, a ratio of 5.51 and a 24 % difference for Cs⁺. The
  difference is tens of percent on a single number in an already
  qualitative leg.
* **The Pauli term still uses the heuristic overlap.** `qm_overlap()`
  dispatches to the exact kernel, but the charge-responsive overlap
  `qm_overlap_q()` (the `exp(−ζR)·ang` estimate) is what the Pauli
  repulsion is built from, and the EHT layer does not use `qm_overlap`
  at all.

---

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
