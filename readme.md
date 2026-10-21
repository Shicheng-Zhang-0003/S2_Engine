# S2 Engine — Nano Chemistry & Biological Simulator

**v9R4** (internal) · **V0.9RC4** (external) · C11 · record path: no dependencies beyond libm
Repository: https://github.com/Shicheng-Zhang-0003/S2_Engine

S2 is a self-contained computational chemistry and biology engine. It places
molecules atom by atom with real condensation chemistry, runs classical
molecular dynamics and a quantum-enhancement layer on them, simulates a
Hodgkin–Huxley neuron, models the KcsA potassium-channel selectivity filter
from deposited coordinates, and seals its results into a hash-verified
record. Three tracks exist — nucleic acids, proteins, electrophysiology —
each validated independently; they are not yet connected to each other.

One rule governs everything here: **every number is sourced, every
approximation is flagged, and every negative result is reported as one.**
Where a demo's answer is "this reduced model cannot decide that", the output
says so in the same sentence as the number.

---

## 1. Build, run, test

Requirements: a C11 compiler (gcc or clang), `make`, libm. Nothing else
for the record. The fullscreen `screen` surface optionally uses ncursesw
(detected via `ncursesw6-config`; without it `s2tui` still builds and
every line command works, only the alternate-screen reports unavailable).

| Command | Effect |
|---|---|
| `make` / `make -j$(nproc)` | build everything: `carbonsim`, `s2tui`, and all eight test binaries (serial record build, byte-deterministic) |
| `make perf` | build everything with **OpenMP multicore** (`-fopenmp -DS2_OPENMP -march=native`) — cell-list O(N) pair search + parallel O(N²) kernels across all cores. The perf build has its **own object directory (`build-perf/`) and its own binaries (`carbonsim-perf`, `s2tui-perf`)**, so it can never contaminate the record build. Performance runs only; **never for the record** (thread scheduling changes FP summation order, so digests differ from the serial record). |
| `./carbonsim` | run the 20-demo record sequence (**~7 minutes**; the polar-WHAM leg dominates); stdout is the byte-deterministic record |
| `./run` | clean build behind a warning gate, run to `runs/<timestamp>.txt`, print old/new digests side by side |
| `./run --accept` | promote that run to `output.txt` + `CURRENT_BASELINE_SHA.txt` + `kcsa.cvmds` (staged atomically: each artefact is written to a temporary name and renamed into place, and the ledger is written last); refuses to promote unless stderr is empty, `make test` is green, **and** the staged `kcsa.cvmds` seal verifies |
| `make test` | build and run the eight unit gates only (19/22/7/202/153/60/31/43) — CI must run `make test && ./verify_scripts.sh && ./s01_verify_record.sh` (plus `make audit` for evidence) |
| `make tui` | build just the interactive terminal (`./s2tui`) |
| `./verify_scripts.sh` | read-only record + spec verifier: baseline digest, `kcsa.cvmds` seal, key/unit compliance, every suite |
| `./s01_verify_record.sh` | builds normally and under ASan/UBSan, asserts both reproduce the recorded digest with empty stderr |
| `make audit` / `make audit-revert` | run the audit's independent oracles and C harnesses; `--revert` links the same regression suite against the pre-fix engine |
| `make clean` | remove `build/`, both binaries, `audit/build` |
| `s2tui` | interactive terminal; `help` prints the command map |

The record build uses `-O3 -g -Wall -Wextra -Werror -ffp-contract=off -std=c11 -Iinclude` and never
`-march=native`: `-ffp-contract=off` pins FP codegen (the `-march` ban alone does not stop FMA
fusion under `-O3` on FMA hosts), and `-Werror` makes warnings fail the build.
Set `SOURCE_DATE_EPOCH` for byte-stable side artifacts (`./run` defaults it to the tree's own
commit time so `kcsa.cvmds` seals identically across repeated runs).

**Current gate status** (all green):

| Gate | Checks | What it pins |
|---|---:|---|
| `test_datastream` | 19 | FIPS 180-4 SHA-256 known-answer vectors, round-trip, tamper rejection, seal-is-final |
| `test_forces` | 22 | analytic dihedrals vs finite-difference oracle, net-zero forces, collinear guard |
| `test_fire` | 7 | FIRE and steepest descent reach the same minima |
| `test_regression` | 202 | independent in-repo oracles: quadrature, finite differences, hand-derived values, QEq ordering guard, QM minimum image, zero-strain bonds, KcsA shell strain |
| `test_external` | 150 | values from outside this repository (NIST CODATA, parm99.dat, FIPS, PDB 1K4C, Latimer–Pitzer–Slansky Born radii, Pauling free-ion polarizabilities) + HO-eps unit pin |
| `test_loop` | 60 | bio/QC/QM loop closure + best-layer: Nernst, kT, Marcus/Shannon/LJ/Coulomb, BJ/penetration, Langevin, flat-bottom, tethers, H-complete ions-hold |
| `test_eht` | 31 | exact overlaps vs 3D-grid STO oracles, Jacobi vs 2x2 forms, H2/H2O SCF laws, fail-loud gates |
| `test_audit4` | 43 | the fourth audit: reference PCG64 vs the upstream expected-output file, EHT mirror symmetry, textbook overlap closed forms, Born induction units, one-electron integral limits, SCF invariants, dipole-coupling scale, cell-list origin, datastream non-finite refusal, L'Hôpital limits, KcsA hydroxyl angle |

plus a warning-free clean build (serial and `make perf`), empty stderr,
byte-identical stdout across repeated runs (serial record only), an intact
`kcsa.cvmds` seal, and the OpenMP build passing all eight gates with
numerically consistent results (within FP reassociation tolerance).

**537 checks in eight gates.** `make test` fails if any suite's exit code
is non-zero **or** if any suite prints a `FAIL` line — both conditions are
checked for every suite, not just some of them.

**CI.** `.github/workflows/ci.yml` runs the whole chain on every push and
pull request: the warning-free build, `make test`, `./verify_scripts.sh`,
`make audit`, `./s01_verify_record.sh`, the OpenMP build and its eight
gates, and finally

```
git diff --exit-code -- output.txt CURRENT_BASELINE_SHA.txt kcsa.cvmds
```

That last step is the one that matters. The October-2026 commits that added
demos 13–16 all carry the message *"all 491 gates green"*, and every one of
them ran `make test` only: the **record** gate was never run, so the
tracked `output.txt` drifted away from HEAD while every commit message
claimed a green tree. The staleness check is now a CI step, so a record
that does not match its source fails the build instead of being
rediscovered months later.

---

## 2. The record and its artifacts

The record is the tracked stdout of `./carbonsim`, and it is the authority
this document describes — not the other way around.

| Artifact | What it is |
|---|---|
| `output.txt` | the verbatim 20-demo record (~1250 lines) |
| `CURRENT_BASELINE_SHA.txt` | its SHA-256 (single source of truth — do not duplicate digests in prose) |
| `kcsa.cvmds` | a schema-1 datastream written by Demo 12: 147 worked-example claims with provenance tags and a SHA-256 payload seal (`DATASTREAM_SPEC.md`). Covers every KcsA leg including Demo 12c's Born layer and contact mismatch |
| `runs/` | candidate runs; untracked. Promotion is explicit (`./run --accept`) precisely so an accidental overwrite cannot masquerade as a record update |
| `audit/` | the audit's own evidence tree: independent math oracles, C harnesses, external references, pre-fix engine copies |

Never write `output.txt` by hand. The stdout/stderr split is deliberate:
stdout is DATA (the record — never wall time, never progress), stderr carries
human display (timings, heartbeats, leg progress) and only when stderr is a
TTY; `CARBON_QUIET` silences it. Pipes and redirects therefore stay honest.
See `include/display.h`.

---

## 3. Architecture

Unit system: **Å, fs, eV, AMU, elementary charge.** All fundamental constants
are 2019 CODATA values, and every conversion factor is derived in-line from
primaries (`include/constants.h`) rather than hand-typed, so reciprocals are
exact by construction.

Current tree: **34 498 lines** (src 26 220 · include 3 847 · tests 4 431).

| File | Lines | Role |
|---|---:|---|
| `src/tui.c` | 7380 | interactive terminal: full POSIX shell + `screen` fullscreen twin, chemistry lab, gas/barostat (not in the record) |
| `src/main.c` | 5223 | the 20-demo record program and datastream consumer (Demo 7 A-T + monomer-subtracted legs; fail-closed exit status) |
| `src/qm.c` | 3291 | QEq-lite (topological hard core, not literal Rappe shielding) + induced dipoles (first-order + coupled SCF), Pauli, TT/BJ dispersion, penetration-damped Coulomb, overlap, SCF driver, **exact STO overlaps (ss/sps/pps/ppi), two-electron integrals (Hubbard-screened Coulomb + Wolfsberg-Helmholtz exchange), Roothaan–Hartree–Fock SCF solver + FD forces** |
| `src/nucleobases.c` | 1279 | five bases, deoxyribose, T-p-A dinucleotide, pairing/geometry helpers |
| `src/forces.c` | 1365 | non-bonded pairs (upper-bound OOB-checked), bonded terms, harmonic/flat-bottom restraints, analytic dihedral gradients, energy breakdown; QM caps fail loud; **OpenMP cell-list pair search + private-force accumulators** |
| `src/sim.c` | 1089 | lifecycle, molecule constructors, topology rebuild (P6-guarded both twins), ion/restraint plumbing (32 restraints), **cell-list build for O(N) non-bonded pairs** |
| `src/integrator.c` | 924 | Velocity Verlet, **reference PCG64** (O'Neill 2014, XSL-RR 128/64) / Knuth-MMIX LCG, Maxwell–Boltzmann, Berendsen/Andersen/exact-OU-Langevin, steepest descent, FIRE, **standard C11 64x64->128 multiply** |
| `src/quantum.c` | 547 | Slater screening and orbital energies, hydrogenic radial profiles, Clementi–Raimondi exponents, real spherical harmonics |
| `src/kcsa_filter.c` | 964 | real 1K4C TVGYG filter, ion sizing, sites, binding/dehydration legs, C-alpha tethers, H-completion (26/subunit), **backbone + side-chain dihedral torsions** |
| `src/qm_born.c` | 129 | Born ionic self-energy, Latimer–Pitzer–Slansky radii, Pauling free-ion polarizabilities |
| `src/qm_fluctuation.c` | 775 | variance/entropy machinery, KcsA shell strain, site-integrity gate, **ion axis constraint for fluctuation sampling** |
| `src/aminoacids.c` | 513 | glycine/alanine/dipeptide/polyalanine builders |
| `src/datastream.c` | 445 | schema-1 writer, self-contained FIPS 180-4 SHA-256, seal verifier (fseeko/ftello, record-tree source-hash, \r/bracket hygiene) |
| `src/periodic_table.c` | 323 | H–Kr element data, UFF ε/σ, Madelung electron configurations |
| `src/tui_view.c` | 383 | viewport cell model shared by the ANSI frame and the curses screen |
| `src/tui_screen.c` | 171 | ncursesw alternate-screen monitor (optional dep, TTY-gated, piped-safe) |
| `src/neuron.c` | 173 | Hodgkin–Huxley 1952 (squid giant axon) |
| `src/loop.c` | 84 | bio/QC/QM loop closure: Nernst, kT scale, LJ/Coulomb from same primaries |
| `tests/test_regression.c` | 1831 | the 202-check in-repo oracle suite |
| `tests/test_external.c` | 1087 | the 150-check outside-reference suite |
| `tests/test_loop.c` | 273 | the 60-check loop-closure + best-layer suite |
| `src/qm_eht_overlap.c` | 467 | exact STO overlaps: prolate A/B machinery, GL quadrature for p-pi, **and `eht_overlap_ao()` — the single shared global-axis AO overlap kernel** |
| `src/qm_eht_scf.c` | 697 | Hoffmann EHT + SCC (Jacobi, Lowdin, Mulliken, Mayer, trust mixing) |
| `tests/test_eht.c` | 331 | the 31-check EHT oracle suite |
| `tests/test_audit4.c` | 519 | the 43-check fourth-audit suite (one check per defect found and fixed) |
| `include/qm_eht.h` | 128 | EHT API + E1-E11 abstraction inventory |
| `include/loop.h` | 55 | loop bridge contract + abstraction inventory |
| `include/` | 3847 | types, constants, per-module contracts; `amber_lj.h` and `display.h` carry shared policy |
| `tests/` (other) | 378 | datastream, forces, fire suites |

### 3.1 Multicore and cell lists (performance build)

**The engine was single-threaded through v9R3.** `make perf` (v9R4+) adds OpenMP multicore
parallelisation plus cell-list spatial partitioning:

| Component | Serial (record) | `make perf` | Scaling |
|-----------|----------------|-------------|---------|
| Non-bonded pairs (LJ+Coulomb) | O(N²), one core | **O(N) cell list + OpenMP over cells** | Near-linear to core count |
| Induction (polarisation) | O(N²), one core | **OpenMP over atoms** | Near-linear |
| Pauli repulsion (FD forces) | O(N²)·6N evals, one core | **OpenMP over pairs** | Near-linear |
| Dispersion (TT-damped C6/r⁶) | O(N²), one core | **OpenMP over pairs** | Near-linear |
| QEq charge equilibration | O(n³) Gaussian elimination, ≤128 atoms | Serial (small n, memory bound) | N/A |

**Cell lists** (`sim_rebuild_cell_list`): a spatial grid with cell size ≥ cutoff. Build is O(N);
pair enumeration visits only the same and 26 neighbouring cells, reducing the pair loop from
O(N²) to O(N + N_neighbours). The list is rebuilt automatically by `forces_calculate` when
`sim->use_cell_list` is set and the atom count exceeds 512 (PBC boxes use the serial path;
cell-list PBC minimum-image is future work).

**The perf build is isolated from the record build.** It has its own object
directory (`build-perf/`), its own binaries (`carbonsim-perf`, `s2tui-perf`)
and its own build-identity strings, which are what the sealed datastream
records. Before the fourth audit the perf target shared `build/*.o` and the
record binary names, which produced three separate failures: `make perf`
after `make` was a **silent no-op** (nothing was recompiled with `-fopenmp`,
and the "multicore" binary contained no OpenMP code at all), `make clean &&
make perf && make` **failed to link** the record binary, and a perf run
**overwrote the tracked record binary**. The record build and the
performance build can no longer touch each other.

**Thread safety:** every OpenMP region uses per-thread private force accumulators that are
reduced after the parallel region, so there is no race on `atoms[].force`. Energy sums use
OpenMP reduction (floating-point reassociation — same total, different order), which is why
the perf build is **not byte-identical** to the serial record and must never be used for
record production.

**Record path:** `make` (default) builds without `-fopenmp`, so `#pragma omp` is a no-op,
`S2_OPENMP` is 0, and the record is bit-identical to before. The OpenMP code is compiled in
but dormant.

---

## 4. What the engine models

### Quantum orbital layer
Slater screening (`Z_eff = Z − S`, Slater's own 1930 iron example reproduced),
orbital energies from Slater `n*`, hydrogenic radial profiles `P(r) = r²|R_nl|²`
with exact expectation values, Clementi–Raimondi SCF exponents for spatial
ranges, and all 49 real spherical harmonics (orthonormal through l = 3).
Demo 1 includes a per-element computed-vs-NIST ledger instead of a blanket
accuracy claim.

### Exact two-center STO integrals and Roothaan-Hartree-Fock SCF

**One overlap kernel, grid-validated.** `eht_overlap_ao()` in
`src/qm_eht_overlap.c` computes the overlap of two Slater-type AOs given
in **global Cartesian axes** - `(n, zeta, l, m)` per AO, with `m = 1, 2, 3`
for p_x, p_y, p_z - from the exact prolate-spheroidal integrals (Rosen
1931; Roothaan 1951). It is the single place where the Slater-Koster
direction cosines and the s-p sigma sign are applied, and both the EHT
solver and the Roothaan-Hartree-Fock solver assemble their overlap
matrices with it. Before the fourth audit they each had their own copy,
and the two copies **disagreed**:

* the EHT copy carried a spurious minus sign in its "p on the first atom"
  branch, so the overlap matrix of a symmetric molecule was not symmetric
  under its own mirror. For a bent H-C-H it produced
  `S(s_H1, p_x) == S(s_H2, p_x)` where the physics demands
  `S(s_H1, p_x) == -S(s_H2, p_x)`. The EHT SCC then **failed to converge
  for H-C-H at all**, and N2 came out with the wrong MO ordering, a
  HOMO-LUMO gap of 0.0000 eV and a Mayer bond order of 2.053 instead of
  the textbook 3.000;
* the RHF copy resolved each ATOM's valence shell instead of each AO's own
  `(n, l)`, and for C/N/O that shell is 2p - so the s-channel wrapper
  refused every O-H, C-H and N-H pair and the **2s-1s sigma-bonding
  channel of water and methane was exactly zero** in the overlap matrix;
* the RHF copy also never populated the principal quantum number per AO,
  so the VSIP core diagonal always fell back to the first ionisation
  potential and the O 2s diagonal came out at -13.598 eV instead of
  -32.3 eV.

Channels: s-s, s-p sigma, p-p sigma, p-p pi (2p only). The textbook
closed forms (e^{-zR}(1+zR+(zR)^2/3) for 1s-1s, e^{-p}(1+p+2p^2/5+p^3/15)
for 2p-2p pi) are pinned in `tests/test_audit4.c` to 1e-12 and 1e-10
respectively, and the kernel is validated against brute-force 3D Cartesian
grid integration in `tests/test_eht.c`.

**Two-electron integrals.** The two-centre Coulomb kernel is a
**Hubbard-screened Coulomb**,

```
gamma_ab(R) = COULOMB_MD / sqrt( R^2 + (28.8 / Ubar_ab)^2 ),  Ubar_ab = (U_a + U_b)/2
```

with `U = IE - EA` (the chemical hardness). It is continuous with the
on-site term (`gamma_aa -> U_a`) and tends to `COULOMB_MD/R` at large R.
The kernel is **externally validated**: the true two-centre Coulomb integral
<2p_O | 1/r | 1s_H> for the engine's own valence densities at the water
O-H distance (0.9572 A) is **5.7156 eV** (grid integration of the actual
charge densities in atomic units), and the kernel returns 5.7720 eV - a
**1.0 % deviation**. The bare 1/R Coulomb would be 15.0 eV (+162 %).

> The kernel used to be labelled "Mataga-Nishimoto", and it was not one:
> Mataga & Nishimoto (1957) and its DFTB descendant interpolate as
> `1/(R + 0.5/(k_A+k_B))`, and the measured deviations from that form are
> -35 % at 0.5 A and -13 % at 1 A. Substituting the textbook form would
> make the model *worse* (+19 % against the oracle), so the kernel was
> kept and the **label was corrected**. The tree now has ONE screened
> Coulomb kernel and it is the one with an external check.

The exchange is the Wolfsberg-Helmholtz resonance integral
`beta_ab = -K_WH * S * (H_aa + H_bb)/2` with `K_WH = 1.75` (Hoffmann
1963), and `qm_eri_zzdo()` builds the RHF Fock contribution `J - 0.5*K`
from a diagonal density matrix under ZDO.

**The Roothaan-Hartree-Fock SCF solver** (`qm_hf_scf()`) assembles S from
the shared kernel and H_core from the Hoffmann VSIP diagonal, solves
F.C = S.C.eps by Loewdin orthogonalization + a Jacobi eigensolver, fills by
Aufbau, and iterates the density with damped mixing. Fail codes are loud:
-1 bad input, -2 singular overlap, -3 non-convergence, and **-4 the
solution is not bound**. Minimal valence basis: H 1s; C/N/O 2s+2px+2py+2pz.

**Honesty block H1-H5** (in `include/qm.h`): this module is *internally
consistent and numerically sound* - the overlap matrix is exact, the density
is normalised, and the finite-difference gradient of a smooth energy
surface is accurate to 1e-5 eV/A - but it is a **qualitative**
electronic-structure layer, not a predictive one:

* the core Hamiltonian is the EHT diagonal `H_ii = -VSIP`; the exact
  one-electron integrals are implemented and grid-validated but cover 1s-1s
  only and are **not** what this solver uses;
* the off-diagonal core element is the Wolfsberg-Helmholtz *resonance*
  integral with `K = 1.75`, a constant calibrated for the EHT
  band-energy formulation, not for an RHF trace expression;
* the two-electron term is the ZDO approximation;
* the consequence is that the H2 curve has its minimum at **R = 1.20 A**
  against the experimental 0.7414 A, and water in this minimal basis does
  not produce a bound solution at all;
* the validated electronic-structure module in this tree is the **EHT one**
  (Demo 18), which converges, gives the correct H2O Mulliken charges
  (O -0.54 e) and reports its own negative result rather than hiding it.

The demos say all of this out loud, and the not-bound gate refuses to
report a solution whose lowest orbital is unbound.

### Wiring quantum SCF into classical MD (`use_qm_scf`)
`qm_hf_scf_forces()` computes the SCF electronic energy and its
finite-difference gradient with respect to atomic positions (same FD pattern
as `qm_pauli_forces()`: 6*N extra SCF solves per force call). It is exposed
as an opt-in force term in `forces_calculate()` via `sim->use_qm_scf`:

* **Default off (0):** the record path is completely classical - every
  existing demo and the byte-deterministic record are unaffected.
* **Opt-in (1):** the SCF energy is added to `potential_energy` and reported
  as `E_qm_scf_total`; finite-difference forces are accumulated onto each
  atom alongside LJ, Coulomb, bonded, angle, dihedral, Pauli, dispersion,
  and induction terms. Independent of `use_scf` (QEq charge loop) - both
  can be set simultaneously.
* Supports H/C/N/O atoms within the minimal valence basis (<=128 atoms);
  larger or unsupported systems return rc=-1 and U=0.
* **Fail-closed:** a displaced SCF solve that fails leaves that force
  component at zero and says so. It used to leave the energy at 0.0 - a
  *finite* number - which the finite difference consumed as a real energy
  and turned into a force of order E/2h, about 1.5e6 eV/A. That was an
  integrator-explosion generator, and it was silent.

Because the SCF is O(nao^4) per solve and forces need 6*N solves per step,
this is intended for small quantum-MD test systems, not production-scale
dynamics.

**What the quantum-MD demos actually show.** Demo 13 (H2) works: the H2
relaxes towards the SCF minimum and the force falls from 19 to 1.6 eV/A.
Demo 14 (H2O) does **not**: the SCF is not bound for water in this basis,
so the demo says so, reports the quantum-MD claim as UNRESOLVED, and runs
the classical trajectory as the reference. (Before the audit the same demo
ran the classical trajectory while *printing* a quantum-MD table, and the
water exploded: the identical trajectory without the quantum term is stable
at 250-330 K for 500 steps, and with it the temperature reached the 2000 K
cap at step 300.)

### Classical molecular dynamics
* Velocity Verlet (symplectic) with **Berendsen** weak-coupling (steering;
  legacy record paths) and **Andersen** stochastic collisions (rigorously
  canonical; the WHAM free-energy legs).
* Lennard-Jones (Lorentz–Berthelot) + Coulomb with an optional dielectric
  divisor; per-atom LJ overrides so jointly-parameterized models (TIP3P water,
  AMBER ff99 biomolecular atoms) stay internally consistent; harmonic
  restraints enter the force loop like every other term; per-axis PBC and a
  CHARMM-style cutoff switch are implemented (the switch is opt-in).
* Bonded terms: harmonic stretch, harmonic bend, and **dihedrals with exact
  analytic chain-rule gradients**. The finite-difference oracle that proved
  them and a 22-check agreement suite remain in-tree.
* Temperature uses the correct `3N − 3` degrees of freedom, with further
  constrained-DOF accounting. Two RNGs: the historical Knuth-MMIX LCG
  (default — keeps the record byte-stable) and the **reference PCG64**
  (O'Neill 2014, XSL-RR 128/64 output) — see §6.1.
* Minimization: steepest descent (monotonic, clash relief) and FIRE
  (smooth-basin polishing), with a per-atom displacement cap and a
  size-relative divergence floor.

### Quantum-enhancement layer
Every term is opt-in, reported as its own ledger entry, and never folded into
the base charges: induced dipoles (first-order, and self-consistent
coupled-dipole solves with Thole-damped coupling and Hellmann–Feynman forces),
overlap Pauli repulsion, Slater–Kirkwood dispersion with Tang–Toennies
damping, and a QEq + dipole SCF charge-equilibration loop (1-pin and 2-pin).
Hard bounds are documented and fail closed: coupled dipole solve ≤ 64 atoms,
QEq ≤ 128, first-order polar ≤ 256.

### Semiempirical SCF layer (EHT)
Hoffmann-1963 Extended Hückel with self-consistent charge (Demo 18,
`qm_eht.h`): exact prolate-spheroidal STO overlaps (grid-validated, not
the `qm_overlap` heuristic), VSIP Hamiltonian with Wolfsberg-Helmholtz
coupling, Löwdin-orthogonalized SCF, **Hubbard-screened Coulomb gammas**
(the U values come from the in-tree Mulliken J; the screening constant and
the hardness were both wrong before the fourth audit — the old kernel
screened with the mean ionisation potential and returned 10.09 eV for the
water O–H pair against a grid oracle of 5.7156 eV, a 76 % error; the
current kernel returns 5.7720 eV, a 1.0 % error), Mulliken charges, Mayer
bond orders, HOMO–LUMO gaps computed from the **occupations the solver
actually used** (the previous index-based formula reported the splitting
*inside* a degenerate pair, and N₂ came out with a gap of 0.0000 eV).
Single-zeta minimal basis (H/C/N/O), ≤ 64 AOs, fail closed. An opt-in BO
switch scales harmonic bonds by an overlap-gated Mayer interpolant, so
bonds weaken and let go instead of pulling to infinity (E6–E12 inventory
in `qm_eht.h`).

### Biopolymer condensation chemistry
Real chemistry, not decoration: leaving groups are genuinely removed following
the PDB Chemical Component Dictionary's own flags. Nucleic acids: the five
bases from PDB CCD ideal coordinates with Aduri et al. 2007 RESP charges and
AMBER ff99 LJ types, tautomer-corrected cytosine, deoxyribose, and a T-p-A
dinucleotide with a real phosphodiester bridge (net −1 e). Proteins: glycine,
alanine, a Gly–Ala dipeptide with a real 1.33 Å peptide bond, a poly-alanine
chain builder, and an α-helix whose i,i+4 backbone hydrogen bond emerges from
Coulomb+LJ alone.

### Electrophysiology
The 4-variable Hodgkin–Huxley 1952 model, RK4-integrated (Euler would need an
impractically small step during the sodium upstroke), with the removable rate
singularities handled by Taylor expansions. Resting gating values reproduce a
published worked example to four decimals.

### KcsA ion-channel program
Two tracks, deliberately:

* **Demo 12 — legacy constructed cage** (kept as the honest vacuum baseline):
  four carbonyl oxygens in 4-fold symmetry, single ion on-axis, fixed-radius
  sites, a radius scan, and a two-ring antiprism with the ring z-separation
  sourced from 1K4C (3.084 Å).
* **Demo 12b — the real filter**: the deposited PDB 1K4C chain C TVGYG
  structure — 164 atoms, net charge 0.000 e, C4-symmetric, CN = 8 at all four
  sites, Thr75 OG1 as the innermost gate, GLY79's carbonyl correctly excluded
  (4.82 Å off-axis). It is rigid by construction; the module documents why a
  flexible-filter calculation collapses without a real protein force field and
  a solvation model.

Each ion is pushed through every leg and the results are printed side by side
in a forensic table rather than summarized into one number. Current record,
in that table's sign convention (`+ = Na⁺ favored`):

| Leg | Value | Direction |
|---|---|---|
| vacuum site legs (point/JC/SCF/v2/relaxed/stiff) | +0.21 … +0.75 eV | Na⁺ |
| computed exchange, rigid filter | −0.5762 eV (deterministic) | K⁺ |
| computed exchange, relaxed filter | −0.1557 ± 0.9041 eV | UNDECIDED (within spread) |
| two-ion exchange | −0.2706 eV | K⁺/KK |
| SCF-polar U(z) barrier, d(K−Na) | +2.5912 eV | Na⁺ pays less |
| polar-WHAM barrier gap | +0.1155 ± 0.0886 eV | UNDECIDED (within spread) |
| fixed-charge WHAM gap | +0.1558 eV | Na⁺ |
| knock-on pair / conductive / landscape | +1.54 / +4.76 / +9.99 eV | Na⁺ |
| measured dehydration (Marcus 1991 TATB) | K⁺ enters 0.726 eV cheaper | K⁺ |

The legs disagree, and that disagreement is the result. **No quantitative
selectivity free energy is claimed.** Only the measured dehydration free
energy is a robust K⁺ advantage in this model; the two-leg sum is compared
with the experimental −0.179 eV reference as a scale check, never as a
prediction. The umbrella sampler prints its measured autocorrelation time τ
and effective sample size with every barrier.

#### Demo 12c — what the ionic layer actually shows

Demo 12c exists because the forensic table above is not the whole story, and
because two questions the tree previously could not ask at all now have
answers.

**Where the desolvation leg comes from.** Until now the hydration free energy
was a transcribed table with nothing able to check it or vary it. There is now
a Born self-energy to compare it against — and it disagrees with the table by
+18.7 kJ/mol, which is *supposed* to happen: the Born expression is a
linear-response continuum result with no room for the cavity and dispersion
terms the real ion–water interaction has, and for monovalent cations it
under-predicts the differential by about that much.

| | ΔG(Na⁺−K⁺), kJ/mol |
|---|---|
| Born, r_Born 1.80 / 2.08 Å | −51.3 |
| measured, Marcus 1991 TATB | −70.0 |
| Born's deficit | **+18.7** |

The Born radius is not unique, and that is reported rather than hidden:
refitting it to the measured hydration free energies shifts the differential
by a further −20.4 kJ/mol. That is the honest size of the uncertainty
attributable to one unsolvable choice.

**Whether the filter actually excludes Na⁺.** This is a geometric question
with a geometric answer, and it is the sharpest thing this engine knows about
KcsA. On the pore axis the mean ion–oxygen distance over the eight
coordinating oxygens is **2.777 Å**, and that is the *shortest the deposited
geometry permits* — an ion cannot get closer than the oxygens are. Measured
against each ion's own preferred contact (Shannon VIII radius + 1.40 Å
carbonyl oxygen):

| ion | prefers | achieves | mismatch |
|---|---|---|---|
| K⁺ | 2.910 Å | 2.777 Å | **−0.133 Å** (squeezed inside its preference) |
| Na⁺ | 2.580 Å | 2.777 Å | **+0.197 Å** (held outside it) |

The cage is built for K⁺, so the larger ion is squeezed 0.13 Å inside what it
wants while the smaller one is held 0.20 Å outside what it wants and cannot
get closer anywhere on the axis. **That asymmetry is the exclusion mechanism,
and it is a property of the deposited coordinates alone** — no force field, no
dielectric, no fitted parameter.

**Why that still does not produce a selectivity number here.** Both ions relax
to the *same* on-axis minimum at the same 2.777 Å, and their ion-induced shell
strains agree to **1.2%** — a ratio of 1.0122 that is exactly invariant across
a 25× stiffness scan, because strain goes as 1/k and the ratio does not. The
reason is quantitative: the ion–oxygen Lennard-Jones wall here is 0.05
kcal/mol, roughly four orders of magnitude below the Coulomb pull at 2.8 Å, so
nothing in this potential rewards a cation for having the right radius.
Turning on the exchange (overlap-Pauli) and dispersion terms moves the
difference by under 0.1%.

So the honest conclusion is a sharp negative one: **in this model the entire
K⁺/Na⁺ discrimination is the tabulated hydration free energy.** The geometric
mismatch is real and measurable, and it still does not convert into a
binding-energy difference.

**The entropic term, and why it is reported unresolved.** Dixit & Asthagiri
(*J. Chem. Theory Comput.* **5**, 1037, 2009; with Asthagiri, Roux & Cross,
*Biophys. J.* **95**, 1269, 2008) show that the *mean* binding energy slightly
**favours Na⁺** — Na⁺ is the better fit for an eight-carbonyl cage — and that
KcsA selectivity comes from the **fluctuation** term, `Δ(βσ²/2) ≈ 5–6
kcal/mol`, driven by carbonyl–carbonyl repulsion strain of about 7.3 kcal/mol.
That is the whole mechanism, and it is a term a rigid-cage model cannot
produce.

This release adds the machinery (`qm_fluctuation`): it runs the deposited
filter as a breathing peptide with the ion inside and measures the variance of
the ion's binding energy over an explicit Langevin NVT trajectory, with
streaming (Welford) moments, standard errors on every variance, a
stiffness-independent linear-elastic strain model, and a **site-integrity
gate** that fails closed if the filter does not hold its 8-fold site.

**It does not yet converge, and the reason was structural.** The engine had no
side-chain torsion terms, and four of the eight filter oxygens belong to
Thr75's side chain. Without torsions the hydroxyl leaves the shell, the
coordination number falls from 8 to 3–5, and the O–O diagonals close from the
crystallographic ~4.5–5.5 Å to ~3.5 Å — measured at every temperature from
50 K to 300 K and every Cα restraint stiffness from 0.5 to 16 eV/Å². No
thermostat and no restraint recovered it, because the missing term was the
physics and not the parameters.

This release adds the missing protein physics: `kcsa_add_dihedrals()` builds
backbone (phi/psi/omega) and side-chain (chi1) torsions restrained toward the
deposited geometry, and the fluctuation measurement now holds the ion on the
pore axis at the deposited site (the symmetric-cage equilibrium is unstable —
any radial displacement draws the cation toward the nearest oxygen because the
ion-O LJ wall of 0.05 kcal/mol is four orders of magnitude too weak to balance
the Coulomb 1/r pull). With torsions + axis constraint, the filter holds its
8-fold coordination and the binding-energy variance can be measured. Full
convergence of the fluctuation term (CN = 8 held for the full sampling window,
stable variance across temperatures and stiffnesses) requires parameter tuning
of the dihedral force constants, tether stiffness, and sampling length — that
tuning is the remaining step before the fluctuation term can produce a
quantitative selectivity number. The leg is reported **UNRESOLVED** until that
tuning is complete.

### Datastream
A text, self-describing, line-oriented format (`DATASTREAM_SPEC.md`, schema 1)
with a SHA-256 integrity seal over the payload, written by Demo 12 as
`kcsa.cvmds`. Claims carry provenance tags (`computed`, `computed-jc2008-params`,
`Marcus1991-TATB`, `1K4C-LINK`, …) and a validator can check external-citation
claims against reference values. The writer fails closed: truncated or
unsealed output is never written.

Two properties of the writer are worth stating because they shaped how
Demo 12c reached the file at all. `ds_open` buffers the whole payload and
truncates at `ds_close` — there is no append path and no reopen — and Demo 12
runs *before* Demos 12b/12c, so the writer cannot reach a value Demo 12c has
not computed yet. Rather than hoist Demo 12's hundred-odd locals to file
scope, Demo 12c's measurement is factored into `kcsa_ionic_measure`, a pure
deterministic function of the deposited coordinates called by both the
printer and the writer. The printed table and the sealed claims therefore
come from one code path and cannot drift apart in principle; 22 checks in
`test_external` pin the values so a loss of determinism breaks a gate rather
than producing two plausible, disagreeing numbers.

A claim that is *not* computable is sealed explicitly rather than omitted:
`kcsa.entropic.computed` is recorded as `0 kcal/mol` tagged
`UNRESOLVED-no-sidechain-torsions`, so the gap sits in the sealed file
beside its published target (`kcsa.entropic.target`, 5.5 kcal/mol from
Dixit & Asthagiri) rather than only in prose.

### Interactive terminal — `s2tui`
A separate binary so the record path can never leak wall time or keystrokes
into stdout. It is a **fully POSIX command surface** over one live
`Simulation` plus one HH neuron: every command name is a POSIX command and
the shell grammar is POSIX (pipes, redirection, `&&`/`||`, `$( )`, globs,
`fc -l`, `tput`, `sh -c '...'` to reach the host), while operands are S2
words — species, atom indices, scopes. Non-POSIX names are refused with
status 127.

The centrepiece is the **pure simulator mode**: `dd if=petri of=world`
builds a 48 Å periodic dish at 310 K holding ~780 atoms of interacting
matter (waters, Na/K/Cl ions, organic monomers) plus an HH neuron, and two
keybind surfaces edit that same world:

* `ps` — live control monitor (space run/pause, `s` step, `+/-` speed,
  `ijkl` rotate, arrows strafe across the environment, and
  create/augment/adapt keys: `w` water, `N` Na⁺, `K` K⁺,
  `C` Cl⁻, `u` base, `a` alanine, `g` glycine, `d` sugar, `H`/`L`
  heat/cool, `f` freeze, `r` replace, `p` clone, `x` delete, `m`
  minimise, `n` neuron, `e` catalysis). Piped, it prints one snapshot.
* `vi world` — modal editor (hjkl pan, `w`/`b` select, `i` insert,
  `r` replace, `x` delete, `p` clone, `u` undo, `/` search, `:` command
  mode running any shell command against the live world).
* `screen` — fullscreen twin of `ps` on the alternate screen (ncursesw
  build, real terminal required): same world, same keybinds, terminfo
  colors, resize handling, status footer. Both surfaces consume the
  same viewport cells, so the two pictures cannot disagree. Piped,
  dumb-`TERM`, missing-lib or failed init degrades to a typed line
  with status 0; `S2TUI_SCREEN=0` forces line mode. Esc quits like `q`.

The POSIX mapping for the rest of the world: `touch` creates matter,
`ln` links, `ln -s` restrains, `unlink` clears, `rm atom <i>` deletes,
`fsck` detects topology, `sleep <n>` advances world time, `kill -SIG`
freezes / stimulates / heats, `nice` minimises, `df`/`du` report state,
`cp <species> x.mol` exports a template and `dd if=x.mol` instantiates it,
`make x.rxn` fires a reaction rule, `sync` saves, `man` documents (30 topics incl. test/help/tput/fc/ps-keys/vi-keys; `ps`/`vi` live keybinds also in this section). Demos
remain live system tests (`test all` = 15/15). TRUST: `!`/`sh -c`/`source` run with YOUR
privileges (host exec + filesystem write) — never source untrusted files;
set `S2_NO_HOST=1` to disable the host escape. `save`/`load` (`S2SAVE1`) is LOSSY:
dihedrals, explicit angle overrides, box/PBC, step/time, QM flags and RNG state are NOT
restored (a loud WARNING prints on every load); `S2SAVE2` with full state is future work.
The orthographic ASCII grid
has depth shading, bonds, restraint anchors and a z-slab cutaway. Stream
contract: stdout is DATA, stderr is diagnostics.

---

## 5. The demos (20)

Run `./carbonsim` for the record, or `test demo <n>` inside `s2tui` for the
fast live smoke check.

12b is the real-filter twin of 12, 17 is the duplex, 18 is EHT, and 13–16
are the quantum-layer demos added in v9R4.

| # | Demo | Current headline |
|---|---|---|
| 1 | Quantum orbital structure | H/C/N/O Slater orbitals, exact H-like expectations, Clementi–Raimondi exponents, per-element NIST ledger |
| 2 | H–H bond curve | covalent and van der Waals regimes as two different physics |
| 3 | H₂O molecular dynamics | Berendsen NVT at 300 K, stable O–H geometry |
| 4 | Cyclic (H₂O)₃ | hydrogen-bonded ring emerges from Coulomb+LJ alone — HELD |
| 5 | CH₄ | tetrahedral 109.47°, zero-strain placement |
| 6 | Five nucleobases | PDB CCD geometry, RESP charges, planarity, tautomer-corrected cytosine |
| 7 | Watson–Crick pairing | G–C > A–U, both bound; magnitudes qualitative (see limitations) |
| 8 | T-p-A dinucleotide | real sugar–phosphate backbone, net charge −1.00 e |
| 9 | Hodgkin–Huxley neuron | action potential, +40 mV peak, 4 spikes / 50 ms |
| 10 | Gly–Ala dipeptide | genuine 1.33 Å peptide bond |
| 11 | α-helix | i,i+4 H-bond emerges: H···O 2.1204 Å, N···O 3.0982 Å |
| 12 | KcsA legacy cage | full selectivity leg program; vacuum sites favor Na⁺, exchange favors K⁺ |
| 12b | Real KcsA filter | deposited TVGYG, 164 atoms, C4 symmetry, CN = 8 at every site |
| 12c | KcsA ionic layer | Born self-energy vs the measured table; filter contact mismatch; entropic term's status |
| 13 | H₂ quantum MD | the Roothaan–Hartree–Fock layer drives an H₂ trajectory; the force falls 19 → 1.6 eV/Å as the molecule relaxes to the SCF minimum |
| 14 | H₂O quantum MD | **UNRESOLVED** — the SCF is not bound for water in this minimal basis, so the demo says so and runs the classical reference instead of printing a quantum-MD table it cannot support |
| 15 | H₂ dissociation curve | the SCF curve **has** an interior minimum (R = 1.20 Å against the experimental 0.7414 Å) after the missing nuclear–nuclear repulsion was added; the deviation is printed, not hidden |
| 16 | H₂ geometry optimization | the SCF-force optimiser converges to the SCF minimum and the deviation from the experimental bond length is printed |
| 17 | DNA duplex | G–C / A–T stack with B-DNA rise and twist; G–C BREATHING / A–T HELD |
| 18 | EHT semiempirical SCF | H2 band dissociation, H2O Mulliken + gap, BO switch dissolves the wall |

---

## 6. Provenance and verification

Two standards, both structural rather than aspirational: parameters are
checked against primary sources, and claims are checked by independent
oracles. `make audit` runs the evidence.

| Subsystem | Source | How it was checked |
|---|---|---|
| Fundamental constants | 2019 CODATA | derived in-line; reciprocals exact by construction |
| Periodic-table LJ | UFF / Rappé et al. 1992 | all 36 elements, ε and σ, cell-by-cell; σ = x1/2^(1/6) |
| Biomolecular LJ | AMBER ff99 (`parm99.dat`) | atom-by-atom mapping; published σ reproduced (N 3.24979 vs 3.2500, …) |
| Nucleobase charges | Aduri et al. 2007 RESP | neutrality and atom-name mapping cross-checks |
| Biomolecular geometry | RCSB PDB CCD ideal coordinates | re-measured from placed Cartesians, not read back |
| KcsA structure | PDB 1K4C chain C | deposited coordinates; C4 symmetry and CN = 8 verified |
| Ion hydration ΔG | Marcus 1991, TATB absolute scale | primary table (−295.3/−365.3 kJ/mol K⁺/Na⁺) |
| Born radii | Latimer, Pitzer & Slansky 1936 | corroborated for Na⁺ by Hummer/Pratt/García arXiv:physics/9709042 (1.80 Å) |
| Free-ion α | Pauling compilation | cross-tabulated; K⁺/Na⁺ ratio 4.6 gated |
| KcsA entropic term | Dixit & Asthagiri, JCTC 5, 1037 (2009) | mechanism, sign and ~5–6 kcal/mol scale recorded as the target |
| Ion LJ sizes | Joung–Cheatham 2008 | parameters + derived conversion |
| Water | TIP3P | geometry and σ against literature |
| Ionic radii | Shannon 1976 (VIII) | contact distances against deposited geometry |
| Neuron | Hodgkin–Huxley 1952 | rate functions and a published worked example, 4 decimals |
| SHA-256 | FIPS 180-4 | NIST known-answer vectors |
| RNG | PCG64 (O'Neill 2014) | the REFERENCE pcg64 (XSL-RR 128/64), checked against the upstream `pcg-c` expected-output file — see §6.2 |

The second-audit pass verified the whole tree against these sources, including
things that were wrong and are now fixed: the KcsA dehydration leg had been
carrying hydration *enthalpies* under a free-energy label (corrected to the
Marcus 1991 TATB free energies, K⁺ advantage 1.368 → 0.726 eV), a glycine
charge table carried a sign-inverted neutrality correction, and several
constants had duplicated spellings free to drift. Each correction is pinned
by a regression check that fails on the old code. The row-by-row history
lives in `release_note_v9R4.md` and the commit log, not here.

### 6.1 Third audit pass (2026-10-10) — independent oracle sweep

A full mathematical, programming and operational re-audit ran an
**independent verification harness** built outside the tree (`/tmp`, linked
against the shipped engine objects, 78 checks) rather than reusing the
project's own oracles, so a shared mistaken assumption could not pass twice.
Its findings:

* **Constants — clean.** Every SI literal in `include/constants.h` matches
  CODATA 2019 exactly (h, c, e, k, N_A by redefinition; mₑ, m_p, m_n, ε₀,
  a₀, E_h, amu to the published uncertainty), and every derived conversion
  (`COULOMB_MD` = 14.3996454784 eV·Å, `MD_FORCE_CONV` = 9.648533212e-3,
  `KCAL_MOL_TO_EV` = 0.04336410414, `AMU_AFS2_TO_EV` = 103.642696) was
  re-derived in Python from its primaries and agrees.
* **Periodic table — clean, 36/36.** All UFF σ and ε were re-verified
  cell-by-cell against an independent machine-readable transcription of
  Rappé Table II (openbabel `UFF.prm`, cross-checked against cp2k and
  PorousMaterials.jl). Two values are commonly mis-transcribed and the tree
  has both **correct**: Si σ = 3.82641 Å (x₁ = 4.295, not the widely
  circulated 4.195) and Ar σ = 3.44600 Å (x₁ = 3.868, not 3.912). Every σ is
  UFF x₁ / 2^(1/6), never the x₁ itself.
* **AMBER ff99 — clean, 13/13.** All LJ types in `include/amber_lj.h` match
  `audit/external/parm99.dat`'s `MOD4 RE` section exactly (H 0.6000/0.0157 …
  P 2.1000/0.2000), and `sigma = 2·R*/2^(1/6)` reproduces AMBER's published
  sigmas. (One documentation typo fixed: the N amide σ is 3.25000 Å, not
  the "3.24979" the header comment claimed.)
* **Analytic forces — clean, to machine precision.** Every force term was
  checked against a central finite difference of its own energy at
  h = 10⁻⁶ Å on systems the shipped suites do not cover: CHARMM cutoff
  switching, flat-bottom restraints, PBC minimum-image pairs, the coupled
  3N×3N dipole solve, Pauli and dispersion. Worst deviation
  **5.2 × 10⁻¹⁰ eV/Å**.
* **Integrator — symplectic.** A harmonic oscillator conserves energy to
  0 eV over 20 000 steps; the water-dimer energy error is bounded, not
  secular (late-window max |ΔE| / early-window = 0.98); Maxwell–Boltzmann,
  Andersen and Langevin all reproduce the correct stationary ⟨v²⟩ = 3kT/m
  and Berendsen relaxes 1000 K → 300 K.
* **Hodgkin–Huxley — clean.** Resting gating values reproduce the published
  1952 worked example to four decimals (m 0.052932, h 0.596121, n 0.317677);
  both removable singularities take their L'Hôpital limits exactly; the
  action potential peaks at +40.3 mV.
* **KcsA — clean.** Marcus 1991 TATB hydration free energies, Shannon 1976
  VIII radii, the kJ/mol-per-eV conversion, and CN = 8 at all four
  deposited 1K4C K⁺ sites all verify.
* **SHA-256 / datastream — one real gap, fixed.** The FIPS 180-4 vectors pass,
  and a tampered payload byte is rejected — but content appended *after* the
  seal verified clean (Q3, now rejected and pinned).
* **QEq — one real physics defect, now refused.** See §7.

The four defects fixed in this pass are listed in `release_note_v9R4.md`
under the third-audit section, and each is pinned by a regression check that
fails on the pre-fix code.

### 6.2 Fourth audit pass (2026-10-21) — full mathematical, programming and operational re-audit

A second full audit, run over the whole tree at once: four independent
reviewers worked in parallel over the quantum layer, the EHT/quantum
modules, the record program and the KcsA/chemistry code, and the
operational pipeline (build, gates, verifier, promotion, CI). Every
finding was then **re-derived from a primary source before it was acted
on** — and that discipline changed the outcome in both directions:

* The claim that the EHT SCC gamma was "not the Ohno interpolation" was
  **correct**, and the label was wrong — but replacing the kernel with the
  textbook Ohno form would have made the model *worse* (+19 % against a
  grid oracle, where the existing kernel is within 1.0 %). The kernel was
  kept, the label was corrected, and the accuracy was pinned.
* The claim that the reference PCG64 was "not reference pcg64" was
  **correct** — and the fix was to implement the real thing and check it
  against the upstream `pcg-c` expected-output file rather than against
  itself.

**What the audit found and what was done about it.** The defects are
listed in full in `release_note_v9R4.md` under the fourth-audit section.
The headline items:

| Area | Defect | Consequence | Status |
|---|---|---|---|
| EHT overlap | sign error in the s–pσ branch | the overlap matrix was not symmetric under the molecular mirror; the SCC did not converge for H–C–H; N₂ had a 0.0000 eV gap and a 2.05 Mayer bond order | fixed; one shared kernel now serves both SCF modules |
| RHF overlap | per-atom valence shells instead of per-AO | the 2s–1s channel of every O–H / C–H / N–H pair was exactly zero | fixed |
| RHF core | `ao_n[]` never populated | the VSIP diagonal silently fell back to the first ionisation potential (O 2s: −13.6 eV instead of −32.3 eV) | fixed |
| RHF energy | the density was passed to `qm_eri_zzdo()` with the wrong stride | the entire two-electron Fock term was built on garbage populations, and the Fock matrix came out asymmetric for a symmetric molecule | fixed |
| RHF energy | the eigenvalue sort overwrote its own source | the LUMO was replaced by a duplicate of the HOMO | fixed |
| RHF energy | the nuclear–nuclear repulsion was missing | the H₂ curve had no minimum at all | fixed; the minimum is now at R = 1.20 Å |
| RHF forces | a failed displaced solve left the energy at 0.0 (a *finite* number) | a force of order E/2h ≈ 1.5 × 10⁶ eV/Å, silently | fixed; fail-closed now |
| Induction | the coupled-dipole block carried a spurious 14.4× Coulomb factor | dipoles 3.6× too large, able to anti-align with the field, so the induced-dipole energy could come out **positive** | fixed |
| Born | the induction energy was missing the 1/(4πε₀) factor, and the "cage field" was the *potential* (5.18 V) rather than the *field* (1.863 V/Å) | the printed induction energy was **111× too large** | fixed; the value is now −0.006253 eV (K⁺) |
| Coulomb kernel | screened with the mean ionisation potential | 10.09 eV for the water O–H pair against a grid oracle of 5.7156 eV (76 % error) | fixed; the kernel now agrees to 1.0 % |
| Fluctuation | the paired energy difference was taken at two different geometries | 5–18 % of the reported variance was an artefact, and larger than the K⁺-vs-Na⁺ difference the argument rested on | fixed |
| KcsA | the hydroxyl cone was built around the C→O axis | ∠C–O–H came out at **70.5°** instead of 109.5° on the innermost filter ligand | fixed |
| Forces | the cell-list path allocated its work buffers twice | every `forces_calculate()` call in the perf path leaked ~432 kB | fixed |
| Forces | the cell-list path ignored the degenerate-cutoff rule and the angle-list exclusion | the two pair paths could disagree about the same system | fixed |
| Sim | the cell-list grid was built with unchecked `int` products | a user-settable small cutoff could overflow into a **heap out-of-bounds write** | fixed; the count is computed in `double`/`size_t` and refused |
| Sim | the grid origin was discarded | `sim_cell_index_atom()` returned −1 for every atom of any system not near the origin | fixed |
| Integrator | the temperature cap wrote to stderr in batch mode | 69 lines of stderr, so the empty-stderr gate — and `./run`'s promotion — could never succeed | fixed; the cap is counted and reported as data |
| RNG | "PCG64" was an invention with a 64-bit increment and a made-up output permutation | no published vector could check it, and the readme contradicted itself | fixed; the reference pcg64 is implemented and pinned against the upstream expected-output file |
| Datastream | a NaN or inf could be sealed into a signed file | a non-finite claim in a file that verified "clean" | fixed; the writer fails closed |
| Datastream | `ds_verify_file(NULL)` called `fopen(NULL)` | undefined behaviour (crashes on musl and macOS) | fixed |
| Pipeline | `verify_scripts.sh` tested `fail` **before** the TUI checks | a red TUI check still printed "VERIFY PASSED" and exited 0 | fixed |
| Pipeline | `./run --accept` copied an **unsealed** `kcsa.cvmds` after warning about it | the tracked, sealed artefact could be replaced by a tampered one | fixed; the promotion refuses |
| Pipeline | `./run --accept` used `cp` in place | an interruption mid-copy left a truncated record | fixed; every artefact is written to a temporary name and renamed |
| Pipeline | `make perf` shared the record's objects and binary names | `make perf` was a silent no-op, or broke the record link, or overwrote the tracked binary | fixed; the perf build is fully isolated |
| Pipeline | three of the eight gates had no `FAIL`-line fallback | a suite could print FAIL and the gate would still pass | fixed |
| Pipeline | there was no CI at all | the stale record reached HEAD with every commit message claiming a green tree | fixed; `.github/workflows/ci.yml` now gates the record |
| Audit | `audit/thirdpass_uff_amber.py` was never executed and compared the engine against the LAMMPS table | 20 red findings that were all artefacts of the checker | replaced by `audit/math/a06_uff_amber.py`, now part of `make audit` |
| Docs | the readme claimed a reference PCG64 in one place and denied it in another | — | fixed |

**New evidence, not just new checks.** The audit added external oracles
that did not exist before:

* the upstream `pcg-c` expected-output file (published PCG64 vectors);
* a brute-force 3D Cartesian grid integration of the *actual* valence
  charge densities, which is what the Born-induction and screened-Coulomb
  numbers are now checked against;
* the textbook closed forms for the STO overlaps and the L'Hôpital limits
  of the Hodgkin–Huxley rate functions;
* the mirror symmetry of a bent H–C–H, used as an oracle for the overlap
  matrix;
* the Open Babel UFF parameter file (which cites Rappé et al., *JACS* **114**
  (1992) 10024) as the reference for all 36 Lennard-Jones rows.

All of these are in `tests/test_audit4.c` (the eighth gate) and
`audit/math/a06_uff_amber.py` (now part of `make audit`).

---

## 7. Known limitations

Stated plainly, because the numbers above cannot be read correctly without them.

* **No bulk solvent.** Everything runs in vacuum or with a relative-permittivity
  divisor; the ion program adds explicit octahedral 6-water clusters but no
  periodic water box, Ewald/PME, or membrane potential.
* **No quantitative KcsA selectivity free energy.** The real filter is rigid,
  the legs disagree, and the umbrella `±` is a three-seed spread, not a
  standard error. Read every K⁺/Na⁺ number as "in this reduced model".
* **The KcsA entropic term is implemented but does not yet converge.** The
  machinery to measure `Δ(βσ²/2)` exists and is gated, and side-chain + backbone
  dihedral torsions (`kcsa_add_dihedrals`) plus an ion axis constraint have been
  added to hold the filter's 8-fold site under Langevin dynamics. Full
  convergence — CN = 8 held for the full sampling window with a stable variance
  across temperatures and stiffnesses — requires parameter tuning of the
  dihedral force constants, tether stiffness, and sampling length. Demo 12c
  reports this as UNRESOLVED with the evidence. It is *not* filled with a
  plausible number.
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
* **`qm_overlap` is a heuristic, not an overlap integral.** It drops the
  Slater polynomial prefactor and maximizes `m` per atom; `qm_bond_order` and
  `qm_pauli` inherit that. The EHT layer does NOT use it: `qm_eht_overlap.c`
  integrates exact STO overlaps (prolate machinery + documented p-pi
  quadrature), validated against brute-force 3D grids by `test_eht`.
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
  The RNG labelled PCG64 is S2-PCG64-like (O'Neill-inspired, NOT reference
  pcg-random.org pcg64 — 64-bit increment, non-reference output permutation);
  the LCG remains the record default. QEq is QEq-lite (topological 1-2/1-3
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
  The H₂ curve therefore has its minimum at 1.20 Å against the
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

## 8. Interpreting the output

* **Sign conventions differ per block and are printed where they matter.**
  The KcsA forensic table uses `dU = E(K⁺) − E(Na⁺)`, so `+ = Na⁺ favored`;
  exchange lines use `negative = K⁺ selective`.
* **`±` values are spreads, not errors.** The WHAM `±` is the sample standard
  deviation across 3 independent seed repeats — with n = 3 its own uncertainty
  is ~52%, and the output says so. Within-run sampling error is characterized
  separately by τ and N_eff, also printed. Rule: when |value| < spread the
  verdict is UNDECIDED (the record now prints this); 4-decimal barriers with
  N_eff ~15 are not determined to 1e-4 — read 2 decimals.
* **How to read τ/N_eff.** τ = integrated autocorrelation time (samples, Geyer
  IPS, stopped before the first non-positive autocovariance); N_eff ≈ retained
  samples / τ per window. The record prints mean-τ and N_eff per repeat plus
  retained/dropped bin counts. With N_eff ~15–30 per window, barriers move ~0.1 eV
  across seeds (observed: polar-WHAM gap −0.04 → +0.12 across the HO-eps fix,
  both within spread) — which is why no selectivity free energy is claimed.
* **Loop closure (one paragraph).** The three tracks never share dynamics, but
  they share primaries: `src/loop.c` recomputes RT/F, kT, Nernst reversals
  (squid 440/50 Na⁺, 20/400 K⁺ at 279.45 K → +52.3/−72.1 mV), the 1000:1
  selectivity scale (−kT ln 1000 ≈ −0.179 eV at 300 K), LJ minima
  (2^(1/6)·σ) and Coulomb from the same CODATA constants, and enumerates the
  11 irreducible reduced-model choices (`LOOP_ABSTRACTIONS`) so no new
  assumption slips in unlisted — adding one breaks `test_loop` by design.
* **Verdicts** ("HELD", "BREATHING", "EMERGED") are the model's own
  thresholds, calibrated to the model's baseline behavior — not experimental
  pass/fail.
* **The record is stdout and only stdout.** stderr is display; if you need to
  quote a number, quote it from `output.txt` and check the digest.
* **One distance convention everywhere (Q2).** Every pair term — classical
  LJ/Coulomb, QEq, induced-dipole field and tensor, Pauli, dispersion —
  now measures a pair through the same `qm_sep`/minimum-image helper, so
  the two halves of a single potential energy can never disagree about a
  pair's distance. `qm.c` previously contained zero calls to
  `vec3_pbc_box` while `forces.c` contained three, which made every
  quantum-enhancement term silently non-periodic in a periodic box. Every
  demo is vacuum, so no recorded number changed.
* **A sealed datastream is sealed to the end of the file (Q3).** `ds_verify_file`
  now requires the `payload-sha256:` line to be the last content: exactly
  one optional trailing newline, then EOF. Content appended after the seal
  is covered by no digest, and previously verified clean. A claim no digest
  covers is not a claim this format can make.

---

## 9. Repository map

```
carbonsim, s2tui        built binaries (carbonsim tracked, s2tui untracked)
output.txt              the record (stdout of ./carbonsim)
CURRENT_BASELINE_SHA.txt  its SHA-256
kcsa.cvmds              sealed datastream written by Demo 12 (promoted atomically with output.txt)
runs/                   staged candidates (ignored; ring buffer, never committed — 6 unpromoted
                        diffs beside 2 matching pairs is normal pre-promotion, not failure)
run, verify_scripts.sh  record workflow and read-only verifier
                        (CI chain: make test && ./verify_scripts.sh && ./s01_verify_record.sh && make audit)
s01_verify_record.sh    ASan/UBSan record reproducibility harness
makefile                one command builds every executable (-Werror, -ffp-contract=off;
                        `make run` execs ./run — there is exactly one record path)
NOTE: `build/*.o` + `carbonsim` are TRACKED per repo custom, so every verifier
(`make clean` first) leaves `git status` showing rebuilt objects. Review only
`git diff output.txt CURRENT_BASELINE_SHA.txt kcsa.cvmds` for record moves.
The tracked objects embed the absolute build path in their DWARF, so they
are not byte-reproducible across machines even though the RECORD is.
readme.md               this file
release_note_v9R4.md    the release's corrections and rationale (history)
DATASTREAM_SPEC.md      the .cvmds format, schema 1
src/, include/, tests/  the engine and its suites
audit/                  independent oracles, C harnesses, external references,
                        pre-fix engine copies (make audit / make audit-revert)
.github/workflows/ci.yml   the CI chain: build, gates, verifier, audit, ASan
                         record reproducibility, the OpenMP build, and the
                         record-staleness check
472-MGB/                S2->MPE display bridge. NOT part of this repo's history:
                        it is versioned independently at
                        github.com/Shicheng-Zhang-0003/MPE-Magi-Bridge
                        and is gitignored here (clone it alongside). Developed
                        in this tree, copied into the 475-MPE tree as
                        v15S/src/ecosystem/mgb/, which is the live integration.
                        Display-only: it spawns s2tui, configures it through the
                        shell S2 already has, and reads S2SAVE1 frames from `sync`.
                        Positions flow S2 -> MPE only; never back.
```

### Consuming the record from another engine (472-MGB)

`472-MGB` drives the interactive `s2tui` shell over pipes — `stdout` is data,
`stderr` carries the prompt — and reads a strict `S2SAVE1` frame via `sync`
(see `DATASTREAM_SPEC.md` for the on-disk record format, which is a different,
sealed artifact). It never touches `output.txt` or the batch `./carbonsim`
path: one record path, one writer.

```bash
make -C 472-MGB test     # 33/33 parser + mapper + spawn-naming gates
make -C 472-MGB probe    # live petri dish through ../s2tui
```

The petri dish currently stages as a **single static snapshot**: bodies are
static and bonds carry `k = 0`, so it renders but does not move or vibrate.
That is Phase 0 by design; follow mode is open work tracked in
`475-MPE/v15S/REMAINING_WORK.md`.
