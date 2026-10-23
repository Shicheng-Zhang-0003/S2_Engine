# S2 Engine — Nano Chemistry & Biological Simulator

**v9R5** (internal) · **V0.9RC5** (external) · C11 · record path: no dependencies beyond libm
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

## Documentation map

The long-form documentation is split by subject. This page is the front
door: what the engine is, how to build and run it, how it is put
together, and what may not be cited.

| Document | What is in it |
|---|---|
| [`readme.md`](readme.md) | this page: build, run, architecture, the record in brief, the headline limitations, the repository map |
| [`docs/physics.md`](docs/physics.md) | every physics layer: the orbital layer, the exact STO integrals and the Roothaan–Hartree–Fock solver, the quantum–MD wiring, classical MD, the quantum-enhancement layer, the EHT layer, biopolymer condensation, electrophysiology, and the datastream |
| [`docs/kcsA.md`](docs/kcsA.md) | the potassium-channel program in full: the deposited 1K4C filter, the ionic layer (Born, induction, contact mismatch), the shell strain, and the entropic term — including the measured fact that the flexible site does not survive |
| [`docs/demos.md`](docs/demos.md) | the 20 demos and what each one actually claims |
| [`docs/verification.md`](docs/verification.md) | provenance: where every parameter comes from, the eight gates, and the third and fourth audit passes |
| [`docs/operations.md`](docs/operations.md) | the record and its artifacts, how to interpret the output, the promotion pipeline, and CI |
| [`docs/limitations.md`](docs/limitations.md) | the full list of known limitations, stated plainly |
| [`docs/terminal.md`](docs/terminal.md) | the interactive terminal `s2tui`: a POSIX shell over one live simulation |
| [`DATASTREAM_SPEC.md`](DATASTREAM_SPEC.md) | the `.cvmds` format: schema 1, the SHA-256 seal, the claim and provenance vocabulary |
| [`release_note_v9R4.md`](release_note_v9R4.md) | the release's correction history: every defect that was found, what it was, and what replaced it |
| [`audit/README.md`](audit/README.md) | the audit's own working material: independent oracles, C harnesses, external references, pre-fix engine copies |

## 1. Build, run, test

Requirements: a C11 compiler (gcc or clang), `make`, libm. Nothing else
for the record. The fullscreen `screen` surface optionally uses ncursesw
(detected via `ncursesw6-config`; without it `s2tui` still builds and
every line command works, only the alternate-screen reports unavailable).

| Command | Effect |
|---|---|
| `make` / `make -j$(nproc)` | build everything: `carbonsim`, `s2tui`, and all eight test binaries (serial record build, byte-deterministic) |
| `make perf` | build everything with **OpenMP multicore** (`-fopenmp -DS2_OPENMP -march=native`) — cell-list O(N) pair search + parallel O(N²) kernels across all cores. The perf build has its **own object directory (`build-perf/`) and its own binaries (`carbonsim-perf`, `s2tui-perf`)**, so it can never contaminate the record build. Performance runs only; **never for the record** (thread scheduling changes FP summation order, so digests differ from the serial record). |
| `./carbonsim` | run the 20-demo record sequence (**~11 minutes**; the polar-WHAM leg dominates); stdout is the byte-deterministic record |
| `./run` | clean build behind a warning gate, run to `runs/<timestamp>.txt`, print old/new digests side by side |
| `./run --accept` | promote that run to `output.txt` + `CURRENT_BASELINE_SHA.txt` + `kcsa.cvmds` (staged atomically: each artefact is written to a temporary name and renamed into place, and the ledger is written last); refuses to promote unless stderr is empty, `make test` is green, **and** the staged `kcsa.cvmds` seal verifies |
| `make test` | build and run the eight unit gates only (19/22/7/202/153/60/31/80) — CI must run `make test && ./verify_scripts.sh && ./s01_verify_record.sh` (plus `make audit` for evidence) |
| `make tui` | build just the interactive terminal (`./s2tui`) |
| `./verify_scripts.sh` | read-only record + spec verifier: baseline digest, `kcsa.cvmds` seal, key/unit compliance, every suite |
| `./s01_verify_record.sh` | builds normally and under ASan/UBSan, asserts both reproduce the recorded digest with empty stderr |
| `make audit` / `make audit-revert` | run the audit's independent oracles and C harnesses; `--revert` links the same regression suite against the pre-fix engine |
| `make clean` | remove `build/`, `build-perf/`, all four binaries, and `audit/build` |
| `s2tui` | interactive terminal; `help` prints the command map |

The record build uses `-O3 -g -Wall -Wextra -Werror -ffp-contract=off -std=c11 -Iinclude` and never
`-march=native`: `-ffp-contract=off` pins FP codegen (the `-march` ban alone does not stop FMA
fusion under `-O3` on FMA hosts), and `-Werror` makes warnings fail the build.
Set `SOURCE_DATE_EPOCH` for byte-stable side artifacts (`./run` defaults it to the tree's own
commit time so `kcsa.cvmds` seals identically across repeated runs).

**Current gate status: all green.** The eight gates, what each one pins,
and the 574 checks they total are in §5 below. On top of them:
a warning-free clean build (serial and `make perf`), an empty stderr on the
record, byte-identical stdout across repeated runs (serial record only),
an intact `kcsa.cvmds` seal, and the OpenMP build passing all eight gates
with numerically consistent results (within FP reassociation tolerance).

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

## 2. Architecture

Unit system: **Å, fs, eV, AMU, elementary charge.** All fundamental constants
are 2019 CODATA values, and every conversion factor is derived in-line from
primaries (`include/constants.h`) rather than hand-typed, so reciprocals are
exact by construction.

Current tree: **34 818 lines** (src 26 337 · include 3 864 · tests 4 617).

| File | Lines | Role |
|---|---:|---|
| `src/tui.c` | 7430 | interactive terminal: full POSIX shell + `screen` fullscreen twin, chemistry lab, gas/barostat (not in the record) |
| `src/main.c` | 5233 | the 20-demo record program and datastream consumer (Demo 7 A-T + monomer-subtracted legs; fail-closed exit status) |
| `src/qm.c` | 3291 | QEq-lite (topological hard core, not literal Rappe shielding) + induced dipoles (first-order + coupled SCF), Pauli, TT/BJ dispersion, penetration-damped Coulomb, overlap, SCF driver, **exact STO overlaps (ss/sps/pps/ppi), two-electron integrals (Hubbard-screened Coulomb + Wolfsberg-Helmholtz exchange), Roothaan–Hartree–Fock SCF solver + FD forces** |
| `src/nucleobases.c` | 1279 | five bases, deoxyribose, T-p-A dinucleotide, pairing/geometry helpers |
| `src/forces.c` | 1365 | non-bonded pairs (upper-bound OOB-checked), bonded terms, harmonic/flat-bottom restraints, analytic dihedral gradients, energy breakdown; QM caps fail loud; **OpenMP cell-list pair search + private-force accumulators** |
| `src/sim.c` | 1089 | lifecycle, molecule constructors, topology rebuild (P6-guarded both twins), ion/restraint plumbing (32 restraints), **cell-list build for O(N) non-bonded pairs** |
| `src/integrator.c` | 924 | Velocity Verlet, **reference PCG64** (O'Neill 2014, XSL-RR 128/64) / Knuth-MMIX LCG, Maxwell–Boltzmann, Berendsen/Andersen/exact-OU-Langevin, steepest descent, FIRE, **standard C11 64x64->128 multiply** |
| `src/quantum.c` | 547 | Slater screening and orbital energies, hydrogenic radial profiles, Clementi–Raimondi exponents, real spherical harmonics |
| `src/kcsa_filter.c` | 980 | real 1K4C TVGYG filter, ion sizing, sites, binding/dehydration legs, C-alpha tethers, H-completion (26/subunit), **backbone + side-chain dihedral torsions** |
| `src/qm_born.c` | 128 | Born ionic self-energy, Latimer–Pitzer–Slansky radii, Pauling free-ion polarizabilities |
| `src/qm_fluctuation.c` | 815 | variance/entropy machinery, KcsA shell strain, site-integrity gate, **ion axis constraint for fluctuation sampling** |
| `src/aminoacids.c` | 513 | glycine/alanine/dipeptide/polyalanine builders |
| `src/datastream.c` | 445 | schema-1 writer, self-contained FIPS 180-4 SHA-256, seal verifier (fseeko/ftello, record-tree source-hash, \r/bracket hygiene) |
| `src/periodic_table.c` | 323 | H–Kr element data, UFF ε/σ, Madelung electron configurations |
| `src/tui_view.c` | 383 | viewport cell model shared by the ANSI frame and the curses screen |
| `src/tui_screen.c` | 171 | ncursesw alternate-screen monitor (optional dep, TTY-gated, piped-safe) |
| `src/neuron.c` | 173 | Hodgkin–Huxley 1952 (squid giant axon) |
| `src/loop.c` | 84 | bio/QC/QM loop closure: Nernst, kT scale, LJ/Coulomb from same primaries |
| `tests/test_regression.c` | 1831 | the 202-check in-repo oracle suite |
| `tests/test_external.c` | 1095 | the 150-check outside-reference suite |
| `tests/test_loop.c` | 273 | the 60-check loop-closure + best-layer suite |
| `src/qm_eht_overlap.c` | 467 | exact STO overlaps: prolate A/B machinery, GL quadrature for p-pi, **and `eht_overlap_ao()` — the single shared global-axis AO overlap kernel** |
| `src/qm_eht_scf.c` | 697 | Hoffmann EHT + SCC (Jacobi, Lowdin, Mulliken, Mayer, trust mixing) |
| `tests/test_eht.c` | 331 | the 31-check EHT oracle suite |
| `tests/test_audit4.c` | 696 | the 76-check fourth-audit suite (one check per defect found and fixed) |
| `include/qm_eht.h` | 128 | EHT API + E1-E11 abstraction inventory |
| `include/loop.h` | 55 | loop bridge contract + abstraction inventory |
| `include/` | 3864 | types, constants, per-module contracts; `amber_lj.h` and `display.h` carry shared policy |
| `tests/` (other) | 378 | datastream, forces, fire suites |

### 2.1 Multicore and cell lists (performance build)

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

## 3. What the engine models

Three tracks exist — nucleic acids, proteins, electrophysiology — each
validated independently; they are not yet connected to each other. Every
layer below is documented in full in
[`docs/physics.md`](docs/physics.md), with its honesty block.

| Layer | What it is | Documented in |
|---|---|---|
| Quantum orbital layer | Slater screening, Slater `n*` orbital energies, hydrogenic radial profiles with exact expectation values, Clementi–Raimondi exponents, all 16 real spherical harmonics (l = 0…3) | [physics](docs/physics.md#quantum-orbital-layer) |
| Exact STO integrals | two-centre Slater-type overlap integrals in prolate-spheroidal coordinates, grid-validated; the **shared** global-axis kernel `eht_overlap_ao()` both SCF modules use | [physics](docs/physics.md#exact-two-center-sto-integrals-and-roothaan-hartree-fock-scf) |
| Roothaan–Hartree–Fock SCF | Löwdin orthogonalization + Jacobi eigensolver, Aufbau filling, damped density mixing, finite-difference forces. **Qualitative layer, not a predictive one** — honesty block H1–H5 in `include/qm.h` | [physics](docs/physics.md#exact-two-center-sto-integrals-and-roothaan-hartree-fock-scf) |
| Semiempirical SCF (EHT) | Hoffmann 1963 with self-consistent charge: exact overlaps, VSIP Hamiltonian, Wolfsberg–Helmholtz coupling, Mulliken charges, Mayer bond orders, HOMO–LUMO gaps | [physics](docs/physics.md#semiempirical-scf-layer-eht) |
| Classical MD | Velocity Verlet, Berendsen / Andersen / exact-OU Langevin, harmonic bonds/angles/analytic dihedrals, flat-bottom restraints, steepest descent and FIRE, the **reference** PCG64 | [physics](docs/physics.md#classical-molecular-dynamics) |
| Quantum enhancement | induced dipoles (first-order and coupled SCF with Thole damping), overlap Pauli, Slater–Kirkwood dispersion with Tang–Toennies damping, QEq-lite charge equilibration — every term opt-in and separately reported | [physics](docs/physics.md#quantum-enhancement-layer) |
| Biopolymer chemistry | leaving groups removed per the PDB Chemical Component Dictionary, CCD ideal coordinates, Aduri 2007 RESP charges, AMBER ff99 LJ types, a real T-p-A phosphodiester, a Gly–Ala peptide bond, an α-helix whose i,i+4 H-bond emerges from Coulomb+LJ alone | [physics](docs/physics.md#biopolymer-condensation-chemistry) |
| Electrophysiology | the 4-variable Hodgkin–Huxley 1952 squid giant axon, RK4-integrated, removable singularities by Taylor expansion | [physics](docs/physics.md#electrophysiology) |
| KcsA ion channel | the deposited 1K4C chain-C TVGYG filter, the ionic layer, the contact-mismatch result, and the entropic term's honest status | [kcsA](docs/kcsA.md) |

## 4. The demos

The record is a 20-demo sequence. [`docs/demos.md`](docs/demos.md) has
the table of what each one claims; every headline in it is printed by
the demo itself, and the record — not the table — is the authority.

A note on the two that report negative results, because they are the
interesting ones:

* **Demo 14 (H₂O quantum MD) is UNRESOLVED.** The Roothaan–Hartree–Fock
  layer does not produce a *bound* solution for water in this minimal
  valence basis, so the solver refuses (status −4) and the demo says so
  and runs the classical reference instead of printing a quantum-MD table
  it cannot support.
* **Demo 15/16 (H₂ dissociation and geometry optimisation)** work, and
  the H₂ minimum sits at **R = 1.30 Å** against the experimental
  0.7414 Å. The deviation is printed, not hidden.

## 5. Verification at a glance

**574 checks in eight gates**, plus 20 independent audit harnesses. Every
gate is built on an oracle *outside* the thing it checks.

| Gate | Checks | What it pins |
|---|---:|---|
| `test_datastream` | 19 | FIPS 180-4 SHA-256 known-answer vectors, round-trip, tamper rejection, seal-is-final |
| `test_forces` | 22 | analytic dihedrals vs a finite-difference oracle, net-zero forces, collinear guard |
| `test_fire` | 7 | FIRE and steepest descent reach the same minima |
| `test_regression` | 202 | in-repo oracles: quadrature, finite differences, hand-derived values, QEq ordering guard, QM minimum image, zero-strain bonds, KcsA shell strain |
| `test_external` | 153 | values from outside this repository: NIST CODATA, `parm99.dat`, FIPS, PDB 1K4C, Latimer–Pitzer–Slansky Born radii, Pauling polarizabilities, HO-eps unit pin |
| `test_loop` | 60 | bio/QC/QM loop closure: Nernst, kT, Marcus/Shannon/LJ/Coulomb, BJ/penetration, Langevin, flat-bottom, tethers, H-complete ions-hold |
| `test_eht` | 31 | exact overlaps vs 3D-grid STO oracles, Jacobi vs 2×2 forms, H₂/H₂O SCF laws, fail-loud gates |
| `test_audit4` | 80 | the fourth audit: the reference PCG64 vs the upstream expected-output file, EHT mirror symmetry, textbook overlap closed forms, Born induction units, one-electron integral limits, SCF invariants, the coupled-dipole operator, the cell-list origin and overflow refusal, the datastream's non-finite refusal, the Hodgkin–Huxley L'Hôpital limits, **the whole KcsA chain re-derived from the deposited geometry**, and the site-integrity gate, and the measured collapse of the flexible site |

`make test` fails if any suite's exit code is non-zero **or** if any suite
prints a `FAIL` line — both are checked for every suite, not just some of
them. The full provenance tables, the two audit passes, and the external
oracles are in [`docs/verification.md`](docs/verification.md).

**This release is V0.9RC5.** It follows V0.9RC4 (the deep audit that
produced the datastream, the KcsA ionic layer, the exact STO overlaps,
the Roothaan-Hartree-Fock layer, and the multicore build) and V0.9RC3
(the first audit pass). V0.9RC5's headline is the fourth audit: a full
mathematical, programming and operational re-audit that found and fixed
30 physics defects, 23 programming defects and 38 operational defects,
added 183 new checks (80 in the new eighth gate), split the 1061-line
readme into a documentation set, and corrected the record itself (which
had silently drifted away from HEAD for several months). The complete
correction history, with the measurement behind each one, is in
[`release_note_v9R4.md`](release_note_v9R4.md) under section 11.

The record also reproduces **byte-for-byte under AddressSanitizer and
UndefinedBehaviorSanitizer** with an empty stderr, which is what
`./s01_verify_record.sh` asserts.

## 6. The limitations that matter most

The full list is in [`docs/limitations.md`](docs/limitations.md). These
are the ones that change how a number may be read:

* **No quantitative KcsA selectivity free energy is claimed.** The legs
  disagree, the umbrella `±` is a three-seed spread rather than a standard
  error, and only the measured dehydration free energy is a robust K⁺
  advantage. The two-leg sum is compared with the experimental −0.179 eV
  as a scale check, never as a prediction.
* **The KcsA entropic term is implemented but does not converge, and the
  flexible site does not survive long enough to sample.** The 8-fold site
  holds for ~50 fs and then the ion's Coulomb pull closes the shell (CN 8.0
  → 6.2 at 100 fs → 4.0 at 200 fs, O–O 4.5–5.5 Å → 3.28 Å) at every
  temperature, equilibration length and tether stiffness tested. The
  site-integrity gate now **refuses** to return a variance for a collapsed
  site, so this cannot be papered over with a plausible number.
* **The Roothaan–Hartree–Fock layer is qualitative.** Its overlap matrix
  is exact and its forces are a clean gradient, but the core Hamiltonian is
  the EHT VSIP diagonal and the two-electron term is a ZDO
  approximation. The H₂ minimum sits at 1.30 Å against 0.7414 Å, and
  water is not bound at all in this basis.
* **No bulk solvent.** Everything runs in vacuum or with a
  relative-permittivity divisor; there is no periodic water box, Ewald/PME
  or membrane potential.
* **Base-pairing energetics are qualitative.** The G–C > A–T ordering is
  right; the magnitudes are not (monomer-subtracted interactions underbind
  gas-phase ab initio references, and no single dielectric fixes both
  pairs at once).
* **Harmonic bonds cannot break** unless the EHT switch is on.
* **The three tracks are unconnected.** Gating is not derived from protein
  structure.
* **QEq-lite is unscreened** and refuses solutions that invert
  electronegativity; the screened Coulomb it would need is not
  implemented, and inventing the screening radii would be an unsourced
  parameter.



## 7. Repository map

```
carbonsim, s2tui        built binaries (carbonsim tracked, s2tui untracked)
carbonsim-perf,         the OpenMP performance build's PRIVATE binaries and
s2tui-perf              object directory (build-perf/) - untracked, and
                        never able to contaminate the record build
output.txt              the record (stdout of ./carbonsim)
CURRENT_BASELINE_SHA.txt  its SHA-256
kcsa.cvmds              sealed datastream written by Demo 12 (promoted atomically with output.txt)
runs/                   staged candidates (ignored, never committed). NOT a ring
                        buffer: nothing prunes them, so the directory grows
                        until you remove it yourself. A staged run that
                        differs from the record is normal - it is a
                        candidate, not a failure.
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
docs/                   the long-form documentation, split by subject:
                        physics, kcsA, demos, verification, operations,
                        limitations, terminal (the map is at the top of
                        readme.md)
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
