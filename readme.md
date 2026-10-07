# S2 Engine — Nano Chemistry & Biological Simulator

**v9R4** (internal) · **V0.9RC4** (external) · C11 · no dependencies beyond libm
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

Requirements: a C11 compiler (gcc or clang), `make`, libm. Nothing else.

| Command | Effect |
|---|---|
| `make` / `make -j$(nproc)` | build everything: `carbonsim`, `s2tui`, and all six test binaries |
| `./carbonsim` | run the 14-demo record sequence (~2 minutes); stdout is the byte-deterministic record |
| `./run` | clean build behind a warning gate, run to `runs/<timestamp>.txt`, print old/new digests side by side |
| `./run --accept` | promote that run to `output.txt` + `CURRENT_BASELINE_SHA.txt`; refuses to promote a run whose stderr is non-empty |
| `make test` | build and run every gate — what CI should invoke |
| `make tui` | build just the interactive terminal (`./s2tui`) |
| `./verify_scripts.sh` | read-only record + spec verifier: baseline digest, `kcsa.cvmds` seal, key/unit compliance, every suite |
| `./s01_verify_record.sh` | builds normally and under ASan/UBSan, asserts both reproduce the recorded digest with empty stderr |
| `make audit` / `make audit-revert` | run the audit's independent oracles and C harnesses; `--revert` links the same regression suite against the pre-fix engine |
| `make clean` | remove `build/`, both binaries, `audit/build` |
| `s2tui` | interactive terminal; `help` prints the command map |

The record build uses `-O3 -g -Wall -Wextra -std=c11 -Iinclude` and never
`-march=native`: floating-point contraction must be portable across machines,
and the record is only meaningful if it reproduces byte-for-byte elsewhere.
Set `SOURCE_DATE_EPOCH` for byte-stable side artifacts.

**Current gate status** (all green):

| Gate | Checks | What it pins |
|---|---:|---|
| `test_datastream` | 17 | FIPS 180-4 SHA-256 known-answer vectors, round-trip, tamper rejection |
| `test_forces` | 22 | analytic dihedrals vs finite-difference oracle, net-zero forces, collinear guard |
| `test_fire` | 7 | FIRE and steepest descent reach the same minima |
| `test_regression` | 165 | independent in-repo oracles: quadrature, finite differences, hand-derived values |
| `test_external` | 49 | values from outside this repository (NIST CODATA, parm99.dat, FIPS, PDB 1K4C) |
| `test_loop` | 54 | bio/QC/QM loop closure + best-layer: Nernst, kT, Marcus/Shannon/LJ/Coulomb, BJ/penetration, Langevin, flat-bottom, tethers, H-complete ions-hold |

plus a warning-free clean build, empty stderr, byte-identical stdout across
repeated runs, and an intact `kcsa.cvmds` seal.

---

## 2. The record and its artifacts

The record is the tracked stdout of `./carbonsim`, and it is the authority
this document describes — not the other way around.

| Artifact | What it is |
|---|---|
| `output.txt` | the verbatim 14-demo record (~1080 lines) |
| `CURRENT_BASELINE_SHA.txt` | its SHA-256 (single source of truth — do not duplicate digests in prose) |
| `kcsa.cvmds` | a schema-1 datastream written by Demo 12: worked-example claims with provenance tags and a SHA-256 payload seal (`DATASTREAM_SPEC.md`) |
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

Current tree: **25 137 lines** (src 19 770 · include 2 774 · tests 2 593).

| File | Lines | Role |
|---|---:|---|
| `src/tui.c` | 6880 | interactive terminal: full POSIX shell ($(()), ${}, ``, ~, if/for/while/case), chemistry lab, gas/barostat, ASCII grid (not in the record) |
| `src/main.c` | 4189 | the 14-demo record program and datastream consumer (Demo 7 A-T + monomer-subtracted legs) |
| `src/qm.c` | 1890 | QEq, induced dipoles (first-order + coupled SCF), Pauli, TT/BJ dispersion, penetration-damped Coulomb, overlap, SCF driver |
| `src/nucleobases.c` | 1279 | five bases, deoxyribose, T-p-A dinucleotide, pairing/geometry helpers |
| `src/forces.c` | 927 | non-bonded pairs, bonded terms, harmonic/flat-bottom restraints, analytic dihedral gradients, energy breakdown |
| `src/sim.c` | 847 | lifecycle, molecule constructors, topology rebuild, ion/restraint plumbing (32 restraints) |
| `src/integrator.c` | 780 | Velocity Verlet, PCG64/LCG, Maxwell–Boltzmann, Berendsen/Andersen/Langevin, steepest descent, FIRE |
| `src/quantum.c` | 547 | Slater screening and orbital energies, hydrogenic radial profiles, Clementi–Raimondi exponents, real spherical harmonics |
| `src/kcsa_filter.c` | 767 | real 1K4C TVGYG filter, ion sizing, sites, binding/dehydration legs, C-alpha tethers, H-completion (26/subunit) |
| `src/aminoacids.c` | 513 | glycine/alanine/dipeptide/polyalanine builders |
| `src/datastream.c` | 374 | schema-1 writer, self-contained FIPS 180-4 SHA-256, seal verifier |
| `src/periodic_table.c` | 294 | H–Kr element data, UFF ε/σ, Madelung electron configurations |
| `src/tui_view.c` | 243 | TUI viewport/renderer |
| `src/neuron.c` | 164 | Hodgkin–Huxley 1952 (squid giant axon) |
| `src/loop.c` | 76 | bio/QC/QM loop closure: Nernst, kT scale, LJ/Coulomb from same primaries |
| `tests/test_regression.c` | 1466 | the 165-check in-repo oracle suite |
| `tests/test_external.c` | 536 | the 49-check outside-reference suite |
| `tests/test_loop.c` | 256 | the 54-check loop-closure + best-layer suite |
| `include/loop.h` | 55 | loop bridge contract + abstraction inventory |
| `include/` | 2758 | types, constants, per-module contracts; `amber_lj.h` and `display.h` carry shared policy |
| `tests/` (other) | 335 | datastream, forces, fire suites |

---

## 4. What the engine models

### Quantum orbital layer
Slater screening (`Z_eff = Z − S`, Slater's own 1930 iron example reproduced),
orbital energies from Slater `n*`, hydrogenic radial profiles `P(r) = r²|R_nl|²`
with exact expectation values, Clementi–Raimondi SCF exponents for spatial
ranges, and all 49 real spherical harmonics (orthonormal through l = 3).
Demo 1 includes a per-element computed-vs-NIST ledger instead of a blanket
accuracy claim.

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
  (default — keeps the record byte-stable) and PCG64.
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
| computed exchange, relaxed filter | −0.1557 ± 0.9041 eV | K⁺ (within noise) |
| two-ion exchange | −0.2706 eV | K⁺/KK |
| SCF-polar U(z) barrier, d(K−Na) | +2.5912 eV | Na⁺ pays less |
| polar-WHAM barrier gap | −0.0422 ± 0.2332 eV | K⁺ (within its error) |
| fixed-charge WHAM gap | +0.1572 eV | Na⁺ |
| knock-on pair / conductive / landscape | +1.54 / +4.76 / +9.99 eV | Na⁺ |
| measured dehydration (Marcus 1991 TATB) | K⁺ enters 0.726 eV cheaper | K⁺ |

The legs disagree, and that disagreement is the result. **No quantitative
selectivity free energy is claimed.** Only the measured dehydration free
energy is a robust K⁺ advantage in this model; the two-leg sum is compared
with the experimental −0.179 eV reference as a scale check, never as a
prediction. The umbrella sampler prints its measured autocorrelation time τ
and effective sample size with every barrier.

### Datastream
A text, self-describing, line-oriented format (`DATASTREAM_SPEC.md`, schema 1)
with a SHA-256 integrity seal over the payload, written by Demo 12 as
`kcsa.cvmds`. Claims carry provenance tags (`computed`, `computed-jc2008-params`,
`Marcus1991-TATB`, `1K4C-LINK`, …) and a validator can check external-citation
claims against reference values. The writer fails closed: truncated or
unsealed output is never written.

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
  camera keys, and create/augment/adapt keys: `w` water, `i` Na⁺, `K` K⁺,
  `C` Cl⁻, `u` base, `a` alanine, `g` glycine, `d` sugar, `H`/`L`
  heat/cool, `f` freeze, `r` replace, `p` clone, `x` delete, `m`
  minimise, `n` neuron, `e` catalysis). Piped, it prints one snapshot.
* `vi world` — modal editor (hjkl pan, `w`/`b` select, `i` insert,
  `r` replace, `x` delete, `p` clone, `u` undo, `/` search, `:` command
  mode running any shell command against the live world).

The POSIX mapping for the rest of the world: `touch` creates matter,
`ln` links, `ln -s` restrains, `unlink` clears, `rm atom <i>` deletes,
`fsck` detects topology, `sleep <n>` advances world time, `kill -SIG`
freezes / stimulates / heats, `nice` minimises, `df`/`du` report state,
`cp <species> x.mol` exports a template and `dd if=x.mol` instantiates it,
`make x.rxn` fires a reaction rule, `sync` saves, `man` documents. Demos
remain live system tests (`test all` = 14/14). The orthographic ASCII grid
has depth shading, bonds, restraint anchors and a z-slab cutaway. Stream
contract: stdout is DATA, stderr is diagnostics.

---

## 5. The demos (14)

Run `./carbonsim` for the record, or `test demo <n>` inside `s2tui` for the
fast live smoke check.

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
| 11 | α-helix | i,i+4 H-bond emerges: H···O 2.1634 Å, N···O 3.1256 Å |
| 12 | KcsA legacy cage | full selectivity leg program; vacuum sites favor Na⁺, exchange favors K⁺ |
| 12b | Real KcsA filter | deposited TVGYG, 164 atoms, C4 symmetry, CN = 8 at every site |
| 17 | DNA duplex | G–C / A–T stack with B-DNA rise and twist; G–C BREATHING / A–T HELD |

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
| Ion LJ sizes | Joung–Cheatham 2008 | parameters + derived conversion |
| Water | TIP3P | geometry and σ against literature |
| Ionic radii | Shannon 1976 (VIII) | contact distances against deposited geometry |
| Neuron | Hodgkin–Huxley 1952 | rate functions and a published worked example, 4 decimals |
| SHA-256 | FIPS 180-4 | NIST known-answer vectors |
| RNG | PCG64 (O'Neill 2014) | canonical 128-bit multiply + XSL-RR permutation |

The second-audit pass verified the whole tree against these sources, including
things that were wrong and are now fixed: the KcsA dehydration leg had been
carrying hydration *enthalpies* under a free-energy label (corrected to the
Marcus 1991 TATB free energies, K⁺ advantage 1.368 → 0.726 eV), a glycine
charge table carried a sign-inverted neutrality correction, and several
constants had duplicated spellings free to drift. Each correction is pinned
by a regression check that fails on the old code. The row-by-row history
lives in `release_note_v9R4.md` and the commit log, not here.

---

## 7. Known limitations

Stated plainly, because the numbers above cannot be read correctly without them.

* **No bulk solvent.** Everything runs in vacuum or with a relative-permittivity
  divisor; the ion program adds explicit octahedral 6-water clusters but no
  periodic water box, Ewald/PME, or membrane potential.
* **No quantitative KcsA selectivity free energy.** The real filter is rigid,
  the legs disagree, and the umbrella `±` is a three-seed spread, not a
  standard error. Read every K⁺/Na⁺ number as "in this reduced model".
* **The three tracks are unconnected.** Nucleic acids, proteins and
  electrophysiology each run standalone; gating is not derived from protein
  structure.
* **Base-pairing energetics are qualitative.** The G–C > A–T ordering is
  correct; monomer-subtracted interactions underbind gas-phase ab initio
  references (~3.5x GC, ~2.4x AT; ratio 1.5 vs 2.2 real), and no single
  dielectric fixes both pairs at once. Prior 4x-overshoot read compared
  dimer totals to interaction refs; Demo 7 now prints both.
* **Harmonic bonds cannot break.** Full separation costs infinite energy.
* **Non-nucleobase charges are approximations.** Sugar, phosphate and
  amino-acid partial charges are charge-balanced but not verified RESP fits;
  the glycine table is now exactly neutral.
* **`qm_overlap` is a heuristic, not an overlap integral.** It drops the
  Slater polynomial prefactor and maximizes `m` per atom; `qm_bond_order` and
  `qm_pauli` inherit that.
* **1-4 non-bonded scaling is deliberately absent** — AMBER's scaling
  presupposes its co-fitted torsions, which this force field does not use; the
  choice is documented in `forces.c` and was tested against the helix.
* **Cutoff switching is implemented but off by default.** Every current demo
  is gas-phase and wants the plain hard cutoff; enable the switch before any
  condensed-phase system.
* **O(N²) non-bonded loop, no neighbour list.** Fine at the current system
  sizes (≤ hundreds of atoms); pair exclusion is O(coordination).
* **Filter flexibility was tried and collapses** without a real protein force
  field: a restraint weak enough to let the filter respond cannot hold a +1
  ion's Coulomb field. The failure is documented, not tuned away.

---

## 8. Interpreting the output

* **Sign conventions differ per block and are printed where they matter.**
  The KcsA forensic table uses `dU = E(K⁺) − E(Na⁺)`, so `+ = Na⁺ favored`;
  exchange lines use `negative = K⁺ selective`.
* **`±` values are spreads, not errors.** The WHAM `±` is the sample standard
  deviation across 3 independent seed repeats — with n = 3 its own uncertainty
  is ~52%, and the output says so. Within-run sampling error is characterized
  separately by τ and N_eff, also printed.
* **Verdicts** ("HELD", "BREATHING", "EMERGED") are the model's own
  thresholds, calibrated to the model's baseline behavior — not experimental
  pass/fail.
* **The record is stdout and only stdout.** stderr is display; if you need to
  quote a number, quote it from `output.txt` and check the digest.

---

## 9. Repository map

```
carbonsim, s2tui        built binaries (carbonsim tracked, s2tui untracked)
output.txt              the record (stdout of ./carbonsim)
CURRENT_BASELINE_SHA.txt  its SHA-256
kcsa.cvmds              sealed datastream written by Demo 12
run, verify_scripts.sh  record workflow and read-only verifier
s01_verify_record.sh    ASan/UBSan record reproducibility harness
makefile                one command builds every executable
readme.md               this file
release_note_v9R4.md    the release's corrections and rationale (history)
DATASTREAM_SPEC.md      the .cvmds format, schema 1
src/, include/, tests/  the engine and its suites
audit/                  independent oracles, C harnesses, external references,
                        pre-fix engine copies (make audit / make audit-revert)
```
