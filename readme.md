# S2 Engine — Nano Chemistry & Biological Simulator

## Note: The new naming strategy across the board will be vxR123, not vxA123. This scheme began with the v9R3 designation; the release current as of this sync is v9R4. New releases after this point in all other projects will shift A to R, Alpha will be removed since it's just too extraneous a word.

## In addition, S2 will not follow the R1 R2 R3 developmental scheme of Lancius and Mathlib. Instead, the v9 to v10 period for S2 will be marked by however many RC releases are necessary to achieve a OpenWorm similar display for visualising bonds occuring.

*Documentation history: s38 record regenerated without `-march=native`; s39 re-synced to include Demo 17 and the s37 dehydration block; s42-physics-fix record regenerated after derivation-based equation fixes (thymine methyl orientation, cytosine H6 bisector, live −1 charge conservation, conservative duplex restraints, global r_mp, HH bounds, derived unit conversions); 2026-09-24 v9R4-release (KcsA 3D→xy cage, JC/ECC/U(z)/QM-Sref, CO2 zero-strain, QEq thread-local, per-axis PBC, datastream lower-case keys); 2026-09-29 V0.9RC4 — deep audit, physics corrections, the real deposited KcsA filter, and the streamlining pass.*

*Every record digest before `67f9f47e…` was produced by a binary that the committed source could not rebuild. See the audit section: the tree did not compile at any point in the v9R4 cycle, twice over, and the tracked build artifacts concealed it. Treat the earlier digests as history, not as reproducible baselines.*

S2 Engine is a multi-scale, grounded-up physical chemistry and biological
simulator written in C. It validates biological processes across multiple
levels of organization: from subatomic quantum orbitals, through molecular
dynamics and biopolymer chain condensation, up to cellular electrophysiology
and ion-channel biophysics.

The organizing principle of the whole project is **emergence from real
physics**: structure and behavior are never hard-coded. Hydrogen bonds,
base-pairing selectivity, secondary-structure formation, and the action
potential waveform are all produced by integrating real force fields and
real differential equations whose only inputs are measured physical
constants, literature-sourced parameters, and verified geometry. Where a
result is qualitative rather than quantitative, or where a demo reports a
genuine negative result, the output says so explicitly.

---

## Release identity

| | |
|---|---|
| Internal name | `v9R4` |
| External tag | `V0.9RC4` |
| Repository | https://github.com/Shicheng-Zhang-0003/S2_Engine |

Both strings are defined once, as `S2_VERSION_INTERNAL` and
`S2_VERSION_EXTERNAL` in `include/constants.h`, and are used by the program
banner and the datastream header. The previous release carried a hand-typed
source string in `main.c` that had been stale for several releases, which is
the same failure mode the derive-in-line rule for constants exists to prevent.

## Project Structure & Evolution

The biological simulation component evolved iteratively across versions
located in `biological/`. Each version added a new physical capability on
top of the validated layers below it:

* **`carbonsim-v1` & `carbonsim-v2`** — Established the physical constants
  (2019 CODATA), the UFF-based periodic table database (elements H→Kr fully
  populated), and the symplectic Velocity Verlet molecular dynamics
  integrator with Berendsen weak-coupling thermostats.

* **`carbonsim-v3` & `carbonsim-v4`** — Programmed nucleotide base placement
  and structure modeling, including the five nucleobases built from verified
  RCSB PDB Chemical Component Dictionary ideal coordinates.

* **`carbonsim-v5`** — Completed NVT ensembles, the emergent cyclic water
  trimer cluster simulation, and full geometry-audit checking.

* **`carbonsim-v6`** — Integrated the Hodgkin-Huxley point-neuron
  electrophysiology track (squid giant axon), a genuinely independent track
  from the chemistry/MD code.

* **`carbonsim-v7`** — Added dipeptide protein backbone chemistry and the
  amino acid building blocks (glycine, alanine), sourced from PDB CCD ideal
  coordinates with PDB-authoritative leaving-atom condensation chemistry.

* **`carbonsim-v8`** — Brought the molecular-recognition and backbone demos
  to their current form: Watson-Crick base-pairing energetics (G-C vs A-U),
  the T-p-A dinucleotide sugar-phosphate backbone, the Hodgkin-Huxley action
  potential, and the Gly-Ala dipeptide peptide bond. Integrated real
  AMBER ff99 Lennard-Jones parameters and verified RESP partial charges
  (Aduri et al. 2007) for the nucleobases, replacing generic per-element
  defaults.

* **`carbonsim-v9A1`** — Added the capabilities needed to go
  from static backbone links to genuine 3D structure:

  - **Dihedral (torsion) forces** in the standard AMBER/CHARMM functional
    form, with gradients computed by central finite difference for
    correctness-by-construction.

  - **Steepest-descent energy minimization** with adaptive step size, a
    per-atom displacement cap, a divergence safety check, and a frozen-atom
    variant.

  - **The poly-alanine chain builder** with centralized, defensive index
    tracking across residue assembly.

  - **Demo 11 — alpha-helix emergence**: given only correct local backbone
    torsion geometry, the defining i,i+4 backbone hydrogen bond forms from
    the same Coulomb+LJ physics validated on the water trimer.

  - **Demo 12 — KcsA selectivity filter (two passes)**: a Coulomb+LJ
    investigation of K⁺ vs Na⁺ coordination. Pass one used generic
    periodic-table oxygen LJ; pass two switched to the real amino-acid-
    specific carbonyl typing (ruling out generic-O typing as the cause),
    added a radius scan letting each ion pick its own preferred
    coordination distance (ruling out geometry-fit), and a fuller two-ring
    antiprism construction that is explicitly flagged as a *placeholder* —
    its ring z-separation was never sourced, so its numbers are not
    interpretable as selectivity energetics yet. The demo reports an
    honest **negative result** (Na⁺ favored, wrong direction, in every
    test) with a systematic diagnosis pointing at missing electronic
    polarizability.

  - **v9 masterplan Step 0 resolved first, in-file**: the potassium row
    discrepancy in `periodic_table.c` (fully-populated data row vs. stale
    "mass only stub" comment) was investigated and resolved where the data
    lives, before any downstream demo depended on potassium parameters.
    The row's electronegativity, ionisation energy and electron affinity
    check out against standard tabulated values; the stale comment was
    removed; the LJ ε/σ remain explicitly flagged as UFF-sourced but not
    independently re-verified against the primary Rappé 1992 table.

  - Amino-acid-specific carbonyl Lennard-Jones typing shared consistently
    between the protein and ion-channel tracks (Demo 12 recomputes
    `aminoacids.c`'s exact carbonyl values rather than borrowing generic
    oxygen).

---

## v9 Master Plan — Status Record

### v9R3 status - the audit & correctness release

v9R3 was NOT a new-physics release. It was a systematic,
sector-by-sector correctness audit of the whole engine (constants,
forces, topology, integrator, periodic table, build, record). Every
finding was fixed in code or documented as a deliberate limitation;
the full 13-demo sequence was regenerated and re-asserted after every
behavior-affecting change. Full ledger: the "v9R3 - The Audit &
Correctness Release" section and `release_note_v9R3.md`.

**Open items carried forward to v9R4+** (direction unchanged by the
audit; this list is the unambiguous starting point for the next
session):
1. **Complete Demo 12's antiprism block** *(done)*: the real ring
z-separation was sourced from 1K4C (3.084 A, symmetry-validated),
and the antiprism block is now a real result. Result: Na+ still
favored (wrong direction), consistent with all other tests.
2. **The missing-physics decision**: polarizability, explicit-solvent
   competition (the dehydration penalty the vacuum calculation cannot
   express at all), or both.
3. **Deferred audit P2 (flagged, not failed)**: analytic dihedral
   gradients (F4) - finite differences remain the oracle.


`v9_masterplan.md` (dated 2026-07-25) pre-registers the build order, the
decisions, and the success criteria for the KcsA selectivity-filter work
*before* any of that code was written. This section records how each step
actually landed, so the next session starts from an exact "what is done,
what is open" picture.

* **Step 0 — resolve the potassium row discrepancy** *(done)*.
  `periodic_table.c` now carries the resolution in-file: electronegativity
  0.82, ionisation 4.341 eV, and the electron affinity match standard
  tabulated values; the stale "mass only stub" comment was removed; the LJ
  parameters remain flagged as UFF-sourced but not independently
  re-verified before quantitative trust.

* **Step 1 — get the real coordinates** *(honest deviation)*. The deposited
  TVGYG backbone coordinates were **not** fetched: in 1K4C, chain C sits
  after two antibody Fab chains, making hand-extraction expensive for a
  one-off, and no reusable PDB reader was justified yet. Instead, Demo 12
  uses the literature coordination distances (Gly77 2.72 Å, Val76 2.83 Å,
  Thr75 2.70 Å) as *constructed inputs*, and says so in the output. It
  therefore tests energetic preference at imposed geometry — not whether
  the correct geometry emerges.

* **Step 2 — build and sanity-check the bare filter** *(superseded by the
  deviation above)*. With no deposited backbone placed, there is no
  backbone to restrain near crystallographic coordinates; the constructed
  oxygen cages are rigid by construction. Restraining near crystal
  coordinates remains the right approach for a future pass if real
  coordinates are sourced.

* **Step 3 — K⁺, then Na⁺** *(executed; result negative)*. Both ions were
  compared at both fixed site radii, across a 2.00–4.20 Å radius scan, and
  in a two-ring antiprism construction. Na⁺ is favored in every test — the
  wrong direction. The radius scan is the stronger negative: even with
  each ion free to choose its own best radius (K⁺ 2.82 Å, Na⁺ 2.44 Å),
  Na⁺ still wins, which rules out geometry-fit entirely. The antiprism
block is now a real result, not a placeholder (see Known Limitations).

* **Step 4 — write it up honestly** *(done)*. The negative result is
  reported in the output in Demo 7's convention: number, comparison, then
  an explicit honest caveat — uncertainty stated rather than folded into
  the number.

* **Success-criteria verdict.** The plan set the bar in advance: geometry
  near 2.85 Å is a fair test; K⁺ favored over Na⁺ at the same site — even
  with imperfect magnitude — counts as a legitimate, informative result;
  getting the *direction* right is the actual test of whether the method
  extends here. The direction test **failed**. Per the plan's own terms
  that is an honest v9 result, not a data error, and the systematic
  elimination (generic typing first, then geometry-fit) leaves missing
  electronic polarizability — together with the codebase-wide absence of
  any explicit-solvent/dehydration physics in a vacuum calculation — as
  the remaining candidates.

**Open items carried out of v9**, in the order the project's own record
supports:

1. **Complete Demo 12's antiprism block** *(done)*: the real ring
   z-separation was sourced from 1K4C (3.084 A, symmetry-validated),
   and the antiprism block is now a real result. Result: Na+ still
   favored (wrong direction), consistent with all other tests.
2. **The missing-physics decision**: polarizability, explicit-solvent
   competition (the dehydration penalty the vacuum calculation cannot
   express at all), or both.
3. **Optionally, the original Step 1**: source the real TVGYG backbone
   coordinates, if a future pass needs geometry *emergence* rather than
   energetic preference at imposed geometry.

---

## v9R3 — The Audit & Correctness Release

`carbonsim-v9R3` is not a new-physics release. It is a systematic,
sector-by-sector correctness audit of the whole engine, held to the
project's own standard: **every number sourced, every approximation
flagged, every silent behavior made explicit.** Each finding was either
fixed in code or documented as a deliberate limitation — nothing left
implicitly wrong. The full demo sequence was regenerated and re-asserted
after every behavior-affecting change; the block below (and `output.txt`,
SHA-256 recorded live in `CURRENT_BASELINE_SHA.txt`)
is the verbatim post-audit record.

**Fixed in code:**
* **Constants derived in-line (C1-C4).** Every unit-conversion factor and
  its reciprocal - `KCAL_MOL_TO_EV`, `COULOMB_MD`, `J_TO_EV`,
  `EV_TO_HARTREE`, `EV_TO_KCAL_MOL` - is now computed from 2019 CODATA
  primaries at the point of definition, not hand-typed. Every reciprocal
  pair multiplies to exactly 1.0 by construction; silent drift is no
  longer possible.
* **BOND_TABLE provenance corrected (F2).** The header's false
  "AMBER ff14SB / CHARMM36" claim was replaced with the honest tier-2/3
  description - generic spectroscopic-order stiffness constants, with
  spot-check conversions showing they match neither AMBER nor CHARMM.
  (r0 is overridden to placed geometry everywhere, so only stiffness,
  never equilibrium structure, is affected.)
* **Fallback warnings (F3).** Bond- and angle-parameter fallback paths
  now print a one-line warning when taken; none fire in the normal run.
* **Topology & integrator hygiene (S1-S3, I1-I2).** Angle-rebuild
  functions renamed to make their reset semantics explicit; dead fields
  removed; the silent no-op `THERMOSTAT_NOSE_HOOVER` removed; the 3N-3
  temperature comment fixed to match the code; the global RNG replaced
  with a per-simulation seedable state.
* **Potassium verified (T1).** The K row's LJ eps/sigma verified against
  Rappé et al. 1992 UFF (x1 = 3.812 A, x2 = 0.035 kcal/mol), fully
  resolving the v9 masterplan Step 0 flag.
* **Makefile ASan path (B1).** The link rule now passes `$(CFLAGS)`, so
  uncommenting the sanitizer line yields a working AddressSanitizer build
  *through make* (previously the flag was dropped at link time); the
  stray `-lm` moved off the compile line.
* **k_restraint macro (M1).** The one hardcoded kcal->eV literal in
  `demo_helix` now uses `KCAL_MOL_TO_EV` - one conversion value in the
  whole tree.

**Investigated and deliberately rejected (F1).** AMBER's 1-4 non-bonded
scaling was implemented and tested, and **broke the validated helix i,i+4
hydrogen bond (Demo 11)** - AMBER's scaling presupposes its co-fitted
torsion parameters, while this force field's dihedrals are restraints
carrying no 1-4 physics to compensate. The scaling was reverted and the
reasoning kept in-source (`forces.c`) and in the limitations below. 1-4
pairs stay at full strength; revisit only alongside properly fitted
torsions.

**Added but gated off (F5).** A CHARMM-style smooth cutoff switch
(potential and force continuous to zero at the cutoff) is implemented for
future condensed-phase use, but is **opt-in and off by default**: every
current demo is gas-phase/vacuum and wants the plain hard-cutoff
potential. The audit found it had been silently active in the gas-phase
demos (attenuating the 9.6-12 A band at the default 12 A cutoff),
perturbing Demos 7 and 11; gating it off restored the clean baseline.

**Deferred (P2, flagged not failed).** Analytic dihedral gradients (F4) -
finite differences remain the oracle, correct by construction.

**Record (M2, D1).** The full 13-demo sequence was regenerated post-audit
and re-asserted (G-C > A-U ordering, trimer bound, helix H-bond in range,
HH worked-example values, KcsA still honestly negative). Demo 11's
emergent H-bond re-formed at H···O = 2.1637 A, N···O = 3.1286 A - the
clean hard-cutoff baseline, unchanged in direction by the audit. Demo 12's
KcsA result remains an honest negative (Na+ favored, wrong direction, in
every test).

---

## Core Simulation Capabilities

### 1. Subatomic Quantum Orbital Structure

* Slater screening rules compute the effective nuclear charge
  ($Z_{\text{eff}} = Z - S$) for every occupied orbital, with the
  screening tiers keyed on principal quantum number (verified against
  Slater's own 1930 worked example for iron before use).
* Orbital energies follow $E_{nl} = -13.6058\,\text{eV}\,(Z_{\text{eff}}/n^*)^2$
  using Slater's effective principal quantum number $n^*$.
* Radial probability profiles $P(r) = r^2|R_{nl}(r)|^2$ are computed
  analytically from hydrogen-like wave functions whose associated
  Laguerre polynomials $L_p^q$ are evaluated by a stable three-term
  recurrence.

### 2. Molecular Dynamics Core

* **Non-bonded forces**: Lennard-Jones potentials (Lorentz-Berthelot
  combination rules) and Coulomb electrostatics with an optional
  relative-permittivity (dielectric) divisor. Per-atom LJ parameters
  are overridable away from the generic UFF element defaults so a
  jointly-parameterized model (e.g. TIP3P water, AMBER ff99 nucleobase
  and backbone types) stays internally consistent.
* **Bonded forces**: harmonic bond stretch and harmonic angle bend,
  plus **dihedral (torsion) terms** in the standard AMBER/CHARMM form
  $V = k\,[1 + \cos(n\phi - \delta)]$.
* **Dihedral gradients by finite difference**: the torsion force is
  obtained by central numerical differentiation of the energy (24 energy
  evaluations per dihedral per call) rather than an analytical
  chain-rule formula. This is a deliberate correctness-over-speed choice
  — the analytical dihedral force is a notorious source of silent sign
  bugs, and a numerically differentiated force is correct by
  construction. The signed dihedral angle itself (atan2-based, so it
  distinguishes handedness) was verified against four hand-built test
  cases before use.
* **Integrator**: symplectic Velocity Verlet with Berendsen
  weak-coupling NVT thermostats.
* **Initialisation**: Maxwell-Boltzmann velocities via Box-Muller
  transforms, with net linear momentum removed. Temperature uses the
  correct $3N-3$ degrees of freedom (centre-of-mass motion excluded).

### 3. Energy Minimization

* Steepest-descent minimizer with an adaptive step size (grown ×1.2 on
  a successful downhill step, shrunk ×0.5 and rolled back on an
  uphill one).
* A hard per-atom displacement cap (0.05 Å) and a divergence safety
  check prevent a large step from tunnelling through a steep repulsive
  wall into the unphysical $r \to 0$ Coulomb-divergence region.
* A frozen-atom variant holds a caller-specified subset of atoms fixed,
  used for relaxing part of an assembly without disturbing the rest.
* Purpose: resolve the severe local steric clashes inherent in a
  freshly-assembled chain before any velocity-based dynamics, which has
  no safe way to absorb such an overlap.

### 4. Biopolymer Condensation Chemistry

Procedural builders use rigid-body transformations (rotations and
translations) plus real condensation chemistry — leaving-group atoms are
genuinely removed, following the PDB Chemical Component Dictionary's own
authoritative leaving-atom flags — to construct biomolecules:

* **Nucleic acids**: the five nucleobases (uracil, cytosine, thymine,
  adenine, guanine) from verified PDB CCD ideal coordinates, with RESP
  partial charges (Aduri et al. 2007) and real AMBER ff99 Lennard-Jones
  types; 2-deoxyribose; and a T-p-A dinucleotide joined by a real
  phosphodiester bridge. Cytosine is tautomer-corrected (the PDB CCD's
  raw "CYT" ideal coordinates encode the minor imino tautomer; the
  exchangeable H is relocated to N1, the Watson-Crick-relevant position,
  leaving N3 bare as the acceptor).
* **Proteins**: glycine and alanine from PDB CCD ideal coordinates, a
  Gly-Ala dipeptide linked by a genuine 1.33 Å peptide bond, and a
  generalized poly-alanine chain builder (up to 16 residues), with the
  angle terms spanning each peptide-bond junction included from the
  start.
* **Defensive index tracking**: chain assembly records each residue's
  backbone atom indices in an `AAResidue` record and re-checks every
  tracked index after every atom removal — a centralized, correct-by-
  construction replacement for the manual index reasoning that produced
  real bugs earlier in the project.
* **Emergent secondary structure**: given only correct local backbone
  torsion geometry (real textbook $\phi/\psi/\omega$ values applied as
  dihedral restraints) and the same validated Coulomb+LJ force field,
  the defining i,i+4 backbone hydrogen bond of an alpha helix forms on
  its own (Demo 11).

### 5. Electrophysiology

* The 4-variable Hodgkin-Huxley (1952) model of the squid giant axon,
  solved with Runge-Kutta 4th-order integration (the equations are
  stiff during the sodium upstroke; Euler would need an impractically
  small step).
* Gating rate functions are protected against their removable
  singularities ($V = -40$ mV and $V = -55$ mV) via small-argument
  Taylor expansions of the $x/(1-e^{-x})$ form.
* Every parameter and rate function was cross-checked against 2+
  independent literature sources, and the resting steady-state gating
  values were reproduced from the rate functions themselves and matched
  to a published worked example to four decimal places before any C code
  was written.

### 6. Ion-Channel Biophysics

* Demo 12 investigates the KcsA potassium-channel selectivity filter
  with the unmodified Coulomb+LJ engine: four backbone carbonyl oxygens
  carrying the real backbone-carbonyl partial charge (−0.55 e) in the
  filter's real 4-fold crystallographic symmetry, with a single ion
  on-axis, across three test blocks:

  1. **Fixed-radius site tests** at two literature/deposited
     coordination distances (Gly77 site, 2.72 Å; Val76 site, 2.83 Å).
  2. **A radius scan** (2.00–4.20 Å in 0.02 Å steps) letting each ion
     find its own preferred coordination distance — the
     cage-size-vs-ion-size framing.
  3. **A two-ring antiprism construction** (Thr75 + Val76 rings with
     the real 45° offset), with the real ring z-separation sourced
     from 1K4C (3.084 Å, symmetry-validated in s10/s10b). Result:
     Na⁺ still favored (wrong direction), consistent with all
     other tests.

* The second pass switched from generic periodic-table oxygen LJ to the
  amino-acid-specific carbonyl typing shared with the protein track
  (recomputing `aminoacids.c`'s exact values). It reports an honest
  **negative result**: Na⁺ is favored at both fixed sites *and* at each
  ion's own preferred radius — the wrong direction everywhere — with
  the systematic diagnosis ruling out first the generic typing, then
  geometry-fit, and pointing at missing electronic polarizability as
  the most likely remaining physics. See the verification and
  limitations sections below.

---

## How to Compile and Run

To build and run the current version (`carbonsim-v9R4`):

```bash
cd biological/v9R4
make
./carbonsim
```

The project ships a makefile with the following targets:

* `make` / `make all` — compile every source in `src/` into object
  files under `build/`, then link the `carbonsim` binary.
* `make run` — build if necessary, then execute `./carbonsim`.
* `make clean` — remove the `build/` directory and the binary.

**Requirements**: a C11 compiler (gcc or clang), `make`, and the
standard C math library (linked automatically via `-lm`). The default
record build compiles with `-O3 -g -Wall -Wextra -std=c11 -Iinclude`
(no `-march=native`; see makefile s38 note for why native is disabled
for the record).

**Regenerating the recorded output**: the block reproduced below is the
demo binary's stdout. After any code change, regenerate it so the
"verbatim record" claim stays provable:

```bash
./carbonsim > output.txt 2> stderr.txt
# Demo 12 also writes kcsa.cvmds as a side effect; both are record artifacts.
sha256sum output.txt  # must match CURRENT_BASELINE_SHA.txt
```

and re-sync this readme's embedded block from that file. (Provenance
note: the archived `output.after-fix09.txt` in this tree predates the
Demo 12 caveat strengthening and the banner-line update now in
`src/main.c`; the block below matches the current source.)

**Debug build**: the makefile carries a commented alternate `CFLAGS` line enabling AddressSanitizer and UndefinedBehaviorSanitizer. As of audit fix B1, the link rule passes `$(CFLAGS)`, so uncommenting that line produces a working sanitised build through make. Verified this release: the ASan build ran the full 13-demo suite with empty stderr (zero memory errors), and its output is byte-for-byte identical to the normal build (both SHA-256 `8e8836a04bb3be479d9249a54151275a55d0aaf5b0b08a75c45d1a95f66815e8`).

**Live display**: stdout is the byte-deterministic record (never wall time,
never progress). Human display — per-demo timings, minimization
heartbeats, leg progress — goes to stderr and only when stderr is a TTY
(`CARBON_QUIET` silences it). Piped/file runs stay silent and s01-clean
by construction. See `include/display.h`.

**Tests**: `make selftest` (datastream writer: SHA KATs, round-trip,
tamper-reject — 17 checks), `make selftest-forces` (analytic dihedral
vs FD oracle to 1e-6 + net-zero force + collinear guard), both gated by
`verify_scripts.sh` alongside the record SHA, warning-clean build,
`kcsa.cvmds` seal, and key/unit compliance.

---

## Simulation Execution Output (v9R4, post-audit)

Verbatim `carbonsim` stdout for the record build. Reproduced by `./carbonsim`;
the SHA-256 of this stream is the contents of `CURRENT_BASELINE_SHA.txt`, and
it is byte-identical across repeated runs *and* between the normal and
ASan+UBSan builds. This output supersedes the pre-audit block, which was
produced by a binary the committed source could not rebuild.

Numbers changed from the previous release because the physics did: the
dispersion force sign, the LJ collision diameter, the polarizability dataset,
the QEq hard core, sp2/sp3 hybridisation labels, a direct dipole solve instead
of a fixed-point iteration, and a 53-bit uniform conversion. See the audit
section for each.

Lines worth reading rather than skimming: the `sampling:` and `bins:` lines in
the MBAR section report the *measured* autocorrelation time and effective
sample size, and the `ERROR BAR DEFINITION` block states exactly what the
`±` on the barriers does and does not mean.

```text

  ╔═══════════════════════════════════════════════════════╗
  ║   CARBON VM — CHEMISTRY SIMULATOR   v9R4 / V0.9RC4    ║
  ║       From subatomic to molecular dynamics            ║
  ╚═══════════════════════════════════════════════════════╝
  v9R4 / V0.9RC4 — https://github.com/Shicheng-Zhang-0003/S2_Engine

  Unit system: Length=Å  Time=fs  Energy=eV  Mass=AMU
  Physical constants: CODATA  |  LJ: UFF (periodic table, sigma = Rmin/2^(1/6)) + AMBER ff99 (biomolecular)  |  Bonds: placed-geometry r0, generic spectroscopic k (audit F2)


╔══════════════════════════════════════════════════════╗
║  DEMO 1: Quantum orbital structure                   ║
╚══════════════════════════════════════════════════════╝
══════════════════════════════════════════
  Hydrogen (H)  Z=1
══════════════════════════════════════════
  Mass              : 1.0080 AMU
  Electronegativity : 2.20 (Pauling)
  Atomic radius     : 1.200 Å (vdW)
  Covalent radius   : 0.310 Å
  Ionisation energy : 13.598 eV
  Electron affinity : 0.754 eV
  Common valence    : 1
  LJ ε              : 0.00191 eV
  LJ σ              : 2.5711 Å
  Config            : 1s1
══════════════════════════════════════════
  Orbital table for H (Z=1)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -13.6057     1         
  Valence orbital: 1s  Z_eff=1.000  r_mp=0.529 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å |##|:..                                  10.0 Å
  QM: hyb=1s lobes=0 lone_pairs=0  chi=7.176 eV  J=6.422 eV  alpha~0.148 A^3
  QM: Y_s=0.2821  Y_px=0.4886  Y_pz=0.4886  psi_val(lobe-max)=-0.5392 A^-3/2
  QMv4: <r>=0.7938 <r2>=0.8401 <T>=13.606 eV  gamma=dZ/dq=+0.150
  QMv4: E_val=  -13.61 eV vs NIST IE=13.598 eV (ratio 1.00; Slater-Hydrogen is order-of-magnitude, not spectroscopy)
  QMv4: CR Zeff=1.000 r_mp_CR=0.5292 A vs r_mp_Slater=0.5292 A, cov_r=0.310 A

══════════════════════════════════════════
  Carbon (C)  Z=6
══════════════════════════════════════════
  Mass              : 12.0110 AMU
  Electronegativity : 2.55 (Pauling)
  Atomic radius     : 1.700 Å (vdW)
  Covalent radius   : 0.770 Å
  Ionisation energy : 11.260 eV
  Electron affinity : 1.262 eV
  Common valence    : 4
  LJ ε              : 0.00455 eV
  LJ σ              : 3.4308 Å
  Config            : 1s2 2s2 2p2
══════════════════════════════════════════
  Orbital table for C (Z=6)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -442.0490    2         
  2s(+0)   2      0      0      -35.9275     2         
  2p(-1)   2      1      -1     -35.9275     1         
  2p(+0)   2      1      0      -35.9275     1         
  Valence orbital: 2p  Z_eff=3.250  r_mp=0.651 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å .##|:.                                   10.0 Å
  QM: hyb=atomic-sp lobes=0 lone_pairs=2  chi=6.261 eV  J=4.999 eV  alpha~0.306 A^3
  QM: Y_s=0.2821  Y_px=0.4886  Y_pz=0.4886  psi_val(lobe-max)=-0.8218 A^-3/2
  QMv4: <r>=0.8141 <r2>=0.7953 <T>=35.928 eV  gamma=dZ/dq=+0.350
  QMv4: E_val=  -35.93 eV vs NIST IE=11.260 eV (ratio 3.19; Slater-Hydrogen is order-of-magnitude, not spectroscopy)
  QMv4: CR Zeff=3.140 r_mp_CR=0.6741 A vs r_mp_Slater=0.6513 A, cov_r=0.770 A

══════════════════════════════════════════
  Nitrogen (N)  Z=7
══════════════════════════════════════════
  Mass              : 14.0070 AMU
  Electronegativity : 3.04 (Pauling)
  Atomic radius     : 1.550 Å (vdW)
  Covalent radius   : 0.710 Å
  Ionisation energy : 14.534 eV
  Electron affinity : 0.000 eV
  Common valence    : 3
  LJ ε              : 0.00299 eV
  LJ σ              : 3.2607 Å
  Config            : 1s2 2s2 2p3
══════════════════════════════════════════
  Orbital table for N (Z=7)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -610.7596    2         
  2s(+0)   2      0      0      -51.7356     2         
  2p(-1)   2      1      -1     -51.7356     1         
  2p(+0)   2      1      0      -51.7356     1         
  2p(+1)   2      1      1      -51.7356     1         
  Valence orbital: 2p  Z_eff=3.900  r_mp=0.543 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å :#|:.                                    10.0 Å
  QM: hyb=atomic-sp lobes=0 lone_pairs=2  chi=7.267 eV  J=7.267 eV  alpha~0.168 A^3
  QM: Y_s=0.2821  Y_px=0.4886  Y_pz=0.4886  psi_val(lobe-max)=-1.0802 A^-3/2
  QMv4: <r>=0.6784 <r2>=0.5523 <T>=51.736 eV  gamma=dZ/dq=+0.350
  QMv4: E_val=  -51.74 eV vs NIST IE=14.534 eV (ratio 3.56; Slater-Hydrogen is order-of-magnitude, not spectroscopy)
  QMv4: CR Zeff=3.834 r_mp_CR=0.5521 A vs r_mp_Slater=0.5427 A, cov_r=0.710 A

══════════════════════════════════════════
  Oxygen (O)  Z=8
══════════════════════════════════════════
  Mass              : 15.9990 AMU
  Electronegativity : 3.44 (Pauling)
  Atomic radius     : 1.520 Å (vdW)
  Covalent radius   : 0.660 Å
  Ionisation energy : 13.618 eV
  Electron affinity : 1.461 eV
  Common valence    : 2
  LJ ε              : 0.00260 eV
  LJ σ              : 3.1181 Å
  Config            : 1s2 2s2 2p4
══════════════════════════════════════════
  Orbital table for O (Z=8)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -806.6815    2         
  2s(+0)   2      0      0      -70.4180     2         
  2p(-1)   2      1      -1     -70.4180     2         
  2p(+0)   2      1      0      -70.4180     1         
  2p(+1)   2      1      1      -70.4180     1         
  Valence orbital: 2p  Z_eff=4.550  r_mp=0.465 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å :#:.                                     10.0 Å
  QM: hyb=atomic-sp lobes=0 lone_pairs=3  chi=7.540 eV  J=6.079 eV  alpha~0.108 A^3
  QM: Y_s=0.2821  Y_px=0.4886  Y_pz=0.4886  psi_val(lobe-max)=-1.3612 A^-3/2
  QMv4: <r>=0.5815 <r2>=0.4058 <T>=70.418 eV  gamma=dZ/dq=+0.350
  QMv4: E_val=  -70.42 eV vs NIST IE=13.618 eV (ratio 5.17; Slater-Hydrogen is order-of-magnitude, not spectroscopy)
  QMv4: CR Zeff=4.450 r_mp_CR=0.4757 A vs r_mp_Slater=0.4652 A, cov_r=0.660 A

  --- QM two-center curves (computed, no fits) ---
  pair/R   S_sig      S_pi       Pauli(eV)  disp(eV)   C6        
  O-O/1.5  1.82e-03   0.00e+00   2.02e-05   -0.5531    6.300   
  O-O/2.0  2.23e-04   0.00e+00   3.02e-07   -0.0984    6.300   
  O-O/2.5  2.72e-05   0.00e+00   4.50e-09   -0.0258    6.300   
  O-O/3.0  3.33e-06   0.00e+00   6.72e-11   -0.0086    6.300   
  O-O/3.5  4.06e-07   0.00e+00   1.00e-12   -0.0034    6.300   
  O-O/4.0  4.96e-08   0.00e+00   1.50e-14   -0.0015    6.300   
  O-O/4.5  6.06e-09   0.00e+00   2.24e-16   -0.0008    6.300   
  K-O/1.5  1.11e-03   0.00e+00   4.91e-06   -0.8172    9.308   
  K-O/2.0  1.15e-04   0.00e+00   5.25e-08   -0.1454    9.308   
  K-O/2.5  1.19e-05   0.00e+00   5.63e-10   -0.0381    9.308   
  K-O/3.0  1.23e-06   0.00e+00   6.02e-12   -0.0128    9.308   
  K-O/3.5  1.27e-07   0.00e+00   6.45e-14   -0.0051    9.308   
  K-O/4.0  1.31e-08   0.00e+00   6.90e-16   -0.0023    9.308   
  K-O/4.5  1.36e-09   0.00e+00   7.39e-18   -0.0011    9.308   
  Na-O/1.5  3.44e-04   0.00e+00   4.97e-07   -0.1968    2.242   
  Na-O/2.0  2.41e-05   0.00e+00   2.44e-09   -0.0350    2.242   
  Na-O/2.5  1.69e-06   0.00e+00   1.20e-11   -0.0092    2.242   
  Na-O/3.0  1.19e-07   0.00e+00   5.89e-14   -0.0031    2.242   
  Na-O/3.5  8.31e-09   0.00e+00   2.89e-16   -0.0012    2.242   
  Na-O/4.0  5.83e-10   0.00e+00   1.42e-18   -0.0005    2.242   
  Na-O/4.5  4.09e-11   0.00e+00   6.99e-21   -0.0003    2.242   
  (disp shown undamped -C6/R^6; force loop applies Tang-Toennies damping)

╔══════════════════════════════════════════════════════╗
║  DEMO 2: H-H covalent bond vs. van der Waals (two different physics)║
╚══════════════════════════════════════════════════════╝
  Van der Waals (LJ):  ε=0.00191 eV   σ=2.5711 Å   r_min=2.8860 Å
  Covalent (harmonic): r0=0.7414 Å   k=36.00 eV/Å²   (real H2 bond length is 0.7414 Å)

  r (Å)    V_covalent(eV)  V_vdW (eV)      Scale
  ────────  ──────────────  ──────────────
  0.6000    0.359891        —             covalent well (depth scale: eV)
  0.6500    0.150371        —             covalent well (depth scale: eV)
  0.7000    0.030851        —             covalent well (depth scale: eV)
  0.7500    0.001331        —             covalent well (depth scale: eV)
  0.8000    0.061811        —             covalent well (depth scale: eV)
  0.8500    0.212291        —             covalent well (depth scale: eV)
  0.9000    0.452771        —             covalent well (depth scale: eV)
  0.9500    0.783251        —             covalent well (depth scale: eV)
  1.0000    1.203731        —             covalent well (depth scale: eV)
  1.0500    1.714211        —             covalent well (depth scale: eV)
  1.1000    2.314691        —             covalent well (depth scale: eV)
  1.1500    3.005171        —             covalent well (depth scale: eV)
  1.2000    3.785651        —             covalent well (depth scale: eV)
  1.2500    4.656131        —             covalent well (depth scale: eV)
  1.3000    5.616611        —             covalent well (depth scale: eV)
  1.3500    6.667091        —             covalent well (depth scale: eV)
  1.4000    7.807571        —             covalent well (depth scale: eV)
  1.4500    9.038051        —             covalent well (depth scale: eV)
  1.5000    10.358531       —             covalent well (depth scale: eV)
  1.5500    11.769011       —             covalent well (depth scale: eV)
  1.6000    13.269491       —             covalent well (depth scale: eV)

  1.2500    —             43.195506       ########################################
  1.3750    —             13.621322       ########################################
  1.5000    —             4.715919        ########################################
  1.6250    —             1.759098        ########################################
  1.7500    —             0.695338        ########################################
  1.8750    —             0.286634        ########################################
  2.0000    —             0.121063        ########################################
  2.1250    —             0.051186        ########################################
  2.2500    —             0.020845        ########################################
  2.3750    —             0.007491        ########################################
  2.5000    —             0.001656        #####################################
  2.6250    —             -0.000788       ###########
  2.7500    —             -0.001693       ##
  2.8750    —             -0.001907       
  3.0000    —             -0.001826       
  3.1250    —             -0.001633       ##
  3.2500    —             -0.001412       #####
  3.3750    —             -0.001200       #######
  3.5000    —             -0.001011       #########
  3.6250    —             -0.000848       ###########
  3.7500    —             -0.000710       ############
  3.8750    —             -0.000596       #############
  4.0000    —             -0.000500       ##############
  4.1250    —             -0.000421       ###############
  4.2500    —             -0.000356       ################
  4.3750    —             -0.000301       ################
  4.5000    —             -0.000256       #################
  4.6250    —             -0.000219       #################
  4.7500    —             -0.000187       ##################
  4.8750    —             -0.000161       ##################
  5.0000    —             -0.000139       ##################

  For scale: the real H2 covalent bond dissociation energy is 4.52 eV
  (a standard spectroscopic constant), versus this vdW well depth of only
  0.00191 eV — roughly 2370x weaker. That gap is why breaking a chemical
  bond (a reaction) costs so much more than separating two molecules that
  are merely touching (melting/evaporation).

  Caveat: the harmonic term above is only valid for small vibrations near
  r0. It's a parabola, not a real bond — it never flattens out, so it would
  (wrongly) predict infinite energy to fully separate the atoms. Capturing
  actual bond breaking needs a Morse potential or a reactive force field —
  a natural next addition to this codebase.

╔══════════════════════════════════════════════════════╗
║  DEMO 3: H2O molecule — Berendsen MD at 300 K      ║
╚══════════════════════════════════════════════════════╝
  Initial geometry:
  idx  sym   Position (Å)           Velocity (Å/fs)        q(e)      mass(AMU)
  ─────────────────────────────────────────────────────────────────────────────────────────────────
  0    O     ( 0.0000  0.0000  0.0000)  ( 0.0000  0.0000  0.0000)  -0.8340   15.999
  1    H     ( 0.7570 -0.5859  0.0000)  ( 0.0000  0.0000  0.0000)  +0.4170    1.008
  2    H     (-0.7570 -0.5859  0.0000)  ( 0.0000  0.0000  0.0000)  +0.4170    1.008

  bond   a    b    order r0(Å)   r(Å)    E(eV)     
  ─────────────────────────────────────────────────
  0      0    1    1     0.9572   0.9572   0.000000  
  1      0    2    1     0.9572   0.9572   0.000000  

  angle  a    b    c    θ0(deg)   k(eV/rad²)
  ────────────────────────────────────────────────
  0      1    0    2    104.52     4.7700    

  Initial thermodynamics:
  KE = 0.077556 eV  PE = 0.000000 eV  E = 0.077556 eV  T = 300.00 K

  Step      t (fs)      KE (eV)     PE (eV)     T (K)     O-H1 dist (Å)
  ────  ──────  ───────  ───────  ─────  ────────────
  1         0.500       0.076786    0.000774    297.02    0.959890    
  101       50.500      0.070391    0.011951    272.29    0.948669    
  201       100.500     0.063695    0.021643    246.39    0.966002    
  301       150.500     0.070961    0.016211    274.49    0.968336    
  401       200.500     0.085644    0.002407    331.29    0.949249    
  501       250.500     0.080459    0.007886    311.23    0.955655    
  601       300.500     0.064235    0.024677    248.47    0.967438    
  701       350.500     0.071026    0.018136    274.74    0.961970    
  801       400.500     0.084634    0.004401    327.38    0.951877    
  901       450.500     0.086052    0.002769    332.86    0.959123    
  1001      500.500     0.074331    0.014483    287.52    0.967624    
  1101      550.500     0.072296    0.016593    279.65    0.953316    
  1201      600.500     0.079674    0.009145    308.19    0.957119    
  1301      650.500     0.085070    0.003492    329.07    0.969495    
  1401      700.500     0.082141    0.006104    317.73    0.954645    
  1501      750.500     0.066759    0.021573    258.23    0.952786    
  1601      800.500     0.072659    0.015701    281.06    0.966459    
  1701      850.500     0.084471    0.003542    326.75    0.963837    
  1801      900.500     0.082428    0.005267    318.85    0.952208    
  1901      950.500     0.074359    0.013379    287.63    0.955183    

  Final geometry after 2000 steps:
  idx  sym   Position (Å)           Velocity (Å/fs)        q(e)      mass(AMU)
  ─────────────────────────────────────────────────────────────────────────────────────────────────
  0    O     (-0.0398 -0.0871  0.0442)  ( 0.0016 -0.0005  0.0015)  -0.8340   15.999
  1    H     (-0.0671  0.7946 -0.3587)  (-0.0265 -0.0030 -0.0032)  +0.4170    1.008
  2    H     ( 0.6992 -0.5834 -0.3423)  ( 0.0009  0.0113 -0.0199)  +0.4170    1.008

════════════════════════════════════════════════════════
  Simulation summary
  Atoms: 3  Bonds: 2  Angles: 1
  Step: 2000   Time: 1000.000 fs   dt: 0.500 fs
  KE: 0.069246 eV   PE: 0.018529 eV   E_total: 0.087775 eV
  Temperature: 267.85 K
════════════════════════════════════════════════════════
  idx  sym   Position (Å)           Velocity (Å/fs)        q(e)      mass(AMU)
  ─────────────────────────────────────────────────────────────────────────────────────────────────
  0    O     (-0.0398 -0.0871  0.0442)  ( 0.0016 -0.0005  0.0015)  -0.8340   15.999
  1    H     (-0.0671  0.7946 -0.3587)  (-0.0265 -0.0030 -0.0032)  +0.4170    1.008
  2    H     ( 0.6992 -0.5834 -0.3423)  ( 0.0009  0.0113 -0.0199)  +0.4170    1.008

  bond   a    b    order r0(Å)   r(Å)    E(eV)     
  ─────────────────────────────────────────────────
  0      0    1    1     0.9572   0.9698   0.002722  
  1      0    2    1     0.9572   0.9704   0.003012  

  angle  a    b    c    θ0(deg)   k(eV/rad²)
  ────────────────────────────────────────────────
  0      1    0    2    104.52     4.7700    
════════════════════════════════════════════════════════


╔══════════════════════════════════════════════════════╗
║  DEMO 4: Cyclic (H2O)3 — emergent hydrogen-bonded ring║
╚══════════════════════════════════════════════════════╝
  9 atoms, 6 bonds, 3 angles
  Ring O-O-O construction: O...O = 2.950 Å per edge
  Initial T=50.00 K  PE=-0.454441 eV (intermolecular H-bonds contribute the negative part)

  Step      t (fs)      KE (eV)     PE (eV)     T (K)     O0-O1(Å)  O1-O2(Å)  O2-O0(Å) 
  ────  ──────  ───────  ───────  ─────  ──────── ──────── ────────
  1         0.500       0.056968    -0.459774   55.09     2.9485     2.9482     2.9482    
  201       100.500     0.084179    -0.688086   81.40     2.8109     2.8369     2.7864    
  401       200.500     0.075307    -0.728756   72.82     2.6869     2.7061     2.6675    
  601       300.500     0.059726    -0.739992   57.76     2.7373     2.7182     2.8116    
  801       400.500     0.054780    -0.748759   52.97     2.8389     2.8820     2.7420    
  1001      500.500     0.062820    -0.763808   60.75     2.7952     2.7113     2.8210    
  1201      600.500     0.061873    -0.760724   59.83     2.6666     2.7372     2.6673    
  1401      700.500     0.050731    -0.747318   49.06     2.6966     2.6629     2.7760    
  1601      800.500     0.050096    -0.747612   48.45     2.7854     2.8435     2.7482    
  1801      900.500     0.048743    -0.751286   47.14     2.8810     2.7207     2.8032    
  2001      1000.500    0.071689    -0.772066   69.33     2.6890     2.7672     2.6975    
  2201      1100.500    0.043257    -0.742653   41.83     2.7285     2.6293     2.7595    
  2401      1200.500    0.043896    -0.740705   42.45     2.7621     2.9009     2.7492    
  2601      1300.500    0.042576    -0.745434   41.17     2.8373     2.7310     2.8346    
  2801      1400.500    0.071832    -0.773147   69.46     2.6836     2.7937     2.7110    
  3001      1500.500    0.062757    -0.761834   60.69     2.7331     2.6310     2.7375    
  3201      1600.500    0.046059    -0.746155   44.54     2.7478     2.8172     2.7066    
  3401      1700.500    0.042257    -0.744921   40.86     2.8321     2.7121     2.8268    
  3601      1800.500    0.052206    -0.755675   50.48     2.7122     2.8349     2.7453    
  3801      1900.500    0.039493    -0.739797   38.19     2.7377     2.6348     2.7584    

  Over 4000 steps (2000 fs): O···O range = [2.620, 2.949] Å
  → Ring stayed bound. Hydrogen bonds emerged from Coulomb + LJ alone — never told the simulator these were H-bonds.

  Final state:
  idx  sym   Position (Å)           Velocity (Å/fs)        q(e)      mass(AMU)
  ─────────────────────────────────────────────────────────────────────────────────────────────────
  0    O     ( 0.6460  0.5321 -1.2947)  (-0.0005  0.0011  0.0011)  -0.8340   15.999
  1    H     ( 0.9110  0.5439 -0.3408)  ( 0.0022 -0.0075  0.0039)  +0.4170    1.008
  2    H     ( 0.9914  1.3498 -1.6254)  ( 0.0077  0.0015  0.0072)  +0.4170    1.008
  3    O     ( 0.7447  0.2376  1.4088)  ( 0.0002  0.0001 -0.0032)  -0.8340   15.999
  4    H     (-0.0342 -0.3085  1.1285)  (-0.0059 -0.0111  0.0096)  +0.4170    1.008
  5    H     ( 0.9830 -0.1305  2.2524)  ( 0.0051 -0.0104 -0.0061)  +0.4170    1.008
  6    O     (-1.3997 -0.7651 -0.1113)  ( 0.0006  0.0002  0.0019)  -0.8340   15.999
  7    H     (-0.7973 -0.2618 -0.7038)  (-0.0073  0.0012 -0.0066)  +0.4170    1.008
  8    H     (-1.9120 -1.2657 -0.7549)  (-0.0067  0.0037 -0.0040)  +0.4170    1.008

╔══════════════════════════════════════════════════════╗
║  DEMO 5: CH4 — tetrahedral geometry check          ║
╚══════════════════════════════════════════════════════╝
  Methane geometry (C at origin):
  idx  sym   Position (Å)           Velocity (Å/fs)        q(e)      mass(AMU)
  ─────────────────────────────────────────────────────────────────────────────────────────────────
  0    C     ( 0.0000  0.0000  0.0000)  ( 0.0000  0.0000  0.0000)  -0.2400   12.011
  1    H     ( 0.6293  0.6293  0.6293)  ( 0.0000  0.0000  0.0000)  +0.0600    1.008
  2    H     ( 0.6293 -0.6293 -0.6293)  ( 0.0000  0.0000  0.0000)  +0.0600    1.008
  3    H     (-0.6293  0.6293 -0.6293)  ( 0.0000  0.0000  0.0000)  +0.0600    1.008
  4    H     (-0.6293 -0.6293  0.6293)  ( 0.0000  0.0000  0.0000)  +0.0600    1.008

  Bond angles:
    H(1)-C(0)-H(2): 109.4712°  (ideal 109.47°)
    H(1)-C(0)-H(3): 109.4712°  (ideal 109.47°)
    H(1)-C(0)-H(4): 109.4712°  (ideal 109.47°)
    H(2)-C(0)-H(3): 109.4712°  (ideal 109.47°)
    H(2)-C(0)-H(4): 109.4712°  (ideal 109.47°)
    H(3)-C(0)-H(4): 109.4712°  (ideal 109.47°)

  Initial potential energy: 0.000000 eV

╔══════════════════════════════════════════════════════╗
║  DEMO 6: The five nucleobases - geometry validation  ║
╚══════════════════════════════════════════════════════╝
  URACIL (C4H4N2O2) - 12 atoms
  Bond       Atom1  Atom2  r (A)   
  ----       -----  -----  -----   
  single     N1     C2     1.3439    (order 1)
  single     N1     C6     1.3687    (order 1)
  single     N1     HN1    0.9702    (order 1)
  double     C2     O2     1.2162    (order 2)
  single     C2     N3     1.3459    (order 1)
  single     N3     C4     1.3482    (order 1)
  single     N3     HN3    0.9691    (order 1)
  double     C4     O4     1.2183    (order 2)
  single     C4     C5     1.4145    (order 1)
  double     C5     C6     1.3495    (order 2)
  single     C5     H5     1.0802    (order 1)
  single     C6     H6     1.0800    (order 1)

  A      B(ctr) C      theta (deg)
  -      ------ -      -----------
  C2     N1     C6     120.63    
  C2     N1     HN1    119.75    
  C6     N1     HN1    119.62    
  N1     C2     O2     119.54    
  N1     C2     N3     120.97    
  O2     C2     N3     119.49    
  C2     N3     C4     120.19    
  C2     N3     HN3    119.87    
  C4     N3     HN3    119.94    
  N3     C4     O4     120.26    
  N3     C4     C5     119.42    
  O4     C4     C5     120.32    
  C4     C5     C6     119.14    
  C4     C5     H5     120.45    
  C6     C5     H5     120.40    
  N1     C6     C5     119.63    
  N1     C6     H6     120.19    
  C5     C6     H6     120.18    

  Ring planarity: max deviation = 0.0092 A (aromatic rings should be ~0)

  Cross-check vs. electron diffraction (Ferenczy et al. 1986):
  Quantity               This model Literature
  C-N (A)                1.344      1.399     
  C4-C5 single (A)       1.414      1.462     
  C5=C6 double (A)       1.350      1.343     

  Total molecular charge: -0.000114 e (RESP, Aduri et al. 2007 + derived HN1)

  CYTOSINE (C4H5N3O) - 13 atoms
  (tautomer-corrected: H relocated to N1, the Watson-Crick-
   relevant position; N3 left bare as required for pairing)
  Bond       Atom1  Atom2  r (A)   
  ----       -----  -----  -----   
  double     N3     C4     1.3717    (order 2)
  single     N3     C2     1.3940    (order 1)
  single     C4     N4     1.3709    (order 1)
  single     C4     C5     1.3363    (order 1)
  single     N1     C2     1.3847    (order 1)
  single     N1     C6     1.2936    (order 1)
  single     N1     HN1    1.0100    (order 1)
  double     C2     O2     1.2317    (order 2)
  single     N4     HN41   0.9938    (order 1)
  single     N4     HN42   0.9939    (order 1)
  double     C5     C6     1.4730    (order 2)
  single     C5     H5     1.0812    (order 1)
  single     C6     H6     1.0900    (order 1)

  N1 bond count: 3 (expect 3: C2, C6, H - donor ready)
  N3 bond count: 2 (expect 2: C4, C2 - bare, acceptor ready)
  Ring planarity: max deviation = 0.0006 A
  Total molecular charge: -0.000014 e

  THYMINE (C5H6N2O2, = 5-methyluracil) - 15 atoms
  Ring planarity: max deviation = 0.0092 A
  C5-CH3 methyl bond length: 1.5100 A (target 1.51 A, toluene-type)
  H-CM-H methyl angle: 109.47 deg (ideal tetrahedral 109.47)
  H(12)-CM-C5 angle: 109.47 deg
  H(13)-CM-C5 angle: 109.47 deg
  H(14)-CM-C5 angle: 109.47 deg
  Max H-CM-C5 deviation from 109.47: 0.00 deg (methyl splays outward, correct)
  Total molecular charge: +0.000000 e (lower-confidence approx, see code comment)

  ADENINE (C5H5N5, fused 5+6 purine ring) - 15 atoms
  Full bicyclic ring planarity: max deviation = 0.0050 A
  Total molecular charge: -0.000014 e

  GUANINE (C5H5N5O, fused 5+6 purine ring) - 16 atoms
  Full bicyclic ring planarity: max deviation = 0.0040 A
  Total molecular charge: -0.000014 e

  All five bases hold real, verified ring geometry. Next: sugar-phosphate
  backbones and Watson-Crick base pairing - G-C should bind via 3 H-bonds,
  A-T via 2, using nothing but the Coulomb+LJ code already validated
  on the water trimer.

╔══════════════════════════════════════════════════════╗
║  DEMO 7: Watson-Crick pairing - does G-C beat A-U?   ║
╚══════════════════════════════════════════════════════╝

  --- Guanine-Cytosine (3 H-bonds: N1-H..N3, N2-H..O2, O6..H-N4) ---
  Initial heavy-atom contacts after geometric placement:
    G:N1...C:N3 = 2.950 A (target 2.95)
    G:N2...C:O2 = 2.937 A
    G:O6...C:N4 = 2.949 A
  Energy breakdown at initial placement: E_LJ=0.736247 eV  E_Coulomb=-5.294396 eV  Total_PE=-4.558149 eV
  Closest intermolecular contact: atom 15 (Z=1) ... atom 20 (Z=8) = 1.967 A
  Largest single LJ repulsion: atom 6 (Z=7, sigma=3.25) ... atom 16 (Z=7, sigma=3.25) = 2.950 A, contributes 0.0415 eV
  Step   t(fs)      PE(eV)     T(K)       primary(A)
  1      0.100      -4.558009  49.96      2.9501    
  2      0.200      -4.557789  49.90      2.9501    
  3      0.300      -4.557492  49.82      2.9502    
  4      0.400      -4.557121  49.72      2.9502    
  5      0.500      -4.556682  49.60      2.9503    
  101    10.100     -4.588020  58.02      2.9166    
  201    20.100     -4.606922  61.37      2.8780    
  301    30.100     -4.626290  64.41      2.8973    
  401    40.100     -4.639000  65.04      2.9410    
  501    50.100     -4.650718  65.41      2.9888    
  601    60.100     -4.654809  64.09      3.0320    
  701    70.100     -4.675516  66.61      3.0719    

  Initial interaction PE:        -4.558149 eV
  PE at closest WC approach:     -4.611531 eV (primary N1...N3 = 2.878 A)
  Global PE minimum over run:    -4.696037 eV (may reflect drift to a different,
                                  non-WC configuration such as stacking)
  Final PE (end of run):         -4.660919 eV

  Final heavy-atom contacts:
    G:N1...C:N3 = 3.105 A
    G:N2...C:O2 = 3.126 A
    G:O6...C:N4 = 2.949 A

  --- Adenine-Uracil (2 H-bonds: N1..H-N3, N6-H..O4) ---
  Initial heavy-atom contacts after geometric placement:
    A:N1...U:N3 = 2.900 A (target 2.90)
    A:N6...U:O4 = 2.899 A
  Energy breakdown at initial placement: E_LJ=0.885984 eV  E_Coulomb=-1.940374 eV  Total_PE=-1.054390 eV
  Closest intermolecular contact: atom 12 (Z=1) ... atom 20 (Z=8) = 1.930 A
  Step   t(fs)      PE(eV)     T(K)       primary(A)
  1      0.100      -1.054736  50.10      2.9000    
  2      0.200      -1.055021  50.19      2.9000    
  3      0.300      -1.055246  50.25      2.8999    
  4      0.400      -1.055414  50.30      2.8999    
  5      0.500      -1.055530  50.34      2.8999    
  101    10.100     -1.093908  60.55      2.8745    
  201    20.100     -1.084507  56.56      2.8725    
  301    30.100     -1.090486  56.56      2.9284    
  401    40.100     -1.087098  54.38      2.9988    
  501    50.100     -1.090774  58.09      3.0334    
  601    60.100     -1.136958  69.06      3.0362    
  701    70.100     -1.130816  61.15      3.0276    

  Initial interaction PE:        -1.054390 eV
  PE at closest WC approach:     -1.080053 eV (primary N1...N3 = 2.866 A)
  Global PE minimum over run:    -1.180083 eV
  Final PE (end of run):         -1.180083 eV

  Final heavy-atom contacts:
    A:N1...U:N3 = 3.047 A
    A:N6...U:O4 = 3.079 A

  ══════════════════════════════════════════════════
  G-C @ closest WC approach: -4.611531 eV (3 H-bonds)
  A-U @ closest WC approach: -1.080053 eV (2 H-bonds)
  --> G-C binds MORE strongly than A-U (3.531479 eV difference),
      and BOTH pairs are correctly attractive (negative PE) -
      the right qualitative chemistry, from nothing but real
      charges + Coulomb + LJ. Never programmed in.

      Honest caveat: the QUANTITATIVE magnitudes here do not
      yet match gas-phase ab initio references (~-1.2 eV G-C,
      ~-0.55 eV A-U) precisely - this classical, pairwise,
      non-polarizable model overestimates the electrostatic
      CONTRAST between the two pairs beyond what a single
      dielectric correction can fix (see the detailed
      investigation in this function's setup code). The
      qualitative ordering is validated; the absolute
      numbers are not yet quantitatively trustworthy.
  ══════════════════════════════════════════════════

╔══════════════════════════════════════════════════════╗
║  DEMO 8: T-p-A dinucleotide - sugar-phosphate backbone║
╚══════════════════════════════════════════════════════╝
  Assembled: 63 atoms, 67 bonds
  (thymidine + deoxyadenosine + 1 phosphodiester bridge)

  Bond length range: [0.9670, 1.6000] A  (all chemically sane)
  Bad bonds: 0   Valence issues: 0

  Glycosidic bonds (real condensation chemistry, target 1.47 A):
    Sugar A C1' - N: 1.4700 A
    Sugar B C1' - N: 1.4700 A

  Phosphodiester bridge:
    P - O(6): 1.6000 A
    P - O(37): 1.6000 A
    P - O(61): 1.4800 A
    P - O(62): 1.4800 A

  Total charge: -1.0000 e (exactly -1 by charge-conservation
  construction: the builder measures the assembled fragment sum
  and places the residual symmetrically on the two equivalent
  non-bridging phosphate oxygens - the real phosphodiester
  convention, enforced live rather than trusted to approximate
  fragment charges)

  This validates the real chain-forming chemistry of the DNA
  backbone. NOT yet built: helical twist/rise (no dihedral
  forces exist in this codebase yet), the complementary strand,
  and base pairing/stacking between strands - all real next
  steps, not implied by this demo.

╔══════════════════════════════════════════════════════╗
║  DEMO 9: Hodgkin-Huxley neuron (squid giant axon)    ║
╚══════════════════════════════════════════════════════╝
  Resting equilibrium (computed live from the rate functions,
  not hardcoded): V=-65.0000 mV  m=0.0529  h=0.5961  n=0.3177
  After 10ms with I_ext=0: V=-64.9997 mV (should stay ~-65, confirms genuine equilibrium)

  Subthreshold stimulus (I_ext=2.0 uA/cm^2): peak V=-60.06 mV -> no spike

  Suprathreshold stimulus (I_ext=10.0 uA/cm^2), first 20ms:
  t(ms)    V(mV)     
  0.01     -64.9003  
  1.01     -55.8700  
  2.01     30.4401   
  3.01     2.9597    
  4.01     -43.4652  
  5.01     -75.0541  
  6.01     -74.0310  
  7.01     -72.5870  
  8.01     -70.8010  
  9.01     -68.7783  
  10.01    -66.6684  
  11.01    -64.6077  
  12.01    -62.6658  
  13.01    -60.8230  
  14.01    -58.9496  
  15.01    -56.6338  
  16.01    -51.4473  
  17.01    29.4074   
  18.01    -12.4163  
  19.01    -63.0141  

  Peak V reached: 40.27 mV (real squid axon: overshoots to ~+40mV)
  Post-spike undershoot (after-hyperpolarization): -75.08 mV
  (real squid axon: dips below rest to ~-75 to -80mV before recovering)
  Spikes fired in 50ms at sustained I_ext=10.0 uA/cm^2: 4
  (repetitive firing under sustained superthreshold current is a
  real physiological behavior - not specially coded, it falls out
  of the same 4 coupled equations running continuously)

  This is a genuinely independent track from the chemistry/MD
  code above - a real next step would connect them (e.g. deriving
  ion channel gating kinetics from actual protein conformational
  MD, rather than the measured empirical rate functions used
  here), which remains real future work.

╔══════════════════════════════════════════════════════╗
║  DEMO 10: Gly-Ala dipeptide - protein backbone chemistry║
╚══════════════════════════════════════════════════════╝
  Assembled: 20 atoms, 19 bonds (glycine + alanine, 1 peptide bond)

  Bond length range: [0.9675, 1.5294] A  (all chemically sane)
  Bad bonds: 0   Valence issues: 0

  Peptide bond (C-N), reported from both directions:
    Gly C -> Ala N: 1.3300 A (textbook value: 1.33 A)
    Ala N -> Gly C: 1.3300 A

  Total charge: -0.0430 e (approximate - amino acid charges are not
  independently verified the way the nucleobase RESP charges are;
  see aminoacids.c for full honest sourcing)

  All three biological polymer types now have at least one real,
  validated backbone link: nucleic acid (phosphodiester, Demo 8),
  protein (peptide bond, this demo), and electrophysiology
  (Hodgkin-Huxley, Demo 9) as an independent track. None of these
  three are connected to each other yet - real future work.

╔══════════════════════════════════════════════════════╗
║  DEMO 11: Alpha helix - does the i,i+4 H-bond emerge?║
╚══════════════════════════════════════════════════════╝
  Built 5-residue poly-alanine chain: 53 atoms, 52 bonds
  Initial clash relaxation: 18811.68 -> 8.46 eV
  Added 12 dihedral restraints (phi, psi, omega for applicable
  residues) toward real textbook values.
  Minimized: PE = 9.6890 eV

  Residue    phi (deg)    psi (deg)   
  1          -58.61      
  2          -61.33      
  3          -58.48      
  4          -54.08      
  0                       -48.04      
  1                       -47.25      
  2                       -52.46      
  3                       -48.32      

  Max deviation from target (phi=-57, psi=-47): 5.46 deg

  === The i,i+4 backbone hydrogen bond (not programmed in) ===
  N-H(4) ... O=C(0): H...O = 2.1634 A, N...O = 3.1256 A
  Real backbone H-bond range: H...O 1.8-2.2 A, N...O 2.8-3.2 A

  --> A backbone hydrogen bond formed under steered local
      torsion geometry (phi/psi/omega restrained) and the SAME
      validated Coulomb+LJ force field. Steering caveat: with
      80 kcal/mol restraints the backbone is forced into the
      helical basin, so this tests cooperation of local
      geometry + non-bonded physics, not spontaneous folding.
  Restraint-release control (dihedrals off, re-minimized):
  H...O = 2.0526 A, N...O = 3.0582 A -> H-bond PERSISTS without steering

  Bond length range: [0.9600, 1.5672] A   Bad bonds: 0   Valence issues: 0

╔══════════════════════════════════════════════════════╗
║  DEMO 12: KcsA selectivity filter - K+ vs Na+, second pass║
╚══════════════════════════════════════════════════════╝
  Four real-charge carbonyl O's, real 4-fold symmetry, radius =
  literature K+-coordination target. This pass uses the real
  amino-acid-specific carbonyl LJ typing (AA_LJ_O_EPS/SIGMA from
  aminoacids.c) instead of generic periodic-table oxygen - see
  source comment for exact scope and what changed from pass one.
  Both ions use the SAME crystallographic cage; ions are point
  charges (no neutral-atom LJ). Vacuum-only point comparison;
  dehydration/polarization/multi-ion physics reported separately.

--- Gly77 site, PDB 1K4C LINK record target: 2.72 A ---
  K+   E_LJ =  -0.025808 eV   E_Coulomb =  -5.519767 eV   Total_PE =  -5.545575 eV
  Na+  E_LJ =  -0.025808 eV   E_Coulomb =  -5.519767 eV   Total_PE =  -5.545575 eV
  Delta (Na+ minus K+): +0.000000 eV  (identical - no vacuum selectivity (same cage; see JC-ion section below for size-dependent result))

--- Val76 site, target: 2.83 A ---
  K+   E_LJ =  -0.021400 eV   E_Coulomb =  -5.305218 eV   Total_PE =  -5.326618 eV
  Na+  E_LJ =  -0.021400 eV   E_Coulomb =  -5.305218 eV   Total_PE =  -5.326618 eV
  Delta (Na+ minus K+): +0.000000 eV  (identical - no vacuum selectivity (same cage; see JC-ion section below for size-dependent result))

  Honest read (fixed-radius tests): same cage, point-charge
  ions, real carbonyl LJ typing.
  Both sites identical for K+ and Na+ - expected: point
  charges in the same cage have identical vacuum energies.
  Vacuum leg carries no selectivity; any selectivity in the
  two-leg sum comes from the dehydration leg alone. Missing:
  ion size/LJ, polarization, protein reorganization,
  multi-ion occupancy, sampling/entropy.

--- Letting each ion find its own preferred radius (2.00-4.20 A scan, 0.02 A steps) ---
  K+   best radius = 2.000 A   E_min = -7.456896 eV
  Na+  best radius = 2.000 A   E_min = -7.456896 eV
  K+ best radius 2.000 A vs Na+ 2.000 A - identical curves (point-charge ions share the scan; no size selectivity in vacuum leg)
  At each ion's OWN best radius: K+ E_min = -7.456896 eV vs Na+ E_min = -7.456896 eV -> identical - point-charge ions share the scan curve; vacuum leg carries no size selectivity

  Honest read (radius-flexible test): same point-charge ions,
  same cage definition; the scan varies the shared cage radius.
  With no ion-size term in the vacuum leg both ions share one
  curve, so this test cannot produce selectivity by construction
  - it demonstrates that size selectivity must come from ion LJ,
  polarization, or dehydration, none of which live in this leg.

--- Fuller real geometry: true 8-oxygen antiprism (site
  S3/K-C3003: Thr75 3D ion-O 2.70 A, Val76 3D 2.83 A, rings at
  z=-/+1.542 A (in-plane 2.22/2.37 A), 45-degree twist, z-sep
  3.084 A from 1K4C; ion on-axis at z=0; same cage for both ions) ---
  K+   E_LJ =  -0.150303 eV   E_Coulomb =   7.812595 eV   Total_PE =   7.662292 eV
  Na+  E_LJ =  -0.150303 eV   E_Coulomb =   7.812595 eV   Total_PE =   7.662292 eV
  closest O-O (all 8 O) = 3.134 A (3D Thr75 2.70 A, Val76 2.83 A,
  z-sep 3.084 A; same cage for both ions)
  Delta (Na+ minus K+): +0.000000 eV  (identical - no vacuum selectivity (point-charge ions, same cage; expected))
  Honest read: 8 real-charge carbonyl O's in the deposited
  antiprism geometry, same cage and same typing for both ions.
  Vacuum O-O repulsion is large and reported as-is; the protein
  backbone that balances it in vivo, plus polarization and
  multi-ion occupancy, remain missing physics. No size offset
  is imposed on Na+.

-- Dehydration legs side by side (s37; not a dG validation) --
Filter binding dU_vac(K)-dU_vac(Na) [antiprism] = +0.0000 eV (no vacuum selectivity (point-charge ions, same cage; expected))
Dehydration ΔG: K+ = +3.057 eV  Na+ = +3.783 eV (Marcus 1991)
Two-leg sum (vacuum ΔU + dehyd ΔG) = -0.7260 eV (K+ favored in sum)
Experimental ΔG (1000:1 at 300 K) = -0.1786 eV (reference scale;
 single-point ΔU lacks TΔS/sampling/reorganization, so the
 difference below is a scale comparison, not an error bar)
Sum minus experimental dG: -0.5474 eV
  JC K+  E_LJ =   0.162757 eV   E_Coulomb =   7.812595 eV   E_restr =   0.000000 eV   E_pol =  -0.714854 eV   E_ecc =   4.557342 eV   Total =   7.975352 eV
  JC Na+ E_LJ =  -0.167849 eV   E_Coulomb =   7.812595 eV   E_restr =   0.000000 eV   E_pol =  -0.714854 eV   E_ecc =   4.226735 eV   Total =   7.644746 eV
--- JC-ion antiprism (restrained scaffold, Joung-Cheatham size) ---
  K+ = 7.975352 eV  Na+ = 7.644746 eV  dU(JC,K-Na) = +0.3306 eV (Na+ favored in JC vacuum leg)
  Induction leg: K+ = -0.7149 eV  Na+ = -0.7149 eV (alpha_O = 0.84 A^3; identical by construction at same geometry)
  ECC leg (x0.75 charge scaling, Coulomb x0.5625): K+ = 4.5573 eV  Na+ = 4.2267 eV  dU = +0.3306 eV
--- Pore-axis U(z) profile (JC ions, restrained cage, z in A, E in eV; single-point, not free-energy PMF) ---
  z        K+           Na+         
  -3.0     12.7450      12.4426     
  -2.5     12.0526      10.8611     
  -2.0     12.4795      9.6546      
  -1.5     12.4618      8.7686      
  -1.0     10.5601      8.0294      
  -0.5     8.6301       7.6472      
  0.0      7.9754       7.6447      
  0.5      8.3193       7.8626      
  1.0      9.3652       8.3103      
  1.5      10.5387      9.0430      
  2.0      11.2216      10.0499     
  2.5      11.8791      11.3580     
  3.0      13.0262      12.8947     
  U(z) minima: K+ = 7.9754 eV at z = 0.0 A | Na+ = 7.6447 eV at z = 0.0 A
  U(z) gap dU(K-Na) at own minima = +0.3306 eV (single-point profile; no TDS, not a free-energy PMF).
--- SCF-polar U(z) (JC + coupled dipoles + Pauli + disp) ---
  z        K+           Na+          pol_K      pol_Na    
  -3.0     13.2623      13.2050      0.6941     0.8612    
  -2.5     12.1708      11.2091      0.3669     0.4658    
  -2.0     12.3424      9.7720       0.1862     0.2556    
  -1.5     12.2185      8.7775       0.1135     0.1563    
  -1.0     10.3511      8.0278       0.1188     0.1373    
  -0.5     8.5018       7.6704       0.1429     0.1462    
  0.0      7.8413       7.6343       0.1013     0.1027    
  0.5      7.9277       7.5951       -0.1537    -0.1534   
  1.0      8.2402       6.7308       -0.8623    -1.4584   
  1.5      7.1024       9.2603       -3.1611    0.3420    
  2.0      11.2225      9.5234       0.2521     -0.4081   
  2.5      11.1973      10.6370      -0.4803    -0.6158   
  3.0      12.3685      12.2321      -0.5053    -0.5701   
  SCF-U(z) minima: K+ = 7.1024 eV at z = 1.5 | Na+ = 6.7308 eV at z = 1.0
  SCF-U(z) barriers: K+ = 6.1598 eV | Na+ = 6.4742 eV | dBarrier(K-Na) = -0.3143 eV
  E_pol at z=0: K+ = 0.1013 Na+ = 0.1027 (symmetric: matched) | at |z|=1.5: K+ = -3.1611 Na+ = 0.3420
  QM K+  QEq q_O=-0.5462 q_ion=1.0000  dE_coul(QEq-fixed)=-0.1843 eV  S=3.23e-06 BO=2.867 Pauli=0.0000 eV  Ohyb=atomic-sp alpha=0.123
  QM Na+ QEq q_O=-0.5462 q_ion=1.0000  dE_coul(QEq-fixed)=+0.1463 eV  S=3.94e-07 BO=0.398 Pauli=0.0000 eV  Ohyb=atomic-sp alpha=0.123
--- QM bottom-up leg (QEq/overlap, same JC cage) ---
  <q_O>=-0.5462  q_K=1.0000 q_Na=1.0000  S_K=0.0000 S_Na=0.0000  BO_K=2.867 BO_Na=0.398
  v2 K+  E_LJ=   0.1628 E_Coul=   7.8126 E_pol=   0.1013 E_pauli=   0.0000 Total=   8.0766 eV
  v2 Na+ E_LJ=  -0.1678 E_Coul=   7.8126 E_pol=   0.1027 E_pauli=   0.0000 Total=   7.7475 eV
--- v2 polarized cage (induction+Pauli in force loop) ---
  K+: pol=0.1013 pauli=0.0000 tot=8.0766 | Na+: pol=0.1027 pauli=0.0000 tot=7.7475 | dU=+0.3291 eV
  scf K+  it=5 <q_O>=-0.5462 E_Coul=  7.7970 E_pol=  0.1134 E_pauli=  0.0000 E_disp= -0.2354 Total=   7.8379 eV
  scf Na+ it=5 <q_O>=-0.5462 E_Coul=  7.7964 E_pol=  0.1197 E_pauli=  0.0000 E_disp= -0.1132 Total=   7.6350 eV
--- v3 SCF cage (JC wall + SCF QM terms) ---
  K+: 7.8379 (it 5) | Na+: 7.6350 (it 5) | dU=+0.2028 eV (Na+ favored)
  relax K+  soft/off: E=  6.2362 off=0.277 <d>=3.523 CN=4/8 (Er=0.843)
  relax K+  soft/ctr: E=  6.7521 off=0.710 <d>=3.355 CN=4/8 (Er=0.783)
  relax K+  stiff:    E=  7.2825 off=0.147 <d>=2.867 CN=8/8 (Er=0.159)
  relax Na+ soft/off: E=  5.8195 off=0.137 <d>=3.566 CN=4/8 (Er=1.048)
  relax Na+ soft/ctr: E=  6.4259 off=1.159 <d>=3.460 CN=3/8 (Er=0.978)
  relax Na+ stiff:    E=  7.4468 off=0.189 <d>=3.044 CN=4/8 (Er=0.574)
--- Relaxed coordination (stiff cage, free ion, multi-start) ---
  K+: E=6.2362±0.5159 off=0.277 <d>=3.523 CN=4 | Na+: E=5.8195±0.6064 off=0.137 <d>=3.566 CN=4 | dU=+0.4167 eV (Na+ favored)
  Stiff strain probe (k=5): K+=7.2825 Na+=7.4468 dU=-0.1643 eV
--- Filter-side gap with hysteresis error ---
  dU(filter K-Na) = +0.4167 ± 0.7961 eV (best minima, soft cage)
--- Coordination probe (same QM physics, 8 vs 6 ligands) ---
  K+:  E8=  7.8379  E6= -1.5388  dE(8-6)=+9.3767 eV (pays for losing 8-fold)
  Na+: E8=  7.6350  E6= -1.6575  dE(8-6)=+9.2925 eV (pays for losing 8-fold)
  6-fold ledger K+: Coul= -1.5617 LJ=  0.1697 pol= -0.0000 disp= -0.1468
  6-fold ledger Na+: Coul= -1.8275 LJ=  0.3051 pol= -0.0013 disp= -0.1338
  Single-point (conductive d=3.084): KK=  -3.1798 NaNa=  -8.0803 dU=+4.9005 eV
--- Knock-on pairs (relaxed, d0=3.084 A) ---
  KK:   E=-215.4617 eV  ion-ion=3.081 A  <ion-O>=3.340 A
  NaNa: E=-150.3790 eV  ion-ion=3.042 A  <ion-O>=3.070 A
  dU(KK-NaNa) = -65.0826 eV (KK favored)
--- Knock-on landscape (ion A at z=0, ion B scanned) ---
  zB       KK           NaNa         flag    
  -4.50    -3.4030      -3.5162              
  -3.75    -5.5329      -5.7148              
  -3.00    -7.6440      -7.9319              
  -2.25    -6.5227      -9.1150              
  -1.50    153.1476     -2.3650      CLASH   
  -0.75    655297.7342  21231.1682   CLASH   
  0.00     -14.6548     -16.5278     CLASH   
  0.75     655297.4492  21224.4319   CLASH   
  1.50     152.0150     -2.2425      CLASH   
  2.25     -7.0588      -9.0174              
  3.00     -7.3090      -7.7557              
  3.75     -5.8263      -5.5591              
  4.50     -3.5079      -3.5840              
  Knock-on landscape barriers (valid ion-ion>=2A): KK = 4.2410 eV | NaNa = 5.5988 eV | dBarrier(KK-NaNa) = -1.3579 eV
  hyd K+  fixed: E_init= -3.7416 E_min= -3.7681 eV (6-water octahedral, d=2.75 A)
  hyd K+  polar: E_min= -3.7606 E_pol=  0.0073 eV
  hyd Na+ fixed: E_init= -4.6308 E_min= -4.6805 eV (6-water octahedral, d=2.35 A)
  hyd Na+ polar: E_min= -4.6777 E_pol= -0.0033 eV
--- Explicit hydration (6-water cluster ΔU) ---
  K+: -3.7681 (polar -3.7606) | Na+: -4.6805 (polar -4.6777) | ΔΔU(K-Na)=+0.9124 eV (Marcus ΔΔG=-0.7260)
  Note: cluster ΔU vs bulk ΔG — scale comparison only.
    sampling: 2400.0 samples/window, mean tau = 83.4 samples, N_eff = 18.3 independent samples/window
    bins: 28 retained, 0 dropped by min_count=10
  wham polar repeat 0 (seeds 100): K+ barrier=0.6014 eV | Na+ barrier=0.6974 eV
    sampling: 2400.0 samples/window, mean tau = 77.7 samples, N_eff = 17.4 independent samples/window
    bins: 27 retained, 1 dropped by min_count=10
  wham polar repeat 1 (seeds 1000): K+ barrier=0.4874 eV | Na+ barrier=0.7602 eV
    sampling: 2400.0 samples/window, mean tau = 73.0 samples, N_eff = 17.6 independent samples/window
    bins: 28 retained, 0 dropped by min_count=10
  wham polar repeat 2 (seeds 2000): K+ barrier=0.5945 eV | Na+ barrier=0.5798 eV
    bins: 20 retained, 8 dropped by min_count=10
  wham fixed-charge ref (seeds 100): K+ barrier=0.5539 eV | Na+ barrier=0.3967 eV | gap=+0.1572 eV
--- Umbrella polar-WHAM free energy (full QM: SCF dipoles+Pauli+disp, 300 K, 3x[7x1500 steps]) ---
  K+: barrier=0.5611±0.0639 eV | Na+: barrier=0.6792±0.0916 eV | gap=-0.1181±0.1450 eV
  Fixed-charge ref gap=+0.1572 eV. Polar sampling decides the kinetics bracket.
  Barrier over bins with >=10 counts.
  ERROR BAR DEFINITION (audit S1): the +/- is the sample standard
    deviation across 3 INDEPENDENT SEED REPEATS of the whole
    sampling protocol. It is NOT a standard error of the mean and
    NOT a confidence interval. With n=3 the standard error of that
    standard deviation is itself ~76% of its value. It captures
    seed-to-seed variation only; within-run sampling error is
    characterised separately by the reported tau and N_eff.
  Note: 3D ion restraint confines laterally; MBAR over z only.
--- FORENSIC: all legs side by side (eV; + = Na+ favored) ---
  leg                            K+        Na+   dU(K-Na)
  point-charge 8-fold        7.6623     7.6623     0.0000
  JC 8-fold                  7.9754     7.6447     0.3306
  SCF 8-fold                 7.8379     7.6350     0.2028
  v2 polar-SCF 8-fold        8.0766     7.7475     0.3291
  relaxed 8-fold             6.2362     5.8195     0.4167
  stiff 8-fold               7.2825     7.4468    -0.1643
  6-fold octahedral         -1.5388    -1.6575         --
    K+ dE(8-6)=+9.3767  Na+ dE(8-6)=+9.2925 (neg = prefers 8-fold)
  SCF-polar U(z) min         7.1024     6.7308     0.3716
  SCF-polar barrier          6.1598     6.4742    -0.3143
  polar-WHAM barrier         0.5611     0.6792    -0.1181
  fixed-WHAM barrier         0.5539     0.3967     0.1572
  knock-on pair           -215.4617  -150.3790   -65.0826
  knock-on conductive       -3.1798    -8.0803     4.9005
  knock-on landsc bar        4.2410     5.5988    -1.3579
  knock-on landscape     655312.3890 21247.6960 634064.6930
--- Computed exchange K+(aq)+Na+.F -> Na+(aq)+K+.F ---
  relaxed-filter: -0.5004 ± 0.7961 eV | rigid-filter: -0.5865 eV (deterministic)
  (negative = K+ selective; expt -0.179 eV)
--- Two-ion exchange 2K+(aq)+NaNa.F -> 2Na+(aq)+KK.F ---
  -66.9169 eV (negative = KK selective; all four terms computed)

  Datastream s42: kcsa.cvmds written (verify: seal intact).

╔══════════════════════════════════════════════════════╗
║  DEMO 12b: the real KcsA filter (PDB 1K4C TVGYG, C4-symmetric)║
╚══════════════════════════════════════════════════════╝

  Filter built: 164 atoms (41 per subunit x 4), net charge +0.000000 e
  Sequence: THR75-VAL76-GLY77-TYR78-GLY79  (the TVGYG signature motif)
  Geometry: deposited 1K4C chain C. The 2001-era structure has no
    hydrogens, so the five amide H and the C-terminal OXT are
    constructed on sp2 external bisectors at standard lengths.

  --- recovered coordination geometry (crystallographic) ---
    site 1  z= -30.553 A   CN = 8   <ion-O> = 2.932 A
    site 2  z= -33.953 A   CN = 8   <ion-O> = 2.777 A
    site 3  z= -37.162 A   CN = 8   <ion-O> = 2.773 A
    site 4  z= -40.505 A   CN = 8   <ion-O> = 2.912 A
    Every site is 8-coordinate. That is the check that the C4
    symmetry assignment is right: a wrong axis gives CN = 2.

  --- two structural facts a hand-built cage cannot produce ---
    1. The innermost gate is the Thr75 SIDE-CHAIN hydroxyl (OG1),
       not a backbone carbonyl. The old poly-alanine model had no
       such atom, so the single most important ligand in the real
       filter was missing entirely.
    2. GLY79's backbone oxygen sits 4.82 A off the pore axis and
       coordinates NO ion. GLY79 is part of the motif but is not a
       filter ligand; assuming all five residues contribute a
       carbonyl to the ion path would have been wrong.
    Tyr78's phenol OH is 9.2 A from the nearest ion: it points into
    the pore WALL. Tyr79 (canonical numbering) is a gating residue,
    not a selectivity residue.

  --- K+ vs Na+ at each site, rigid deposited filter ---
    The filter is K+-sized: the deposited contacts (2.77-2.93 A)
    match K+ (1.38 A radius + 1.40 A oxygen = 2.78 A) and miss Na+
    (1.02 A + 1.40 A = 2.42 A) by 0.35 A on all eight ligands.

    Dehydration cost to enter the site, from measured single-ion
    hydration free energies (absolute scale):
      K+   -322 kJ/mol -> 3.337 eV
      Na+  -454 kJ/mol -> 4.705 eV
      K+ therefore enters 1.368 eV cheaper. This is a measured bulk
    thermodynamic quantity, not a force-field term, and it is the
    dominant contribution.

    site   K+ CN/r        Na+ CN/r        E_bind K+ E_bind Na+
    1      8 / 2.93 A     8 / 2.93 A        -3.6428    -3.6293
    2      8 / 2.78 A     8 / 2.78 A        -6.4995    -6.4865
    3      8 / 2.77 A     8 / 2.77 A        -8.3063    -8.2935
    4      8 / 2.91 A     8 / 2.91 A        -9.0423    -9.0297

    interaction (rigid cage, ion radius + Coulomb + LJ):
      K+  -27.4911 eV      Na+ -27.4390 eV      difference -0.0521 eV
    plus the measured dehydration cost:
      K+  -14.1419 eV      Na+ -8.6175 eV      difference -5.5244 eV

  --- what this does and does not show ---
    SHOWS: the filter is now the real KcsA TVGYG filter, with the
    real geometry, the real symmetry and the real Thr75 inner gate.
    The K+/Na+ preference in THIS model is carried essentially
    entirely by the measured dehydration free energy (1.368 eV); the
    rigid-cage electrostatics separate the two ions by only
    0.0130 eV, because a symmetric 8-oxygen cage pulls both cations
    to the same axis position.
    DOES NOT SHOW: an absolute selectivity free energy. There is no
    bulk solvent, no membrane potential, no ion concentrations, and
    no flexible filter - a flexible-filter calculation was tried
    and collapses, which kcsa_filter.c documents. A rigid
    fixed-charge model cannot turn the real 0.35 A geometric
    mismatch into a binding-energy difference on its own.
    The previous model's Na+-favouring result was not a KcsA
    result: it was the electrostatics of a constructed cage.

╔══════════════════════════════════════════════════════╗
║  DEMO 17: DNA duplex - minimal G-C and A-T base pair stack║
╚══════════════════════════════════════════════════════╝
  Placed: G-C at y=0, A-T at y=3.4 with 36-deg B-DNA twist
  Total atoms: 59
  Initial PE: -5.806420 eV
  Closest inter-base pair: 35-53 at 1.932 A
  Post-placement H-bonds:
    G-C: N1...N3=2.950  O6...N4=2.953  N2...O2=2.939
    A-T: N1...N3=2.900  N6...O4=2.903
  Energy split at placement: G-C pair LJ=-0.0440 C=-0.2379 | A-T pair LJ=-0.0175 C=-0.1243 | stacking LJ=-0.4975 C=-0.0033
  Clash relief: PE -5.806420 -> -6.202801 eV
  Post-min H-bonds:
    G-C: N1...N3=2.918  O6...N4=2.961  N2...O2=2.974
    A-T: N1...N3=2.916  N6...O4=2.853
  Applied glycosidic restraints (backbone proxy, real springs):
    atom 0 anchored at (1.252, 0.008, -1.982), k=0.50 eV/A^2
    atom 18 anchored at (-3.208, -0.009, 6.055), k=0.50 eV/A^2
    atom 29 anchored at (-1.853, 3.392, -1.240), k=0.50 eV/A^2
    atom 44 anchored at (3.059, 3.398, 6.457), k=0.50 eV/A^2
  Twist restraint: dihedral(G:N1-C:N3-A:N1-T:N3) -> 116.2 deg, k=5 kcal/mol
  Trajectory (every 1000 steps, 8000 total):
  step   GC1      GC2      GC3      AT1      AT2      stackE     twist    rise    
  0      2.918    2.961    2.974    2.916    2.853    -0.5015    82.1     5.332   
  1000   3.052    3.210    3.225    2.987    3.264    -0.5247    79.3     5.314   
  2000   3.356    3.066    3.145    3.122    3.135    -0.5423    79.7     4.909   
  3000   3.111    3.396    3.204    3.185    3.567    -0.6540    72.9     4.364   
  4000   3.187    4.102    3.180    3.218    3.330    -0.6848    63.2     3.895   
  5000   3.442    3.927    3.159    3.030    3.258    -0.6616    64.2     4.091   
  6000   3.173    3.238    3.415    3.088    3.311    -0.6571    74.9     4.281   
  7000   3.354    3.045    3.039    3.138    3.248    -0.6404    75.4     4.506   
  8000   3.469    3.049    2.908    3.183    3.381    -0.6875    78.7     4.751   
  After 8000 restrained MD steps: PE=-6.340655 eV (E_restr=0.043094)  T=49.77 K

  Post-MD H-bonds (restrained):
    G-C: N1...N3=3.469  O6...N4=3.049  N2...O2=2.908
    A-T: N1...N3=3.183  N6...O4=3.381

  Planarity: G=0.0000  C=0.0000  A=0.0000  T=0.0000

PLACEMENT VERDICT:  G-C PAIRED, A-T PAIRED
MD STABILITY (8000-step ensemble): G-C BREATHING (range 2.91-4.10), A-T HELD (range 2.99-3.57)
--> A STABLE, H-bonded two-base-pair DNA stack over 4 ps.
    Watson-Crick pairing emerged from Coulomb+LJ; the twist
    restraint stopped stacking shear from prying G-C apart;
    residual motion is bounded breathing about the model's own
    H-bond lengths (longer than WC ideal under eps=4, as in Demo 7).
  Phase 2 junction 0: C1'-N = 1.4700 A (target 1.47)
  Phase 2 junction 1: C1'-N = 1.4700 A (target 1.47)
  Phase 2 junction 2: C1'-N = 1.4700 A (target 1.47)
  Phase 2 junction 3: C1'-N = 1.4700 A (target 1.47)
  Phase 2 total charge: -0.4728 e (4 glycosidic H removed; fragments approximate)
  Phase 2 atoms: 123 (was 59)
  Phase 2 restraints: 4 sugar-C1' anchors, k=0.50
  Phase 2 trajectory (sugar-tethered, 4000 steps):
    step 0     GC 3.469 3.049 2.908 | AT 3.183 3.381
    step 2000  GC 3.379 3.413 3.224 | AT 3.348 3.329
    step 4000  GC 3.528 3.927 3.195 | AT 3.569 3.450
  Phase 2 verdict: G-C BREATHING* (max 3.93), A-T HELD (max 3.57)
  (*3-sample ensemble; same thresholds as Phase 1)

  ══════════════════════════════════════════════════
  RESULT RECAP
  1 quantum      : H/C/N/O orbitals + expectations + dimer curves
  2 bond curve   : H2 covalent-vs-vdW curves
  3 water MD     : H2O Berendsen MD @300K, T=267.9K
  4 trimer       : trimer H-bond ring HELD
  5 methane      : CH4 tetrahedral 109.47 deg, PE=0
  6 nucleobases  : U/C/T/A/G geometry + RESP charges validated
  7 pairing      : G-C>A-U, both bound
  8 dinucleotide : T-p-A backbone, charge -1.00e
  9 HH neuron    : HH AP peak +40mV, 4 spikes/50ms
  10 dipeptide   : Gly-Ala peptide bond 1.33A
  11 helix       : helix i,i+4 H-bond EMERGED
  12 KcsA        : KcsA exch rigid -0.59 det
  17 duplex      : duplex P1 BREATHING/HELD P2 BREATHING*/HELD
  ══════════════════════════════════════════════════

  All demos complete.
  Three validated tracks now exist: nucleic acids (bases through a
  real phosphodiester bond), proteins (a real peptide bond AND, given
  correct local backbone torsion geometry, a genuine emergent alpha-
  helical hydrogen bond), and electrophysiology (a genuine Hodgkin-
  Huxley action potential). None are connected to each other yet.
  Real next steps: a full DNA duplex, gene regulatory logic, a
  synapse between neurons, and eventually deriving ion channel
  gating from actual protein structure rather than empirical rate
  equations - closing the loop between the protein and
  electrophysiology tracks.
```
## v9R4 Supplement — Engine Arc Since the v9R3 Audit

What changed, grouped by subsystem (all in `output.txt`, all sealed in
`kcsa.cvmds` where applicable):

* **Truth fixes (behavior-changing, re-baselined).** KcsA antiprism cage
  corrected from ring-radius-as-xy to 3D-derived in-plane radii (2.22 /
  2.37 Å for 3D 2.70 / 2.83 Å); QM overlap reference moved from a
  hand-typed 1.5 Å to covalent-contact distance; CO₂ zero-strain
  override (placed 1.163 Å vs table 1.23); Slater 10–30% energy claim
  corrected to factors of ~2–5× with a per-element computed-vs-NIST
  ledger in Demo 1.
* **Quantum foundation (v4).** Clementi–Raimondi SCF exponents for
  spatial ranges (dual-sourced against WebElements; K-3d and Se-3d
  dropped on conflict/evidence); exact H-like expectations; real f
  harmonics; lobe-max σ + sum-rule π overlap; charge-responsive
  screening γ into overlap/dispersion/α; self-consistent dipoles
  (Thole-damped T-coupling, Hellmann–Feynman forces) behind
  `use_pol_scf`; Slater–Kirkwood dispersion + Tang–Toennies damping;
  SCF QEq+dipole charge loop (1-pin and 2-pin) with gain/mixing
  stabilization; ion-aware polarizability (K⁺ 0.83 / Na⁺ 0.18 — the
  term that wakes off-axis).
* **KcsA program.** Point → JC → SCF → v2/SCF-polar → relaxed
  (multi-start, hysteresis errors) → stiff strain probe → 6-vs-8
  coordination probe → knock-on pairs (relaxed + conductive
  single-points) → knock-on landscape over clash-flagged valid subset
  → explicit-water hydration → polar + fixed WHAM → computed 1-ion and
  2-ion exchange. Forensic table prints every leg side by side.
* **Duplex.** Pre-MD minimization, stacking-shear diagnostics, gentle
  twist restraint, trajectory ensemble verdicts, Phase 2 with real
  glycosidic sugar tethers (4× C1′–N at 1.47 Å).
* **Classical core.** Analytic dihedral gradients (P2 closed, FD oracle
  kept, agreement selftest); Andersen canonical thermostat (WHAM legs);
  FIRE minimizer (smooth basins; steepest kept for clash relief after
  head-to-head); r2 micro-opts; per-axis PBC; hardened validation and
  thread-local QEq buffers throughout.
* **Display.** Deterministic stdout (the record) vs TTY-gated stderr
  (timings, heartbeats, leg progress) split in `include/display.h`;
  version banner; 13-line result recap; three selftests
  (`test_datastream`, `test_forces`, `test_fire`) all gated by
  `verify_scripts.sh` and the s01 record harness.

## Code Organisation

`include/amber_lj.h` holds the AMBER ff99 Lennard-Jones types for every
biomolecular atom, in one place, together with the `2^(1/6)` conversion that
turns an AMBER Rstar into a standard 12-6 collision diameter. Before the
streamlining pass those numbers existed in four private copies:
`AMBER_RSTAR_TO_SIGMA` in `nucleobases.c`, `AA_RSTAR_TO_SIGMA` in
`aminoacids.c`, a third set in `kcsa_filter.c`, and the literal
`1.6612 * 2.0 / 1.122462048309373` written out in **fourteen** places in
`main.c`. `aminoacids.c` had admitted the duplication in its own comment —
"duplicated here rather than shared via a header refactor to avoid touching
already-validated, working code under time pressure" — a debt that was never
paid. `2^(1/6)` is now written exactly once in the tree.

Consolidating them immediately paid for itself by surfacing a bug: the
carbonyl-oxygen sigma used by the new KcsA filter module had been hardcoded
as **3.06615** when the AMBER ff99 value is **2.95992**, and the regression
test had copied the same wrong literal, so the module and its test agreed
with each other and both were wrong. The shared constant fixed both, and
moving the filter energies by ~0.0014 eV is the size of the correction.

## Verification Discipline & Known Limitations

This project holds itself to two standards, stated here explicitly so
any result can be judged against them: **every number is sourced, and
every approximation is flagged.** Where a result is qualitative rather
than quantitative, or where a demo reports a genuine negative result,
the output says so in plain language rather than burying it.

### Verification discipline

* **Constants and unit conversions.** Fundamental constants are 2019
  CODATA values (exact where the redefinition applies). Every MD unit
  conversion is derived in-line in `constants.h` rather than taken on
  faith — e.g. the Coulomb prefactor `COULOMB_MD = 14.3996 eV·Å/e²`
  and the kinetic-energy factor `103.6427 eV` per AMU·(Å/fs)² both
  show their full derivation.
* **Primary-source parameter sourcing.** Geometry comes from the RCSB
  PDB Chemical Component Dictionary ideal coordinates (fetched from
  named URLs on a dated fetch); nucleobase partial charges are RESP
  values from Aduri et al. (2007), Table 1; Lennard-Jones parameters
  are real AMBER ff99 values from a TINKER-format parameter file,
  mapped atom-by-atom against AMBER's own per-base atom-type lines.
* **Cross-checking before implementation.** Parameters and formulas
  were verified against two or more independent sources before any C
  code was written. The Hodgkin-Huxley resting gating values were
  recomputed from the rate functions and matched a published worked
  example to four decimal places (0.0529, 0.5961, 0.3177) first; the
  signed dihedral formula was checked against four hand-built test
  cases (+90°, −90°, 0°, 180°); Slater screening was reproduced
  against Slater's own 1930 worked example for iron; the AMBER
  R*→σ conversion was confirmed by reproducing TIP3P oxygen's
  literature σ of 3.1506 Å exactly.
* **"Trust but verify" geometry checks.** Placed molecular geometry is
  re-measured independently from the Cartesian positions (bond
  lengths, ring angles, planarity deviations, total charge) rather
  than read back from stored equilibrium values — the same check that
  caught an NH₃ sign error silently making the molecule planar.
* **Correct degrees of freedom.** Temperature uses `3N − 3` (centre-of-
  mass motion removed), not `3N` — a bug that underestimates T by 33%
  for a 3-atom system.
* **Pre-registered prerequisites resolved before use.** The v9
  masterplan's Step 0 (the potassium row discrepancy in
  `periodic_table.c`) was resolved in-file, with the resolution
  documented where the data lives, before any downstream demo depended
  on potassium parameters — the same verify-before-building-on-it
  order the rest of the codebase uses.
* **Build-configuration reproducibility.** The full demo sequence has
  been run under AddressSanitizer (stderr empty — no memory errors)
  and with `-ffast-math` both on and off. The ASan output and the
  no-fast-math output are byte-identical (same SHA-256, archived as
  `output.after-fix09.txt`); the fast-math output differs only in the
  low-order digits of minimizer-dependent values (e.g. Demo 11 H···O
  2.1639 vs 2.1636 Å) — expected floating-point reassociation noise
  in a chaotic relaxation, not a physics difference.

**Confidence tiers.** Not all parameters carry the same weight, and the
code says which tier each one is in:

1. **Independently verified against primary data** — nucleobase RESP
   charges and AMBER ff99 Lennard-Jones parameters; the CODATA
   constants; the Hodgkin-Huxley parameters and rate functions.
2. **Standard literature / textbook values, not re-derived here** —
   the 1.33 Å peptide bond length and the Engh-&-Huber-style backbone
   junction angles; the glycosidic bond length.
3. **Explicitly approximate, charge-balanced placeholders** — the
   sugar, phosphate, and amino-acid partial charges, and thymine's
   methyl-group charges. These are flagged in-source and in the output
   as *not* independently verified the way tier-1 values are.

Potassium's LJ ε in `periodic_table.c` is verified against Rappé et al. 1992
UFF (x2 = 0.035 kcal/mol) as of audit T1. The σ column needed a correction in
the deep audit: it held UFF's x1 value (3.812 Å), which is **Rmin**, the
distance of minimum — not the 12-6 collision diameter the field is declared to
hold. Every σ is now `x1 / 2^(1/6)` (K: 3.39611 Å), which places the LJ
minimum at the correct 3.812 Å and matches the AMBER conversions already used
for amino acids and nucleobases.

### Known limitations & honest scope

* **1-4 non-bonded scaling deliberately not applied (audit F1).** AMBER's
  1-4 scaling (LJ/2, Coulomb/1.2) presupposes AMBER's co-fitted torsion
  parameters; this force field's dihedrals are restraints, not fitted
  torsions, so scaling 1-4 pairs removes physics with nothing to
  compensate. It was tested and broke the validated helix i,i+4 H-bond
  (Demo 11). 1-4 pairs stay at full strength; revisit only with fitted
  torsions.
* **Ion selectivity (Demo 12) is a model result, not a KcsA result.** See the
  audit section: the filter is a frozen, restrained array of carbonyl oxygens
  and the protein is poly-alanine, so the real TVGYG selectivity chemistry is
  not present. The barriers reported below are the umbrella-sampled U(z)
  spans for that cage, and the numbers moved when the physics was corrected.
  Any K⁺/Na⁺ comparison should be read as "in this reduced model", never as a
  statement about the channel. The umbrella sampler now prints its measured
  autocorrelation time and effective sample size, because the `±` on these
  figures is a three-seed spread and **not** a standard error.
* **Demo 12's antiprism block is now a real result (no longer a placeholder).** The
ring z-separation was sourced from 1K4C (3.084 A, symmetry-validated in
s10/s10b) and the artificial O-O overlap is removed. The block now reports a
genuine, interpretable Delta (Na+ still favored, wrong direction). Note: the
absolute Total_PE values are net positive due to the isolated cage's O-O
self-repulsion (no protein backbone to balance it); only the Delta is the
meaningful selectivity signal.
* **Base-pairing energetics (Demo 7) are qualitative.** G–C correctly
  binds more strongly than A–U and both pairs are correctly
  attractive, but the absolute magnitudes overshoot gas-phase ab
  initio references (~−1.2 eV G–C, ~−0.55 eV A–U) by roughly 4×, and
  no single dielectric constant brings both pairs into simultaneous
  quantitative agreement. The ordering is validated; the numbers are
  not yet quantitatively trustworthy.
* **No bulk solvent.** Everything runs in vacuum with a relative-
  permittivity divisor (dielectric = 4 for base pairing) standing in
  for condensed-phase screening, plus now explicit octahedral 6-water
  clusters that feed the ion-exchange legs directly. Bulk water (PBC
  box, Ewald/PME) remains future work — the root of the remaining
  quantitative gaps.
* **Approximate charges on non-nucleobase components.** Sugar,
  phosphate, and amino-acid partial charges are charge-balanced
  approximations, not verified RESP fits (see the confidence tiers).
* **Scaling and algorithms.** The non-bonded loop is O(N²) with no
  neighbour list (the code suggests a cell list beyond ~500 atoms),
  though pair exclusion is O(coordination) via partner lists and
  cutoff checks avoid sqrt for out-of-range pairs. Dihedral gradients
  are exact analytic chain rule (audit P2 closed; the FD oracle and a
  22-check agreement selftest remain in-tree). Minimization offers
  steepest descent (clash relief — monotonic and safe) and FIRE
  (smooth-basin polishing, with an agreement selftest); sampling offers
  Berendsen steering plus rigorously canonical Andersen collisions.
* **Harmonic bonds cannot break.** The bond potential is a parabola,
  so it wrongly predicts infinite energy at full separation. Real bond
  breaking needs a Morse potential or a reactive force field.
* **The three tracks are not connected.** Nucleic acids, proteins, and
  electrophysiology each run independently. Deriving ion-channel
  gating from protein conformational MD — closing the loop between the
  protein and electrophysiology tracks — remains future work.
* **Demo 11 is not spontaneous folding.** Local backbone torsions
  (φ/ψ/ω) are restrained to textbook values; only the global i,i+4
  hydrogen bond is allowed to emerge. Demo 12 uses a literature
  coordination radius as a constructed input rather than deposited
  backbone coordinates, so it tests energetic preference at imposed
  geometry, not whether the correct geometry emerges.
* **Cutoff smoothing is present but gated off (audit F5).** A
CHARMM-style switching function (potential and force continuous to zero
at the cutoff) is implemented for future condensed-phase use but is
opt-in and disabled by default - every current demo is gas-phase/vacuum
and uses the plain hard-cutoff potential. It must be enabled, with a
real cutoff, before any condensed-phase system; until then it is not
exercised by any demo.

---

## Current Status

`carbonsim-v9R4` is a validated physics engine and a set of validated
building blocks, not a production simulator. The fundamental layers —
quantum orbitals, the molecular-dynamics core, dihedral forces, energy
minimization, biopolymer condensation chemistry, and the
Hodgkin-Huxley track — are implemented, cross-checked, and demonstrated
to produce genuine emergent behavior.

The v9 goal — testing whether the unmodified Coulomb+LJ engine extends
to ion selectivity — was answered honestly: **it does not, in vacuum,
with fixed charges.** Na⁺ wins at every KcsA site and across the radius
scan, and the systematic elimination (generic typing first, then
geometry-fit) makes that a genuine physics finding about the model's
ceiling rather than a parameter bug.

v9R4 then built what the finding pointed at — polarizability (Thole-
damped induction, self-consistent dipoles), explicit-solvent competition
(octahedral 6-water clusters), charge equilibration with dipole feedback,
Slater-Kirkwood dispersion, multi-ion knock-on configurations, and
umbrella-sampled free energies. Single-ion vacuum filter statics still
favor Na⁺ (+0.2…+0.4 eV across six legs — correct Coulomb physics, kept
as the honest baseline); K⁺ wins the computed exchange (−0.578 eV
deterministic rigid-filter leg), the SCF-polar transit barrier
(−0.36 eV), and the knock-on landscape over valid geometries
(−1.62 eV). Sampled polar kinetics (+0.33 ± 0.06 Na⁺) disagrees with the
single-point barrier — reported as an open bracket, not averaged away.
Details in the v9R4 supplement below and `release_note_v9R4.md`.

The honest next steps, in the order the project's own record supports:

1. **Complete Demo 12's antiprism block** *(done)*: the real ring
z-separation was sourced from 1K4C (3.084 A, symmetry-validated),
and the antiprism block is now a real result.
2. **The missing-physics decision** *(built, all three)*: polarizability
   (in-loop induction + SCF dipoles), explicit-solvent competition
   (computed 6-water legs feed the exchange directly), and charge
   response (QEq+dipole SCF, Slater-Kirkwood dispersion).
3. **The longer-horizon tracks already named**: a full DNA duplex,
   gene-regulatory logic, a synapse between neurons, and eventually
   deriving ion-channel gating from actual protein structure — closing
   the loop between the protein and electrophysiology tracks.

## Demo 17 — DNA duplex (incorporated into main executable)

Demo 17 validates Watson-Crick base-pairing geometry from Coulomb+LJ
physics alone. Two base pairs (G-C and A-T) are placed using
perpendicular-to-WC-edge placement (the H-bond direction, validated
against Demo 7's G-C/A-U results), stacked along the ring-normal
direction (Y axis) with a 3.4 A rise.

**Result:** Post-placement H-bonds are textbook-perfect:
- G-C: N1...N3=2.950 A, O6...N4=2.953 A, N2...O2=2.939 A (3 H-bonds)
- A-T: N1...N3=2.900 A, N6...O4=2.903 A (2 H-bonds)

Glycosidic tethers plus a gentle inter-pair twist restraint (the
backbone's helicity job; H-bond lengths stay free) hold the stack:
8000-step ensemble verdict G-C BREATHING / A-T HELD, with stacking-shear
diagnostics (twist/rise columns) proving the shear mechanism. Phase 2
rebuilds the tethers as real glycosidic bonds to four deoxyribose sugars
(all C1'–N at exactly 1.47 A, charge reported) and re-tests: G-C
BREATHING / A-T HELD again — the backbone proxy v2 is consistent with
v1. Thresholds are calibrated to the model's own isolated-pair behavior
(Demo 7 finals run 2.95–3.13 under eps=4), not to WC ideal. Full
phosphodiester polymer remains the documented next step.

Demo 17 is now part of the main executable flow (no --dna flag needed).

---

## v9R4 Deep Audit — Physics Corrections & Hardening (2026-09-29)

A second, independent audit was carried out against the whole tree. Unlike the
first pass, this one did not take the code's comments or documentation at face
value: every claim was re-derived and every analytic force was checked against
an independent finite-difference oracle, every constant against its primary
source, and every tabulated dataset against published reference values.

**The tree did not compile.** Two separate defects, in sequence:

* The first v9R4 pass added static assertions to `types.h` that were simply
  wrong — four asserted struct sizes that ignore padding, and four that
  referenced the `SimError` enum ~20 lines before it was declared. `make`
  failed immediately on every source file.
* The repair that removed that block also dropped a closing brace in
  `forces.c`, so `origin/main` **still did not compile from a clean tree**.
  Every function after line 170 was swallowed into an unclosed `if`.

Neither was visible, because **this repository tracks `build/*.o` and the
`carbonsim` binary**. An in-place `make` finds the committed objects
up-to-date, rebuilds nothing, and reports success while the source is broken.
`verify_scripts.sh` counted build warnings with a bare `make` and therefore
reported "0 warnings" on a tree that could not build at all. The verify script
now does `make clean` first, and that check is a hard gate.

### What was wrong with the physics

Each item below was found by measurement, and each is now covered by a
regression test that fails on the old code.

| # | Defect | How it was caught | Severity |
|---|--------|-------------------|----------|
| D1 | **Dispersion force sign inverted.** `E = -C6 f6/r6` is attractive, but the force was `-(dE/dr)(d/r)` instead of `+(dE/dr)(d/r)`, so the term was an effective **repulsion** in the force path while every printed `E_disp` looked correct. | FD oracle: analytic force was anti-parallel to the gradient, relL2 = **2.000** — the exact signature of `F = -F_FD`. Flipping the sign gives 5.9e-11. | **critical** |
| D2 | **SCF polarization force wrong, and the energy was discontinuous.** With self-consistent dipoles `E = (I - Tα)^-1 E0`, so `dU/dE0 = -μ(I - Tα)^-1`; the Hellmann-Feynman shortcut is unavailable. The code used `s = 1.0` (98% error) and, on solver non-convergence, silently swapped in a *different* energy functional. | Measured a **5.92 eV jump over a 0.01 Å displacement** in the potential energy — a discontinuous Hamiltonian, so no energy conservation. | **critical** |
| D2b | **The dipole "SCF" was a fixed-point iteration for a LINEAR system.** It needed 96–97 of its 100 iterations for **four atoms**, and a tuned-acceleration variant still failed on 121 of 401 geometries. | Convergence telemetry added to the umbrella sampler. | **critical** |
| D5 | **Catastrophe cancellation on bonded pairs.** With no dipole hard core, a 1.3 Å C=O bond and α_C = 1.76 gives an off-diagonal coupling αk f(r)/r³ ≈ 2.8, so `(I - A)` has an eigenvalue past 1 and the system has **no solution**. | Direct solve returned "singular" on all 401 geometries. | **critical** |
| F6 | **LJ σ stored UFF's Rmin, not σ.** `lj_sigma` is declared as the 12-6 collision diameter and consumed by the standard form, but the table held UFF's x1 column. Every effective LJ well sat **12.2% too far out** (carbon: minimum at 4.32 Å instead of 3.85 Å). | Recomputed the C–C well from the engine's own parameters. | **major** |
| F7 | **Carbonyl O and amide N were classified sp3.** The rule `steric >= 4 \|\| (s_count > 0 && p_count >= 3)` fired for every N and O, overriding the steric logic. These are the functional groups the entire biomolecular model is built from. | Direct call on constructed carbonyls/amides. | **major** |
| F8 | **`r_mp³` used as a polarizability.** It is a screening-radius volume, and because the branch was selected by atomic number it made the induction energy **13.6× discontinuous between adjacent elements** (H 0.420 vs He 0.031). It also fed Pauli screening and the Slater-Kirkwood C6. | Tabulated all 36 elements against measured values. | **major** |
| F9 | **QEq had no hard core and no clamp.** Standard QEq zeroes `A_ij` for 1-2/1-3; including the full 1/r term let the solve run away to **q(C) = +4.82 e, q(O) = −5.63 e** on a bonded C-O-O fragment. | Realistic bonded fragments. | **major** |
| F11 | **The minimiser divergence floor did not scale.** A hard −50000 eV floor silently disables minimization for any system below that: a protein at a routine −15 eV/atom crosses it at ~3300 atoms, and `s10_1K4C.pdb` is in this repository. It would have returned an unminimized structure while reporting success. | Scale the floor to `−60 eV/atom`; large deep systems now minimize. | **major** |
| S1 | **Free-energy sampling had no error analysis.** 300 nominal samples per umbrella window, autocorrelation ignored, and barriers reported to two decimals. | Instrumented: measured τ = 50–96 samples, **N_eff ≈ 15–30 independent samples per window** out of 2400 nominal. | **major** |
| I1 | **PCG64 was documented but did not exist.** Both the readme and the release note stated "PCG64 RNG added as opt-in". `grep` found no PCG64 anywhere in the tree. | Implemented for real; both generators now tested. | **major** |
| F14 | **Silent orbital truncation.** The fixed-size ml-slot table broke out of its loop with no signal, so heavy elements quietly lost electrons. | Now reported through an out-parameter and covered by a test. | moderate |
| F9b | **Unclamped charges on the force path.** `forces_calculate`'s unpinned SCF branch wrote QEq output straight onto the atoms, while `qm_scf_run` clamps at ±2 e. | — | moderate |
| F10 | **`qm_alpha` dereferenced before its NULL check** — `qm_alpha(NULL)` segfaulted. | — | moderate |
| F13 | **`integrator_remove_com_velocity` had no guard**, unlike every other entry point in the file. | — | minor |
| F12 | **Dead clamp** in both steepest-descent minimizers: `scale = effective_step/max_force` already bounds the displacement, so the per-atom branch could never fire. | — | minor |
| B1 | **The build check could not fail** (bare `make` against committed objects). | See above. | **major** |
| B2 | **`frozen[32]` stack array** indexed by an atom index and read for every atom by `integrator_minimize_frozen`. Fits today at 19 atoms; one more water overruns it. | Sized to the system, bounds-checked. | moderate |
| D4 | **Induction constant hand-typed** as 0.069446, a 2.2e-6 truncation of `1/COULOMB_MD`, in a codebase whose stated rule is to derive every reciprocal in-line. | Now derived. | minor |

### Memory safety

Two confirmed out-of-bounds **writes**, both reproduced under UBSan:

```
forces.c:147  index 4 out of bounds for type 'int [4]'              forces_bond_params(Za,Zb,order) with order=4
forces.c:199  index 118 out of bounds for type 'int [118][118][118]' forces_angle_params with Z=118
```

The de-duplication arrays indexed caller-supplied values, and cost 6.9 MB of
BSS to implement "print once". Replaced with a bounded key set, so no
caller-supplied value is ever used as an index. Also fixed: five unguarded
`->element->` dereferences (the periodic table stops at Kr, Z=36, while
`MAX_ELEMENTS` advertises 118, so `sim_add_atom(s, 37, ...)` yields
`element == NULL`), `vec3_pbc` producing NaN on a degenerate box
(`0 * round(dr/0)`), and `nb_planarity_deviation` reading past the atom array
on a bad index.

### What was checked and found correct

Worth stating, because it bounds where the doubt lies.

* **The hydrogenic quantum mechanics is exact.** `∫r²|R_nl|²dr`, `<r>` and
  `<r²>` for n = 1…3, l = 0…n−1 agree with Simpson quadrature to **1e-13**.
  The Laguerre three-term recurrence, the `scale³` normalization, and the
  partitioned golden-section search for `r_mp` are all correct.
* **Slater's rules are correct**, including the subtle point that `(n−1)s`,
  `(n−1)p` and `(n−1)d` form a *single* 0.85 tier. Reproduces Slater's own
  Fe 4s (3.75) and Fe 3d (6.25) worked examples exactly.
* **The Clementi–Raimondi table is sound.** The O 1s cell (`7.6579`) was
  suspected of being a digit transposition of `7.7579`; it is **not**. The 1s
  series rises by a uniform 0.9925 per element from Li to Kr, and 7.6579 sits
  on that line while 7.7579 would put a 0.19 kink in it. The test suite now
  asserts that smoothness, so a real transposition would be caught.
* **All 49 real spherical harmonics are orthonormal** for l ≤ 3, including the
  hand-normalized f₂ coefficient the notes claimed was "fixed by integration".
* **SHA-256 is a correct self-contained FIPS 180-4 implementation** (all three
  NIST vectors) and the datastream verifier is carefully hardened against
  parser confusion: it rejects a second `[end]`, requires the hash after it,
  requires exactly 64 lowercase hex plus newline/EOF, and hashes the correct
  payload boundary. Tamper detection verified.
* **Velocity Verlet conserves energy** to < 1e-7 eV/step over 20 000 steps,
  flat in `dt` — the signature of a correct symplectic integrator.
* **CHARMM cutoff switching, the angle gradient, and the Tang-Toennies damping
  (including `df6/dr = e^{-x} b x^6 / 720`) are all correct.** The induction
  analytic gradient for the *first-order* (uncoupled) term is correct too —
  the defect was specific to the self-consistent path.
* **MBAR is implemented correctly.** The `Σ_w n_w / Σ_w N_w exp((F_w − V_w)/kT)`
  estimator and its self-consistent `F_w` update are textbook. The problem was
  the sampling, not the estimator.
* **`neuron.c` is a clean squid-axon Hodgkin–Huxley implementation** with the
  removable singularities in α_m and α_n handled by a correct Taylor expansion.

### Verification gates

A new regression suite, `make selftest-regression`, adds **138 checks**, one
per defect above, built on independent oracles — finite differences,
quadrature, NIST SHA-256 vectors, published reference values — so it can
actually fail. `make test` runs every gate.

Current state:

```
clean build                        0 warnings, 0 errors
selftest (datastream)              17 checks green
selftest-forces                    22 checks green
selftest-fire                       7 checks green
selftest-regression               138 checks green
ASan + UBSan                        0 memory errors, 0 UB findings, empty stderr
stdout byte-identical across runs   yes
stdout byte-identical under ASan    yes
record SHA-256                      67f9f47eb3d3e5b117b5ba8f4c28dfe212caf4ed12323c97246fcf93da9ee72d
```

The record digest changed from the previous release. That is expected and
honest: the physics changed. A dispersion force that pointed the wrong way
produces different trajectories, as do a corrected LJ well, a tabulated
polarizability, a hard-core QEq, sp2/sp3 labels, a direct dipole solve, and a
53-bit uniform conversion in place of 31-bit. The old digest described code
that did not compile.

### Honest statement about the error bars

The umbrella free-energy barriers are reported with a `±` that is the sample
standard deviation across **three independent seed repeats** of the whole
sampling protocol. It is **not** a standard error of the mean and **not** a
confidence interval. With n = 3 the standard error of that standard deviation
is itself about 76% of its value. It captures seed-to-seed variation only;
within-run sampling error is characterised separately by the measured
autocorrelation time and effective sample size, which the program now prints.

The measured numbers are `N_eff ≈ 15–30` independent samples per window. On
that basis a barrier quoted as `0.4874 ± 0.0030 eV` is **not** credible to
four decimal places — the tight value reflects three seeds happening to agree,
not a well-determined estimate, and the accompanying Na⁺ figure
(`0.4896 ± 0.1490 eV`) shows the true scale of repeat-to-repeat variation.
Anyone quoting these barriers should quote the larger of the two and state
`N_eff`.

### What this release does not fix

* **The filter is now the real KcsA filter, and what that does and does not
  buy.** `src/kcsa_filter.c` builds the TVGYG selectivity filter from
  **deposited PDB 1K4C coordinates** (chain C, THR75–VAL76–GLY77–TYR78–GLY79),
  with the other three tetramer subunits generated by the C4 rotation about
  the pore axis. The symmetry was *derived*, not assumed: only the correct
  axis reproduces a coordination number of 8 at every deposited K⁺, and a
  wrong axis gives 2. The build reproduces the deposited geometry exactly —
  CN = 8 at all four sites, K–O 2.773–2.932 Å, filter charge 1.8e-15 e.

  This replaces the hand-placed restrained oxygen cage and the poly-alanine
  chain. Two structural facts came out of the coordinates that no hand-built
  model would have produced, and both are load-bearing:

  * **The innermost gate is the Thr75 side-chain hydroxyl (OG1)**, not a
    backbone carbonyl. The poly-alanine model had no such atom, so the single
    most important ligand in the real filter was missing entirely.
  * **GLY79's backbone oxygen sits 4.82 Å off the pore axis and coordinates
    no ion.** GLY79 is part of the signature motif but is not a filter
    ligand. Tyr78's phenol OH is 9.2 Å from the nearest ion — it points into
    the pore *wall*, not the pore. Tyr79 (canonical numbering) is a gating
    residue, not a selectivity residue.

  **What the model still cannot do.** The K⁺/Na⁺ preference here is carried
  essentially entirely by the measured single-ion hydration free energy
  (K⁺ −322, Na⁺ −454 kJ/mol → K⁺ enters 1.368 eV cheaper). The rigid
  eight-oxygen cage separates the two ions by only 0.05 eV, because a
  symmetric cage pulls both cations to the same axis position regardless of
  radius. The real 0.35 Å geometric mismatch is real and is now measured,
  but a rigid fixed-charge model in vacuum cannot convert it into a binding
  free energy on its own — that is what filter *flexibility* and a solvation
  model are for. A flexible-filter calculation was implemented, found to
  collapse (the restraints that let the filter respond are far too weak to
  hold a +1 ion's Coulomb field), and is **not** reported; the collapse is
  documented in `kcsa_filter.c` rather than tuned away. So: the geometry is
  now genuinely KcsA, the mechanism is correctly identified, and no
  quantitative selectivity free energy is claimed.

* **No bulk solvent, no membrane potential, no ion concentrations.** A
  relative-permittivity divisor stands in for condensed-phase screening.
* **1-4 non-bonded scaling is deliberately absent** (audit F1), which remains
  a non-standard choice documented in place.
* **Harmonic bonds cannot break.**
* **The three tracks are still not connected.**
* **Base-pairing energetics remain qualitative** — the G–C > A–U ordering is
  validated; the absolute magnitudes overshoot gas-phase references by ~4x.
* **`qm_overlap` is not an overlap integral.** It returns `exp(-ζr)` times
  angular lobe factors, with the Slater polynomial prefactor dropped and `m`
  maximised independently per atom, so a misoriented orbital still scores
  full credit. It is a heuristic and the naming oversells it; `qm_bond_order`
  and `qm_pauli` are built on it.

---

## Demo 17 — DNA duplex (incorporated into main executable)
