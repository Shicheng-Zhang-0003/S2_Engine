# S2 Engine — Nano Chemistry & Biological Simulator

## Note: The new naming strategy across the board will be vxR123, not vxA123. This scheme began with the v9R3 designation; the release current as of this sync is v9R4. New releases after this point in all other projects will shift A to R, Alpha will be removed since it's just too extraneous a word.

## In addition, S2 will not follow the R1 R2 R3 developmental scheme of Lancius and Mathlib. Instead, the v9 to v10 period for S2 will be marked by however many RC releases are necessary to achieve a OpenWorm similar display for visualising bonds occuring.

*Documentation synced to the current tree state post-v9R4 (hygiene closure: quantum.c normalization fix, periodic_table.c explicit initializers, makefile suppression removal; readme designation and references updated), 2026-08-23. The embedded output block is unchanged and byte-identical to the record (SHA-256 875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935).*

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
the full 12-demo sequence was regenerated and re-asserted after every
behavior-affecting change. Full ledger: the "v9R3 - The Audit &
Correctness Release" section and `release_note_v9R3.md`.

**Open items carried forward to v9R4+** (direction unchanged by the
audit; this list is the unambiguous starting point for the next
session):
1. **Complete Demo 12's antiprism block**: source the real ring
   z-separation from 1K4C so the placeholder becomes a real result or
   a real exclusion. Until then its numbers must not be cited as
   selectivity energetics.
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
  block is a placeholder, not a result (see Known Limitations).

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

1. **Complete Demo 12's antiprism block**: source the real ring
   z-separation from 1K4C so the placeholder becomes either a real result
   or a real exclusion. Until then its numbers must not be cited as
   selectivity energetics.
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
SHA-256 `875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935`)
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

**Record (M2, D1).** The full 12-demo sequence was regenerated post-audit
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
     the real 45° offset), which is explicitly flagged as a
     *placeholder*: the real ring z-separation was not sourced this
     pass, both rings sit at z=0, and the resulting O–O overlap
     dominates the energies.

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
build compiles with `-O3 -Wall -Wextra -std=c11 -march=native`.

**Regenerating the recorded output**: the block reproduced below is the
demo binary's stdout. After any code change, regenerate it so the
"verbatim record" claim stays provable:

```bash
./carbonsim > output.txt
```

and re-sync this readme's embedded block from that file. (Provenance
note: the archived `output.after-fix09.txt` in this tree predates the
Demo 12 caveat strengthening and the banner-line update now in
`src/main.c`; the block below matches the current source.)

**Debug build**: the makefile carries a commented alternate `CFLAGS` line enabling AddressSanitizer and UndefinedBehaviorSanitizer. As of audit fix B1, the link rule passes `$(CFLAGS)`, so uncommenting that line produces a working sanitised build through make. Verified this release: the ASan build ran the full 12-demo suite with empty stderr (zero memory errors), and its output is byte-for-byte identical to the normal build (both SHA-256 `875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935`).

---

## Simulation Execution Output (v9R4)

Below is the complete, unedited output of the `carbonsim-v9R4`
demonstration run — all twelve demos, from the quantum orbital tables
through the Hodgkin-Huxley action potential and the KcsA selectivity
filter. It is reproduced verbatim from `output.txt`.

Note the output itself here: Demo 7 flags that its base-pairing
magnitudes are qualitatively but not yet quantitatively trustworthy;
Demo 11 states explicitly that it is not a spontaneous-folding test;
and Demo 12 reports a genuine negative result with a systematic
diagnosis, and labels its antiprism block a placeholder rather than a
result. These caveats are part of the record, not footnotes.

```text

  ╔═══════════════════════════════════════════════════════╗
  ║       CARBON VM — CHEMISTRY SIMULATOR                 ║
  ║       From subatomic to molecular dynamics            ║
  ╚═══════════════════════════════════════════════════════╝

  Unit system: Length=Å  Time=fs  Energy=eV  Mass=AMU
  Physical constants: 2019 CODATA  |  LJ: UFF defaults + AMBER ff99 overrides  |  Bonds: placed-geometry r0, generic spectroscopic k (audit F2)


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
  LJ σ              : 2.8860 Å
  Config            : 1s1
══════════════════════════════════════════
  Orbital table for H (Z=1)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -13.6058     1         
  Valence orbital: 1s  Z_eff=1.000  r_mp=0.529 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å |##|:..                                  10.0 Å

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
  LJ σ              : 3.8510 Å
  Config            : 1s2 2s2 2p2
══════════════════════════════════════════
  Orbital table for C (Z=6)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -442.0524    2         
  2s(+0)   2      0      0      -35.9278     2         
  2p(-1)   2      1      -1     -35.9278     1         
  2p(+0)   2      1      0      -35.9278     1         
  Valence orbital: 2s  Z_eff=3.250  r_mp=0.853 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å  :##|:.                                  10.0 Å

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
  LJ σ              : 3.6600 Å
  Config            : 1s2 2s2 2p3
══════════════════════════════════════════
  Orbital table for N (Z=7)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -610.7644    2         
  2s(+0)   2      0      0      -51.7361     2         
  2p(-1)   2      1      -1     -51.7361     1         
  2p(+0)   2      1      0      -51.7361     1         
  2p(+1)   2      1      1      -51.7361     1         
  Valence orbital: 2s  Z_eff=3.900  r_mp=0.710 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å  |#|..                                   10.0 Å

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
  LJ σ              : 3.5000 Å
  Config            : 1s2 2s2 2p4
══════════════════════════════════════════
  Orbital table for O (Z=8)
  Orbital  n      l      ml     Energy(eV)   Occ       
  -------  --     --     --     ----------   ---       
  1s(+0)   1      0      0      -806.6879    2         
  2s(+0)   2      0      0      -70.4185     2         
  2p(-1)   2      1      -1     -70.4185     2         
  2p(+0)   2      1      0      -70.4185     1         
  2p(+1)   2      1      1      -70.4185     1         
  Valence orbital: 2s  Z_eff=4.550  r_mp=0.609 Å
  Radial probability P(r) = r²|R_nl(r)|²:
  0 Å  ##:.                                    10.0 Å


╔══════════════════════════════════════════════════════╗
║  DEMO 2: H-H covalent bond vs. van der Waals (two different physics)║
╚══════════════════════════════════════════════════════╝
  Van der Waals (LJ):  ε=0.00191 eV   σ=2.8860 Å   r_min=3.2394 Å
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

  1.2500    —             173.941035      ########################################
  1.3750    —             55.138778       ########################################
  1.5000    —             19.251154       ########################################
  1.6250    —             7.276015        ########################################
  1.7500    —             2.934933        ########################################
  1.8750    —             1.248046        ########################################
  2.0000    —             0.553166        ########################################
  2.1250    —             0.252640        ########################################
  2.2500    —             0.117371        ########################################
  2.3750    —             0.054539        ########################################
  2.5000    —             0.024686        ########################################
  2.6250    —             0.010325        ########################################
  2.7500    —             0.003425        ########################################
  2.8750    —             0.000181        #####################
  3.0000    —             -0.001255       ######
  3.1250    —             -0.001797       #
  3.2500    —             -0.001907       
  3.3750    —             -0.001817       
  3.5000    —             -0.001645       ##
  3.6250    —             -0.001449       ####
  3.7500    —             -0.001256       ######
  3.8750    —             -0.001080       ########
  4.0000    —             -0.000925       ##########
  4.1250    —             -0.000790       ###########
  4.2500    —             -0.000675       ############
  4.3750    —             -0.000577       #############
  4.5000    —             -0.000494       ##############
  4.6250    —             -0.000424       ###############
  4.7500    —             -0.000365       ################
  4.8750    —             -0.000314       ################
  5.0000    —             -0.000272       #################

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
║  DEMO 3: H2O molecule — NVT MD at 300 K            ║
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
  single     C6     H6     1.1016    (order 1)

  N1 bond count: 3 (expect 3: C2, C6, H - donor ready)
  N3 bond count: 2 (expect 2: C4, C2 - bare, acceptor ready)
  Ring planarity: max deviation = 0.0006 A
  Total molecular charge: -0.000014 e

  THYMINE (C5H6N2O2, = 5-methyluracil) - 15 atoms
  Ring planarity: max deviation = 0.0092 A
  C5-CH3 methyl bond length: 1.5100 A (target 1.51 A, toluene-type)
  H-CM-H methyl angle: 109.47 deg (ideal tetrahedral 109.47)
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
  Energy breakdown at initial placement: E_LJ=0.736544 eV  E_Coulomb=-5.294565 eV  Total_PE=-4.558021 eV
  Closest intermolecular contact: atom 15 (Z=1) ... atom 20 (Z=8) = 1.967 A
  Largest single LJ repulsion: atom 6 (Z=7, sigma=3.25) ... atom 16 (Z=7, sigma=3.25) = 2.950 A, contributes 0.0415 eV
  Step   t(fs)      PE(eV)     T(K)       primary(A)
  1      0.100      -4.557881  49.96      2.9501    
  2      0.200      -4.557660  49.90      2.9501    
  3      0.300      -4.557362  49.82      2.9502    
  4      0.400      -4.556990  49.72      2.9502    
  5      0.500      -4.556548  49.60      2.9503    
  101    10.100     -4.588062  58.06      2.9166    
  201    20.100     -4.606500  61.28      2.8780    
  301    30.100     -4.626525  64.52      2.8972    
  401    40.100     -4.638784  65.00      2.9410    
  501    50.100     -4.650537  65.40      2.9889    
  601    60.100     -4.655535  64.29      3.0321    
  701    70.100     -4.674488  66.36      3.0718    

  Initial interaction PE:        -4.558021 eV
  PE at closest WC approach:     -4.611601 eV (primary N1...N3 = 2.878 A)
  Global PE minimum over run:    -4.695777 eV (may reflect drift to a different,
                                  non-WC configuration such as stacking)
  Final PE (end of run):         -4.662334 eV

  Final heavy-atom contacts:
    G:N1...C:N3 = 3.105 A
    G:N2...C:O2 = 3.126 A
    G:O6...C:N4 = 2.948 A

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
  G-C @ closest WC approach: -4.611601 eV (3 H-bonds)
  A-U @ closest WC approach: -1.080053 eV (2 H-bonds)
  --> G-C binds MORE strongly than A-U (3.531549 eV difference),
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

  Total charge: -1.0550 e (real backbone convention: -1 per
  phosphodiester; approximate here since the sugar and phosphate
  charges are not independently verified the way the nucleobase
  RESP charges are - see nucleobases.c for full honest sourcing)

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

  Total charge: -0.3044 e (approximate - amino acid charges are not
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
  Initial clash relaxation: 18811.75 -> 8.50 eV
  Added 12 dihedral restraints (phi, psi, omega for applicable
  residues) toward real textbook values.
  Minimized: PE = 9.7447 eV

  Residue    phi (deg)    psi (deg)   
  1          -58.75      
  2          -61.51      
  3          -58.29      
  4          -53.96      
  0                       -48.38      
  1                       -47.55      
  2                       -52.66      
  3                       -48.26      

  Max deviation from target (phi=-57, psi=-47): 5.66 deg

  === The i,i+4 backbone hydrogen bond (not programmed in) ===
  N-H(4) ... O=C(0): H...O = 2.1637 A, N...O = 3.1286 A
  Real backbone H-bond range: H...O 1.8-2.2 A, N...O 2.8-3.2 A

  --> A real backbone hydrogen bond formed. Given only the
      correct LOCAL torsion geometry (phi/psi/omega, all real,
      textbook values) and the SAME validated Coulomb+LJ force
      field already proven on the water trimer and G-C/A-U
      pairing, the defining GLOBAL structural feature of an
      alpha helix emerged on its own.

  Bond length range: [0.9602, 1.5669] A   Bad bonds: 0   Valence issues: 0

╔══════════════════════════════════════════════════════╗
║  DEMO 12: KcsA selectivity filter - K+ vs Na+, second pass║
╚══════════════════════════════════════════════════════╝
  Four real-charge carbonyl O's, real 4-fold symmetry, radius =
  literature K+-coordination target. This pass uses the real
  amino-acid-specific carbonyl LJ typing (AA_LJ_O_EPS/SIGMA from
  aminoacids.c) instead of generic periodic-table oxygen - see
  source comment for exact scope and what changed from pass one.

--- Gly77 site, PDB 1K4C LINK record target: 2.72 A ---
  K+   E_LJ =   0.576511 eV   E_Coulomb =  -5.515810 eV   Total_PE =  -4.939299 eV
  Na+  E_LJ =   0.039705 eV   E_Coulomb =  -5.515810 eV   Total_PE =  -5.476105 eV
  Delta (Na+ minus K+): -0.536806 eV  (Na+ favored, wrong direction)

--- Val76 site, target: 2.83 A ---
  K+   E_LJ =   0.315958 eV   E_Coulomb =  -5.301414 eV   Total_PE =  -4.985456 eV
  Na+  E_LJ =   0.003689 eV   E_Coulomb =  -5.301414 eV   Total_PE =  -5.297725 eV
  Delta (Na+ minus K+): -0.312269 eV  (Na+ favored, wrong direction)

  Honest read (fixed-radius tests): switching from generic-O to
  amino-acid-specific carbonyl LJ typing is the only change from
  the first pass.
  Still Na+-favored at both sites even with correct LJ typing.
  That rules out generic-O typing as the (sole) cause and
  points more toward the remaining candidates: fixed (not
  relaxed) geometry, and the lack of polarizability.

--- Letting each ion find its own preferred radius (2.00-4.20 A scan, 0.02 A steps) ---
  K+   natural radius = 2.820 A   E_min = -4.986164 eV
  Na+  natural radius = 2.440 A   E_min = -5.780797 eV
  K+ natural radius 2.820 A vs Na+ 2.440 A - K+ prefers the larger cage, as expected for the bigger ion
  At each ion's OWN best radius: K+ E_min = -4.986164 eV vs Na+ E_min = -5.780797 eV -> Na+ favored, wrong direction

  Honest read (radius-flexible test): this is the same real
  charges and the same real carbonyl LJ typing as above - the
  only thing that changed is letting each ion pick its own
  distance instead of forcing both to the crystal's real but
  shared 2.72/2.83 A. Na+ still winning even with each ion free to pick its own
  radius is the stronger negative result of the two - it
  says this isn't a geometry-fit problem at all, and
  polarizability moves from 'the last remaining candidate'
  to 'the most likely one.'

--- Fuller real geometry: Thr75-O + Val76-O together (the real
  antiprism site S3/K-C3003 actually shows, not a single ring) ---
  K+   E_LJ =  14.808098 eV   E_Coulomb =  12.402345 eV   Total_PE =  27.210443 eV
  Na+  E_LJ =  13.903269 eV   E_Coulomb =  12.402345 eV   Total_PE =  26.305614 eV
  Delta (Na+ minus K+): -0.904829 eV  (Na+ favored, wrong direction)
  Honest read: both real target distances (Thr75 2.70 A, Val76
  2.83 A) at once, real 45-degree antiprism offset, same real
  charges and typing as every test above - the only new variable
  is twice the coordinating oxygens. CAVEAT - stronger than a
  geometry note: because the real ring z-separation was not
  sourced this pass, both rings sit at z=0, which puts oxygens
  from the two rings ~2.1 A apart, deep inside their mutual LJ
  repulsion. The large positive energies above are dominated by
  that artificial O-O overlap, NOT by ion coordination, so the
  numbers in this block are NOT interpretable as selectivity
  energetics until a sourced z-separation is added. Treat this
  block as a placeholder, not a result.

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

---

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

Potassium's LJ ε/σ in `periodic_table.c` are verified against Rappé et
al. 1992 UFF (K: x1=3.812 Å, x2=0.035 kcal/mol) as of audit T1,
resolving the masterplan Step 0 flag.

### Known limitations & honest scope

* **1-4 non-bonded scaling deliberately not applied (audit F1).** AMBER's
  1-4 scaling (LJ/2, Coulomb/1.2) presupposes AMBER's co-fitted torsion
  parameters; this force field's dihedrals are restraints, not fitted
  torsions, so scaling 1-4 pairs removes physics with nothing to
  compensate. It was tested and broke the validated helix i,i+4 H-bond
  (Demo 11). 1-4 pairs stay at full strength; revisit only with fitted
  torsions.
* **Ion selectivity (Demo 12) is a negative result.** The Coulomb+LJ
  model favors Na⁺ over K⁺ at every KcsA filter site and across a
  radius scan — the wrong direction. The second pass ruled out generic
  oxygen typing (switching to the real amino-acid-specific carbonyl LJ
  changed nothing) and the radius scan ruled out geometry-fit (Na⁺
  still wins at each ion's own preferred radius), pointing at missing
  electronic **polarizability** as the most likely missing physics.
  This is reported as a finding, not hidden.
* **Demo 12's antiprism block is a placeholder, not a result.** The
  real ring z-separation was not sourced this pass, so both rings sit
  at z=0 and oxygens from the two rings overlap inside their mutual LJ
  repulsion; the block's energies are dominated by that artifact. It
  must not be read as selectivity energetics until a sourced
  z-separation (1K4C) is added.
* **Base-pairing energetics (Demo 7) are qualitative.** G–C correctly
  binds more strongly than A–U and both pairs are correctly
  attractive, but the absolute magnitudes overshoot gas-phase ab
  initio references (~−1.2 eV G–C, ~−0.55 eV A–U) by roughly 4×, and
  no single dielectric constant brings both pairs into simultaneous
  quantitative agreement. The ordering is validated; the numbers are
  not yet quantitatively trustworthy.
* **No explicit solvent.** Everything runs in vacuum with a relative-
  permittivity divisor (dielectric = 4 for base pairing) standing in
  for condensed-phase screening. This is the root of the quantitative
  base-pairing error above and a major limitation for any real
  biomolecular energetics — including the ion-selectivity question,
  where the dehydration penalty the real channel balances against
  coordination cannot be expressed at all in a vacuum calculation.
* **Approximate charges on non-nucleobase components.** Sugar,
  phosphate, and amino-acid partial charges are charge-balanced
  approximations, not verified RESP fits (see the confidence tiers).
* **Scaling and algorithms.** The non-bonded loop is O(N²) with no
  neighbour list (the code suggests a cell list beyond ~500 atoms);
  dihedral forces use 24 finite-difference energy evaluations per
  dihedral per call (correct-by-construction, but expensive); and the
  minimizer is steepest descent only — no conjugate gradient or
  L-BFGS.
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

The honest next steps, in the order the project's own record supports:

1. **Complete Demo 12's antiprism block**: source the real ring
   z-separation from 1K4C so the placeholder becomes a real result or
   a real exclusion.
2. **The missing-physics decision**: polarizability, explicit solvent
   competition, or both — the vacuum cage calculation cannot express
   the dehydration term that real selectivity depends on.
3. **The longer-horizon tracks already named**: a full DNA duplex,
   gene-regulatory logic, a synapse between neurons, and eventually
   deriving ion-channel gating from actual protein structure — closing
   the loop between the protein and electrophysiology tracks.
