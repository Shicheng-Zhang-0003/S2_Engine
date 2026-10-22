# Provenance and verification

*S2 Engine v9R4 · where every number comes from, and how each claim
is checked against something outside this repository.*

Two standards, both structural rather than aspirational: parameters
are checked against primary sources, and claims are checked by
independent oracles.

<a id="6-provenance-and-verification"></a>
## 6. Provenance and verification

Two standards, both structural rather than aspirational: parameters are
checked against primary sources, and claims are checked by independent
oracles. `make audit` runs the evidence.

| Subsystem | Source | How it was checked |
|---|---|---|
| Fundamental constants | 2019 CODATA | derived in-line; reciprocals exact by construction |
| Periodic-table LJ | UFF / Rappé et al. 1992 | all 36 elements, ε and σ, cell-by-cell; σ = x1/2^(1/6) |
| Biomolecular LJ | AMBER ff99 (`parm99.dat`) | atom-by-atom mapping; published σ reproduced (N 3.25000, C 3.39967, …) |
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

<a id="61-third-audit-pass-2026-10-10-—-independent-oracle-sweep"></a>
### 6.1 Third audit pass (2026-10-10) — independent oracle sweep

A full mathematical, programming and operational re-audit ran an
**independent verification harness** built outside the tree (`/tmp`, linked
against the shipped engine objects, 78 checks) rather than reusing the
project's own oracles, so a shared mistaken assumption could not pass twice.
Its findings:

* **Constants — clean.** Every SI literal in `include/constants.h` matches
  CODATA 2019 exactly (h, c, e, k, N_A by redefinition; mₑ, m_p, m_n, ε₀,
  a₀, E_h, amu to the published uncertainty), and every derived conversion
  (`COULOMB_MD` = 14.3996454785 eV·Å, `MD_FORCE_CONV` = 9.6485332157e-3,
  `KCAL_MOL_TO_EV` = 0.0433641042, `AMU_AFS2_TO_EV` = 103.6426965) was
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

<a id="62-fourth-audit-pass-2026-10-21-—-full-mathematical-programming-and-operational-re-audit"></a>
### 6.2 Fourth audit pass (2026-10-21) — full mathematical, programming and operational re-audit

The eighth gate, `tests/test_audit4.c` (80 checks), exists because of this
pass: one check per defect, each built on an oracle *outside* the thing it
checks — the upstream `pcg-c` expected-output file, textbook overlap
closed forms, the mirror symmetry of a bent H–C–H, a 3D grid integration
of the actual valence charge densities, the L'Hôpital limits of the
Hodgkin–Huxley rate functions, and the KcsA chain re-derived from the
deposited geometry. The audit's own working material is in
[`../audit/README.md`](../audit/README.md).

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
| RHF energy | the nuclear–nuclear repulsion was missing | the H₂ curve had no minimum at all | fixed; the minimum is now at R = 1.300 Å |
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

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
