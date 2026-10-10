# Carbon VM Structured Datastream — Schema v1

**Status:** specification (implementation lands in s41; first consumer in s42).
**Companion:** `MASTERPLAN.md` Part 5. This file is the single source of
truth for the simulator's machine-readable output. The printf demo log
remains for humans; this datastream exists for verification.

---

## 1. Purpose
Every scientific claim the simulator makes must be written here as a
structured record — computed once, stored, and checkable by an external
validator against a cited reference. This replaces "parse a printf line"
with "read a typed field."

## 2. File format
- Plain text, UTF-8, LF line endings.
- `#` begins a comment (whole line). Blank lines are ignored.
- `[section]` starts a named section.
- Data rows are whitespace-separated fields. A section may begin with a
  `# field field ...` line naming its columns.
- File extension `.cvmds`. One file per run, e.g. `kcsa.cvmds` in the tree
  root (historical spec said `run/kcsa.cvmds`, but `run` is an executable
  script in this tree, not a directory; root path is the documented
  deviation, see Demo 12 consumer comment).
- The **payload** is every line from the first `[header]` line through the
  line immediately before `[end]`. `[end]` carries its SHA-256.
- `wallclock` is live by default; set `SOURCE_DATE_EPOCH` (seconds since
  epoch) for a byte-stable record artifact.

## 3. Sections

### 3.1 `[header]` — required
Key/value pairs, one per line, `key: value`.
- `schema:` integer schema version (this document = `1`)
- `demo:` identifier, e.g. `kcsa`, `helix`, `water_md`, `basepairing`
- `build-compiler:` e.g. `gcc 13.3.0` (full version; major-only `gcc 13`
  is a legacy coarse value, still accepted by verifiers)
- `build-flags:` exact CFLAGS used (a *record* build must NOT carry `-march=native`)
- `source-hash:` SHA-256 of the source tree or git commit when VCS is
  available; otherwise a `record-tree-*` token describing the tree state
  (no-VCS fallback, e.g. `record-tree-v9R4-unversioned`; validators check
  presence and the `record-tree-` prefix, not hash equality, in that mode).
  Extra free-form suffixes after `;` are FORBIDDEN (full-audit O10).
- `rng-seed:` integer seed for any stochastic initialisation
- `wallclock:` ISO-8601 timestamp of the run

Additional `kebab-case: value` header lines beyond the seven above (e.g.
`version-internal`, `cage-geometry`) are PERMITTED, are covered by the seal,
and MUST be ignored by verifiers (full-audit O12).

### 3.2 `[atoms]` — optional
Static per-atom identity, one row per atom.
Columns: `idx  Z  symbol  mass_amu  charge_e  lj_eps_ev  lj_sigma_a`

`lj_sigma_a` is the 12-6 Lennard-Jones **collision diameter** — the distance at
which the potential crosses zero — and is the quantity consumed by the
`4 eps [(sigma/r)^12 - (sigma/r)^6]` form in `pair_nonbonded_core()`. It is NOT
UFF's `x1` column, which is Rmin (the distance of minimum) and exceeds sigma by
2^(1/6). The periodic table was corrected during the v9R4 deep audit to store
`x1 / 2^(1/6)`; amino acid and nucleobase sites use the equivalent AMBER
`Rstar -> sigma` conversion. The minimum of any pair potential reconstructed
from a record's `lj_eps_ev`/`lj_sigma_a` columns should therefore come out at
`2^(1/6) * sigma`, and can be used as a validation check on a consumer.

### 3.3 `[steps]` — optional
Per-timestep scalar observables, one row per step.
Columns: `step  time_fs  KE_ev  PE_ev  E_ev  T_k`

### 3.4 `[claims]` — required for any demo asserting a scientific result
One row per claim.
Columns: `key  value  unit  provenance`
- `key`: dot-separated namespace, `<demo>.<subsystem>.<quantity>`,
  lower-case except documented element/energy shorthands (`e_k/e_na`
  for K+/Na+ ions, `de_*`, `s_*`, `bo_*`, `q_*`, `alpha_o`; legacy files
  may carry `E_K/E_Na/ddG/ddU/q_O/dE/S/BO/alpha_O` — verifiers accept both
  but new writers must use lower-case)
- `value`: signed decimal number, or a short token string
- `unit`: `eV`, `A`, `A^3`, `K`, `fs`, `e` (elementary charge),
  `dimensionless`, or `-`
- `provenance`: `computed` or `computed-<leg>(-<detail>)*` where `<leg>` is one
  of `jc2008-params`, `ecc-scaled`, `v2`, `v2-induction`, `v2-pauli`, `scf-qm`,
  `relaxed-scf`, `stiff-scaffold`, `6fold`, `exchange`, `knockon`, `scf-uz`,
  `explicit-6water`, `explicit-6water-polar`, `polar-wham`, `polar-wham-300k`,
  `wham-300k`, `scf-uz`, or a leg documented in §6 (full-audit O13; verifiers
  MUST whitelist, not substring-match), `computed-jc2008-params` (JC ion sizes),
  `computed-ecc-scaled`, or a citation tag — `Marcus1991-TATB`
  (TATB-based absolute single-ion hydration FREE energies, Marcus 1991
  Faraday Trans. 87, 2995; full-audit N1 corrects the earlier
  `Marcus1997-TATB` tag, whose -322/-454 values were old absolute-scale
  hydration ENTHALPIES; legacy files may carry `Marcus1997-TATB` —
  verifiers accept both but new writers must use `Marcus1991-TATB`),
  `1K4C-LINK`, `Aduri2007`, `expt-1000:1@300K`, etc. `JC2008` alone is
  deprecated: JC numbers are computed with JC parameters, not measured by JC.

### 3.5 `[end]` — required
- `payload-sha256:` SHA-256 over the payload (defined in §2).

**FULL-AUDIT Q3 (2026-10-10) — the seal must terminate the file.** The
`payload-sha256:` line is the last content in the file: exactly one optional
trailing newline, then EOF. Nothing may follow it. Content after the seal is
outside the hashed payload *and* outside every structural check, so it would
otherwise verify clean while being covered by no digest — a claim this
format cannot honestly make. `ds_verify_file` rejects such a file, and
`test_datastream` pins the rejection. (Measured before the fix: appending a
claim line to a sealed file left `ds_verify_file() == 0`, i.e. accepted.)

## 4. Naming & sign conventions
- Claim keys are lower-case hierarchical (`<demo>.<subsystem>.<quantity>`)
  with the §3.4 element/energy exceptions; new keys must be lower-case.
- Energies in eV, distances in Å (volumes `A^3`), charges in `e`,
  temperatures in K unless noted.
- `ddg/ddu` are differences; sign convention stated ONCE in the Demo 12
  forensic header + readme §8 (`dU = E(K+) − E(Na+)`, `+ = Na+ favored`;
  exchange lines `negative = K+ selective`; FULL-AUDIT O14: per-claim
  trailing `#` comments are NOT emitted by ds_add_claim and are NOT
  required — the demo header carries the convention; legacy `ddG/ddU`
  accepted).
- Pore-axis single-point `U(z)` profiles are NOT free-energy PMFs;
  label them `U(z)`, never `PMF`, unless sampling/entropy is included.

## 5. Validation contract
A validator reads a `.cvmds` file and, for each claim whose provenance is
an external citation, compares `value` against the masterplan Part 6
reference within a stated tolerance. Claims with provenance `computed` are
checked only for internal consistency (e.g.
`ddG_corrected == ddG_vacuum + dehyd.K - dehyd.Na`).

## 6. KcsA worked example (first consumer, wired in s42)
The Demo 12 antiprism + dehydration run must emit at least (current
lower-case key form; legacy upper-case files accepted):

```
[claims]
kcsa.antiprism.e_k         <eV>    eV  computed
kcsa.antiprism.e_na        <eV>    eV  computed
kcsa.antiprism.ddg_vacuum  <eV>    eV  computed      # e_k - e_na; + means Na+ favored
kcsa.dehyd.k               +3.061  eV  Marcus1991-TATB
kcsa.dehyd.na              +3.786  eV  Marcus1991-TATB
kcsa.ddg_corrected         <eV>    eV  computed      # ddg_vacuum + dehyd.k - dehyd.na
kcsa.ddg_experimental      -0.179  eV  expt-1000:1@300K
kcsa.ddg_deviation         <eV>    eV  computed      # corrected - experimental
```
Plus JC (`computed-jc2008-params`), ECC (`computed-ecc-scaled`), `U(z)`
(`kcsa.pmf.*`, single-point, not free energy), induction (`kcsa.pol.*`),
v2 in-loop polarization/Pauli (`kcsa.v2.*`: `computed-v2-induction`,
`computed-v2-pauli`, `computed-v2`), v3 SCF (`kcsa.scf.*`:
`computed-scf-qm`, hybrid JC wall + SCF charges + polar/Pauli/disp),
relaxed coordination (`kcsa.relax.*`: `computed-relaxed-scf`, free ion
multi-start with hysteresis; `kcsa.strain.*`: `computed-stiff-scaffold`;
`kcsa.cn.*`: `computed-6fold` octahedral probe;
`kcsa.exchange.*`: `computed-exchange`, rigid deterministic +
relaxed with error, 2-ion `exchange2`; `kcsa.pair.*` + `kcsa.knock.*`:
`computed-knockon` pair minima, conductive single-points, and valid-
subset landscape barriers; `kcsa.scfuz.*`: `computed-scf-uz` SCF-polar pore
profile with K+-favored single-point barrier), explicit 6-water hydration
(`kcsa.hyd.*`: `computed-explicit-6water[-polar]`), WHAM free energy
(`kcsa.fe.*`: `computed-polar-wham-300k` 3-repeat polar sampling,
`computed-wham-300k` fixed-charge reference), and QM
(`kcsa.qm.*`, charges in `e`, `alpha_o` in `A^3`, `S_ref` at
covalent-contact distance) legs as emitted by Demo 12.
