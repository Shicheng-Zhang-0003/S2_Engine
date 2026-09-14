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
- File extension `.cvmds`. One file per run, e.g. `run/kcsa.cvmds`.
- The **payload** is every line from the first `[header]` line through the
  line immediately before `[end]`. `[end]` carries its SHA-256.

## 3. Sections

### 3.1 `[header]` — required
Key/value pairs, one per line, `key: value`.
- `schema:` integer schema version (this document = `1`)
- `demo:` identifier, e.g. `kcsa`, `helix`, `water_md`, `basepairing`
- `build-compiler:` e.g. `gcc 13.2.0`
- `build-flags:` exact CFLAGS used (a *record* build must NOT carry `-march=native`)
- `source-hash:` SHA-256 of the source tree, or git commit if available
- `rng-seed:` integer seed for any stochastic initialisation
- `wallclock:` ISO-8601 timestamp of the run

### 3.2 `[atoms]` — optional
Static per-atom identity, one row per atom.
Columns: `idx  Z  symbol  mass_amu  charge_e  lj_eps_ev  lj_sigma_a`

### 3.3 `[steps]` — optional
Per-timestep scalar observables, one row per step.
Columns: `step  time_fs  KE_ev  PE_ev  E_ev  T_k`

### 3.4 `[claims]` — required for any demo asserting a scientific result
One row per claim.
Columns: `key  value  unit  provenance`
- `key`: dot-separated namespace, `<demo>.<subsystem>.<quantity>`
- `value`: signed decimal number, or a short token string
- `unit`: `eV`, `A`, `K`, `fs`, `dimensionless`, or `-`
- `provenance`: `computed`, or a citation tag — `Marcus1991`, `1K4C-LINK`,
  `Aduri2007`, `expt-1000:1@300K`, etc.

### 3.5 `[end]` — required
- `payload-sha256:` SHA-256 over the payload (defined in §2).

## 4. Naming & sign conventions
- Claim keys are lower-case and hierarchical: `<demo>.<subsystem>.<quantity>`.
- Energies in eV, distances in Å, temperatures in K unless noted.
- `ddG` is a free-energy difference; its sign convention is stated in the
  claim's trailing comment or provenance.

## 5. Validation contract
A validator reads a `.cvmds` file and, for each claim whose provenance is
an external citation, compares `value` against the masterplan Part 6
reference within a stated tolerance. Claims with provenance `computed` are
checked only for internal consistency (e.g.
`ddG_corrected == ddG_vacuum + dehyd.K - dehyd.Na`).

## 6. KcsA worked example (first consumer, wired in s42)
The Demo 12 antiprism + dehydration run must emit at least:

```
[claims]
kcsa.antiprism.E_K         <eV>    eV  computed
kcsa.antiprism.E_Na        <eV>    eV  computed
kcsa.antiprism.ddG_vacuum  <eV>    eV  computed      # E_K - E_Na; + means Na+ favored
kcsa.dehyd.K               +3.057  eV  Marcus1991
kcsa.dehyd.Na              +3.783  eV  Marcus1991
kcsa.ddG_corrected         <eV>    eV  computed      # ddG_vacuum + dehyd.K - dehyd.Na
kcsa.ddG_experimental      -0.179  eV  expt-1000:1@300K
kcsa.ddG_deviation         <eV>    eV  computed      # corrected - experimental
```
