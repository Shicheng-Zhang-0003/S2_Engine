# Operations — the record, the artifacts, and the pipeline

*S2 Engine v9R4 · how the record is produced, verified, and
promoted, and how to read it.*

<a id="the-record-and-its-artifacts"></a>
## The record and its artifacts

The record is the tracked stdout of `./carbonsim`, and it is the authority
this document describes — not the other way around.

| Artifact | What it is |
|---|---|
| `output.txt` | the verbatim 20-demo record (~1250 lines) |
| `CURRENT_BASELINE_SHA.txt` | its SHA-256 (single source of truth — do not duplicate digests in prose) |
| `kcsa.cvmds` | a schema-1 datastream written by Demo 12: 147 worked-example claims with provenance tags and a SHA-256 payload seal (`DATASTREAM_SPEC.md`). Covers every KcsA leg including Demo 12c's Born layer and contact mismatch |
| `runs/` | candidate runs; untracked. Promotion is explicit (`./run --accept`) precisely so an accidental overwrite cannot masquerade as a record update. **Not** a ring buffer — nothing prunes it, so it grows until you remove it |
| `audit/` | the audit's own evidence tree: independent math oracles, C harnesses, external references, pre-fix engine copies |
| `docs/` | this documentation set (see the map in `readme.md`) |
| `build-perf/`, `carbonsim-perf`, `s2tui-perf` | the OpenMP performance build's PRIVATE output. It has its own object directory and its own binary names so that it can never contaminate the record build. Untracked. |

Never write `output.txt` by hand. The stdout/stderr split is deliberate:
stdout is DATA (the record — never wall time, never progress), stderr carries
human display (timings, heartbeats, leg progress) and only when stderr is a
TTY; `CARBON_QUIET` silences it. Pipes and redirects therefore stay honest.
See `include/display.h`.

**The record's stderr must be empty.** That is a gate, not a convention:
`verify_scripts.sh` and `./s01_verify_record.sh` both assert it, and
`./run --accept` refuses to promote a run whose stderr is not empty. The
fourth audit found the temperature cap writing 69 lines to stderr in batch
mode, which made the gate — and therefore the promotion — impossible to
satisfy. A cap event is now counted in `sim->temp_cap_events` and reported
as DATA.

### The promotion pipeline

| Step | Command | What it guarantees |
|---|---|---|
| build | `make clean && make` | a warning-free build from the tracked sources |
| stage | `./run` | writes `runs/<stamp>.{txt,stderr.txt,cvmds}`; prints the recorded and the new digest side by side; **never** touches `output.txt` |
| promote | `./run --accept` | refuses unless stderr is empty, `make test` is green, and the staged `kcsa.cvmds` seal verifies; then writes each artefact to a temporary name and renames it into place, the ledger last |
| verify | `./verify_scripts.sh` | the baseline digest, the seal (twice — through the engine and through an independent Python re-computation), a clean build, all eight gates with their pinned counts, the verdict-slot structure, and the TUI contract |
| reproduce | `./s01_verify_record.sh` | the record reproduces **byte-for-byte** under AddressSanitizer and UndefinedBehaviorSanitizer with an empty stderr, in both builds |
| audit | `make audit` | the independent oracles and C harnesses, against a clean-first build |
| CI | `.github/workflows/ci.yml` | all of the above, ending in `git diff --exit-code` on the three record artefacts |

**Why the last CI step exists.** The October-2026 commits that added demos
13–16 all carry the message *"all 491 gates green"* and all ran `make test`
only. The **record** gate was never run, so the tracked `output.txt`
drifted away from HEAD while every commit message claimed a green tree.
The staleness check is now a CI step, so a record that does not match its
source fails the build instead of being rediscovered months later.

**The OpenMP performance build is isolated.** `make perf` builds
`build-perf/` with `carbonsim-perf` and `s2tui-perf`. It used to share the
record build's objects and binary names, which produced three separate
failures: `make perf` after `make` was a **silent no-op** with no OpenMP
code in the "multicore" binary, `make clean && make perf && make` failed
to link the record, and a perf run **overwrote the tracked record
binary**. Never promote a perf run to the record: thread scheduling
changes the floating-point summation order, so the digests differ.

---
<a id="interpreting-the-output"></a>
## Interpreting the output

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

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
