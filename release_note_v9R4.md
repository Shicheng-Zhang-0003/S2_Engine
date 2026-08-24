# carbonsim v9R4 — the hygiene closure release

Not a new-physics release. Completion of the correctness and hygiene
debt carried out of the v9R3 audit, executed as seven single-purpose,
individually-gated scripts (s01–s07) run from the folder containing
this tree. Every code change was proven behavior-neutral by a
byte-identical record; every documentation change was gated on the
record staying intact.

## Highlights
- **The record never moved.** Every regeneration across the whole arc
  is byte-for-byte identical to the v9R3 post-audit record — the same
  SHA-256 at every step, asserted after every step.
- **Both build configurations compile warning-clean under -Wall and
  -Wextra, with nothing suppressed.** The `-Wno-missing-braces` flag
  was removed, not extended, once the flood it masked was eliminated
  at the source.

## What the arc fixed
- **Q1 — quantum normalization.** `quantum_radial_wavefunction`
  computed the hydrogen-like normalization denominator with (n+l)!
  CUBED; the correct formula uses (n+l)! to the first power. The old
  wavefunctions were under-normalized by a factor of (n+l)! — 2x for
  any 2s orbital, 6x for 2p; H 1s happened to be exact because
  1! = 1. Fixed. Proven behavior-neutral: every consumer is either
  shape-normalized (the ASCII plot divides by its own P_max) or an
  argmax (quantum_most_probable_radius), so the record could not move
  — and did not.
- **B1–B3 — initializer braces.** The PERIODIC_TABLE `ground_config`
  initializers `{0}` triggered `-Wmissing-braces` under the ASan
  build line (which carries `-Wall` without the suppression). Three
  passes: (1) GCC's own fixit `{0}` → `{{0}}` on all 36 element rows
  — insufficient, because ElectronConfig nests
  `int config[MAX_SHELLS][4]` and GCC descends a level; (2) `{{{0}}}`
  — silenced `-Wmissing-braces` but lost the `{0}` exemption from
  `-Wmissing-field-initializers` under `-Wextra`; (3) fully explicit
  `{{{0}}, 0, 0}` on all 36 rows — names every ElectronConfig field,
  relies on no compiler exemption, runtime-identical (the fields were
  zero before and are zero now). Sentinel row untouched by design.
- **M1 — dead suppression removed.** `-Wno-missing-braces` removed
  from the makefile's normal CFLAGS line, with the reasoning kept
  in-source above the surviving `CFLAGS += -march=native`. Verified,
  not assumed: the normal build (now unsuppressed) and an ASan build
  with `-Wall -Wextra` were both gated at zero warnings of any class.
- **N1 — naming hygiene.** The readme's operational references now
  point at v9R4: build path, embedded-output header, current-version
  line, Current Status. Historical entries (the v9A1 version-history
  line, the v9R3 audit section body) left untouched by design.

## Gate discipline — including the misses the gates caught
- s06's first version carried a post-fix check that counted ANY line
  containing the flag string — including the script's own explanatory
  comment — as "flag still present." Its own gate caught the false
  failure before the script could pass; the fixed check greps active
  (non-comment) makefile lines only.
- s04's first version died on a sed expression error; caught by the
  script's own post-fix verification before any gate ran, and the
  tree was left in a known, asserted state.
- Every script re-asserted the full record through the s01 harness
  (normal + ASan builds, output SHAs, both program stderrs). The
  record SHA never changed across the whole arc.

## Scripts in this arc
s01 record-verification harness · s02 quantum normalization fix ·
s03–s05 initializer braces (three passes) · s06 dead-suppression
removal · s07 naming hygiene · s08 this release record.

## Open items (unchanged by this release)
v9R4 closed hygiene debt, not physics. The v9R3 open items carry
forward untouched:
1. **Complete Demo 12's antiprism block**: source the real ring
   z-separation from 1K4C so the placeholder becomes a real result
   or a real exclusion. Until then its numbers must not be cited as
   selectivity energetics.
2. **The missing-physics decision**: polarizability, explicit-solvent
   competition (the dehydration penalty the vacuum calculation cannot
   express at all), or both.
3. **Deferred audit P2 (flagged, not failed)**: analytic dihedral
   gradients (F4) — finite differences remain the oracle.

## Record (SHAs)
- Normal build output (`output.txt`):
  `875a2c0cf30ccf4fc6ebff8b0b64547063c7a80c5250d1b32c72af3048bc6935`
- ASan build output (`output.asan.txt`): same SHA — byte-for-byte
  identical to the normal build. ASan stderr: empty (no memory
  errors). The readme's embedded output block is unchanged and
  matches this record.

## Build
```bash
cd biological/v9R4
make
./carbonsim > output.txt
```
Requires a C11 compiler and `make`. Both the default build and the
sanitizer build compile warning-clean (-Wall -Wextra) with no
suppression flags.
