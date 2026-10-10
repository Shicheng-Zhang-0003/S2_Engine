# audit/ — everything that lives outside the engine's own source

This directory holds the working material of the v9R4 audit. It was
previously split across `biological/temp/audit/` and `biological/temp_audit/`,
one and two directories above the engine tree, with absolute paths baked into
every script. Moving it inside `v9R4/` means every path here is relative to
the tree root, so the audit is reproducible from a fresh clone rather than
only on the machine that produced it.

**None of this compiles into the engine.** `src/`, `include/` and `tests/`
are the product; this is the evidence and the harnesses that produced it. The
makefile does not reference this directory and `make test` does not run
anything here.

## Layout

| path | what it is |
|---|---|
| `notes/FINDINGS.md` | the audit's running findings log |
| `math/` | independent mathematical oracles in Python (symbolic gradient checks, constants, pair and dihedral derivatives) |
| `harness/` | C harnesses linked against the engine objects — the dipole solver, QEq conservation, minimizer guards, KcsA geometry |
| `external/` | checks against sources outside the repository, plus the fetched reference files |
| `probe_thirdpass.c` | the 2026-10-10 third-pass verification harness (78 independent checks: FD forces for every term, NVE conservation, thermostat ensembles, HH closed forms, QEq ordering, KcsA, SHA-256) |
| `thirdpass_uff_amber.py` | third-pass table verification: all 36 UFF σ/ε and 13 AMBER ff99 types against primary transcriptions |
| `external/uff_reference.prm` | machine-readable transcription of UFF Table II (Rappé 1992), vendored so the table check runs offline |
| `orig/` | verbatim copies of the **pre-fix** `qm.c`, `integrator.c` and `constants.h`, used to prove a regression test actually goes red on the defect it claims to catch |
| `history/` | one-off fix and verify scripts from the pre-v9R4 pass, kept for provenance |
| `stale_artifacts/` | artifacts found misplaced, kept rather than deleted |
| `scratch/` | derived data: build logs and captured runs. **Not tracked** |

## Building and running the third-pass harness

The third-pass evidence is reproducible like the rest of the audit. It links
the *shipped* engine objects, so build first, then:

```bash
make                                  # build the engine
gcc -O2 -std=c11 -Iinclude -o /tmp/probe_thirdpass \
    audit/probe_thirdpass.c \
    build/{aminoacids,datastream,forces,integrator,kcsa_filter,loop,neuron,nucleobases,periodic_table,qm,qm_eht_overlap,qm_eht_scf,quantum,sim}.o -lm
/tmp/probe_thirdpass                  # prints per-check PASS/FAIL, exits non-zero on any
python3 audit/thirdpass_uff_amber.py  # UFF + AMBER table verification (offline)
```

The harness deliberately uses **independent** oracles — central finite
differences, closed-form hydrogenic results, statistical-ensemble
identities — rather than the project's own, so a shared mistaken assumption
cannot pass twice.

## `external/` reference files

Fetched once, on 2026-10-01, and kept so the external checks are reproducible
without a network. Hashes:

```
parm99.dat    41630 bytes  sha256 d2ec7da8eeb08a4bdc0f49a57798e2a7...
1K4C.cif.gz  135170 bytes  sha256 5af16a4ad3a1b93d8662a3b9e3e248da...
```

Sources:

* `parm99.dat` — the original AMBER ff99 distribution,
  `http://archive.ambermd.org/200609/att-0158/parm99.dat`. The NONBONDED
  section is what pins `R*` and `epsilon` for `O`, `OH`, `OS`, `N`, `C`, `CT`.
* `1K4C.cif.gz` — RCSB PDB entry 1K4C (Zhou et al., *Nature* 1998), the
  2.0 Å KcsA selectivity-filter structure. `gunzip -k 1K4C.cif.gz` to read it;
  `ext_check.py` parses the decompressed file.

The CODATA values are **not** kept as files. NIST refuses scripted requests
to `physics.nist.gov/cgi-bin/cuu` (every request returns "ill-formed request"),
so the constants in `tests/test_external.c` were transcribed by hand from the
rendered NIST pages and are cross-checked against the derivations that exist
for them — `h/2pi`, `E_h/e`, `4184/N_A/e`. Where no derivation exists, the
transcription stands on the citation in the source comment.

## `orig/` is deliberately not path-fixed

`orig/qm_orig.c` and `orig/integrator.c` are byte-identical to their commits
(`src/qm.c` at `7283356`). Their `#include "../include/..."` lines refer to
that original location, one level under the tree root. Rewriting them to
`../../include` would make them compile here but would stop them being a
faithful copy of the code under audit, which is the only reason they exist.
So they are left alone and `run_audit.sh` builds them from a scratch
directory at the right depth.

## Running it

```bash
./audit/run_audit.sh            # every harness and oracle
./audit/run_audit.sh --revert   # also link the suite against the pre-fix engine
```

`--revert` is the interesting half. It compiles `orig/qm_orig.c` and
`orig/integrator.c`, links `tests/test_regression.c` against *those* instead of
the fixed engine objects, and reports which checks go red. Four do — the
coupled-dipole residual, dipole rotation covariance, QEq conservation under
the bound, and SCF charge conservation. If a future change makes that number
zero, the regression suite has stopped being evidence of anything and the
mutation in `--revert` needs re-deriving.