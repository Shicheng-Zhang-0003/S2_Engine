# KcsA engine-path plug-in: what the numbers actually do

*Diagnostic record, not a claim of a working selectivity model. The sealed
engine was not modified. Reproduce with
`side_fixes/probes/run_kcsa_engine_probe.sh`.*

## Why this exists

`side_fix_kcsa_stabilized.c` is a separate sampler. The question this probe
answers is narrower: **if the side-fix parameters are compiled into the
engine's own `qm_kcsa_binding_series()`, does the engine's KcsA leg hold the
site, and do its numbers become usable?**

The probe copies `src/qm_fluctuation.c` to a temp directory, makes three
compile-time overrides (`QM_KCSA_EPS`, `QM_KCSA_DT`, `QM_KCSA_GAMMA`), and
links the copy against the engine's own objects. No sealed file changes; the
temp copy is deleted on exit.

Protocol: real 1K4C TVGYG tetramer + 104 side-chain hydrogens, S2 site,
300 K, K+ vs Na+, 200 equilibration steps + N sampled steps, Langevin,
C-alpha tether 8 eV/A^2. The engine defaults are `eps=1`, `dt=1.0 fs`,
`gamma=0.001`; the side fix uses `eps=4`, `dt=0.5 fs`, `gamma=0.02`.

Sign convention is the engine's (`include/qm_fluctuation.h`):

```
ddG = 0.5 * beta * (var_K - var_Na),   NEGATIVE = K+ favoured.
```

## Baseline: the engine's collapse reproduces

| ion | CN | <O-O> | <Ebind> | var | gate |
|---|---:|---:|---:|---:|---|
| K+ | 5.000 | 4.118 A | -4.2729 eV | 0.4753 eV^2 | refuses |
| Na+ | 4.565 | 4.138 A | -3.7700 eV | 0.8116 eV^2 | refuses |

`qm_kcsa_binding_series()` returns -1 and `ok=0`, so no variance is offered.
This is the engine's documented collapse, reproduced on the engine path.

## No single side-fix parameter holds the site

200 fs, tether 8, engine sampler:

| change | K CN | Na CN | ddG (engine sign) | gate |
|---|---:|---:|---:|---|
| eps=4 only | 6.365 | 5.510 | - | refuses |
| dt=0.5 only | 6.790 | 5.210 | - | refuses |
| gamma=0.02 only | 5.810 | 6.665 | - | refuses |
| eps=4 + dt=0.5 | 7.865 | 7.675 | -0.86 kcal/mol | pass |
| eps=4 + gamma=0.02 | 8.000 | 8.000 | **-5.55 kcal/mol** | pass |
| dt=0.5 + gamma=0.02 | 7.085 | 7.380 | -14.50 kcal/mol | pass |
| full set | 8.000 | 8.000 | -14.26 kcal/mol | pass |

The stabilization is a four-parameter combination, including the timestep
and friction. It is not a derived physical mechanism, and the result is
strongly parameter-dependent: the same 200 fs window gives -0.86, -5.55 and
-14.5 kcal/mol for three different subsets.

## "Held" does not mean crystallographic

Full set, engine sampler:

| window | K CN | Na CN | K <O-O> | Na <O-O> | ddG (engine sign) | gate |
|---:|---:|---:|---:|---:|---:|---|
| 200 fs | 8.000 | 8.000 | 3.965 A | 3.930 A | -14.26 kcal/mol | pass |
| 300 fs | 8.000 | 8.000 | 3.976 A | 4.013 A | -7.66 kcal/mol | pass |
| 500 fs | 7.862 | 7.960 | 4.164 A | 4.147 A | **+14.77 kcal/mol** | pass |
| 900 fs | 6.973 | 7.016 | 4.287 A | 4.243 A | +5.29 kcal/mol | K+ fails |

The deposited O-O diagonals are 4.5-5.5 A. At 200 fs the shell is at
3.93-3.97 A: the coordination number is 8, but the geometry is compressed by
roughly 0.5-1.5 A. The engine's gate checks CN only, so it passes this state.
The side fix's test was renamed from "remains crystallographic" to "stays
above the collapse floor" for the same reason.

## The entropic term is not converged

The full set gives -14.26 kcal/mol at 200 fs, -7.66 at 300 fs, and
**+14.77 at 500 fs** - a sign flip. At 900 fs one site is below the CN >= 7
gate, so the number is not usable. The apparent 200 fs agreement of the
`eps+gamma` subset with the literature 5-6 kcal/mol target is not robust:
at 500 fs both sites fail and the same subset gives +33.9 kcal/mol.

This is a direct, engine-path confirmation of the existing verdict: the
KcsA entropic term is **UNRESOLVED**. No selectivity free energy should be
quoted from these windows.

## Relationship to the side-fix module

The side fix's own 200 fs run reports CN 8.00/8.00 and
`<O-O> = 3.97/3.96 A`, matching the engine-path full set. It also had a
sign bug: it returned `0.5*beta*(var_Na - var_K)`, the opposite of the
engine convention. That is fixed; it now calls `qm_fluctuation_term_eV()`
and the 200 fs value is **-0.4901 +/- 0.0609 eV** (K+ favoured), still
UNRESOLVED. See `side_fix_kcsa_stabilized.h`, K2-K4.

## Reproduce

```sh
# 200 fs matrix (8 configurations)
make -C side_fixes probe-kcsa

# add the 300/500/900 fs windows
./side_fixes/probes/run_kcsa_engine_probe.sh --long
```

The probe is diagnostic only. It is not part of `make test`, not part of the
record, and does not modify `src/`, `include/`, `tests/` or `makefile`.
