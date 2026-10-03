# AUDIT FINDINGS — S2 Engine v9R4 / V0.9RC4

Severity: **CRITICAL** breaks physics silently · **MAJOR** wrong result, bounded ·
**MINOR** hygiene/doc · **DOC** documentation states a falsehood.

---

## C1 — CRITICAL — `qm_solve_dipoles` does not solve its documented equation

`src/qm.c:1002-1123`. The function documents a **3N × 3N** system
`(I − A) mu = alpha (*) E0` and the release note calls it "an exact 128-bit
Gaussian elimination … no iteration count, no convergence test, no failure
mode". It solves neither.

Four independent structural defects, any one of which is fatal:

| id | defect | evidence |
|----|--------|----------|
| D-1 | **Wrong system size.** Elimination loops run `c < N`, `r < N`; `M` is an N×N array. A system of N atoms has 3N dipole unknowns. | `for (int c = 0; c < N; c++)` at qm.c:1077, 1089–1102 |
| D-2 | **Wrong block layout.** A 3×3 tensor block at row stride 3N / column stride 3 needs offsets `{0,1,2, N,N+1,N+2, 2N,2N+1,2N+2}`. The code writes **3 of 9** — (xx),(xy),(xz) — at `+0, +N, +2N`, which are rows i, i+1, i+2 of an N-wide matrix. yy, yz, zz, zx, zy are never written. | qm.c:1070-1072 |
| D-3 | **Uninitialised reads.** Only `M[0 … N*N−1]` is zeroed, but writes reach index `3N² − N − 1`. For **every N>1** the off-diagonal writes accumulate onto stale thread-local memory from a previous call. | qm.c:1029 zeroing vs 1070-1072 writes |
| D-4 | **rhs permutation mismatch.** Filled `rhs[3*i+c]`, read `rhs[comp*N+r]`. Equal only when N=1. | qm.c:1021 vs 1090, 1107 |

Measured with a corrected harness (A06; the engine's own objects linked
against the HEAD version and against the fixed version, same binary
otherwise):

| | worst residual of the documented equation | vs independent correct solve |
|---|---|---|
| **HEAD** | **2.25** | **1.00** |
| after fix | 6.5e-14 | 7.4e-13 |

> **Correction to my own first measurement.** An earlier run of this harness
> reported residuals of 5.4 and 15.6. Those were *my* bug: the accumulation
> loop was wrapped in a spurious `for (b=0;b<3;b++)` that applied each tensor
> term three times. The residual was inflated by exactly that factor. The
> table above is from the corrected harness, and the corrected numbers are the
> ones that count. The defect is real and large — a relative difference of
> **1.0** from a correct solve means the function returned an entirely
> different vector, not a slightly-wrong one — but the "1.4e+01" figure was
> never the engine's.

**Why every shipped test passes.** The FD force test verifies the force is
the gradient of *whatever* `U = −½ Σ mu·E0` the function returns — which is
true by construction. The continuity test passes because stale thread-local
memory is a *deterministic* function of the previous call. "Dipole solver
converges" only checks `rc == 0`, and the N×N matrix actually solved is
nonsingular. No shipped demo uses `use_pol_scf` with N other than 9 or 10.

**Blast radius:** every `use_pol_scf = 1` path — Demo 12's SCF-polar U(z),
v2/SCF-polar 8-fold legs, the polar WHAM windows, the knock-on pair runs.
This is the engine's headline v4 quantum feature.

---

## C2 — MAJOR — the ±2 e clamp silently breaks charge conservation

`src/qm.c:569-570` (HEAD) clamps `out_q[i]` to ±2 e **after** the augmented
solve. The augmented system carries `Σ q = total_q` as one of its rows; a
post-hoc clamp violates it. `qm_scf_run` has the same defect compounded by
0.3/0.7 under-relaxation — and mixing is itself enough to break the
constraint, because `qprev` starts from the atoms' *stored* charges, which
need not sum to `target` at all.

| | `qm_qeq` | `qm_scf_charges` |
|---|---|---|
| **HEAD** worst \|Σq\| | **8.0 e** (233/4000 clusters broken) | **8.3 e** |
| after fix | 1.1e-15 | 8.9e-16 |

The release note claims "total charge still conserved to 1e-9". True only
when no atom reaches the clamp and no iteration runs — which is exactly the
case the shipped suite tests.

**Fix.** The bound is now enforced by contracting the *deviation* from the
constraint-consistent mean, and the constraint is re-imposed after
under-relaxation by a *uniform* shift. Both preserve the group-charge
ordering bit-for-bit (a uniform shift leaves every `q_i − q_j` untouched; a
contraction about the mean is monotone), which is the thing the
under-relaxation exists to protect. The physics claim — beyond ±2 e no partial
charge is meaningful — is unchanged.

**Not a defect, recorded so it is not re-litigated:** QEq's
electronegativity ordering *does* invert for an **unbonded** C/O pair closer
than ~2.3 Å, because the 1/r coupling then exceeds the hardness `J` and the
solve enters the charge-transfer mode. Measured:

```
bonded   r=1.43 A   q(C)=+0.115  q(O)=-0.115   OK
bonded   r=1.23 A   q(C)=+0.115  q(O)=-0.115   OK
unbonded r=2.80 A   q(C)=+1.614  q(O)=-1.614   OK
unbonded r=2.20 A   q(C)=-0.635  q(O)=+0.635   inverted
unbonded r=1.50 A   q(C)=-0.157  q(O)=+0.157   inverted
```

That is a property of the functional on an unphysical input (two atoms 1.5 Å
apart with no bond is not a molecule), not a coding error — which is why the
shipped suite tests the bonded case and says so.

---

## M1 — MINOR — `PLANCK_HBAR` is a 10-digit truncation of an exact value

`include/constants.h:13` literal `1.054571817e-34` vs `h/(2π) =
1.0545718176461563913e-34`. Relative error **6.13e-10**. Unused, but it is
the only remaining hand-typed reciprocal of an exact constant — precisely the
pattern audit fixes C1–C4 removed everywhere else.

## M2 — MINOR — `HARTREE_TO_EV` literal vs derived `EV_TO_HARTREE`

`constants.h:52-53`: one half of the pair is derived, the other typed. Product
deviates from 1 by 1.3e-15. The comment claims agreement "by construction".

## M3 — MINOR — minimiser tolerance 0 → division by zero

`integrator_minimize` / `_frozen`: `scale = effective_step / max_force` with
no guard when `max_force == 0`.

**Characterised precisely (A12/A12b).** With `force_tolerance == 0` and an
already-converged system, `max_force < force_tolerance` is `0 < 0` = **false**,
so the loop body *does* run: `scale = 0.01/0 = inf`, every displacement becomes
`0 * inf = NaN`, and the `isnan(E_new)` rollback branch rescues it one
iteration later. So the **returned energy is correct** and the observable
behaviour is benign — but the function exits after one wasted iteration rather
than recognising convergence, and it does so through the NaN guard rather than
through a test. `integrator_fire` has carried the explicit `max_force <
1e-300` guard all along, which is an asymmetry with no physical justification.
Fixed by giving both steepest-descent variants the same guard.

I am recording this as MINOR rather than MAJOR specifically because I could
not make it produce a wrong answer. The *latent* form (a caller with tol==0 on a
system whose forces are all zero) is reachable, and the save is real, but the
measured outcome is a correct answer by luck of ordering.

## M4 — MINOR — `verify_scripts.sh`'s `ok()` cannot fail

`verify_scripts.sh:15` — `ok(){ printf ...; }` never sets `fail=1`. All 13
"kcsa cage wiring (truth fixes)" checks are decorative: they `grep -Fc` a
pattern and print the count with a hand-written `(>=N)` in the label. Nothing
compares against the threshold.

## M5 — MINOR — `make test` swallows a regression failure

`makefile:97` — `./build/test_regression | tail -3`. The pipeline's exit
status is `tail`'s. Verified: a deliberately failing suite returns 0 through
this recipe. `verify_scripts.sh` catches it separately, so this is redundancy
loss, not a hole.

## M6 — MINOR — `run` overwrites the sealed record with no integrity check

`run` is `./carbonsim > output.txt`. It neither verifies nor regenerates
`CURRENT_BASELINE_SHA.txt`.

## M7 — MINOR — tracked `.d` files contradict the commit that "untracked" them

`f8a3ed1` is titled "untrack the auto-generated .d dependency files" but only
added `.gitignore` lines; no `git rm --cached`. All 11 `.d` files plus the
stale `output.asan.txt` remain tracked. `.gitignore` does not affect tracked
files. Confirmed the stale-`.d` footgun is real: injecting a nonexistent
header name makes `make` fail hard with "No rule to make target" (exit 2).

## M8 — MINOR — stale artifacts that break the harness

* `output.asan.txt` (tracked) is `58bb9f69…` — the **previous pre-audit**
  release. `s01_verify_record.sh` diffs against it, so **the s01 harness fails
  against the current tree by construction.**
* Parent-level `kcsa.cvmds` has the older header schema
  (`source-hash: record-tree-v9R4-fixed-no-vcs`, no version fields).

---

## D1 — DOC — readme states a digest nothing produces

`readme.md:461`: "both SHA-256 `8e8836a04bb3…`". The record is `67f9f47e…`.
No build of this tree produces `8e8836a0…`.

## D2 — DOC — `kcsa_filter.h` ion sigmas are the *pre-fix* wrong values

`kcsa_filter.h:148` states σ = 1.888 Å (K⁺) / 1.246 Å (Na⁺). The code computes
**1.993 / 1.352**. The documented numbers are exactly what the *old, wrong*
oxygen σ 3.06615 produced — the very constant commit `ef1e8d6` corrected. A
stale header comment left by the fix that introduced the shared constant.

## D3 — DOC — `kcsa_filter.c` describes charge neutralisation that does not exist

Lines 36–55 claim residues are neutral "by a UNIFORM offset added to its
carbons, **computed at build time** from the other atoms in the residue". There
is no such code; `q` is read straight from the table. The *numbers* are
balanced (tetramer nets −1e-16 e, exactly 5 ligand O per subunit) but the
described *mechanism* was never implemented.

## D4 — DOC — `aminoacids.c` claims a sum-to-zero assertion that is absent

Line 17: "Charges sum to exactly zero for each neutral free amino acid.
Verified by sum-to-zero assertion below." There is no assertion, and
**glycine sums to +0.257 e** (alanine is −0.0000). Demo 10's printed
−0.0430 e is honest; the header claim is false.

## D5 — DOC — readme/release-note quote a 4× different cage separation

`readme.md:2028` and `release_note_v9R4.md:186` say the rigid cage separates
the ions by "0.05 eV". That is the **sum over 4 sites** (−0.0521). The
program's own output prints **0.0130 eV**, the per-site mean. The prose and
the record disagree by 4×.

## D6 — DOC — `kcsa_filter.h:173` documents a parameter that does not exist

"`frozen[]` must be sized to sim->num_atoms and is allocated internally".
`kcsa_site_binding` has no such parameter — leftover from the pre-rewrite
signature.

## D7 — DOC — "13-demo suite" is 14 demos

`readme.md:131,281,461` and release note say 13; `output.txt` has 14 banners
(1–12, 12b, 17).

---

## V1 — MINOR — a meaningless number is presented as a leg

`output.txt:900` forensic table, `knock-on landscape` row: **655312 eV**. That
is the raw sum including the CLASH-flagged rows the row above it explicitly
excludes.

## V2 — MINOR — `exchange2 = −66.9 eV` has no sanity caveat

Printed next to `expt −0.179 eV` with no scale warning, and it inherits the
−215 eV knock-on collapse artifact that the surrounding text flags elsewhere.

---

## VERIFIED CORRECT (independent, not self-consistency)

* LJ + Coulomb force coefficients == `(dV/dr)/r` **symbolically** (sympy, exact identity).
* LJ well minimum lands on `2^(1/6)·σ` for H, C, K, AMBER-O to 1e-16.
* Dihedral gradient == **sympy analytic derivative** to 7.6e-14 on 20 random geometries — not a finite difference.
* Dihedral energy == symbolic to 4.6e-15.
* Dihedral translation invariance 2.5e-16, **rotational invariance 1.1e-15** (not covered by the shipped suite).
* Dihedral handedness real: 345/345 mirror reflections flip the sign.
* Thole `f` and `df/dr` exact against sympy; f(2.35)=0.803, f(2.75)=0.926.
* `U = −½Σmu·E0` ≡ `−½·C·Σα|E|²` to 4.6e-16.
* QEq: mu equalisation and charge conservation exact **when no atom is clamped**.
* 34/34 deposited KcsA heavy atoms match 1K4C chain C to **0.000000 Å**.
* Four K⁺ z-values match deposited HETATM exactly; C4 gives CN=8 at every site.
* Constants: all 2018-CODATA literals correct to their publication precision.

---

## C1b — CRITICAL — fifth defect in the same function: a sign error in the off-diagonal

Found while validating the rewrite above, and confirmed by the fact that after
fixing the layout the residual came out as the *negative* of the right-hand
side. Even with a correct 3N×3N layout,
the assembly uses `-=` where the documented operator requires `+=`.

Physical operator:      A_ij = kij * [3 d̂d̂ᵀ − I]
Required matrix entry:  (I − A)_ij = −A_ij = +kij * [I − 3 d̂d̂ᵀ]

The code writes:

    M[...] -= kij * (((a==b) ? 1.0 : 0.0) - 3.0*u[a]*u[b]);   /* = −kij[I−3d̂d̂ᵀ] */

which **adds** A_ij rather than subtracting it, so the assembled system is
I + A instead of I − A. The comment directly above the line reads
`-A_ij = kij * (I - 3 dhat dhat^T)`, i.e. it states the correct entry and then
applies the wrong sign to it.

Measured on a 3-atom K⁺/O/O system with the layout fixed but the sign still
wrong (so this is a clean, isolated demonstration of the sign alone):

```
alpha_0*E0_0/C = (+4.5565, +2.6521, -0.1700)
mu_0 − A mu_0  = (-4.5807, -2.7353, +0.1813)
```

i.e. the residual is the *negative* of the right-hand side — the exact
signature of the whole off-diagonal block having the wrong sign.

**Consequence:** the corrected solve still had to be re-signed, and the
original function was doubly wrong: wrong system size, wrong block layout,
uninitialised reads, mismatched rhs permutation, *and* a flipped off-diagonal
sign. Five independent defects in one function.

---

# EFFECT OF FIXING C1 ON THE RECORD

The dipole solver is the engine's v4 quantum feature and every `use_pol_scf`
demo depends on it. Regenerating the record moves 202 lines. The changes that
matter scientifically:

| quantity | HEAD | after fix | comment |
|---|---|---|---|
| knock-on pair ΔU (KK−NaNa) | **−65.08 eV** | **+1.54 eV** | was a collapse artefact |
| `exchange2` ΔG | **−66.92 eV** | **−0.271 eV** | was ~370× the experimental ΔG |
| knock-on pair E (KK) | −215.46 eV | −7.02 eV | −215 eV is not a binding energy |
| polar-WHAM barrier gap | −0.118 ± 0.145 | −0.042 ± 0.233 | same sign, 3× wider |
| relaxed-filter exchange | −0.500 ± 0.796 | −0.156 ± 0.904 | still consistent with expt |
| rigid-filter exchange | −0.5865 | −0.5762 | ~unchanged (no polarization term) |

The two knock-on numbers are the headline. A −215 eV "binding energy" for two
cations in a filter was never physics — it is the induced-dipole catastrophe
the solver's own hard core was only partly masking, expressed through a broken
linear solve. With the solve correct, the pair energy is −7 eV, which is
comparable to the single-ion binding energies the same demo reports, and the
2-ion exchange ΔG lands within a factor of 1.5 of experiment instead of a
factor of 370.

**Honest note on what is still not right.** The knock-on landscape barrier for
KK moves to 15.98 eV, which is *worse* looking than the old 4.24 eV but is not
an artefact of the solver: the landscape scan has one surviving minimum at
−19.63 eV (z = −3.00) that drags the max up. That is a bin-statistics
artefact in a single-point scan, unrelated to C1, and it is now the largest
number in that section — which is exactly why the section needs the caveat it
never had (finding V2).
