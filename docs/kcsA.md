# KcsA — the potassium channel program

*S2 Engine v9R5 · the selectivity filter, the ionic layer, and the
entropic term, and exactly how far each of them is trustworthy.*

The filter is the deposited 1K4C chain-C structure, not a cage.
The ionic layer is the Born self-energy against the measured
hydration scale. The entropic term is implemented and **measured
not to converge** — and the site-integrity gate now refuses to
return a variance for a collapsed site, so the negative result
cannot be papered over.

<a id="kcsa-ion-channel-program"></a>
### KcsA ion-channel program
Two tracks, deliberately:

* **Demo 12 — legacy constructed cage** (kept as the honest vacuum baseline):
  four carbonyl oxygens in 4-fold symmetry, single ion on-axis, fixed-radius
  sites, a radius scan, and a two-ring antiprism with the ring z-separation
  sourced from 1K4C (3.084 Å).
* **Demo 12b — the real filter**: the deposited PDB 1K4C chain C TVGYG
  structure — 164 atoms, net charge 0.000 e, C4-symmetric, CN = 8 at all four
  sites, Thr75 OG1 as the innermost gate, GLY79's carbonyl correctly excluded
  (4.82 Å off-axis). It is rigid by construction; the module documents why a
  flexible-filter calculation collapses without a real protein force field and
  a solvation model.

Each ion is pushed through every leg and the results are printed side by side
in a forensic table rather than summarized into one number. Current record,
in that table's sign convention (`+ = Na⁺ favored`):

| Leg | Value | Direction |
|---|---|---|
| vacuum site legs (point charges, same cage) | +0.0000 eV | identical — a point-charge ion carries no size |
| JC-ion leg (size-dependent LJ) | dU(JC, K−Na) = +0.3306 eV | Na⁺ |
| computed exchange, rigid filter | −0.5772 eV (deterministic) | K⁺ |
| computed exchange, relaxed filter | −0.3215 ± 0.9382 eV | UNDECIDED (within spread) |
| two-ion exchange | −0.2350 eV | KK selective |
| SCF-polar U(z) barrier, d(K−Na) | −0.3671 eV | Na⁺ pays **more** |
| polar-WHAM barrier gap | +0.2017 ± 0.0044 eV | UNDECIDED (within spread) |
| fixed-charge WHAM gap | +0.1558 eV | Na⁺ |
| knock-on pair / conductive / landscape | +1.5807 / +4.8238 / −1.6244 eV | mixed |
| coordination probe, 8 vs 6 ligands | dE(8−6) = +9.1518 eV (K⁺), +9.0697 eV (Na⁺) | losing the 8th site costs ~9 eV either way |
| measured dehydration (Marcus 1991 TATB) | K⁺ enters 0.7255 eV cheaper | K⁺ |

> The `±` on the relaxed-exchange leg is the spread across three
> independent seeds, **not** a standard error, and the record says so at
> the point of use. The two-ion exchange is a sum of rigid-cage and
> six-water-cluster model terms, not a free energy: its sign is a model
> statement, and its magnitude carries no uncertainty bar because there is
> no solvent model to put one on. Compare it with the experimental
> ~−0.18 eV for scale, never as a prediction.

The legs disagree, and that disagreement is the result. **No quantitative
selectivity free energy is claimed.** Only the measured dehydration free
energy is a robust K⁺ advantage in this model; the two-leg sum is compared
with the experimental −0.179 eV reference as a scale check, never as a
prediction.

The disagreement is not a rounding difference: the rigid-filter exchange
leg says K⁺ (−0.5772 eV) and the SCF-polar pore-axis barrier says the
opposite sign (−0.3671 eV, i.e. Na⁺ pays more to cross). Two legs of the
same model disagreeing about the sign of the selectivity is exactly the
result that forbids a selectivity number — and it is reported rather than
averaged away. The umbrella sampler prints its measured autocorrelation time τ
and effective sample size with every barrier.

<a id="demo-12c-—-what-the-ionic-layer-actually-shows"></a>
#### Demo 12c — what the ionic layer actually shows

Demo 12c exists because the forensic table above is not the whole story, and
because two questions the tree previously could not ask at all now have
answers.

**Where the desolvation leg comes from.** Until now the hydration free energy
was a transcribed table with nothing able to check it or vary it. There is now
a Born self-energy to compare it against — and it disagrees with the table by
+18.7 kJ/mol, which is *supposed* to happen: the Born expression is a
linear-response continuum result with no room for the cavity and dispersion
terms the real ion–water interaction has, and for monovalent cations it
under-predicts the differential by about that much.

| | ΔG(Na⁺−K⁺), kJ/mol |
|---|---|
| Born, r_Born 1.80 / 2.08 Å | −51.3 |
| measured, Marcus 1991 TATB | −70.0 |
| Born's deficit | **+18.7** |

The Born radius is not unique, and that is reported rather than hidden:
refitting it to the measured hydration free energies shifts the differential
by a further −20.4 kJ/mol. That is the honest size of the uncertainty
attributable to one unsolvable choice.

**Whether the filter actually excludes Na⁺.** This is a geometric question
with a geometric answer, and it is the sharpest thing this engine knows about
KcsA. On the pore axis the mean ion–oxygen distance over the eight
coordinating oxygens is **2.777 Å**, and that is the *shortest the deposited
geometry permits* — an ion cannot get closer than the oxygens are. Measured
against each ion's own preferred contact (Shannon VIII radius + 1.40 Å
carbonyl oxygen):

| ion | prefers | achieves | mismatch |
|---|---|---|---|
| K⁺ | 2.910 Å | 2.777 Å | **−0.133 Å** (squeezed inside its preference) |
| Na⁺ | 2.580 Å | 2.777 Å | **+0.197 Å** (held outside it) |

The cage is built for K⁺, so the larger ion is squeezed 0.13 Å inside what it
wants while the smaller one is held 0.20 Å outside what it wants and cannot
get closer anywhere on the axis. **That asymmetry is the exclusion mechanism,
and it is a property of the deposited coordinates alone** — no force field, no
dielectric, no fitted parameter.

**Why that still does not produce a selectivity number here.** Both ions relax
to the *same* on-axis minimum at the same 2.777 Å, and their ion-induced shell
strains agree to **1.2%** — a ratio of 1.0122 that is exactly invariant across
the record's 50× stiffness scan (1 → 50 eV/Å²), because strain goes as 1/k
and the ratio does not. The reason is quantitative: the ion–oxygen
Lennard-Jones wall here is 0.05 kcal/mol, roughly four orders of magnitude
below the Coulomb pull at 2.8 Å, so nothing in this potential rewards a
cation for having the right radius.

With **LJ + Coulomb alone** the two ions' strain energies differ by **1%** —
they are literally the same geometry taken twice. Turning on the exchange
(overlap-Pauli) and dispersion terms, which are size-sensitive by
construction, moves that to **1.2%**. The quantum layer changes the
answer; it just cannot change it enough to matter here.

So the honest conclusion is a sharp negative one: **in this model the entire
K⁺/Na⁺ discrimination is the tabulated hydration free energy.** The geometric
mismatch is real and measurable, and it still does not convert into a
binding-energy difference.

**The entropic term, and why it is reported unresolved.** Dixit & Asthagiri
(*J. Chem. Theory Comput.* **5**, 1037, 2009; with Asthagiri, Roux & Cross,
*Biophys. J.* **95**, 1269, 2008) show that the *mean* binding energy slightly
**favours Na⁺** — Na⁺ is the better fit for an eight-carbonyl cage — and that
KcsA selectivity comes from the **fluctuation** term, `Δ(βσ²/2) ≈ 5–6
kcal/mol`, driven by carbonyl–carbonyl repulsion strain of about 7.3 kcal/mol.
That is the whole mechanism, and it is a term a rigid-cage model cannot
produce.

This release adds the machinery (`qm_fluctuation`): it runs the deposited
filter as a breathing peptide with the ion inside and measures the variance of
the ion's binding energy over an explicit Langevin NVT trajectory, with
streaming (Welford) moments, standard errors on every variance, a
stiffness-independent linear-elastic strain model, and a **site-integrity
gate** that fails closed if the filter does not hold its 8-fold site.

**It does not yet converge, and the reason was structural.** The engine had no
side-chain torsion terms, and four of the eight filter oxygens belong to
Thr75's side chain. Without torsions the hydroxyl leaves the shell, the
coordination number falls from 8 to 3–5, and the O–O diagonals close from the
crystallographic ~4.5–5.5 Å to ~3.5 Å — measured at every temperature from
50 K to 300 K and every Cα restraint stiffness from 0.5 to 16 eV/Å². No
thermostat and no restraint recovered it, because the missing term was the
physics and not the parameters.

The machinery for it is in place: `kcsa_add_dihedrals()` builds backbone
(phi/psi/omega) and side-chain (chi1) torsions restrained toward the deposited
geometry, and the fluctuation measurement holds the ion on the pore axis at
the deposited site (the symmetric-cage equilibrium is unstable — any radial
displacement draws the cation toward the nearest oxygen because the ion-O LJ
wall of 0.05 kcal/mol is four orders of magnitude too weak to balance the
Coulomb 1/r pull).

**And the site is still lost.** This was measured, not assumed, and the
fourth audit corrected a claim in this document that said otherwise:

| trajectory time | mean CN | mean O–O |
|---|---|---|
| 10–50 fs | **8.000** | 4.17–4.28 Å |
| 70 fs | 7.600 | 4.28 Å |
| 100 fs | 6.233 | 4.19 Å |
| 200 fs | 4.000 | 3.28 Å |

against a deposited coordination of 8 and crystallographic O–O diagonals of
4.5–5.5 Å. The filter strain energy in the collapsed state is 6–9 eV, so
what the sampler would measure is the variance of a **collapsing cage**, not
the breathing fluctuation the entropic term needs. Stiffening the angle term
does not fix it (CN 3.2 at 1.0 eV/rad², 2.7 at 2.2 eV/rad² — worse), so this
is not a missing-stiffness parameter but the limit of the reduced protein
model: holding a real 8-coordinate selectivity filter needs a real protein
force field and explicit solvent.

The consequence is enforced in code, not just in prose: the **site-integrity
gate now refuses** to return a binding-energy variance for a site that was
not maintained (`qm_kcsa_binding_series()` returns −1 and clears `ok`), so
the variance of a collapsed cage can never be quoted as the entropic term.
The leg is reported **UNRESOLVED**, and the record says which of the two
reasons is in force.

**The table above is reproduced by a gate**, not asserted in prose:
`tests/test_audit4.c` runs the sampler twice — 20 equilibration steps
(the site is still intact, CN = 8.000) and 200 steps (the shell has
closed, CN = 4.000, O–O = 3.279 Å) — and fails if either endpoint moves.
If a future change makes the filter hold, that check goes red and this
table has to be updated with it.

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
