# Physics — what the engine models

*S2 Engine v9R4 · C11 · units: Å, fs, eV, AMU, elementary charge*

The engine's physics layers, one by one: what each one is, what
it is good for, and where it stops being trustworthy. The honesty
blocks are not decoration — they say which claims may be cited.

<a id="quantum-orbital-layer"></a>
### Quantum orbital layer
Slater screening (`Z_eff = Z − S`, Slater's own 1930 iron example reproduced),
orbital energies from Slater `n*`, hydrogenic radial profiles `P(r) = r²|R_nl|²`
with exact expectation values, Clementi–Raimondi SCF exponents for spatial
ranges, and all 16 real spherical harmonics — l = 0…3, i.e. 1 s + 3 p + 5 d +
7 f — orthonormal to quadrature precision (worst |∫Y_iY_jdΩ − δ_ij| = 7.2 × 10⁻⁵
over the full 16 × 16 matrix).
Demo 1 includes a per-element computed-vs-NIST ledger instead of a blanket
accuracy claim.

<a id="exact-two-center-sto-integrals-and-roothaan-hartree-fock-scf"></a>
### Exact two-center STO integrals and Roothaan-Hartree-Fock SCF

**One overlap kernel, grid-validated.** `eht_overlap_ao()` in
`src/qm_eht_overlap.c` computes the overlap of two Slater-type AOs given
in **global Cartesian axes** - `(n, zeta, l, m)` per AO, with `m = 1, 2, 3`
for p_x, p_y, p_z - from the exact prolate-spheroidal integrals (Rosen
1931; Roothaan 1951). It is the single place where the Slater-Koster
direction cosines and the s-p sigma sign are applied, and both the EHT
solver and the Roothaan-Hartree-Fock solver assemble their overlap
matrices with it. Before the fourth audit they each had their own copy,
and the two copies **disagreed**:

* the EHT copy carried a spurious minus sign in its "p on the first atom"
  branch, so the overlap matrix of a symmetric molecule was not symmetric
  under its own mirror. For a bent H-C-H it produced
  `S(s_H1, p_x) == S(s_H2, p_x)` where the physics demands
  `S(s_H1, p_x) == -S(s_H2, p_x)`. The EHT SCC then **failed to converge
  for H-C-H at all**, and N2 came out with the wrong MO ordering, a
  HOMO-LUMO gap of 0.0000 eV and a Mayer bond order of 2.053 instead of
  the textbook 3.000;
* the RHF copy resolved each ATOM's valence shell instead of each AO's own
  `(n, l)`, and for C/N/O that shell is 2p - so the s-channel wrapper
  refused every O-H, C-H and N-H pair and the **2s-1s sigma-bonding
  channel of water and methane was exactly zero** in the overlap matrix;
* the RHF copy also never populated the principal quantum number per AO,
  so the VSIP core diagonal always fell back to the first ionisation
  potential and the O 2s diagonal came out at -13.598 eV instead of
  -32.3 eV.

Channels: s-s, s-p sigma, p-p sigma, p-p pi (2p only). The textbook
closed forms (e^{-zR}(1+zR+(zR)^2/3) for 1s-1s, e^{-p}(1+p+2p^2/5+p^3/15)
for 2p-2p pi) are pinned in `tests/test_audit4.c` to 1e-12 and 1e-10
respectively, and the kernel is validated against brute-force 3D Cartesian
grid integration in `tests/test_eht.c`.

**Two-electron integrals.** The two-centre Coulomb kernel is a
**Hubbard-screened Coulomb**,

```
gamma_ab(R) = COULOMB_MD / sqrt( R^2 + (28.8 / Ubar_ab)^2 ),  Ubar_ab = (U_a + U_b)/2
```

with `U = IE - EA` (the chemical hardness). It is continuous with the
on-site term (`gamma_aa -> U_a`) and tends to `COULOMB_MD/R` at large R.
The kernel is **externally validated**: the true two-centre Coulomb integral
<2p_O | 1/r | 1s_H> for the engine's own valence densities at the water
O-H distance (0.9572 A) is **5.7156 eV** (grid integration of the actual
charge densities in atomic units), and the kernel returns 5.7720 eV - a
**1.0 % deviation**. The bare 1/R Coulomb would be 15.0 eV (+162 %).

> The kernel used to be labelled "Mataga-Nishimoto", and it was not one:
> Mataga & Nishimoto (1957) and its DFTB descendant interpolate as
> `1/(R + 0.5/(k_A+k_B))`, and the measured deviations from that form are
> -35 % at 0.5 A and -13 % at 1 A. Substituting the textbook form would
> make the model *worse* (+19 % against the oracle), so the kernel was
> kept and the **label was corrected**. The tree now has ONE screened
> Coulomb kernel and it is the one with an external check.

The exchange is the Wolfsberg-Helmholtz resonance integral
`beta_ab = -K_WH * S * (H_aa + H_bb)/2` with `K_WH = 1.75` (Hoffmann
1963), and `qm_eri_zzdo()` builds the RHF Fock contribution `J - 0.5*K`
from a diagonal density matrix under ZDO.

**The Roothaan-Hartree-Fock SCF solver** (`qm_hf_scf()`) assembles S from
the shared kernel and H_core from the Hoffmann VSIP diagonal, solves
F.C = S.C.eps by Loewdin orthogonalization + a Jacobi eigensolver, fills by
Aufbau, and iterates the density with damped mixing. Fail codes are loud:
-1 bad input, -2 singular overlap, -3 non-convergence, and **-4 the
solution is not bound**. Minimal valence basis: H 1s; C/N/O 2s+2px+2py+2pz.

**Honesty block H1-H5** (in `include/qm.h`): this module is *internally
consistent and numerically sound* - the overlap matrix is exact, the density
is normalised, and the finite-difference gradient of a smooth energy
surface is accurate to 1e-5 eV/A - but it is a **qualitative**
electronic-structure layer, not a predictive one:

* the core Hamiltonian is the EHT diagonal `H_ii = -VSIP`; the exact
  one-electron integrals are implemented and grid-validated but cover 1s-1s
  only and are **not** what this solver uses;
* the off-diagonal core element is the Wolfsberg-Helmholtz *resonance*
  integral with `K = 1.75`, a constant calibrated for the EHT
  band-energy formulation, not for an RHF trace expression;
* the two-electron term is the ZDO approximation;
* the consequence is that the H2 curve has its minimum at **R = 1.300 A**
  against the experimental 0.7414 A, and water in this minimal basis does
  not produce a bound solution at all;
* the validated electronic-structure module in this tree is the **EHT one**
  (Demo 18), which converges, gives the correct H2O Mulliken charges
  (O -0.54 e) and reports its own negative result rather than hiding it.

The demos say all of this out loud, and the not-bound gate refuses to
report a solution whose lowest orbital is unbound.

<a id="wiring-quantum-scf-into-classical-md-use_qm_scf"></a>
### Wiring quantum SCF into classical MD (`use_qm_scf`)
`qm_hf_scf_forces()` computes the SCF electronic energy and its
finite-difference gradient with respect to atomic positions (same FD pattern
as `qm_pauli_forces()`: 6*N extra SCF solves per force call). It is exposed
as an opt-in force term in `forces_calculate()` via `sim->use_qm_scf`:

* **Default off (0):** the record path is completely classical - every
  existing demo and the byte-deterministic record are unaffected.
* **Opt-in (1):** the SCF energy is added to `potential_energy` and reported
  as `E_qm_scf_total`; finite-difference forces are accumulated onto each
  atom alongside LJ, Coulomb, bonded, angle, dihedral, Pauli, dispersion,
  and induction terms. Independent of `use_scf` (QEq charge loop) - both
  can be set simultaneously.
* Supports H/C/N/O atoms within the minimal valence basis (<=128 atoms);
  larger or unsupported systems return rc=-1 and U=0.
* **Fail-closed:** a displaced SCF solve that fails leaves that force
  component at zero and says so. It used to leave the energy at 0.0 - a
  *finite* number - which the finite difference consumed as a real energy
  and turned into a force of order E/2h, about 1.5e6 eV/A. That was an
  integrator-explosion generator, and it was silent.

Because the SCF is O(nao^4) per solve and forces need 6*N solves per step,
this is intended for small quantum-MD test systems, not production-scale
dynamics.

**What the quantum-MD demos actually show.** Demo 13 (H2) works: the H2
relaxes towards the SCF minimum and the force falls from 19 to 1.6 eV/A.
Demo 14 (H2O) does **not**: the SCF is not bound for water in this basis,
so the demo says so, reports the quantum-MD claim as UNRESOLVED, and runs
the classical trajectory as the reference. (Before the audit the same demo
ran the classical trajectory while *printing* a quantum-MD table, and the
water exploded: the identical trajectory without the quantum term is stable
at 250-330 K for 500 steps, and with it the temperature reached the 2000 K
cap at step 300.)

<a id="classical-molecular-dynamics"></a>
### Classical molecular dynamics
* Velocity Verlet (symplectic) with **Berendsen** weak-coupling (steering;
  legacy record paths) and **Andersen** stochastic collisions (rigorously
  canonical; the WHAM free-energy legs).
* Lennard-Jones (Lorentz–Berthelot) + Coulomb with an optional dielectric
  divisor; per-atom LJ overrides so jointly-parameterized models (TIP3P water,
  AMBER ff99 biomolecular atoms) stay internally consistent; harmonic
  restraints enter the force loop like every other term; per-axis PBC and a
  CHARMM-style cutoff switch are implemented (the switch is opt-in).
* Bonded terms: harmonic stretch, harmonic bend, and **dihedrals with exact
  analytic chain-rule gradients**. The finite-difference oracle that proved
  them and a 22-check agreement suite remain in-tree.
* Temperature uses the correct `3N − 3` degrees of freedom, with further
  constrained-DOF accounting. Two RNGs: the historical Knuth-MMIX LCG
  (default — keeps the record byte-stable) and the **reference PCG64**
  (O'Neill 2014, XSL-RR 128/64 output) — see §6.1.
* Minimization: steepest descent (monotonic, clash relief) and FIRE
  (smooth-basin polishing), with a per-atom displacement cap and a
  size-relative divergence floor.

<a id="quantum-enhancement-layer"></a>
### Quantum-enhancement layer
Every term is opt-in, reported as its own ledger entry, and never folded into
the base charges: induced dipoles (first-order, and self-consistent
coupled-dipole solves with Thole-damped coupling and Hellmann–Feynman forces),
overlap Pauli repulsion, Slater–Kirkwood dispersion with Tang–Toennies
damping, and a QEq + dipole SCF charge-equilibration loop (1-pin and 2-pin).
Hard bounds are documented and fail closed: coupled dipole solve ≤ 64 atoms,
QEq ≤ 128, first-order polar ≤ 256.

<a id="semiempirical-scf-layer-eht"></a>
### Semiempirical SCF layer (EHT)
Hoffmann-1963 Extended Hückel with self-consistent charge (Demo 18,
`qm_eht.h`): exact prolate-spheroidal STO overlaps (grid-validated, not
the `qm_overlap` heuristic), VSIP Hamiltonian with Wolfsberg-Helmholtz
coupling, Löwdin-orthogonalized SCF, **Hubbard-screened Coulomb gammas**
(the U values come from the in-tree Mulliken J; the screening constant and
the hardness were both wrong before the fourth audit — the old kernel
screened with the mean ionisation potential and returned 10.09 eV for the
water O–H pair against a grid oracle of 5.7156 eV, a 76 % error; the
current kernel returns 5.7720 eV, a 1.0 % error), Mulliken charges, Mayer
bond orders, HOMO–LUMO gaps computed from the **occupations the solver
actually used** (the previous index-based formula reported the splitting
*inside* a degenerate pair, and N₂ came out with a gap of 0.0000 eV).
Single-zeta minimal basis (H/C/N/O), ≤ 64 AOs, fail closed. An opt-in BO
switch scales harmonic bonds by an overlap-gated Mayer interpolant, so
bonds weaken and let go instead of pulling to infinity (E6–E12 inventory
in `qm_eht.h`).

<a id="biopolymer-condensation-chemistry"></a>
### Biopolymer condensation chemistry
Real chemistry, not decoration: leaving groups are genuinely removed following
the PDB Chemical Component Dictionary's own flags. Nucleic acids: the five
bases from PDB CCD ideal coordinates with Aduri et al. 2007 RESP charges and
AMBER ff99 LJ types, tautomer-corrected cytosine, deoxyribose, and a T-p-A
dinucleotide with a real phosphodiester bridge (net −1 e). Proteins: glycine,
alanine, a Gly–Ala dipeptide with a real 1.33 Å peptide bond, a poly-alanine
chain builder, and an α-helix whose i,i+4 backbone hydrogen bond emerges from
Coulomb+LJ alone.

<a id="electrophysiology"></a>
### Electrophysiology
The 4-variable Hodgkin–Huxley 1952 model, RK4-integrated (Euler would need an
impractically small step during the sodium upstroke), with the removable rate
singularities handled by Taylor expansions. Resting gating values reproduce a
published worked example to four decimals.
<a id="datastream"></a>
### Datastream
A text, self-describing, line-oriented format (`DATASTREAM_SPEC.md`, schema 1)
with a SHA-256 integrity seal over the payload, written by Demo 12 as
`kcsa.cvmds`. Claims carry provenance tags (`computed`, `computed-jc2008-params`,
`Marcus1991-TATB`, `1K4C-LINK`, …) and a validator can check external-citation
claims against reference values. The writer fails closed: truncated or
unsealed output is never written.

Two properties of the writer are worth stating because they shaped how
Demo 12c reached the file at all. `ds_open` buffers the whole payload and
truncates at `ds_close` — there is no append path and no reopen — and Demo 12
runs *before* Demos 12b/12c, so the writer cannot reach a value Demo 12c has
not computed yet. Rather than hoist Demo 12's hundred-odd locals to file
scope, Demo 12c's measurement is factored into `kcsa_ionic_measure`, a pure
deterministic function of the deposited coordinates called by both the
printer and the writer. The printed table and the sealed claims therefore
come from one code path and cannot drift apart in principle; 22 checks in
`test_external` pin the values so a loss of determinism breaks a gate rather
than producing two plausible, disagreeing numbers.

A claim that is *not* computable is sealed explicitly rather than omitted:
`kcsa.entropic.computed` is recorded as `0 kcal/mol` tagged
`UNRESOLVED-no-sidechain-torsions`, so the gap sits in the sealed file
beside its published target (`kcsa.entropic.target`, 5.5 kcal/mol from
Dixit & Asthagiri) rather than only in prose.

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
