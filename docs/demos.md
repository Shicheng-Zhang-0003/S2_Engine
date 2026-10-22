# The demos — what each one claims

*S2 Engine v9R4 · the 20-demo record sequence.*

Run `./carbonsim` for the record, or `test demo <n>` inside `s2tui`
for the fast live smoke check. Every headline below is printed by
the demo itself, and the record is the authority — not this table.

<a id="the-demos-20"></a>
## The demos (20)

12b is the real-filter twin of 12, 17 is the duplex, 18 is EHT, and 13–16
are the quantum-layer demos added in v9R4.

| # | Demo | Current headline |
|---|---|---|
| 1 | Quantum orbital structure | H/C/N/O Slater orbitals, exact H-like expectations, Clementi–Raimondi exponents, per-element NIST ledger |
| 2 | H–H bond curve | covalent and van der Waals regimes as two different physics |
| 3 | H₂O molecular dynamics | Berendsen NVT at 300 K, stable O–H geometry |
| 4 | Cyclic (H₂O)₃ | hydrogen-bonded ring emerges from Coulomb+LJ alone — HELD |
| 5 | CH₄ | tetrahedral 109.47°, zero-strain placement |
| 6 | Five nucleobases | PDB CCD geometry, RESP charges, planarity, tautomer-corrected cytosine |
| 7 | Watson–Crick pairing | G–C > A–U, both bound; magnitudes qualitative (see limitations) |
| 8 | T-p-A dinucleotide | real sugar–phosphate backbone, net charge −1.00 e |
| 9 | Hodgkin–Huxley neuron | action potential, +40 mV peak, 4 spikes / 50 ms |
| 10 | Gly–Ala dipeptide | genuine 1.33 Å peptide bond |
| 11 | α-helix | i,i+4 H-bond emerges: H···O 2.1204 Å, N···O 3.0982 Å |
| 12 | KcsA legacy cage | full selectivity leg program; vacuum sites favor Na⁺, exchange favors K⁺ |
| 12b | Real KcsA filter | deposited TVGYG, 164 atoms, C4 symmetry, CN = 8 at every site |
| 12c | KcsA ionic layer | Born self-energy vs the measured table; filter contact mismatch; entropic term's status |
| 13 | H₂ quantum MD | the Roothaan–Hartree–Fock layer drives an H₂ trajectory; the force falls 19 → 1.6 eV/Å as the molecule relaxes to the SCF minimum |
| 14 | H₂O quantum MD | **UNRESOLVED** — the SCF is not bound for water in this minimal basis, so the demo says so and runs the classical reference instead of printing a quantum-MD table it cannot support |
| 15 | H₂ dissociation curve | the SCF curve **has** an interior minimum (R = 1.300 Å, E = −12.690 eV, against the experimental 0.7414 Å) after the missing nuclear–nuclear repulsion was added; the deviation is printed, not hidden |
| 16 | H₂ geometry optimization | the SCF-force optimiser converges to the same minimum (R = 1.300 Å) and prints the +0.559 Å deviation from the experimental bond length |
| 17 | DNA duplex | G–C / A–T stack with B-DNA rise and twist; G–C BREATHING / A–T HELD |
| 18 | EHT semiempirical SCF | H2 band dissociation, H2O Mulliken + gap, BO switch dissolves the wall |

---

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
