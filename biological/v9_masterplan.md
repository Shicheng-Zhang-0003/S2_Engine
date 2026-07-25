# v9 Plan — KcsA Selectivity Filter

No code below. This is the build order and the decisions to make before
any of it gets written.

## Where v8 left off
Three tracks validated separately — nucleic acids, protein backbone,
electrophysiology — explicitly not connected. `neuron.c`'s own comments
name the real gap: ion channel gating currently comes from the 1952
empirical rate equations, not from protein structure. Full gating is
out of reach right now — `sim.h` has no membrane, no voltage, nothing
that couples geometry to an applied field, and a real channel is
thousands of atoms undergoing slow conformational change. The
selectivity filter is the piece of that problem that's static geometry
only, and needs none of that missing infrastructure.

## Goal
Build the minimal KcsA selectivity filter (not the whole channel), and
test whether the existing Coulomb+LJ engine — unmodified — gets the K+
coordination geometry right and prefers K+ over Na+ in the correct
direction.

## Step 0 — prerequisite, do this first and separately
`periodic_table.c` has a real discrepancy: potassium's data row (Z=19)
is fully populated with physically plausible numbers (electronegativity
0.82, ionization 4.341 eV, LJ ε=0.035 kcal/mol, σ=3.812 Å), but the
comment above it still says "mass only stub." Resolve which one is
true — either the comment is stale and the row is trustworthy, or the
row is placeholder-but-plausible and needs real sourcing — before
anything downstream depends on it. Sodium (Z=11) has no such flag and
appears to be in order, but hasn't been given the same scrutiny in this
conversation that, say, iron's Slater screening got; worth a glance
while in that file. This is small and self-contained — doesn't block
anything else below, but should land before Step 2.

## Step 1 — get the real coordinates
Source: PDB 1BL8 (Doyle et al. 1998, original KcsA structure) or 1K4C
(Zhou et al. 2001, higher resolution — this is the structure the 2.85 Å
figure comes from). Need the backbone atoms (N, CA, C, O at minimum)
for residues Thr75–Val76–Gly77–Tyr78–Gly79 ("TVGYG"), across all four
identical subunits — roughly 20–30 atoms total, not the whole channel.

Decision to make here: hand-extract those atoms once (fast, matches how
the nucleobase geometry was sourced — a static array with an honest
provenance comment), or write a small reusable PDB-backbone-atom
reader. Worth the reader only if a full DNA duplex or another real
structure is coming soon after; not worth it for a one-off.

## Step 2 — build and sanity-check the bare filter
Same pattern as Demo 6 before Demo 7: place the real backbone geometry,
check bond lengths and planarity are chemically sane, before doing
anything dynamic with it. Open decision: the filter alone has none of
the surrounding protein that holds it rigid in the real channel, so
letting it fully relax risks collapse into something meaningless.
Recommend restraining the backbone near its crystallographic
coordinates — the same "given correct local geometry" move Demo 11
already used for the alpha helix — rather than a free minimization.

## Step 3 — K+, then Na+
Place one bare K+ ion (needs Step 0 resolved first) at a single central
site — not the full multi-ion single-file arrangement real channels
use; that's a separate, harder problem. Run the existing minimizer.
Report the mean carbonyl-O···K+ distance against the literature target:
2.85 Å, measured range 2.70–3.08 Å.

Then swap the same site to Na+, same protocol, same restraints. Report
the same distance, plus the E_LJ / E_Coulomb / total interaction-energy
breakdown for both ions at the same site, same format as Demo 7's G-C
vs. A-U comparison.

## Step 4 — write it up honestly, in the codebase's own voice
Follow Demo 7's exact convention: report the number, compare directly
to the literature value, then an explicit "Honest caveat" paragraph
rather than folding uncertainty into the number itself.

## Success criteria, set now rather than after the fact
- Geometry near 2.85 Å is a fair test, not a stacked deck — QM/MM
  studies of this exact filter get a K+ coordination number of about
  6.6, classical MD gets about 6.7. Classical isn't obviously
  disadvantaged here.
- K+ binding more favorably than Na+ at the same site, even if the
  margin doesn't match published free-energy numbers, is a legitimate,
  informative v9 result — a Demo-7-shaped outcome, not a failure.
  Getting the *magnitude* of that preference right is the documented
  hard part for non-polarizable force fields; getting the *direction*
  right is the actual test of whether the method extends here at all.

## Explicitly out of scope for v9
Full tetramer beyond the filter loop, membrane/lipid, voltage, gating
conformational change, the multi-ion single-file "knock-on" permeation
mechanism. All real, all later, none needed to answer what v9 is
actually asking.
