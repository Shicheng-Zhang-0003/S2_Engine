#ifndef QM_H
#define QM_H

#include "types.h"

/*
 * qm.h — bottom-up quantum mechanics library (Class 2, v2).
 *
 * Purpose: derive chemistry from quantum numbers instead of transcribing
 * it. Every routine below takes (n, l, ml, Zeff, IE, EA, r_mp) or an
 * Atom/Simulation and returns a derived quantity. No per-element
 * hand tables beyond the periodic-table primaries (IE, EA, Zeff) plus
 * the documented Applequist polarizability set for v2 induction.
 *
 * v1 scope (delivered):
 *  - Real spherical harmonics Y_lm for s, p, d AND f (real cubic set,
 *    unit norms verified by direct integration; g+ still returns 0).
 *  - Full psi = R_nl(r) * Y_lm(theta,phi) via quantum_radial_wavefunction.
 *  - Hybridization analyser: lobe count/directions from s/p census.
 *  - Mulliken chi/J from table IE/EA; QEq charge equilibration solver.
 *  - Overlap estimate S (exponential decay x angular factor), bond order,
 *    Pauli energy E = A*S^2, polarizability alpha = r_mp^3 (order of
 *    magnitude; fixed 0.84 legacy value kept in KcsA energetics).
 *  - One SCF helper: qm_refresh_charges() runs QEq and writes
 *    partial_charge. Integer-config Slater screening and fractional QEq
 *    charges coexist as separate layers in v1 (charge-dependent screening
 *    is the documented next step).
 *
 * v2 scope (this release — full bottom-up emergence wiring):
 *  - Isotropic induced-dipole polarization with ANALYTIC forces +
 *    Thole damping (a=2.0 A, f=1-exp(-(r/a)^3)):
 *    U_pol = -0.5 * C * Σ_i alpha_i |E_i|^2, E_i[V/A] damped per pair.
 *    Undamped first-order induction has a polarization catastrophe
 *    (E~1/r^2, U~-1/r^4 → minimizer collapse, observed -28000 eV);
 *    Thole keeps first-shell (2.35-2.75 A, f=0.80-0.93) near-exact.
 *    Per-atom alpha from Applequist (C/N/O) with r_mp^3 fallback,
 *    documented per call. Dipole-dipole coupling omitted (first-order
 *    induction only — stated limitation, avoids iterative SCF cost and
 *    residual catastrophe at short range).
 *  - Overlap-Pauli energy AND finite-difference forces (same
 *    correctness-by-construction pattern as dihedral FD): E=A*S^2 summed
 *    over eligible nonbonded pairs, forces via central differences.
 *  - Charge-dependent screening note: Slater Zeff still integer-config;
 *    QEq fractional charges live in the Coulomb/polar layers. Full
 *    q-dependent Zeff remains future work (stated, not silently mixed).
 *
 * v4 scope (quantum-foundation release):
 *  - Clementi-Raimondi SCF exponents for spatial ranges (r_mp, overlap
 *    zeta, dispersion b) with Slater fallback; Slater kept for energies.
 *  - Exact H-like expectations <r>, <r^2>, <1/r>, <T> (analytic).
 *  - Real f harmonics (unit norms verified by integration).
 *  - Lobe-max sigma projection (rotational upper bound) + sum-rule pi
 *    estimate; charge-responsive screening gamma into overlap/alpha.
 *  - Self-consistent dipoles (Thole-damped T-coupling, Hellmann-Feynman
 *    forces) via use_pol_scf; Demo 1 accuracy ledger (E vs NIST IE,
 *    CR r_mp vs covalent radius) and no-fit dimer QM curves.
 *
 * Units: lengths A, energies eV, charges e, alpha A^3, S dimensionless.
 */

typedef struct {
    char   label[16];   /* "1s", "sp", "sp2", "sp3", "p", "d", "none" */
    int    n_lobes;     /* directional lobe count (0 for bare s) */
    Vec3   lobes[4];    /* unit lobe directions (up to 4) */
    int    n_lone_pairs;
} QmHybrid;

/* ── Angular functions ─────────────────────────────────────────────── */
/* Real spherical harmonics, Condon-Shortley phase, normalized over the
 * sphere. Convention: ml=-1 -> py, ml=0 -> pz (m=0), ml=+1 -> px for p;
 * d follows the standard real set (z2, xz, yz, x2-y2, xy) keyed on ml.
 * theta: polar angle from +z [0,pi]; phi: azimuth in xy from +x. */
double qm_Y_real(int l, int ml, double theta, double phi);

/* Full orbital value psi = R_nl(r) * Y_lm, in A^(-3/2). r in A. */
double qm_psi(int n, int l, int ml, double Zeff, Vec3 pos);

/* ── Per-atom derived quantities ───────────────────────────────────── */
/* Hybridisation from the atom's electron config census (s/p counts) and
 * its sigma-partner count. With no topology available the atom is treated
 * as unconjugated - the conservative sp3 reading. */
QmHybrid qm_hybridization(const Atom *atom);

/*
 * Topology-aware hybridisation, the correct entry point whenever a
 * Simulation is available.
 *
 * AUDIT FIX F7: amide nitrogen and ammonia nitrogen have IDENTICAL local
 * bond orders (all single) and identical valence, so no bare-Atom
 * classifier can tell them apart - yet the amide is sp2 and the amine is
 * sp3, because the amide's lone pair is delocalised into the adjacent
 * carbonyl pi system and therefore occupies a p orbital rather than a
 * hybrid. This entry point inspects the bonded partners' bond orders to
 * establish conjugation, which is what makes the steric number correct.
 *
 * Verified: water O sp3, ammonia N sp3, carbonyl O sp2, amide N sp2,
 * carbonyl C sp2, methylene C sp3.
 */
QmHybrid qm_hybridization_ctx(const Simulation *sim, int atom_idx);

/* Mulliken electronegativity chi = (IE+EA)/2 [eV] and hardness
 * J = (IE-EA)/2 [eV]. EA<=0 entries use EA=0 (table convention). */
void qm_chi_J(const Element *el, double *chi_out, double *J_out);

/* Polarizability volume estimate alpha = r_mp^3 (A^3), r_mp from the
 * valence orbital. Order-of-magnitude sphere estimate, documented.
 * NOTE scale gap: qm_alpha(O) ~0.10-0.13 A^3 while the KcsA induction
 * leg uses Applequist carbonyl-O 0.84 A^3. Both are order-of-magnitude;
 * the fixed 0.84 is kept for energetics continuity, qm_alpha for
 * diagnostics. Do not mix without rescaling. */
double qm_alpha(const Atom *atom);

/* Reference overlap S_ref for qm_bond_order: same-pair overlap at
 * covalent-contact distance (covalent-radii sum). Use this instead of
 * a hand-typed distance so BO is not systematically underestimated. */
double qm_overlap_ref(const Atom *a, const Atom *b, Vec3 dir_a_to_b);

/* ── Two-center estimates ──────────────────────────────────────────── */
/* Overlap estimate between valence Slater orbitals along the bond axis:
 * S = exp(-zeta*R) * ang, zeta = 0.5*(Zeff_a/(n_a*a0) + Zeff_b/(n_b*a0))
 * in 1/A, R in A, ang = Y_a(dir)*Y_b(-dir) normalized product
 * (|ang|<=1). Documented estimate, not an exact two-center integral. */
double qm_overlap(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b);

/* Charge-responsive overlap: Zeff shifted by gamma*q per atom
 * (gamma from integer-config Slater finite differences, q in e).
 * qm_overlap() == qm_overlap_q(...,0,0) exactly. */
double qm_overlap_q(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b,
                    double qa, double qb);

/* Screening slope dZeff/dq for an atom's (n,l) shell, from integer
 * Slater configs (no tables). Positive ~0.35 intra-shell: anions
 * (q<0) get lower Zeff (diffuse), cations higher (contract). */
double qm_gamma_atom(const Atom *atom, int n, int l);

/* Pi-overlap magnitude estimate for p-p pairs: S_pi = exp(-zeta*R) *
 * ang_pi, where ang_pi comes from the l-shell sum rule
 * (Y_sig^2 + Y_pi1^2 + Y_pi2^2 = (2l+1)/4pi per atom): the sigma
 * channel takes its share, the remainder bounds the pi channels.
 * Returns 0 for s-involving pairs (no pi symmetry). Documented
 * estimate for stacking/aromatic analysis, not an integral. */
double qm_overlap_pi(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b);

/* Bond order from overlap: BO = clamp(3*S/S_ref, 0, 3) with S_ref the
 * overlap of the same pair at covalent-contact distance (caller passes
 * S_ref; library also provides qm_overlap_ref via covalent radii sum). */
double qm_bond_order(double S, double S_ref);

/* Pauli repulsion estimate E = A*S^2 with A = (J_a+J_b)/2 (eV). */
double qm_pauli(double S, double J_a, double J_b);

/* ── v2: per-atom polarizability for induction ─────────────────────── */
/* Applequist atom volumes (A^3) for H/C/N/O + r_mp^3 fallback.
 * Returns 0 for unknown/unsupported (no silent default). */
double qm_polarizability(const Atom *atom);

/* Ion-aware (alpha, IE) for dispersion: closed-shell ions use ion data
 * (K+: alpha 0.83 A^3, IE 31.63 eV second IE; Na+: 0.18, 47.29),
 * not neutral-atom values. Neutral atoms use qm_polarizability +
 * table IE. Writes alpha (A^3) and IE (eV). */
void qm_alpha_IE(const Atom *atom, double *alpha_out, double *IE_out);

/* Slater-Kirkwood C6 [eV·A^6] = 1.5*αa*αb*Ia*Ib/(Ia+Ib), α in A^3,
 * I in eV (conversion factors cancel to exactly 1.5 — see qm.c). */
double qm_c6(const Atom *a, const Atom *b);

/* Tang-Toennies damping f6(x)=1-exp(-x)*Σ_{k=0..6} x^k/k!, x=b*r,
 * b = mean Slater zeta (1/A) of the pair. Returns f and via dfdr. */
double qm_tt_f6(double b, double r, double *dfdr_out);

/* ── v2: induction energy + analytic forces ────────────────────────── */
/* U_pol = -0.5*C*Σ_i alpha_i |E_i|^2 with E_i from all other point
 * charges (1-2/1-3 excluded like LJ when called from forces_calculate;
 * direct callers pass an exclusion callback or NULL for all-pairs).
 * Accumulates analytic forces onto sim->atoms[]. Returns U_pol in eV.
 * Dipole-dipole coupling omitted (first-order only, documented). */
double qm_induction_forces(Simulation *sim, double dielectric);

/* Energy-only induction (no force writes) for diagnostics. */
double qm_induction_energy(const Simulation *sim, double dielectric);

/* ── v4: self-consistent dipoles ───────────────────────────────────── */
/*
 * Solve the coupled-dipole equations for the induced dipoles.
 *
 * AUDIT FIX D2: the self-consistency condition is LINEAR in mu,
 *
 *     mu_i = alpha_i ( E0_i + sum_j f_thole(r_ij) [3 d d - I] mu_j
 *                      / r_ij^3 * k / eps_r )
 *
 * so this is a single N x N system (I - A) mu = alpha (*) E0, solved
 * exactly by Gaussian elimination with partial pivoting. It was
 * previously a damped fixed-point iteration, which needed 96-97 of its
 * 100 iterations even for four atoms and silently failed on ordinary
 * geometries - the caller then swapped the energy functional, which
 * made the potential energy discontinuous (audit D3).
 *
 * mu_out receives N entries in e*A. Returns 0 on success, -1 if the
 * system is singular at this geometry, which is the honest signal that
 * the induced-dipole model has no solution (undamped catastrophe
 * cancellation) rather than a value to paper over.
 */
int qm_solve_dipoles(const Simulation *sim, double dielectric, Vec3 *mu_out);

/*
 * SCF induction energy AND its force, in one conservative call.
 *
 * AUDIT FIX D2/D3: the force is the exact central-difference gradient
 * of the returned energy. The Hellmann-Feynman shortcut is unavailable
 * once the dipoles are self-consistent, because E = (I - T alpha)^-1 E0
 * and so dU/dE0 = -mu (I - T alpha)^-1, which the shared pointwise
 * force kernel cannot express. Measured against finite differences, the
 * previous analytic weights were 98% wrong; correcting the weight to the
 * naive 1/2 made it worse (relL2 3.3), which is the signature of a
 * wrong model rather than a wrong factor.
 *
 * The energy functional is the SAME whether or not the dipole system
 * solves, so E stays continuous in geometry - the previous code
 * substituted a different functional on non-convergence, producing a
 * measured 5.92 eV jump over 0.01 A and no energy conservation.
 *
 * rc_out (may be NULL) receives 0 if the dipoles solved, -1 otherwise.
 * Returns U in eV.
 */
double qm_induction_scf_forces(Simulation *sim, double dielectric, int *rc_out);

/* ── v2: Pauli energy + finite-difference forces ───────────────────── */
/* E = Σ_pairs A*S^2 over eligible nonbonded pairs (cutoff + 1-2/1-3
 * excluded, same set as LJ). Forces via central FD (h=1e-5 A), same
 * pattern as dihedral FD. Accumulates onto atoms, returns total eV. */
double qm_pauli_forces(Simulation *sim);

/* Energy-only Pauli for diagnostics. */
double qm_pauli_energy(const Simulation *sim);

/* ── v3: QM-derived dispersion + SCF driver ────────────────────────── */
/* Damped dispersion E = -Σ C6*f6(b,r)/r^6 with analytic forces.
 * Accumulates onto atoms, returns total eV. Same pair eligibility
 * as Pauli/LJ. */
double qm_dispersion_forces(Simulation *sim);

/* Energy-only dispersion. */
double qm_dispersion_energy(const Simulation *sim);

/* SCF driver: iterate pinned-QEq charges (ion fixed at +1) with dipole
 * reaction-field shift dchi_i = -C*Σ_j alpha_j E_j·T_ij (dipole
 * potential feeding back into electronegativity) until charges change
 * < tol (max 5 iters). Then writes partial_charge. Returns iterations
 * used, or -1 on solver failure (charges untouched). */
int qm_scf_charges(Simulation *sim, double total_q, double dielectric,
                   int pinned_idx, double pinned_q);

/* Two-pin variant (e.g. knock-on ion pair): both pinned atoms hold
 * fixed charges, shell equilibrates around them. Same return contract. */
int qm_scf_charges_2pin(Simulation *sim, double total_q, double dielectric,
                        int pinned0, double q0, int pinned1, double q1);
/* ── Charge equilibration (Rappe-Goddard QEq-lite) ─────────────────── */
/* Solve chi_i + J_i*q_i + sum_j C_ij*q_j = mu (equalization) with
 * sum q = total_q, C_ij = COULOMB_MD/(r_ij*dielectric) (i!=j).
 * Writes out_q[0..n-1]. Returns 0 on success, -1 on singular/failure
 * (out_q untouched). n<=128. O(n^3) Gaussian elimination; demo-scale. */
int qm_qeq(const Simulation *sim, double total_q, double dielectric,
           double *out_q);

/* One SCF helper: refresh all partial_charge from QEq with the given
 * total charge and dielectric. Returns 0 ok, -1 if solver failed
 * (charges untouched). */
int qm_refresh_charges(Simulation *sim, double total_q, double dielectric);

/* Pinned variant: atom pinned_idx holds fixed charge pinned_q (e.g. a
 * closed-shell ion whose charge state neutral-atom chi/J cannot
 * describe); all other atoms equilibrate with sum = total_q - pinned_q.
 * Returns 0 ok, -1 on failure (out_q untouched). */
int qm_qeq_pinned(const Simulation *sim, double total_q, double dielectric,
                  int pinned_idx, double pinned_q, double *out_q);

#endif /* QM_H */
