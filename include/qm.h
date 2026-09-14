#ifndef QM_H
#define QM_H

#include "types.h"

/*
 * qm.h — bottom-up quantum mechanics library (Class 2, v1).
 *
 * Purpose: derive chemistry from quantum numbers instead of transcribing
 * it. Every routine below takes (n, l, ml, Zeff, IE, EA, r_mp) or an
 * Atom/Simulation and returns a derived quantity. No per-element
 * hand tables beyond the periodic-table primaries (IE, EA, Zeff).
 *
 * v1 scope (deliberate):
 *  - Real spherical harmonics Y_lm for s, p, d (real form). f returns 0
 *    with documented cutoff (no f chemistry in H..Kr demos).
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
/* Hybridization from the atom's electron config census (s/p counts). */
QmHybrid qm_hybridization(const Atom *atom);

/* Mulliken electronegativity chi = (IE+EA)/2 [eV] and hardness
 * J = (IE-EA)/2 [eV]. EA<=0 entries use EA=0 (table convention). */
void qm_chi_J(const Element *el, double *chi_out, double *J_out);

/* Polarizability volume estimate alpha = r_mp^3 (A^3), r_mp from the
 * valence orbital. Order-of-magnitude sphere estimate, documented. */
double qm_alpha(const Atom *atom);

/* ── Two-center estimates ──────────────────────────────────────────── */
/* Overlap estimate between valence Slater orbitals along the bond axis:
 * S = exp(-zeta*R) * ang, zeta = 0.5*(Zeff_a/(n_a*a0) + Zeff_b/(n_b*a0))
 * in 1/A, R in A, ang = Y_a(dir)*Y_b(-dir) normalized product
 * (|ang|<=1). Documented estimate, not an exact two-center integral. */
double qm_overlap(const Atom *a, const Atom *b, double R_ang, Vec3 dir_a_to_b);

/* Bond order from overlap: BO = clamp(3*S/S_ref, 0, 3) with S_ref the
 * overlap of the same pair at covalent-contact distance (caller passes
 * S_ref; library also provides qm_overlap_ref via covalent radii sum). */
double qm_bond_order(double S, double S_ref);

/* Pauli repulsion estimate E = A*S^2 with A = (J_a+J_b)/2 (eV). */
double qm_pauli(double S, double J_a, double J_b);

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
