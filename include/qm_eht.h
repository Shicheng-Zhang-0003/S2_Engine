/* qm_eht.h — Extended Hückel + self-consistent charge (SCC-EHT-lite).
 *
 * Scope (read before citing): single-zeta STO minimal valence basis
 * (H 1s; C/N/O 2s+2p), Hoffmann 1963 VSIPs, Wolfsberg-Helmholtz K=1.75,
 * no d orbitals, H/C/N/O only (Z 1,6,7,8), Löwdin-orthogonalized SCF
 * with Ohno-interpolated Hubbard gammas (U from in-tree Mulliken J),
 * at most EHT_MAX_AOS basis functions and EHT_MAX_ATOMS atoms — larger
 * systems refuse loudly (rc), never silently.
 *
 * Honesty block (abstraction inventory for this module):
 *  E1 single-zeta STOs (no split-valence, no polarization functions).
 *  E2 minimal valence basis (no core relaxation, no virtuals beyond it).
 *  E3 Wolfsberg-Helmholtz K=1.75 (empirical, not derived).
 *  E4 Ohno gamma with in-tree U (SCC coupling is modeled, not ab initio).
 *  E5 band energy + bare nuclear repulsion: EHT is NON-VARIATIONAL.
 *    Curves are qualitative (dissociation limits correct, minima/depths
 *    approximate). Never quote as predictions; Demo 18 prints both the
 *    EHT curve and the harmonic wall for the same geometry.
 *  E6 BO-scaled bonds treat the Mayer index as a per-step slow variable
 *    (same class as QEq charges): forces derive from scaled k exactly
 *    within the step; the B(R) dependence itself is not differentiated.
 *  E7 minimal-basis EHT gaps run LARGE (H2 ~35 eV at equilibrium: big-S
 *    regime times K=1.75; textbook EHT trait, not a solver error — the 2x2
 *    closed forms agree to 1e-6). Gate bands are generous by design.
 *  E8 restricted closed-shell determinant: homonuclear Mayer BO persists
 *    at dissociation (B~1 at 4 A) — the RHF static-correlation failure,
 *    inherited, not fixed. The FORCE switch below is therefore an
 *    overlap-gated Mayer interpolant (engineered dissociative switch,
 *    labeled as such), while printed BOs stay pure Mayer with this note.
 */
#ifndef QM_EHT_H
#define QM_EHT_H

#include "types.h"

#define EHT_MAX_AOS 64
#define EHT_MAX_ATOMS 128
/* Jacobi eigensolver caps */
#define EHT_JACOBI_SWEEPS 100
#define EHT_SCC_ITERS 50
#define EHT_SCC_TOL 1e-6
/* Fermi smearing for SCC level crossings (E10), eV. Gapped ground states
 * (gaps >> kT) are unaffected; dissociation limits converge fractionally. */
#define EHT_FERMI_KT 0.1
/* SCC step: constant mixing SMIX composed with a DSTEP trust cap
 * (E10; see solver comment). Fixed points untouched. */
#define EHT_SCC_SMIX 0.25
#define EHT_SCC_DSTEP 0.25

/* ---- exact STO overlaps (qm_eht_overlap.c) ---- */

/* STO exponent from Slater Zeff: zeta = Zeff/(n* a0), A^{-1}. -1 bad. */
double eht_zeta (int Z, int n, int l, const ElectronConfig *cfg);
/* s-s (n = 1,2), s-pσ, p-pσ, p-pπ (2p only). Dimensionless. R in A. */
double eht_overlap_ss (int na, double za, int nb, double zb, double R);
double eht_overlap_sps (int ns, double zs, int np, double zp, double R,
                        int p_on_b);
double eht_overlap_pps (double za, double zb, double R);
double eht_overlap_ppi (double za, double zb, double R);

/* ---- SCF (qm_eht_scf.c) ---- */

/* Hoffmann 1963 valence-state ionization potentials (eV, positive).
 * Hii = -VSIP. 0 = unsupported (Z,n,l). */
double eht_vsip (int Z, int n, int l);
/* Valence electron count (H 1, C 4, N 5, O 6). -1 unsupported Z. */
int eht_valence_electrons (int Z);
/* Valence AO count per atom (H 1, C/N/O 4). 0 unsupported. */
int eht_valence_aos (int Z);

typedef struct {
    int n_ao;                           /* basis size */
    int n_atoms;                        /* atoms covered */
    int n_elec;                         /* valence electrons */
    int ao_atom[EHT_MAX_AOS];           /* AO -> atom index */
    int ao_n[EHT_MAX_AOS];
    int ao_l[EHT_MAX_AOS];
    int ao_m[EHT_MAX_AOS];              /* 0 s; xyz flag 1,2,3 for px,py,pz */
    double eps[EHT_MAX_AOS];            /* MO energies, ascending, eV */
    double P[EHT_MAX_AOS][EHT_MAX_AOS]; /* density matrix */
    double S[EHT_MAX_AOS][EHT_MAX_AOS]; /* overlap at solution (for Mayer) */
    double charges[EHT_MAX_ATOMS];      /* Mulliken, e */
    double gap;                         /* LUMO-HOMO, eV */
    double e_band;                      /* Tr(PH), eV */
    double e_nuc;                       /* screened nuclear repulsion, eV */
    double e_total;                     /* e_band + e_nuc, eV (E5: qualitative) */
    int scc_iters;                      /* SCC iterations used */
    double scc_resid;                   /* final max |dq| */
} qm_eht_t;

/* Full SCC solve. 0 ok; -1 unsupported/over cap; -2 singular S;
 * -3 SCC non-convergence. Charges/BO/gap/energies in out. */
int qm_eht_solve (const Simulation *sim, qm_eht_t *out);
/* Mayer bond order from cached solve. -1.0 on bad indices. */
double qm_eht_bond_order (const qm_eht_t *e, int ia, int ib);
/* Reactive hook (E6): forces.c solves once per step into its own static
 * result when sim->use_eht_bo is set, then scales each harmonic bond by
 * qm_eht_bond_factor (BO treated as a per-step slow variable, same class
 * as QEq charges). Over-cap/unsupported warns once and leaves bonds
 * unscaled (C9 pattern). */
/* Overlap-gated Mayer switch (E8): clamp(B * S_cur/S_eq, 0, 1).
 * Needs sim positions; S channels evaluated live, S_eq at the covalent-
 * radii sum. unity at equilibrium, ->0 dissolved. */
double qm_eht_bond_factor (const qm_eht_t *e, const Simulation *sim,
                           int ia, int ib);

#endif /* QM_EHT_H */
