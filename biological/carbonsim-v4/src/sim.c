#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "../include/sim.h"
#include "../include/periodic_table.h"
#include "../include/quantum.h"
#include "../include/forces.h"
#include "../include/constants.h"

/*
 * sim.c
 * Simulation construction, topology management, and molecule builders.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Lifecycle
 * ══════════════════════════════════════════════════════════════════════════ */
Simulation *sim_create(int atom_capacity, int bond_capacity) {
    Simulation *sim = (Simulation *)calloc(1, sizeof(Simulation));
    if (!sim) return NULL;

    sim->atoms = (Atom *)calloc(atom_capacity, sizeof(Atom));
    sim->bonds = (Bond *)calloc(bond_capacity, sizeof(Bond));
    /* Angles: heuristic — up to 3 × bonds */
    int angle_cap = bond_capacity * 3;
    sim->angles   = (Angle *)calloc(angle_cap, sizeof(Angle));

    if (!sim->atoms || !sim->bonds || !sim->angles) {
        free(sim->atoms); free(sim->bonds); free(sim->angles);
        free(sim);
        return NULL;
    }

    sim->capacity_atoms    = atom_capacity;
    sim->capacity_bonds    = bond_capacity;
    sim->capacity_angles   = angle_cap;

    /* Defaults */
    sim->dt      = 1.0;     /* 1 fs timestep                                */
    sim->cutoff  = 12.0;    /* 12 Å non-bonded cutoff                       */
    sim->use_lj       = 1;
    sim->use_coulomb  = 1;
    sim->use_bonds    = 1;
    sim->use_angles   = 1;
    sim->dielectric   = 1.0;  /* no screening by default */
    sim->ff_type = FF_UFF;

    /* Box: large vacuum (no PBC by default) */
    sim->box.dimensions = vec3(1000.0, 1000.0, 1000.0);
    sim->box.periodic[0] = sim->box.periodic[1] = sim->box.periodic[2] = 0;

    /* Berendsen thermostat off by default */
    sim->thermostat.type               = THERMOSTAT_NONE;
    sim->thermostat.target_temperature = 300.0;
    sim->thermostat.tau                = 100.0;  /* fs */

    return sim;
}

void sim_destroy(Simulation *sim) {
    if (!sim) return;
    free(sim->atoms);
    free(sim->bonds);
    free(sim->angles);
    free(sim->dihedrals);
    free(sim);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Atom management
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_add_atom(Simulation *sim, int Z, Vec3 pos, double partial_charge) {
    if (sim->num_atoms >= sim->capacity_atoms) return SIM_ERR_OVERFLOW;

    const Element *el = pt_element(Z);
    if (!el) return SIM_ERR_BADATOM;

    int idx = sim->num_atoms++;
    Atom *a = &sim->atoms[idx];
    memset(a, 0, sizeof(Atom));

    a->id             = idx;
    a->Z              = Z;
    a->element        = el;
    a->position       = pos;
    a->velocity       = vec3_zero();
    a->force          = vec3_zero();
    a->mass           = el->mass;
    a->partial_charge = partial_charge;
    a->formal_charge  = 0;

    /* Default LJ parameters to the generic UFF element values.
     * Molecule constructors (e.g. sim_place_h2o) may override these
     * per atom for force-field-specific contexts via sim_set_atom_lj(). */
    a->lj_epsilon = el->lj_epsilon;
    a->lj_sigma   = el->lj_sigma;

    /* Electron configuration */
    pt_electron_config(Z, &a->electron_config);

    /* Populate orbital table */
    quantum_fill_orbitals(a);

    return idx;
}

int sim_add_atom_sym(Simulation *sim, const char *symbol,
                     Vec3 pos, double q) {
    const Element *el = pt_by_symbol(symbol);
    if (!el) return SIM_ERR_BADATOM;
    return sim_add_atom(sim, el->Z, pos, q);
}

void sim_set_atom_lj(Simulation *sim, int atom_idx, double epsilon, double sigma) {
    if (atom_idx < 0 || atom_idx >= sim->num_atoms) return;
    sim->atoms[atom_idx].lj_epsilon = epsilon;
    sim->atoms[atom_idx].lj_sigma   = sigma;
}

void sim_set_bond_params(Simulation *sim, int bond_idx, double r0, double k) {
    if (bond_idx < 0 || bond_idx >= sim->num_bonds) return;
    sim->bonds[bond_idx].r0 = r0;
    sim->bonds[bond_idx].k  = k;
}

int sim_add_angle_explicit(Simulation *sim, int a, int b, int c,
                            double theta0, double k) {
    if (a < 0 || a >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (b < 0 || b >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (c < 0 || c >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (sim->num_angles >= sim->capacity_angles) return SIM_ERR_OVERFLOW;

    int idx = sim->num_angles++;
    Angle *ang = &sim->angles[idx];
    ang->atom_a = a;
    ang->atom_b = b;
    ang->atom_c = c;
    ang->theta0 = theta0;
    ang->k      = k;
    return idx;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Bond management
 * ══════════════════════════════════════════════════════════════════════════ */
int sim_add_bond(Simulation *sim, int ia, int ib, int order) {
    if (ia < 0 || ia >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (ib < 0 || ib >= sim->num_atoms) return SIM_ERR_BADATOM;
    if (sim->num_bonds >= sim->capacity_bonds) return SIM_ERR_OVERFLOW;

    int idx = sim->num_bonds++;
    Bond *b = &sim->bonds[idx];

    b->atom_a = ia;
    b->atom_b = ib;
    b->order  = order;

    /* Look up equilibrium parameters */
    BondParam bp;
    forces_bond_params(sim->atoms[ia].Z, sim->atoms[ib].Z, order, &bp);
    b->r0 = bp.r0;
    b->k  = bp.k;

    /* Current length */
    b->current_length = vec3_dist(sim->atoms[ia].position,
                                   sim->atoms[ib].position);
    b->energy = 0.0;

    /* Update atom bond lists */
    Atom *a = &sim->atoms[ia];
    Atom *c = &sim->atoms[ib];
    if (a->num_bonds < MAX_BONDS_PER_ATOM) {
        a->bond_partners[a->num_bonds]  = ib;
        a->bond_orders  [a->num_bonds]  = order;
        a->num_bonds++;
    }
    if (c->num_bonds < MAX_BONDS_PER_ATOM) {
        c->bond_partners[c->num_bonds]  = ia;
        c->bond_orders  [c->num_bonds]  = order;
        c->num_bonds++;
    }

    return idx;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Auto-detect bonds from distances
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_detect_bonds(Simulation *sim) {
    int N = sim->num_atoms;

    for (int i = 0; i < N - 1; i++) {
        for (int j = i + 1; j < N; j++) {
            const Element *ei = sim->atoms[i].element;
            const Element *ej = sim->atoms[j].element;

            double r_cov_sum = ei->covalent_radius + ej->covalent_radius;
            double r_max     = BOND_TOLERANCE * r_cov_sum;
            double r         = vec3_dist(sim->atoms[i].position,
                                          sim->atoms[j].position);

            if (r < r_max) {
                /* Estimate bond order from distance ratio */
                int order = 1;
                double ratio = r / r_cov_sum;
                if (ratio < 0.78) order = 3;
                else if (ratio < 0.87) order = 2;

                sim_add_bond(sim, i, j, order);
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Build angle list from bond topology
 * For each atom b, for every pair of bonds (a-b) and (b-c): add angle a-b-c
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_build_angles(Simulation *sim) {
    sim->num_angles = 0;

    for (int b_idx = 0; b_idx < sim->num_atoms; b_idx++) {
        Atom *b = &sim->atoms[b_idx];
        /* Enumerate all pairs of neighbours of b */
        for (int p = 0; p < b->num_bonds - 1; p++) {
            for (int q = p + 1; q < b->num_bonds; q++) {
                if (sim->num_angles >= sim->capacity_angles) return;

                int ia = b->bond_partners[p];
                int ic = b->bond_partners[q];

                Angle *ang = &sim->angles[sim->num_angles++];
                ang->atom_a = ia;
                ang->atom_b = b_idx;
                ang->atom_c = ic;

                /* Look up angle parameters */
                AngleParam ap;
                forces_angle_params(sim->atoms[ia].Z, b->Z,
                                    sim->atoms[ic].Z, &ap);
                ang->theta0 = ap.theta0;
                ang->k      = ap.k;
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Build angle list with theta0 derived from current Cartesian geometry
 * (rather than the generic element-keyed table) - for ring systems where
 * the seed positions already encode the correct, context-specific angle.
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_build_angles_geometric(Simulation *sim, double k_default) {
    sim->num_angles = 0;

    for (int b_idx = 0; b_idx < sim->num_atoms; b_idx++) {
        Atom *b = &sim->atoms[b_idx];
        for (int p = 0; p < b->num_bonds - 1; p++) {
            for (int q = p + 1; q < b->num_bonds; q++) {
                if (sim->num_angles >= sim->capacity_angles) return;

                int ia = b->bond_partners[p];
                int ic = b->bond_partners[q];

                double theta = vec3_angle(
                    vec3_sub(sim->atoms[ia].position, b->position),
                    vec3_sub(sim->atoms[ic].position, b->position));

                Angle *ang = &sim->angles[sim->num_angles++];
                ang->atom_a = ia;
                ang->atom_b = b_idx;
                ang->atom_c = ic;
                ang->theta0 = theta;
                ang->k      = k_default;
            }
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * PBC helpers
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_set_box(Simulation *sim, double lx, double ly, double lz) {
    sim->box.dimensions = vec3(lx, ly, lz);
    sim->box.periodic[0] = sim->box.periodic[1] = sim->box.periodic[2] = 1;
}

void sim_wrap_positions(Simulation *sim) {
    Vec3 L = sim->box.dimensions;
    for (int i = 0; i < sim->num_atoms; i++) {
        Vec3 *p = &sim->atoms[i].position;
        if (sim->box.periodic[0]) {
            p->x -= L.x * floor(p->x / L.x);
        }
        if (sim->box.periodic[1]) {
            p->y -= L.y * floor(p->y / L.y);
        }
        if (sim->box.periodic[2]) {
            p->z -= L.z * floor(p->z / L.z);
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Molecule constructors
 * ══════════════════════════════════════════════════════════════════════════ */

/* ── H₂ ─────────────────────────────────────────────────────────────────── */
int sim_place_h2(Simulation *sim, Vec3 origin) {
    /* Bond length 0.7414 Å, atoms symmetric about origin */
    int first = sim->num_atoms;
    sim_add_atom(sim, 1, vec3(origin.x - 0.3707, origin.y, origin.z), 0.0);
    sim_add_atom(sim, 1, vec3(origin.x + 0.3707, origin.y, origin.z), 0.0);
    sim_add_bond(sim, first, first+1, 1);
    sim_build_angles(sim);
    return first;
}

/* ── H₂O (TIP3P geometry) ───────────────────────────────────────────────── */
/*
 * TIP3P parameters:
 *   r(OH) = 0.9572 Å
 *   ∠HOH  = 104.52°
 *   Oxygen at origin.
 *   H positions:
 *     H1: ( r sinα, -r cosα, 0 )    α = half-angle = 52.26°
 *     H2: (-r sinα, -r cosα, 0 )
 *   Partial charges: q(O) = -0.834e, q(H) = +0.417e
 */
int sim_place_h2o(Simulation *sim, Vec3 origin) {
    double rOH   = 0.9572;
    double alpha = 52.26 * 3.14159265358979323846 / 180.0;
    double hx    = rOH * sin(alpha);
    double hy    = rOH * cos(alpha);

    int first = sim->num_atoms;
    /* O */
    int iO = sim_add_atom(sim, 8,
        vec3(origin.x,      origin.y,      origin.z),      -0.834);
    /* H1 */
    int iH1 = sim_add_atom(sim, 1,
        vec3(origin.x + hx, origin.y - hy, origin.z),      +0.417);
    /* H2 */
    int iH2 = sim_add_atom(sim, 1,
        vec3(origin.x - hx, origin.y - hy, origin.z),      +0.417);

    /*
     * TIP3P (Jorgensen, Chandrasekhar, Madura, Impey & Klein, J. Chem.
     * Phys. 79, 926 (1983)) — verified parameters:
     *   Oxygen: sigma = 3.15061 Å, epsilon = 0.1521 kcal/mol
     *   Hydrogen: NO Lennard-Jones site at all (sigma = epsilon = 0).
     *     Only the oxygen carries a LJ centre in the original model;
     *     the CHARMM-modified TIP3P variant adds small LJ terms to H,
     *     but that is a distinct, separately-named force field.
     *
     * Generic UFF defaults (sigma_O=3.500 Å) were fit for general
     * organic chemistry, not for this specific, jointly-parameterized
     * water model — using them here would put the O-O equilibrium
     * distance in the wrong place relative to where the TIP3P charges
     * expect it, breaking the delicate cancellation that holds liquid
     * water together. Override explicitly.
     */
    sim_set_atom_lj(sim, iO,  0.1521 * KCAL_MOL_TO_EV, 3.15061);
    sim_set_atom_lj(sim, iH1, 0.0,                      0.0);
    sim_set_atom_lj(sim, iH2, 0.0,                      0.0);

    sim_add_bond(sim, first,   first+1, 1);  /* O-H1 */
    sim_add_bond(sim, first,   first+2, 1);  /* O-H2 */
    sim_build_angles(sim);
    return first;
}

/* ── NH₃ ────────────────────────────────────────────────────────────────── */
/*
 * Pyramidal geometry:
 *   r(NH) = 1.012 Å, ∠HNH = 106.67°
 *   N at origin, H atoms placed in trigonal pyramidal arrangement.
 *   Partial charges: q(N) = -1.02e, q(H) = +0.34e
 */
int sim_place_nh3(Simulation *sim, Vec3 origin) {
    double rNH   = 1.012;
    double hnh   = 106.67 * 3.14159265358979323846 / 180.0;
    /* Tetrahedral-like placement:
     * H at 120° azimuthal spacing, cone half-angle θ_cone */
    double cos_hnh = cos(hnh);
    /* cos(∠HNH) = (v1·v2) where vi = Hi-N = rNH*(sin θ cos φi, sin θ sin φi, -cos θ) */
    /* v1·v2 = rNH² (sin²θ cos(2π/3) - cos²θ) ... setting equal to cos(hnh): */
    /* sin²θ cos(120°) - cos²θ = cos_hnh
     * -0.5 sin²θ - cos²θ = cos_hnh
     * -0.5(1-cos²θ) - cos²θ = cos_hnh
     * -0.5 + 0.5cos²θ - cos²θ = cos_hnh
     * -0.5 - 0.5cos²θ = cos_hnh
     * cos²θ = -2(cos_hnh + 0.5) = -(2cos_hnh + 1)  */
    double cos2_theta = -(2.0*cos_hnh + 1.0);
    if (cos2_theta < 0.0) cos2_theta = 0.0;
    double cos_theta  = sqrt(cos2_theta);
    double sin_theta  = sqrt(1.0 - cos2_theta);

    int first = sim->num_atoms;
    sim_add_atom(sim, 7,
        vec3(origin.x, origin.y, origin.z), -1.02);

    double PI23 = 2.0 * 3.14159265358979323846 / 3.0;
    for (int k = 0; k < 3; k++) {
        double phi = k * PI23;
        sim_add_atom(sim, 1,
            vec3(origin.x + rNH * sin_theta * cos(phi),
                 origin.y + rNH * sin_theta * sin(phi),
                 origin.z - rNH * cos_theta),
            +0.34);
        sim_add_bond(sim, first, first+1+k, 1);
    }
    sim_build_angles(sim);
    return first;
}

/* ── CH₄ ────────────────────────────────────────────────────────────────── */
/*
 * Tetrahedral: r(CH) = 1.090 Å, ∠HCH = 109.47°
 * H positions: vertices of a tetrahedron inscribed in a cube.
 * q(C) = -0.24e, q(H) = +0.06e (OPLS-AA)
 */
int sim_place_ch4(Simulation *sim, Vec3 origin) {
    double rCH = 1.090;
    /* Tetrahedral vertices: (±1,±1,±1) normalised, alternating sign for
     * a true tetrahedron (not all-positive corners of cube). */
    double s = rCH / sqrt(3.0);
    double H[4][3] = {
        { s,  s,  s},
        { s, -s, -s},
        {-s,  s, -s},
        {-s, -s,  s}
    };

    int first = sim->num_atoms;
    sim_add_atom(sim, 6, origin, -0.24);   /* C */
    for (int k = 0; k < 4; k++) {
        sim_add_atom(sim, 1,
            vec3(origin.x + H[k][0],
                 origin.y + H[k][1],
                 origin.z + H[k][2]),
            +0.06);
        sim_add_bond(sim, first, first+1+k, 1);
    }
    sim_build_angles(sim);
    return first;
}

/* ── CO₂ ────────────────────────────────────────────────────────────────── */
/*
 * Linear: r(C=O) = 1.163 Å
 * q(C) = +0.70e, q(O) = -0.35e  (OPLS)
 */
int sim_place_co2(Simulation *sim, Vec3 origin) {
    double rCO = 1.163;
    int first = sim->num_atoms;
    sim_add_atom(sim, 8, vec3(origin.x - rCO, origin.y, origin.z), -0.35);
    sim_add_atom(sim, 6, vec3(origin.x,        origin.y, origin.z), +0.70);
    sim_add_atom(sim, 8, vec3(origin.x + rCO, origin.y, origin.z), -0.35);

    sim_add_bond(sim, first,   first+1, 2);  /* O=C */
    sim_add_bond(sim, first+1, first+2, 2);  /* C=O */
    sim_build_angles(sim);
    return first;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Diagnostics
 * ══════════════════════════════════════════════════════════════════════════ */
void sim_print_atoms(const Simulation *sim) {
    printf("  %-4s %-4s  %-22s  %-22s  %-8s  %-8s\n",
           "idx","sym","Position (Å)","Velocity (Å/fs)","q(e)","mass(AMU)");
    printf("  %s\n",
        "──────────────────────────────────────────────────────────────"
        "───────────────────────────────────");
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        printf("  %-4d %-4s  (%7.4f %7.4f %7.4f)  "
               "(%7.4f %7.4f %7.4f)  %+7.4f  %7.3f\n",
               i, a->element->symbol,
               a->position.x, a->position.y, a->position.z,
               a->velocity.x, a->velocity.y, a->velocity.z,
               a->partial_charge, a->mass);
    }
}

void sim_print_bonds(const Simulation *sim) {
    printf("  %-6s %-4s %-4s %-5s %-8s %-8s %-10s\n",
           "bond","a","b","order","r0(Å)","r(Å)","E(eV)");
    printf("  %s\n","─────────────────────────────────────────────────");
    for (int i = 0; i < sim->num_bonds; i++) {
        const Bond *b = &sim->bonds[i];
        double r = vec3_dist(sim->atoms[b->atom_a].position,
                              sim->atoms[b->atom_b].position);
        double e = 0.5 * b->k * (r - b->r0) * (r - b->r0);
        printf("  %-6d %-4d %-4d %-5d %-8.4f %-8.4f %-10.6f\n",
               i, b->atom_a, b->atom_b, b->order, b->r0, r, e);
    }
}

void sim_print_angles(const Simulation *sim) {
    printf("  %-6s %-4s %-4s %-4s %-10s %-10s\n",
           "angle","a","b","c","θ0(deg)","k(eV/rad²)");
    printf("  %s\n","────────────────────────────────────────────────");
    for (int i = 0; i < sim->num_angles; i++) {
        const Angle *a = &sim->angles[i];
        printf("  %-6d %-4d %-4d %-4d %-10.2f %-10.4f\n",
               i, a->atom_a, a->atom_b, a->atom_c,
               a->theta0 * 180.0 / 3.14159265358979323846,
               a->k);
    }
}

void sim_print_summary(const Simulation *sim) {
    printf("\n════════════════════════════════════════════════════════\n");
    printf("  Simulation summary\n");
    printf("  Atoms: %d  Bonds: %d  Angles: %d\n",
           sim->num_atoms, sim->num_bonds, sim->num_angles);
    printf("  Step: %llu   Time: %.3f fs   dt: %.3f fs\n",
           (unsigned long long)sim->step, sim->time, sim->dt);
    printf("  KE: %.6f eV   PE: %.6f eV   E_total: %.6f eV\n",
           sim->kinetic_energy, sim->potential_energy, sim->total_energy);
    printf("  Temperature: %.2f K\n", sim->temperature);
    printf("════════════════════════════════════════════════════════\n");
    sim_print_atoms(sim);
    printf("\n");
    sim_print_bonds(sim);
    printf("\n");
    sim_print_angles(sim);
    printf("════════════════════════════════════════════════════════\n\n");
}
