#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>
#include "constants.h"
#include "vec3.h"

/*
 * types.h
 * All data structures for the chemistry simulator, organised bottom-up:
 *   QuantumNumbers → ElectronConfig → Element (periodic table entry)
 *   → Atom → Bond → Simulation
 *
 * Design principle: every struct is the minimal faithful representation of
 * its physical counterpart. No magic numbers, no opaque indices — every
 * field corresponds to a measurable physical quantity with documented units.
 */

/* ══════════════════════════════════════════════════════════════════════════
 * Quantum layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * The four quantum numbers that uniquely identify an electron's state
 * in a hydrogen-like orbital.
 *
 *  n  : principal quantum number     1, 2, 3, …
 *  l  : azimuthal quantum number     0 … n-1    (0=s, 1=p, 2=d, 3=f)
 *  ml : magnetic quantum number      -l … +l
 *  ms : spin quantum number          +0.5 (↑) or -0.5 (↓)
 */
typedef struct {
    int    n;
    int    l;
    int    ml;
    double ms;
} QuantumNumbers;

/*
 * A single atomic orbital, characterised by quantum numbers, its
 * current energy, and how many electrons occupy it (0, 1, or 2).
 * orbital_energy is in eV, computed via Slater screening.
 */
typedef struct {
    QuantumNumbers qn;
    double         orbital_energy;  /* eV  (negative = bound)               */
    int            occupation;      /* 0, 1, or 2                            */
} Orbital;

/*
 * Electron configuration stored as a 2D array.
 *   config[shell][subshell] = electron count
 *   shell    : 0-indexed principal quantum number n-1  (n=1..7)
 *   subshell : 0=s, 1=p, 2=d, 3=f
 *
 * Example — Carbon (Z=6): 1s² 2s² 2p²
 *   config[0][0] = 2   (1s)
 *   config[1][0] = 2   (2s)
 *   config[1][1] = 2   (2p)
 */
typedef struct {
    int config[MAX_SHELLS][4];
    int total_electrons;
    int valence_electrons;
} ElectronConfig;

/* ══════════════════════════════════════════════════════════════════════════
 * Periodic table layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * Immutable element data loaded once from the periodic table database.
 * All lengths in Ångström, energies in eV, masses in AMU.
 * lj_epsilon and lj_sigma are UFF Lennard-Jones parameters.
 */
typedef struct {
    int    Z;                        /* atomic number                        */
    char   symbol[4];                /* e.g. "C", "Fe"                      */
    char   name[32];                 /* e.g. "Carbon"                        */
    double mass;                     /* AMU                                  */
    double electronegativity;        /* Pauling scale (-1 = unknown)         */
    double atomic_radius;            /* Å (van der Waals)                    */
    double covalent_radius;          /* Å                                    */
    double vdw_radius;               /* Å                                    */
    double ionization_energy;        /* eV (first ionisation)                */
    double electron_affinity;        /* eV (0 if negative / noble gas)       */
    int    valence;                  /* common valence (dominant)            */
    ElectronConfig ground_config;
    double lj_epsilon;               /* eV  (UFF well depth)                 */
    double lj_sigma;                 /* Å   (UFF collision diameter)         */
} Element;

/* ══════════════════════════════════════════════════════════════════════════
 * Molecular dynamics layer
 * ══════════════════════════════════════════════════════════════════════════ */

/*
 * An atom in the simulation.
 * Classical MD representation: position / velocity / force + electronic state.
 *
 * Units:
 *   position : Å
 *   velocity : Å/fs
 *   force    : eV/Å
 *   mass     : AMU
 *   charge   : elementary charge e
 */
typedef struct {
    int            id;
    int            Z;                /* atomic number                        */
    const Element *element;          /* pointer into the periodic table      */

    /* Classical MD state */
    Vec3   position;                 /* Å                                    */
    Vec3   velocity;                 /* Å/fs                                 */
    Vec3   force;                    /* eV/Å  (accumulated each step)        */
    double mass;                     /* AMU                                  */

    /* Electronic state */
    ElectronConfig electron_config;
    double         partial_charge;   /* e   (partial charge for Coulomb)     */
    int            formal_charge;    /* integer formal charge                */

    /*
     * Per-atom LJ parameters. Default to element->lj_epsilon/lj_sigma (UFF)
     * at creation time, but OVERRIDABLE per atom-type-in-context.
     *
     * This matters because real force-field parameters are not pure
     * element properties — they depend on chemical context. TIP3P's
     * water oxygen (σ=3.1507 Å, ε=0.1521 kcal/mol) is a DIFFERENT atom
     * type from generic UFF oxygen (σ=3.500 Å), even though both are
     * element 8. Mixing them inconsistently (TIP3P charges + UFF LJ)
     * produces unphysical results — see sim_place_h2o() for the fix.
     * TIP3P additionally gives hydrogen NO Lennard-Jones term at all
     * (σ=ε=0); only oxygen carries the LJ site in the original 1983 model.
     */
    double         lj_epsilon;       /* eV  (atom-type-specific, not just element default) */
    double         lj_sigma;         /* Å                                    */

    /* Orbital data (populated by quantum module for small systems) */
    Orbital        orbitals[MAX_ORBITALS];
    int            num_orbitals;

    /* Bonding topology */
    int            bond_partners[MAX_BONDS_PER_ATOM];
    int            bond_orders[MAX_BONDS_PER_ATOM];   /* 1,2,3              */
    int            num_bonds;
} Atom;

/*
 * A covalent bond between two atoms.
 * The harmonic potential is: V = 0.5 * k * (r - r0)²
 *   k  : force constant   eV/Å²
 *   r0 : equilibrium bond length   Å
 */
typedef struct {
    int    atom_a;
    int    atom_b;
    int    order;               /* 1=single, 2=double, 3=triple, 0=aromatic */
    double r0;                  /* equilibrium length, Å                     */
    double k;                   /* force constant, eV/Å²                    */
    double current_length;      /* Å  (updated each step)                   */
    double energy;              /* eV (updated each step)                   */
} Bond;

/*
 * A bond angle between three atoms (a–b–c, b is the central atom).
 * V = 0.5 * k * (θ - θ0)²
 */
typedef struct {
    int    atom_a, atom_b, atom_c;
    double theta0;              /* equilibrium angle, radians                */
    double k;                   /* force constant, eV/rad²                  */
} Angle;

/*
 * A dihedral torsion between four atoms (a–b–c–d).
 * V = k * (1 + cos(n*φ - δ))
 */
typedef struct {
    int    atom_a, atom_b, atom_c, atom_d;
    double k;                   /* eV                                        */
    int    n;                   /* multiplicity                              */
    double delta;               /* phase offset, radians                     */
} Dihedral;

/* ══════════════════════════════════════════════════════════════════════════
 * Simulation box
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    Vec3 dimensions;            /* Å — set to a very large value for vacuum  */
    int  periodic[3];           /* 1 = periodic boundary in that axis        */
} SimBox;

/* ══════════════════════════════════════════════════════════════════════════
 * Force-field parameters for non-bonded interactions
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    FF_NONE   = 0,              /* bare Coulomb only                         */
    FF_UFF    = 1,              /* Universal Force Field (LJ params)         */
    FF_OPLS   = 2,              /* OPLS-AA (for organics)                    */
    FF_CUSTOM = 3
} ForceFieldType;

/* ══════════════════════════════════════════════════════════════════════════
 * Thermostat and barostat
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    THERMOSTAT_NONE      = 0,
    THERMOSTAT_BERENDSEN = 1,   /* simple velocity rescaling                 */
    THERMOSTAT_NOSE_HOOVER = 2  /* proper NVT ensemble                       */
} ThermostatType;

typedef struct {
    ThermostatType type;
    double target_temperature; /* K                                          */
    double tau;                /* coupling time constant, fs                 */
    double xi;                 /* Nosé-Hoover friction variable              */
    double Q;                  /* Nosé-Hoover mass parameter                 */
} Thermostat;

/* ══════════════════════════════════════════════════════════════════════════
 * Top-level simulation state
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    /* Particles */
    Atom    *atoms;
    int      num_atoms;
    int      capacity_atoms;

    /* Bonded topology */
    Bond     *bonds;
    int       num_bonds;
    int       capacity_bonds;

    Angle    *angles;
    int       num_angles;
    int       capacity_angles;

    Dihedral *dihedrals;
    int       num_dihedrals;
    int       capacity_dihedrals;

    /* Box */
    SimBox   box;

    /* Time */
    double   time;              /* current time, fs                          */
    double   dt;                /* timestep, fs                              */
    uint64_t step;

    /* Thermodynamics */
    double   temperature;       /* K  (instantaneous)                        */
    double   kinetic_energy;    /* eV                                        */
    double   potential_energy;  /* eV                                        */
    double   total_energy;      /* eV                                        */
    double   virial;            /* eV (for pressure calculation)             */

    /* Force field */
    ForceFieldType ff_type;
    double         cutoff;      /* non-bonded cutoff, Å                      */
    int            use_coulomb;
    int            use_lj;
    int            use_bonds;
    int            use_angles;

    /* Thermostat */
    Thermostat thermostat;
} Simulation;

/* ══════════════════════════════════════════════════════════════════════════
 * Return codes
 * ══════════════════════════════════════════════════════════════════════════ */

typedef enum {
    SIM_OK               =  0,
    SIM_ERR_ALLOC        = -1,
    SIM_ERR_BADATOM      = -2,
    SIM_ERR_BADBOND      = -3,
    SIM_ERR_OVERFLOW     = -4,
    SIM_ERR_BADPARAM     = -5
} SimError;

#endif /* TYPES_H */
