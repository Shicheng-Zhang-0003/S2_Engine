#include <math.h>
#include <stdio.h>
#include "../include/integrator.h"
#include "../include/forces.h"
#include "../include/constants.h"

/*
 * integrator.c
 *
 * Velocity Verlet + Berendsen thermostat.
 * All unit conversions are explicit and commented.
 */

#define KB_EV   8.617333262e-5   /* Boltzmann constant in eV/K */

/* ══════════════════════════════════════════════════════════════════════════
 * Half-step A: velocity kick + position drift
 *
 * Physics:
 *   a_i = F_i / m_i                    [eV/Å / AMU]
 *   a_i [Å/fs²] = a_i [eV/Å/AMU] × MD_FORCE_CONV
 *
 *   v_i(t+dt/2) = v_i(t) + 0.5 × a_i × dt
 *   r_i(t+dt)   = r_i(t) + v_i(t+dt/2) × dt
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_kick_drift(Simulation *sim) {
    double dt   = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];

        /* a [Å/fs²] = F [eV/Å] / m [AMU] × MD_FORCE_CONV */
        double inv_m = MD_FORCE_CONV / a->mass;

        /* Half-kick: v += 0.5 a dt */
        a->velocity.x += 0.5 * a->force.x * inv_m * dt;
        a->velocity.y += 0.5 * a->force.y * inv_m * dt;
        a->velocity.z += 0.5 * a->force.z * inv_m * dt;

        /* Drift: r += v dt  (using half-kicked velocity) */
        a->position.x += a->velocity.x * dt;
        a->position.y += a->velocity.y * dt;
        a->position.z += a->velocity.z * dt;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Half-step B: final velocity kick with new forces
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_kick(Simulation *sim) {
    double dt = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        double inv_m = MD_FORCE_CONV / a->mass;

        a->velocity.x += 0.5 * a->force.x * inv_m * dt;
        a->velocity.y += 0.5 * a->force.y * inv_m * dt;
        a->velocity.z += 0.5 * a->force.z * inv_m * dt;
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Full velocity Verlet step
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_step(Simulation *sim) {
    /* 1. Kick velocities by half-step, drift positions by full step */
    integrator_kick_drift(sim);

    /* 2. Recalculate forces at new positions */
    forces_calculate(sim);

    /* 3. Complete velocity update with new forces */
    integrator_kick(sim);

    /* 4. Apply thermostat if active */
    if (sim->thermostat.type == THERMOSTAT_BERENDSEN)
        integrator_berendsen(sim);

    /* 5. Update thermodynamics */
    sim->kinetic_energy = integrator_kinetic_energy(sim);
    sim->temperature    = integrator_temperature(sim);
    sim->total_energy   = sim->kinetic_energy + sim->potential_energy;

    sim->step++;
    sim->time += sim->dt;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Kinetic energy
 *
 * KE [eV] = 0.5 × Σ m_i [AMU] × |v_i|² [Å²/fs²] × AMU_AFS2_TO_EV
 *
 * AMU_AFS2_TO_EV = 1 AMU × (Å/fs)² in eV
 *   = 1.66053906660e-27 kg × (1e5 m/s)²
 *   = 1.66053906660e-27 × 1e10 J
 *   = 1.66053906660e-17 / 1.602176634e-19 eV
 *   = 103.6427 eV
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_kinetic_energy(const Simulation *sim) {
    double ke = 0.0;
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        double v2 = vec3_norm2(a->velocity);
        ke += 0.5 * a->mass * v2;
    }
    return ke * AMU_AFS2_TO_EV;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Instantaneous temperature from equipartition theorem
 *
 * <KE> = (3N/2) k_B T  →  T = 2 KE / (3N k_B)
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_temperature(const Simulation *sim) {
    if (sim->num_atoms == 0) return 0.0;
    double ke = integrator_kinetic_energy(sim);
    int dof = 3 * sim->num_atoms;   /* degrees of freedom */
    return (2.0 * ke) / ((double)dof * KB_EV);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Berendsen thermostat
 *
 * Rescales velocities: v_new = λ v
 *   λ = sqrt(1 + (dt/τ)(T_0/T − 1))
 *
 * Safe for λ < 0: clamp (1 + ...) to ≥ 0 to avoid sqrt of negative.
 * This happens when T >> T_0; the clamp is unphysical but prevents crash.
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_berendsen(Simulation *sim) {
    double T     = integrator_temperature(sim);
    double T0    = sim->thermostat.target_temperature;
    double tau   = sim->thermostat.tau;
    double dt    = sim->dt;

    if (T < 1.0e-10) return;  /* avoid divide-by-zero at T≈0 */

    double ratio = T0 / T;
    double arg   = 1.0 + (dt / tau) * (ratio - 1.0);
    if (arg < 0.0) arg = 0.0;
    double lambda = sqrt(arg);

    for (int i = 0; i < sim->num_atoms; i++)
        vec3_iscale(&sim->atoms[i].velocity, lambda);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Maxwell-Boltzmann velocity initialisation
 *
 * Box-Muller transform: given U1, U2 uniform in (0,1),
 *   Z0 = sqrt(-2 ln U1) cos(2π U2)  ~  N(0,1)
 *   Z1 = sqrt(-2 ln U1) sin(2π U2)  ~  N(0,1)
 *
 * Thermal speed for component x:
 *   σ_v = sqrt(k_B T / m)  in Å/fs
 *   v_x = Z × σ_v
 *
 * σ_v [Å/fs]: k_B T [eV] / m [AMU] × (1/AMU_AFS2_TO_EV)
 * ══════════════════════════════════════════════════════════════════════════ */
static unsigned long rng_state;

static double rand_uniform(void) {
    /* LCG with Knuth constants — good enough for velocity initialisation */
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(rng_state >> 33) / (double)(1ULL << 31);
}

static double rand_normal(void) {
    double u1, u2;
    do { u1 = rand_uniform(); } while (u1 < 1.0e-10);
    u2 = rand_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
                                   unsigned long seed) {
    rng_state = (seed == 0) ? 12345678901234567ULL : (unsigned long long)seed;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        /* σ_v [Å/fs] = sqrt(k_B T [eV] / (m [AMU] × AMU_AFS2_TO_EV)) */
        double sigma_v = sqrt(KB_EV * T_init / (a->mass * AMU_AFS2_TO_EV));

        a->velocity.x = rand_normal() * sigma_v;
        a->velocity.y = rand_normal() * sigma_v;
        a->velocity.z = rand_normal() * sigma_v;
    }

    /* Remove centre-of-mass drift */
    integrator_remove_com_velocity(sim);

    /* Rescale to exact target temperature */
    double T_actual = integrator_temperature(sim);
    if (T_actual > 1.0e-10) {
        double scale = sqrt(T_init / T_actual);
        for (int i = 0; i < sim->num_atoms; i++)
            vec3_iscale(&sim->atoms[i].velocity, scale);
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Remove centre-of-mass velocity
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_remove_com_velocity(Simulation *sim) {
    double total_mass = 0.0;
    Vec3   p_com      = vec3_zero();

    for (int i = 0; i < sim->num_atoms; i++) {
        double m = sim->atoms[i].mass;
        total_mass += m;
        vec3_iadd(&p_com, vec3_scale(sim->atoms[i].velocity, m));
    }
    if (total_mass < 1.0e-10) return;

    Vec3 v_com = vec3_scale(p_com, 1.0 / total_mass);
    for (int i = 0; i < sim->num_atoms; i++)
        vec3_isub(&sim->atoms[i].velocity, v_com);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Print step summary
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_print_step(const Simulation *sim) {
    printf("  Step %6llu  t=%8.2f fs  "
           "KE=%9.5f  PE=%9.5f  E=%9.5f eV  T=%7.2f K\n",
           (unsigned long long)sim->step,
           sim->time,
           sim->kinetic_energy,
           sim->potential_energy,
           sim->total_energy,
           sim->temperature);
}
