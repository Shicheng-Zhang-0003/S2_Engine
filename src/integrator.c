#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../include/integrator.h"
#include "../include/forces.h"
#include "../include/constants.h"
#include "../include/display.h"

/*
 * integrator.c
 *
 * Velocity Verlet + Berendsen thermostat.
 * All unit conversions are explicit and commented.
 */

#define KB_EV   (BOLTZMANN_K / EV_TO_J)   /* eV/K, derived in-line per C1-C4 rule */


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
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    if (!isfinite(sim->dt)) return;
    double dt   = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        if (!isfinite(a->force.x) || !isfinite(a->force.y) || !isfinite(a->force.z)) continue;

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
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    if (!isfinite(sim->dt)) return;
    double dt = sim->dt;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        if (!isfinite(a->force.x) || !isfinite(a->force.y) || !isfinite(a->force.z)) continue;
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
    else if (sim->thermostat.type == THERMOSTAT_ANDERSEN)
        integrator_andersen(sim);

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
    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    double ke = 0.0;
    for (int i = 0; i < sim->num_atoms; i++) {
        const Atom *a = &sim->atoms[i];
        double v2 = vec3_norm2(a->velocity);
        if (!isfinite(v2) || !isfinite(a->mass)) continue;
        ke += 0.5 * a->mass * v2;
    }
    return ke * AMU_AFS2_TO_EV;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Instantaneous temperature from equipartition theorem
 *
 * <KE> = dof/2 × k_B × T
 *   → T = 2 × KE / (dof × k_B)
 * dof = 3N - 3 - constrained (see integrator.h for what counts as
 * constrained). Guarded to a minimum of 1 so a fully-pinned system
 * reports 0 KE / 1-dof rather than dividing by zero. (History: an
 * earlier version used 3N, underestimating T by 33% for N=3 water -
 * the COM correction is the -3; `constrained` generalizes it.)
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_temperature(const Simulation *sim) {
    if (!sim || sim->num_atoms < 2) return 0.0;
    double ke  = integrator_kinetic_energy(sim);
    if (!isfinite(ke)) return 0.0;
    long dof = 3L * (long)sim->num_atoms - 3L - (long)sim->num_constrained_dof;
    if (dof < 1) dof = 1;
    return (2.0 * ke) / ((double)dof * KB_EV);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Berendsen thermostat
 *
 * Rescales velocities: v_new = λ v
 *   λ = sqrt(1 + (dt/τ)(T_0/T − 1))
 *
 * arg clamp: if 1 + (dt/τ)(T_0/T - 1) < 0 (extreme overshoot), clamp
 * to 0.01 rather than 0.0 — setting λ=0 would stop ALL atomic motion,
 * making the next step's T also 0, then the thermostat stays locked
 * at λ=0 indefinitely. The minimum clamp λ = 0.1 (arg = 0.01) still
 * cools strongly while keeping motion alive so the thermostat recovers.
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_berendsen(Simulation *sim) {
    if (!sim || !sim->atoms) return;
    double T     = integrator_temperature(sim);
    double T0    = sim->thermostat.target_temperature;
    double tau   = sim->thermostat.tau;
    double dt    = sim->dt;

    if (!(T > 1.0e-10) || !isfinite(T)) return;
    if (!(tau > 1e-12) || !isfinite(tau)) return;
    if (!isfinite(T0) || !isfinite(dt)) return;

    double ratio = T0 / T;
    double arg   = 1.0 + (dt / tau) * (ratio - 1.0);
    if (!isfinite(arg)) return;
    if (arg < 0.01) arg = 0.01;          /* clamp: avoids motion lockout */
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

static double rand_uniform(uint64_t *state) {
    /* LCG with Knuth constants — good enough for velocity initialisation */
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*state >> 33) / (double)(1ULL << 31);
}

static double rand_normal(uint64_t *state) {
    double u1, u2;
    do { u1 = rand_uniform(state); } while (u1 < 1.0e-10);
    u2 = rand_uniform(state);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
                                   unsigned long seed) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    if (!isfinite(T_init) || T_init < 0.0) return;
    sim->rng_state = (seed == 0) ? 12345678901234567ULL : (uint64_t)seed;

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        /* σ_v [Å/fs] = sqrt(k_B T [eV] / (m [AMU] × AMU_AFS2_TO_EV)) */
        double sigma_v = sqrt(KB_EV * T_init / (a->mass * AMU_AFS2_TO_EV));
        if (!isfinite(sigma_v)) sigma_v = 0.0;

        a->velocity.x = rand_normal(&sim->rng_state) * sigma_v;
        a->velocity.y = rand_normal(&sim->rng_state) * sigma_v;
        a->velocity.z = rand_normal(&sim->rng_state) * sigma_v;
    }

    /* Remove centre-of-mass drift */
    integrator_remove_com_velocity(sim);

    /* Rescale to exact target temperature */
    double T_actual = integrator_temperature(sim);
    if (T_actual > 1.0e-10 && isfinite(T_actual) && T_init > 0.0) {
        double scale = sqrt(T_init / T_actual);
        if (!isfinite(scale)) return;
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
 * Andersen thermostat (Andersen, JCP 72, 2384 (1980))
 *
 * Each atom independently suffers a stochastic "collision" with the heat
 * bath with probability p = 1 - exp(-nu*dt) per step; collided atoms get
 * fresh Maxwell-Boltzmann velocities at T0 (Box-Muller, per-atom sigma).
 * Unlike Berendsen rescaling this generates the RIGOROUS canonical (NVT)
 * ensemble — the right thermostat for free-energy sampling (WHAM legs).
 * COM is re-zeroed after kicking (keeps the 3N-3 temperature exact for
 * small systems; documented perturbation of the ensemble, negligible).
 * Default nu = 1/tau if nu <= 0 (tau reuse, documented).
 * ══════════════════════════════════════════════════════════════════════════ */
void integrator_andersen(Simulation *sim) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    double T0 = sim->thermostat.target_temperature;
    double dt = sim->dt;
    double nu = sim->thermostat.nu;
    if (!(nu > 0.0) || !isfinite(nu)) {
        if (!(sim->thermostat.tau > 1e-12) || !isfinite(sim->thermostat.tau)) return;
        nu = 1.0 / sim->thermostat.tau;
    }
    if (!isfinite(T0) || T0 < 0.0 || !isfinite(dt) || dt <= 0.0) return;
    double p = 1.0 - exp(-nu * dt);
    if (!(p > 0.0) || !isfinite(p)) return;
    if (p > 1.0) p = 1.0;
    int kicked = 0;
    for (int i = 0; i < sim->num_atoms; i++) {
        if (rand_uniform(&sim->rng_state) >= p) continue;
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        double sigma_v = sqrt(KB_EV * T0 / (a->mass * AMU_AFS2_TO_EV));
        if (!isfinite(sigma_v)) continue;
        a->velocity.x = rand_normal(&sim->rng_state) * sigma_v;
        a->velocity.y = rand_normal(&sim->rng_state) * sigma_v;
        a->velocity.z = rand_normal(&sim->rng_state) * sigma_v;
        kicked = 1;
    }
    if (kicked) integrator_remove_com_velocity(sim);
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

/* ══════════════════════════════════════════════════════════════════════════
 * Steepest-descent energy minimization with adaptive step size
 *
 * Includes a hard per-atom displacement cap and a divergence safety
 * check - both necessary, confirmed by testing: without them, a
 * single large step can "tunnel through" a steep LJ repulsive wall
 * into the unphysical region where opposite charges approach r->0 and
 * Coulomb attraction diverges, and the minimizer will "successfully"
 * descend toward that fake, non-physical minimum (observed directly:
 * two atoms collapsed to exactly 0.0000 A while PE plunged to
 * -396,041 eV before these safeguards were added).
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_minimize(Simulation *sim, int max_iterations,
                            double initial_step, double force_tolerance) {
    const double MAX_DISPLACEMENT = 0.05;
    const double DIVERGENCE_THRESHOLD = -50000.0;

    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    if (max_iterations < 1) return sim->potential_energy;
    if (!isfinite(initial_step) || initial_step <= 0.0) return sim->potential_energy;
    if (!isfinite(force_tolerance) || force_tolerance < 0.0) return sim->potential_energy;
    if ((size_t)sim->num_atoms > (size_t)2147483647 / sizeof(Vec3)) return sim->potential_energy;
    double step_size = initial_step;
    forces_calculate(sim);
    double E_current = sim->potential_energy;

    Vec3 *saved_positions = (Vec3 *)malloc(sizeof(Vec3) * (size_t)sim->num_atoms);
    if (!saved_positions) return E_current;

    for (int iter = 0; iter < max_iterations; iter++) {
        if (progress_mark(iter, max_iterations))
            progress("minimize %d/%d (E=%.3f)", iter, max_iterations, E_current);
        double max_force = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) {
            double f = vec3_norm(sim->atoms[i].force);
            if (f > max_force) max_force = f;
        }
        if (max_force < force_tolerance) break;

        for (int i = 0; i < sim->num_atoms; i++)
            saved_positions[i] = sim->atoms[i].position;

        double effective_step = (step_size < MAX_DISPLACEMENT) ? step_size : MAX_DISPLACEMENT;
        double scale = effective_step / max_force;

        for (int i = 0; i < sim->num_atoms; i++) {
            Vec3 disp = vec3_scale(sim->atoms[i].force, scale);
            double disp_len = vec3_norm(disp);
            if (disp_len > MAX_DISPLACEMENT)
                disp = vec3_scale(disp, MAX_DISPLACEMENT / disp_len);
            sim->atoms[i].position = vec3_add(sim->atoms[i].position, disp);
        }

        forces_calculate(sim);
        double E_new = sim->potential_energy;

        if (E_new < DIVERGENCE_THRESHOLD || isnan(E_new)) {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].position = saved_positions[i];
            forces_calculate(sim);
            free(saved_positions);
            return sim->potential_energy;
        }

        if (E_new < E_current) {
            E_current = E_new;
            step_size *= 1.2;
        } else {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].position = saved_positions[i];
            forces_calculate(sim);
            step_size *= 0.5;
            if (step_size < 1.0e-8) break;
        }
    }

    free(saved_positions);
    return sim->potential_energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Frozen-atom variant - same safety mechanisms, but atoms with
 * frozen[i] != 0 never move.
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_minimize_frozen(Simulation *sim, const int *frozen,
                                   int max_iterations, double initial_step,
                                   double force_tolerance) {
    const double MAX_DISPLACEMENT = 0.05;
    const double DIVERGENCE_THRESHOLD = -50000.0;

    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    if (!frozen) return integrator_minimize(sim, max_iterations, initial_step, force_tolerance);
    if (max_iterations < 1) return sim->potential_energy;
    if (!isfinite(initial_step) || initial_step <= 0.0) return sim->potential_energy;
    if (!isfinite(force_tolerance) || force_tolerance < 0.0) return sim->potential_energy;
    if ((size_t)sim->num_atoms > (size_t)2147483647 / sizeof(Vec3)) return sim->potential_energy;
    double step_size = initial_step;
    forces_calculate(sim);
    double E_current = sim->potential_energy;

    Vec3 *saved_positions = (Vec3 *)malloc(sizeof(Vec3) * (size_t)sim->num_atoms);
    if (!saved_positions) return E_current;

    for (int iter = 0; iter < max_iterations; iter++) {
        if (progress_mark(iter, max_iterations))
            progress("minimize %d/%d (E=%.3f)", iter, max_iterations, E_current);
        double max_force = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) {
            if (frozen[i]) continue;
            double f = vec3_norm(sim->atoms[i].force);
            if (f > max_force) max_force = f;
        }
        if (max_force < force_tolerance) break;

        for (int i = 0; i < sim->num_atoms; i++)
            saved_positions[i] = sim->atoms[i].position;

        double effective_step = (step_size < MAX_DISPLACEMENT) ? step_size : MAX_DISPLACEMENT;
        double scale = effective_step / max_force;

        for (int i = 0; i < sim->num_atoms; i++) {
            if (frozen[i]) continue;
            Vec3 disp = vec3_scale(sim->atoms[i].force, scale);
            double disp_len = vec3_norm(disp);
            if (disp_len > MAX_DISPLACEMENT)
                disp = vec3_scale(disp, MAX_DISPLACEMENT / disp_len);
            sim->atoms[i].position = vec3_add(sim->atoms[i].position, disp);
        }

        forces_calculate(sim);
        double E_new = sim->potential_energy;

        if (E_new < DIVERGENCE_THRESHOLD || isnan(E_new)) {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].position = saved_positions[i];
            forces_calculate(sim);
            free(saved_positions);
            return sim->potential_energy;
        }

        if (E_new < E_current) {
            E_current = E_new;
            step_size *= 1.2;
        } else {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].position = saved_positions[i];
            forces_calculate(sim);
            step_size *= 0.5;
            if (step_size < 1.0e-8) break;
        }
    }

    free(saved_positions);
    return sim->potential_energy;
}

/* ══════════════════════════════════════════════════════════════════════════
 * FIRE minimization (Bitzek-Koskinen-Gähler-Moseler-Parrinello 2006)
 *
 * Inertial MD with velocity mixing: v -> (1-a)v + a|v|Fhat whenever
 * power P = F.v > 0 (accelerate along downhill direction, dt grows),
 * full stop (v = 0, dt shrinks) on uphill. Same safety rails as the
 * steepest minimizers above: per-step displacement cap, -50000 eV
 * divergence floor with rollback, NaN guards. Velocities and dt are
 * saved on entry and restored on exit (FIRE owns them while running).
 * ══════════════════════════════════════════════════════════════════════════ */
double integrator_fire(Simulation *sim, int max_iterations,
                       double dt_start, double force_tolerance) {
    const double MAX_DISPLACEMENT = 0.05;
    const double DIVERGENCE_THRESHOLD = -50000.0;
    const double F_INC = 1.1, F_DEC = 0.5, F_ALPHA = 0.99;
    const double ALPHA_START = 0.1;
    const int N_MIN = 5;

    if (!sim || !sim->atoms || sim->num_atoms < 1) return 0.0;
    if (max_iterations < 1) return sim->potential_energy;
    if (!isfinite(dt_start) || dt_start <= 0.0 || dt_start > 5.0) return sim->potential_energy;
    if (!isfinite(force_tolerance) || force_tolerance < 0.0) return sim->potential_energy;
    if ((size_t)sim->num_atoms > (size_t)2147483647 / sizeof(Vec3)) return sim->potential_energy;

    Vec3 *saved_vel = (Vec3 *)malloc(sizeof(Vec3) * (size_t)sim->num_atoms);
    Vec3 *saved_pos = (Vec3 *)malloc(sizeof(Vec3) * (size_t)sim->num_atoms);
    if (!saved_vel || !saved_pos) { free(saved_vel); free(saved_pos); return sim->potential_energy; }
    for (int i = 0; i < sim->num_atoms; i++) saved_vel[i] = sim->atoms[i].velocity;
    double dt_saved = sim->dt;

    double dt = dt_start;
    double dt_max = 10.0 * dt_start;
    if (dt_max > 2.0) dt_max = 2.0;
    double alpha = ALPHA_START;
    int n_accel = 0;

    forces_calculate(sim);
    for (int i = 0; i < sim->num_atoms; i++)
        sim->atoms[i].velocity = vec3_zero();

    int iter;
    for (iter = 0; iter < max_iterations; iter++) {
        if (progress_mark(iter, max_iterations))
            progress("fire %d/%d (E=%.3f)", iter, max_iterations, sim->potential_energy);
        double max_force = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) {
            double f = vec3_norm(sim->atoms[i].force);
            if (f > max_force) max_force = f;
        }
        if (max_force < force_tolerance) break;
        if (!isfinite(max_force) || max_force < 1e-300) break;

        /* Power P = F.v (current velocities). */
        double P = 0.0;
        for (int i = 0; i < sim->num_atoms; i++)
            P += vec3_dot(sim->atoms[i].force, sim->atoms[i].velocity);
        if (!isfinite(P)) P = -1.0;

        /* Velocity mixing toward forces. */
        double vnorm = 0.0, fnorm = 0.0;
        for (int i = 0; i < sim->num_atoms; i++) {
            vnorm += vec3_norm2(sim->atoms[i].velocity);
            fnorm += vec3_norm2(sim->atoms[i].force);
        }
        vnorm = sqrt(vnorm);
        fnorm = sqrt(fnorm);
        if (isfinite(vnorm) && isfinite(fnorm) && vnorm > 0.0 && fnorm > 0.0) {
            for (int i = 0; i < sim->num_atoms; i++) {
                Vec3 mix = vec3_add(vec3_scale(sim->atoms[i].velocity, 1.0 - alpha),
                                    vec3_scale(sim->atoms[i].force, alpha * vnorm / fnorm));
                sim->atoms[i].velocity = isfinite(mix.x + mix.y + mix.z) ? mix : vec3_zero();
            }
        }

        if (P > 0.0) {
            n_accel++;
            if (n_accel > N_MIN) {
                dt *= F_INC;
                if (dt > dt_max) dt = dt_max;
                alpha *= F_ALPHA;
            }
        } else {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].velocity = vec3_zero();
            dt *= F_DEC;
            alpha = ALPHA_START;
            n_accel = 0;
        }

        for (int i = 0; i < sim->num_atoms; i++)
            saved_pos[i] = sim->atoms[i].position;
        sim->dt = dt;
        integrator_kick_drift(sim);
        /* Displacement cap: rescale the whole step if any atom flew. */
        {
            double worst = 0.0;
            for (int i = 0; i < sim->num_atoms; i++) {
                double d = vec3_dist(sim->atoms[i].position, saved_pos[i]);
                if (d > worst) worst = d;
            }
            if (worst > MAX_DISPLACEMENT && worst > 0.0) {
                double s = MAX_DISPLACEMENT / worst;
                for (int i = 0; i < sim->num_atoms; i++) {
                    Vec3 dd = vec3_sub(sim->atoms[i].position, saved_pos[i]);
                    sim->atoms[i].position = vec3_add(saved_pos[i], vec3_scale(dd, s));
                }
                for (int i = 0; i < sim->num_atoms; i++)
                    vec3_iscale(&sim->atoms[i].velocity, s);
            }
        }
        forces_calculate(sim);
        integrator_kick(sim);

        double E_new = sim->potential_energy;
        if (E_new < DIVERGENCE_THRESHOLD || isnan(E_new) || !isfinite(E_new)) {
            for (int i = 0; i < sim->num_atoms; i++)
                sim->atoms[i].position = saved_pos[i];
            forces_calculate(sim);
            break;
        }
    }

    for (int i = 0; i < sim->num_atoms; i++)
        sim->atoms[i].velocity = saved_vel[i];
    sim->dt = dt_saved;
    forces_calculate(sim);
    free(saved_vel);
    free(saved_pos);
    return sim->potential_energy;
}
