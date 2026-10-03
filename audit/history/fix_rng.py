#!/usr/bin/env python3
"""
Improve the random number generator in integrator.c to use PCG64
instead of the simple LCG.
"""

with open('TREE/src/integrator.c', 'r') as f:
    content = f.read()

# Replace the LCG with PCG64
old_rng = '''static double rand_uniform(uint64_t *state) {
    /* LCG with Knuth constants — good enough for velocity initialisation */
    *state = *state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(*state >> 33) / (double)(1ULL << 31);
}

static double rand_normal(uint64_t *state) {
    double u1, u2;
    do { u1 = rand_uniform(state); } while (u1 < 1.0e-10);
    u2 = rand_uniform(state);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}'''

new_rng = '''/* PCG64 random number generator (O'Neill, 2014) — better statistical
 * properties than LCG, still fast and lightweight. */
static uint64_t pcg64_state = 0;
static uint64_t pcg64_inc = 1;

static uint64_t pcg64_next(void) {
    uint64_t old = pcg64_state;
    pcg64_state = old * 6364136223846793005ULL + (pcg64_inc | 1);
    uint64_t xorshifted = ((old >> 18u) ^ old) >> 27u;
    uint64_t rot = old >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 63));
}

static double rand_uniform(void) {
    return (double)(pcg64_next() >> 11) * (1.0 / 9007199254740992.0);  // 1/2^53
}

static double rand_normal(void) {
    double u1, u2;
    do { u1 = rand_uniform(); } while (u1 < 1.0e-10);
    u2 = rand_uniform();
    return sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159265358979323846 * u2);
}

/* Seed the PCG64 generator */
static void pcg64_seed(uint64_t seed) {
    pcg64_state = 0;
    pcg64_inc = (seed << 1) | 1;
    pcg64_next();
    pcg64_state += seed;
    pcg64_next();
}'''

content = content.replace(old_rng, new_rng)

# Update the integrator_maxwell_boltzmann function to use the new RNG
old_maxwell = '''void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
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
    }'''

new_maxwell = '''void integrator_maxwell_boltzmann(Simulation *sim, double T_init,
                                    unsigned long seed) {
    if (!sim || !sim->atoms || sim->num_atoms < 1) return;
    if (!isfinite(T_init) || T_init < 0.0) return;
    pcg64_seed((seed == 0) ? 12345678901234567ULL : (uint64_t)seed);
    sim->rng_state = pcg64_state;  /* Keep for backward compatibility */

    for (int i = 0; i < sim->num_atoms; i++) {
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        /* σ_v [Å/fs] = sqrt(k_B T [eV] / (m [AMU] × AMU_AFS2_TO_EV)) */
        double sigma_v = sqrt(KB_EV * T_init / (a->mass * AMU_AFS2_TO_EV));
        if (!isfinite(sigma_v)) sigma_v = 0.0;

        a->velocity.x = rand_normal() * sigma_v;
        a->velocity.y = rand_normal() * sigma_v;
        a->velocity.z = rand_normal() * sigma_v;
    }'''

content = content.replace(old_maxwell, new_maxwell)

# Update Andersen thermostat to use new RNG
old_andersen = '''    int kicked = 0;
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
    }'''

new_andersen = '''    int kicked = 0;
    for (int i = 0; i < sim->num_atoms; i++) {
        if (rand_uniform() >= p) continue;
        Atom *a = &sim->atoms[i];
        if (!(a->mass > 1e-12) || !isfinite(a->mass)) continue;
        double sigma_v = sqrt(KB_EV * T0 / (a->mass * AMU_AFS2_TO_EV));
        if (!isfinite(sigma_v)) continue;
        a->velocity.x = rand_normal() * sigma_v;
        a->velocity.y = rand_normal() * sigma_v;
        a->velocity.z = rand_normal() * sigma_v;
        kicked = 1;
    }'''

content = content.replace(old_andersen, new_andersen)

with open('TREE/src/integrator.c', 'w') as f:
    f.write(content)

print("RNG updated to PCG64!")