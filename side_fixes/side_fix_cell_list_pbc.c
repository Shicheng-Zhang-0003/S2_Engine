#include "side_fix_cell_list_pbc.h"
#include "../include/vec3.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

PbcCellList *side_fix_pbc_cell_list_create(int capacity) {
    if (capacity < 16) capacity = 16;
    PbcCellList *cl = (PbcCellList *)calloc(1, sizeof(PbcCellList));
    if (!cl) return NULL;
    cl->capacity = capacity;
    cl->next = (int *)malloc(capacity * sizeof(int));
    if (!cl->next) { free(cl); return NULL; }
    return cl;
}

void side_fix_pbc_cell_list_destroy(PbcCellList *cl) {
    if (!cl) return;
    if (cl->head) free(cl->head);
    if (cl->next) free(cl->next);
    free(cl);
}

static inline int wrap_cell(int c, int n) {
    int w = c % n;
    return (w < 0) ? w + n : w;
}

int side_fix_pbc_cell_list_build(PbcCellList *cl, const Simulation *sim, double cutoff_A) {
    if (!cl || !sim || !sim->atoms || sim->num_atoms < 1) return -1;
    if (cutoff_A <= 0.0) return -1;

    Vec3 box = sim->box.dimensions;
    if (box.x <= 2.0 * cutoff_A || box.y <= 2.0 * cutoff_A || box.z <= 2.0 * cutoff_A) {
        /* Minimum image convention requires box >= 2 * cutoff */
        return -2;
    }

    /* Resize atom capacity if needed */
    if (sim->num_atoms > cl->capacity) {
        int new_cap = sim->num_atoms * 2;
        int *new_next = (int *)realloc(cl->next, new_cap * sizeof(int));
        if (!new_next) return -1;
        cl->next = new_next;
        cl->capacity = new_cap;
    }

    cl->nx = (int)floor(box.x / cutoff_A);
    cl->ny = (int)floor(box.y / cutoff_A);
    cl->nz = (int)floor(box.z / cutoff_A);
    if (cl->nx < 3) cl->nx = 3;
    if (cl->ny < 3) cl->ny = 3;
    if (cl->nz < 3) cl->nz = 3;

    cl->cell_size_x = box.x / (double)cl->nx;
    cl->cell_size_y = box.y / (double)cl->ny;
    cl->cell_size_z = box.z / (double)cl->nz;

    int total_cells = cl->nx * cl->ny * cl->nz;
    if (total_cells != cl->num_cells) {
        if (cl->head) free(cl->head);
        cl->head = (int *)malloc(total_cells * sizeof(int));
        if (!cl->head) return -1;
        cl->num_cells = total_cells;
    }

    for (int i = 0; i < total_cells; i++) cl->head[i] = -1;

    /* Bin atoms into periodic cells */
    for (int i = 0; i < sim->num_atoms; i++) {
        Vec3 pos = sim->atoms[i].position;
        /* Wrap atom coordinate into primary unit cell [0, L) */
        double px = fmod(pos.x, box.x); if (px < 0.0) px += box.x;
        double py = fmod(pos.y, box.y); if (py < 0.0) py += box.y;
        double pz = fmod(pos.z, box.z); if (pz < 0.0) pz += box.z;

        int cx = (int)(px / cl->cell_size_x);
        int cy = (int)(py / cl->cell_size_y);
        int cz = (int)(pz / cl->cell_size_z);
        if (cx >= cl->nx) cx = cl->nx - 1;
        if (cy >= cl->ny) cy = cl->ny - 1;
        if (cz >= cl->nz) cz = cl->nz - 1;

        int cell_idx = (cx * cl->ny + cy) * cl->nz + cz;
        cl->next[i] = cl->head[cell_idx];
        cl->head[cell_idx] = i;
    }

    return 0;
}

int side_fix_pbc_cell_list_find_pairs(const PbcCellList *cl, const Simulation *sim,
                                      double cutoff_A, PbcPair *out_pairs, int max_pairs) {
    if (!cl || !sim || !out_pairs || max_pairs < 1) return -1;
    double cut2 = cutoff_A * cutoff_A;
    Vec3 box = sim->box.dimensions;
    int pair_count = 0;

    /* Half-stencil neighbor offsets to avoid double-counting pairs */
    int stencil[14][3] = {
        {0, 0, 0},
        {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
        {1, 1, 0}, {1, -1, 0}, {1, 0, 1}, {1, 0, -1},
        {0, 1, 1}, {0, 1, -1},
        {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1}
    };

    for (int cx = 0; cx < cl->nx; cx++) {
        for (int cy = 0; cy < cl->ny; cy++) {
            for (int cz = 0; cz < cl->nz; cz++) {
                int cell1 = (cx * cl->ny + cy) * cl->nz + cz;

                for (int s = 0; s < 14; s++) {
                    int ncx = cx + stencil[s][0];
                    int ncy = cy + stencil[s][1];
                    int ncz = cz + stencil[s][2];

                    /* Compute integer wrap shift */
                    double shift_x = 0.0, shift_y = 0.0, shift_z = 0.0;
                    if (ncx >= cl->nx) shift_x = box.x;
                    else if (ncx < 0)  shift_x = -box.x;

                    if (ncy >= cl->ny) shift_y = box.y;
                    else if (ncy < 0)  shift_y = -box.y;

                    if (ncz >= cl->nz) shift_z = box.z;
                    else if (ncz < 0)  shift_z = -box.z;

                    int wcx = wrap_cell(ncx, cl->nx);
                    int wcy = wrap_cell(ncy, cl->ny);
                    int wcz = wrap_cell(ncz, cl->nz);
                    int cell2 = (wcx * cl->ny + wcy) * cl->nz + wcz;

                    int i = cl->head[cell1];
                    while (i >= 0) {
                        int j = (cell1 == cell2) ? cl->next[i] : cl->head[cell2];
                        while (j >= 0) {
                            Vec3 pos_i = sim->atoms[i].position;
                            Vec3 pos_j = sim->atoms[j].position;

                            /* Wrap coordinates into primary unit cell [0, L) */
                            double pix = fmod(pos_i.x, box.x); if (pix < 0.0) pix += box.x;
                            double piy = fmod(pos_i.y, box.y); if (piy < 0.0) piy += box.y;
                            double piz = fmod(pos_i.z, box.z); if (piz < 0.0) piz += box.z;

                            double pjx = fmod(pos_j.x, box.x); if (pjx < 0.0) pjx += box.x;
                            double pjy = fmod(pos_j.y, box.y); if (pjy < 0.0) pjy += box.y;
                            double pjz = fmod(pos_j.z, box.z); if (pjz < 0.0) pjz += box.z;

                            double dx = (pjx + shift_x) - pix;
                            double dy = (pjy + shift_y) - piy;
                            double dz = (pjz + shift_z) - piz;

                            /* Minimum image wrap */
                            dx -= box.x * round(dx / box.x);
                            dy -= box.y * round(dy / box.y);
                            dz -= box.z * round(dz / box.z);

                            double r2 = dx * dx + dy * dy + dz * dz;
                            if (r2 <= cut2 && r2 > 1e-12) {
                                if (pair_count < max_pairs) {
                                    out_pairs[pair_count].atom_i = i;
                                    out_pairs[pair_count].atom_j = j;
                                    out_pairs[pair_count].dist_A = sqrt(r2);
                                    out_pairs[pair_count].r_ij = vec3(dx, dy, dz);
                                    pair_count++;
                                }
                            }
                            j = cl->next[j];
                        }
                        i = cl->next[i];
                    }
                }
            }
        }
    }

    return pair_count;
}
