#ifndef SIDE_FIX_CELL_LIST_PBC_H
#define SIDE_FIX_CELL_LIST_PBC_H

#include "../include/types.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * side_fix_cell_list_pbc.h
 *
 * Computational fix for periodic non-bonded pair searches:
 *
 * Current limitation:
 *   The engine's cell list spatial partitioning (make perf) only supports
 *   non-periodic vacuum boxes; simulations with periodic boundary conditions
 *   (PBC) fall back to serial O(N^2) all-pairs minimum image checks.
 *
 * Fix implemented here:
 *   Periodic toroidal cell list that wraps across periodic boundaries in O(N).
 *   Enables scalable O(N) pair searches for bulk solvent and condensed-phase
 *   periodic simulations without O(N^2) scaling bottlenecks.
 */

typedef struct {
    int atom_i;
    int atom_j;
    double dist_A;
    Vec3 r_ij;
} PbcPair;

typedef struct {
    int nx, ny, nz;
    double cell_size_x, cell_size_y, cell_size_z;
    int num_cells;
    int *head;      /* head atom index in each cell */
    int *next;      /* linked-list of atoms in each cell */
    int capacity;
} PbcCellList;

/* Create and initialize periodic cell list */
PbcCellList *side_fix_pbc_cell_list_create(int capacity);
void side_fix_pbc_cell_list_destroy(PbcCellList *cl);

/* Build cell list from simulation with periodic box */
int side_fix_pbc_cell_list_build(PbcCellList *cl, const Simulation *sim, double cutoff_A);

/* Enumerate all non-bonded pairs within cutoff using periodic cell list */
int side_fix_pbc_cell_list_find_pairs(const PbcCellList *cl, const Simulation *sim,
                                      double cutoff_A, PbcPair *out_pairs, int max_pairs);

#ifdef __cplusplus
}
#endif

#endif /* SIDE_FIX_CELL_LIST_PBC_H */
