#ifndef PERIODIC_TABLE_H
#define PERIODIC_TABLE_H

#include "types.h"

/*
 * periodic_table.h
 * Immutable database of element properties used by the simulator.
 * Elements 1–36 (H → Kr) are fully populated.
 * Elements 37–118 are stubbed with mass only — extend as needed.
 *
 * Data sources:
 *   Pauling electronegativities : IUPAC 2013
 *   Atomic/covalent/vdW radii  : Alvarez 2013 (Dalton Trans.)
 *   Ionisation energies        : NIST Atomic Spectra Database
 *   Electron affinities        : Andersen 1999 (J. Phys. Chem. Ref. Data)
 *   LJ parameters              : Rappé 1992 UFF (J. Am. Chem. Soc.)
 */

/* The global table — indexed by atomic number Z (1-based, element[0] unused) */
extern const Element PERIODIC_TABLE[MAX_ELEMENTS + 1];

/* ── Accessors ───────────────────────────────────────────────────────────── */

/* Returns pointer to element data for atomic number Z, or NULL if invalid. */
const Element *pt_element(int Z);

/* Look up by symbol string (case-sensitive, e.g. "C", "Fe"). */
const Element *pt_by_symbol(const char *symbol);

/* Populate a ground-state ElectronConfig for atomic number Z. */
void pt_electron_config(int Z, ElectronConfig *cfg);

/* Pretty-print electron configuration string, e.g. "1s2 2s2 2p2" */
void pt_config_string(const ElectronConfig *cfg, char *buf, int buflen);

/* Print a summary of one element to stdout */
void pt_print_element(int Z);

#endif /* PERIODIC_TABLE_H */
