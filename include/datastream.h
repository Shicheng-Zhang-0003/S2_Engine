#ifndef DATASTREAM_H
#define DATASTREAM_H

#include <stddef.h>

/*
* datastream.h — structured datastream writer, schema 1.
*
* Format contract: DATASTREAM_SPEC.md. One .cvmds file per run.
* Sections: [header] (key: value lines), [atoms], [steps], [claims],
* and [end] carrying a SHA-256 seal over everything before it.
*
* Call order:
*   ds_open -> ds_set_header* -> ds_add_atom* -> ds_add_step* -> ds_add_claim* -> ds_close
* Header fields must be set before any atoms/steps/claims are added.
*/

typedef struct DSWriter DSWriter;

/* Open a writer for `path` with demo identifier `demo` (e.g. "kcsa").
* Writes the fixed [header] fields (schema, demo, build identity,
* wallclock) immediately. Returns NULL on allocation failure. */
DSWriter *ds_open(const char *path, const char *demo);

/* Append an extra header field, e.g. ("rng-seed", "42"). */
void ds_set_header(DSWriter *w, const char *key, const char *value);

/* Append one row to [atoms]: idx Z symbol mass charge lj_eps lj_sigma. */
void ds_add_atom(DSWriter *w, int idx, int Z, const char *symbol,
                 double mass_amu, double charge_e,
                 double lj_eps_ev, double lj_sigma_a);

/* Append one row to [steps]: step time_fs KE PE E T. */
void ds_add_step(DSWriter *w, long step, double time_fs,
                 double ke_ev, double pe_ev, double e_ev, double t_k);

/* Append one row to [claims]: key value unit provenance. */
void ds_add_claim(DSWriter *w, const char *key, double value,
                  const char *unit, const char *provenance);

/* Seal the payload (SHA-256 into [end]) and write the file.
* Returns 0 on success, -1 on failure. Frees the writer either way. */
int ds_close(DSWriter *w);

/* SHA-256 of len bytes at data, as 64 lowercase hex chars + NUL. */
void ds_sha256_hex(const void *data, size_t len, char out65[65]);

/* Re-read a .cvmds file and verify its payload seal and structure.
* Returns 0 if intact, -1 otherwise. */
int ds_verify_file(const char *path);

#endif /* DATASTREAM_H */
