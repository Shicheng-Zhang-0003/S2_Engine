/* verify_cvmds.c — FULL-AUDIT O16: thin seal checker that calls the real
 * ds_verify_file() instead of re-implementing the payload-slice rule in
 * python (which forks on any future framing change). Exit 0 = sealed. */
#include <stdio.h>
#include "../include/datastream.h"

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: verify_cvmds <file.cvmds>\n"); return 2; }
    int rc = ds_verify_file(argv[1]);
    if (rc == 0) printf("SEALED %s\n", argv[1]);
    else printf("UNSEALED %s (rc=%d)\n", argv[1], rc);
    return rc == 0 ? 0 : 1;
}
