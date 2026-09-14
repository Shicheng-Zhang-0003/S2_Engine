#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/datastream.h"

/*
* test_datastream.c — selftest for the schema-1 datastream writer.
*
* Checks:
*   1. SHA-256 known-answer vectors, including the empty-string hash
*      (identical to the EMPTY_SHA constant the s01 record harness
*      uses — a built-in cross-check against the system sha256sum).
*   2. Full write path across all four sections.
*   3. Structural re-read of the written file.
*   4. Claim value round-trip through the text format.
*   5. Seal verification, plus rejection of a tampered copy.
*/

static int failures = 0;

static void check(int cond, const char *name) {
    if (cond) printf("  PASS  %s\n", name);
    else { printf("  FAIL  %s\n", name); failures++; }
}

int main(void) {
    char hex[65];
    const char *path = "build/selftest.cvmds";
    const char *tpath = "build/selftest_tampered.cvmds";
    DSWriter *w;
    FILE *f;
    char *content;
    long sz;

    printf("datastream selftest (schema 1)\n");

    /* 1. SHA-256 known-answer tests */
    ds_sha256_hex("", 0, hex);
    check(strcmp(hex,
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855") == 0,
        "sha256 of empty input matches the project EMPTY_SHA");
    ds_sha256_hex("abc", 3, hex);
    check(strcmp(hex,
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad") == 0,
        "sha256 of \"abc\" matches the FIPS 180-4 test vector");

    /* 2. Write a small datastream exercising every section */
    w = ds_open(path, "selftest");
    check(w != NULL, "ds_open succeeds");
    if (!w) { printf("selftest: cannot continue\n"); return 1; }
    ds_set_header(w, "rng-seed", "42");
    ds_set_header(w, "source-hash", "deadbeef");
    ds_add_atom(w, 0, 8, "O", 15.9990, -0.834, 0.006803, 3.15061);
    ds_add_atom(w, 1, 1, "H", 1.008, 0.417, 0.0, 0.0);
    ds_add_step(w, 0, 0.0, 0.077556, 0.0, 0.077556, 300.0);
    ds_add_claim(w, "selftest.value", 8.936812, "eV", "computed");
    ds_add_claim(w, "selftest.reference", -0.1786, "eV", "expt-1000:1@300K");
    check(ds_close(w) == 0, "ds_close writes the file");

    /* 3. Seal verification on the untampered file */
    check(ds_verify_file(path) == 0, "ds_verify_file accepts the untampered file");

    /* 4. Structural re-read */
    f = fopen(path, "rb");
    content = NULL;
    sz = 0;
    if (f) {
        fseek(f, 0, SEEK_END);
        sz = ftell(f);
        fseek(f, 0, SEEK_SET);
        content = (char *)malloc((size_t)sz + 1);
        if (content) {
            if (fread(content, 1, (size_t)sz, f) != (size_t)sz) {
                free(content); content = NULL;
            } else content[sz] = '\0';
        }
        fclose(f);
    }
    check(content != NULL, "written file reads back");
    if (content) {
        check(strncmp(content, "[header]\n", 9) == 0,
              "file starts with [header]");
        check(strstr(content, "schema: 1\n") != NULL, "schema: 1 present");
        check(strstr(content, "rng-seed: 42\n") != NULL,
              "custom header field present");
        check(strstr(content, "\n[atoms]\n") != NULL, "[atoms] section present");
        check(strstr(content, "\n[steps]\n") != NULL, "[steps] section present");
        check(strstr(content, "\n[claims]\n") != NULL, "[claims] section present");
        check(strstr(content, "\n[end]\npayload-sha256: ") != NULL,
              "[end] with payload-sha256 present");
        /* Atom row: same bit-exact parsed comparison as the step row
         * below - %.17g output is only substring-stable for values
         * whose decimal form is short AND exactly representable. */
        {
            char *al = strstr(content, "\n[atoms]\n");
            int aidx = -1, az = -1;
            char asym[16] = {0};
            double am = 0, aq = 0, ae = 0, as = 0;
            int n = 0;
            if (al) n = sscanf(al, "\n[atoms]\n# idx Z symbol mass_amu charge_e lj_eps_ev lj_sigma_a\n"
                                   "%d %d %15s %lg %lg %lg %lg",
                               &aidx, &az, asym, &am, &aq, &ae, &as);
            check(n == 7 && aidx == 0 && az == 8 && strcmp(asym, "O") == 0 &&
                  am == 15.999 && aq == -0.834 && ae == 0.006803 && as == 3.15061,
                  "atom row round-trips exactly (bit-identical doubles)");
        }
        /* Step row: parsed field-by-field with EXACT double equality.
         * %.17g guarantees strtod restores the identical binary double
         * the writer held, so == here proves a lossless round trip
         * (a %.10g-era substring check could not distinguish a
         * truncated value from the true one). */
        {
            char *sl = strstr(content, "\n[steps]\n");
            long pst = -1;
            double tm = 0, ke = 0, pe = 0, en = 0, tk = 0;
            int n = 0;
            if (sl) n = sscanf(sl, "\n[steps]\n# step time_fs KE_ev PE_ev E_ev T_k\n"
                                   "%ld %lg %lg %lg %lg %lg",
                               &pst, &tm, &ke, &pe, &en, &tk);
            check(n == 6 && pst == 0 && tm == 0.0 && ke == 0.077556 &&
                  pe == 0.0 && en == 0.077556 && tk == 300.0,
                  "step row round-trips exactly (bit-identical doubles)");
        }
        {
            char *cl = strstr(content, "selftest.value ");
            double v = 0.0;
            if (cl) v = atof(cl + strlen("selftest.value "));
            check(cl != NULL && fabs(v - 8.936812) < 1e-9,
                  "claim value round-trips exactly");
        }

        /* 5. Tamper test: flip one payload byte, seal must reject */
        {
            char *t = (char *)malloc((size_t)sz + 1);
            FILE *tf;
            memcpy(t, content, (size_t)sz + 1);
            t[12] ^= 0x01;   /* inside the header text, before [end] */
            tf = fopen(tpath, "wb");
            if (tf) { fwrite(t, 1, (size_t)sz, tf); fclose(tf); }
            free(t);
            check(ds_verify_file(tpath) != 0,
                  "ds_verify_file rejects a tampered copy");
        }
        free(content);
    }

    if (failures) {
        printf("selftest: %d FAILURE(S)\n", failures);
        return 1;
    }
    printf("selftest: 17 checks, all passed\n");
    return 0;
}
