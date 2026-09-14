#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <time.h>
#include "../include/datastream.h"

/*
* datastream.c — structured datastream writer, schema 1.
*
* Implements DATASTREAM_SPEC.md: a text, line-oriented, self-describing
* record format with a SHA-256 integrity seal over the payload.
*
* Design notes:
*  - The entire payload is buffered in memory and written at ds_close().
*    These files are small (kB scale), and buffering guarantees the
*    hashed bytes and the written bytes are identical by construction.
*  - SHA-256 is implemented here (FIPS 180-4), self-contained, and is
*    verified against known-answer vectors in tests/test_datastream.c —
*    including the empty-string hash, which is the same EMPTY_SHA the
*    s01 record harness uses, giving a cross-check against sha256sum.
*  - Build identity (build-compiler, build-flags) is injected by the
*    makefile via -DDS_COMPILER_ID / -DDS_CFLAGS_ID; builds that
*    override CFLAGS on the command line (ASan, debug) fall back to
*    "unknown" via the #ifndef defaults below.
*  - Section markers ([header], [atoms], [steps], [claims], [end]) are
*    reserved: claim keys and provenance strings must not contain them.
*/

#ifndef DS_COMPILER_ID
#define DS_COMPILER_ID "unknown"
#endif
#ifndef DS_CFLAGS_ID
#define DS_CFLAGS_ID "unknown"
#endif

struct DSWriter {
    char   path[512];
    char   demo[64];
    char  *buf;
    size_t len;
    size_t cap;
    long   n_atoms, n_steps, n_claims;
    int    failed;
};

/* ── SHA-256 (FIPS 180-4), self-contained ─────────────────────────── */

typedef struct {
    uint32_t state[8];
    uint64_t bitlen;
    uint8_t  buf[64];
    size_t   buflen;
} SHA256;

static const uint32_t K256[64] = {
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,
    0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
    0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
    0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
    0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,
    0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,
    0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
    0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
    0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
    0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,
    0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,
    0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
    0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
    0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

static uint32_t rotr32(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}

static void sha256_init(SHA256 *s) {
    s->state[0]=0x6a09e667u; s->state[1]=0xbb67ae85u;
    s->state[2]=0x3c6ef372u; s->state[3]=0xa54ff53au;
    s->state[4]=0x510e527fu; s->state[5]=0x9b05688cu;
    s->state[6]=0x1f83d9abu; s->state[7]=0x5be0cd19u;
    s->bitlen = 0;
    s->buflen = 0;
}

static void sha256_block(SHA256 *s, const uint8_t *p) {
    uint32_t w[64];
    uint32_t a,b,c,d,e,f,g,h;
    int i;
    for (i = 0; i < 16; i++)
        w[i] = ((uint32_t)p[i*4] << 24) | ((uint32_t)p[i*4+1] << 16) |
               ((uint32_t)p[i*4+2] << 8) | (uint32_t)p[i*4+3];
    for (i = 16; i < 64; i++) {
        uint32_t s0 = rotr32(w[i-15],7) ^ rotr32(w[i-15],18) ^ (w[i-15] >> 3);
        uint32_t s1 = rotr32(w[i-2],17) ^ rotr32(w[i-2],19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    a=s->state[0]; b=s->state[1]; c=s->state[2]; d=s->state[3];
    e=s->state[4]; f=s->state[5]; g=s->state[6]; h=s->state[7];
    for (i = 0; i < 64; i++) {
        uint32_t S1  = rotr32(e,6) ^ rotr32(e,11) ^ rotr32(e,25);
        uint32_t ch  = (e & f) ^ ((~e) & g);
        uint32_t t1  = h + S1 + ch + K256[i] + w[i];
        uint32_t S0  = rotr32(a,2) ^ rotr32(a,13) ^ rotr32(a,22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2  = S0 + maj;
        h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    s->state[0]+=a; s->state[1]+=b; s->state[2]+=c; s->state[3]+=d;
    s->state[4]+=e; s->state[5]+=f; s->state[6]+=g; s->state[7]+=h;
}

static void sha256_update(SHA256 *s, const void *data, size_t len) {
    const uint8_t *p = (const uint8_t *)data;
    s->bitlen += (uint64_t)len * 8u;
    while (len > 0) {
        size_t n = 64 - s->buflen;
        if (n > len) n = len;
        memcpy(s->buf + s->buflen, p, n);
        s->buflen += n; p += n; len -= n;
        if (s->buflen == 64) { sha256_block(s, s->buf); s->buflen = 0; }
    }
}

static void sha256_final(SHA256 *s, uint8_t out[32]) {
    uint64_t bitlen = s->bitlen;
    int i;
    s->buf[s->buflen++] = 0x80;
    if (s->buflen > 56) {
        while (s->buflen < 64) s->buf[s->buflen++] = 0;
        sha256_block(s, s->buf);
        s->buflen = 0;
    }
    while (s->buflen < 56) s->buf[s->buflen++] = 0;
    for (i = 7; i >= 0; i--)
        s->buf[s->buflen++] = (uint8_t)(bitlen >> (i * 8));
    sha256_block(s, s->buf);
    for (i = 0; i < 8; i++) {
        out[i*4]   = (uint8_t)(s->state[i] >> 24);
        out[i*4+1] = (uint8_t)(s->state[i] >> 16);
        out[i*4+2] = (uint8_t)(s->state[i] >> 8);
        out[i*4+3] = (uint8_t)(s->state[i]);
    }
}

void ds_sha256_hex(const void *data, size_t len, char out65[65]) {
    SHA256 s;
    uint8_t d[32];
    int i;
    sha256_init(&s);
    sha256_update(&s, data, len);
    sha256_final(&s, d);
    for (i = 0; i < 32; i++) sprintf(out65 + 2*i, "%02x", d[i]);
    out65[64] = '\0';
}

/* ── payload buffer ────────────────────────────────────────────────── */

static void ds_append(DSWriter *w, const char *fmt, ...) {
    va_list ap, ap2;
    int n;
    va_start(ap, fmt);
    va_copy(ap2, ap);
    n = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (n < 0) { va_end(ap2); w->failed = 1; return; }
    if (w->len + (size_t)n + 1 > w->cap) {
        size_t ncap = w->cap ? w->cap : 4096;
        char *nb;
        while (ncap < w->len + (size_t)n + 1) ncap *= 2;
        nb = (char *)realloc(w->buf, ncap);
        if (!nb) { va_end(ap2); w->failed = 1; return; }
        w->buf = nb;
        w->cap = ncap;
    }
    vsnprintf(w->buf + w->len, (size_t)n + 1, fmt, ap2);
    va_end(ap2);
    w->len += (size_t)n;
}

/* ── public API ────────────────────────────────────────────────────── */

DSWriter *ds_open(const char *path, const char *demo) {
    DSWriter *w = (DSWriter *)calloc(1, sizeof(DSWriter));
    time_t now;
    struct tm *tmv;
    char stamp[32];
    if (!w) return NULL;
    if (snprintf(w->path, sizeof w->path, "%s", path) >= (int)sizeof w->path ||
        snprintf(w->demo, sizeof w->demo, "%s", demo) >= (int)sizeof w->demo) {
        free(w);
        return NULL; /* path/demo must fit - truncation would fork identity */
    }
    now = time(NULL);
    tmv = gmtime(&now);
    if (!tmv) { free(w); return NULL; }
    strftime(stamp, sizeof stamp, "%Y-%m-%dT%H:%M:%SZ", tmv);
    ds_append(w, "[header]\n");
    ds_append(w, "schema: 1\n");
    ds_append(w, "demo: %s\n", w->demo);
    ds_append(w, "build-compiler: %s\n", DS_COMPILER_ID);
    ds_append(w, "build-flags: %s\n", DS_CFLAGS_ID);
    ds_append(w, "wallclock: %s\n", stamp);
    /* NOTE (spec compliance): DATASTREAM_SPEC.md also lists source-hash
     * and rng-seed as header fields. They are run-specific, so the
     * CALLER provides them via ds_set_header() (the s42 consumer does);
     * ds_open writes only what it can know itself. */
    if (w->failed) { free(w); return NULL; }
    return w;
}

/* Section-marker hygiene: keys, symbols, units and provenance tags must
 * not contain anything the file parser treats structurally ('\n', '[',
 * ']'). An smuggled "\n[end]\n" inside a claim key would otherwise fork
 * the payload the verifier hashes (fail-closed today, but fragile).
 * Rejected inputs set the failed flag so ds_close() writes nothing. */
static int ds_text_ok(const char *s) {
    if (!s || !*s) return 0;
    for (; *s; s++)
        if (*s == '\n' || *s == '\r' || *s == '[' || *s == ']') return 0;
    return 1;
}

void ds_set_header(DSWriter *w, const char *key, const char *value) {
    if (!w || w->failed) return;
    if (!ds_text_ok(key) || !value || strchr(value, '\n')) { w->failed = 1; return; }
    ds_append(w, "%s: %s\n", key, value);
}

void ds_add_atom(DSWriter *w, int idx, int Z, const char *symbol,
                 double mass_amu, double charge_e,
                 double lj_eps_ev, double lj_sigma_a) {
    if (!w || w->failed) return;
    if (!ds_text_ok(symbol)) { w->failed = 1; return; }
    if (w->n_atoms++ == 0)
        ds_append(w, "\n[atoms]\n"
                     "# idx Z symbol mass_amu charge_e lj_eps_ev lj_sigma_a\n");
    /* %.17g: round-trip precision. %.10g silently truncated doubles
     * (demonstrated: pi -> 3.141592654), forcing every future validator
     * to budget ~1e-10 relative slop that is pure serialization loss,
     * not physics. 17 significant digits round-trip exactly. */
    ds_append(w, "%d %d %s %.17g %.17g %.17g %.17g\n",
              idx, Z, symbol, mass_amu, charge_e, lj_eps_ev, lj_sigma_a);
}

void ds_add_step(DSWriter *w, long step, double time_fs,
                 double ke_ev, double pe_ev, double e_ev, double t_k) {
    if (!w || w->failed) return;
    if (w->n_steps++ == 0)
        ds_append(w, "\n[steps]\n"
                     "# step time_fs KE_ev PE_ev E_ev T_k\n");
    ds_append(w, "%ld %.17g %.17g %.17g %.17g %.17g\n",
              step, time_fs, ke_ev, pe_ev, e_ev, t_k);
}

void ds_add_claim(DSWriter *w, const char *key, double value,
                  const char *unit, const char *provenance) {
    if (!w || w->failed) return;
    if (!ds_text_ok(key) || !ds_text_ok(unit) || !ds_text_ok(provenance)) {
        w->failed = 1;
        return;
    }
    if (w->n_claims++ == 0)
        ds_append(w, "\n[claims]\n"
                     "# key value unit provenance\n");
    ds_append(w, "%s %.17g %s %s\n", key, value, unit, provenance);
}

int ds_close(DSWriter *w) {
    char hex[65];
    FILE *f;
    int rc = 0;
    if (!w) return -1;
    if (!w->failed) {
        ds_sha256_hex(w->buf, w->len, hex);
        ds_append(w, "\n[end]\npayload-sha256: %s\n", hex);
    }
    if (w->failed) rc = -1;
    else {
        f = fopen(w->path, "wb");
        if (!f || fwrite(w->buf, 1, w->len, f) != w->len) rc = -1;
        if (f) fclose(f);
    }
    free(w->buf);
    free(w);
    return rc;
}

int ds_verify_file(const char *path) {
    FILE *f = fopen(path, "rb");
    long sz;
    char *buf;
    char *endline, *endline2, *hashfield, *hexstart;
    char hex[65];
    int rc = -1;
    if (!f) return -1;
    fseek(f, 0, SEEK_END);
    sz = ftell(f);
    fseek(f, 0, SEEK_SET);
    if (sz < 90) { fclose(f); return -1; }
    buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return -1; }
    if (fread(buf, 1, (size_t)sz, f) != (size_t)sz) {
        free(buf); fclose(f); return -1;
    }
    fclose(f);
    buf[sz] = '\0';
    if (strncmp(buf, "[header]\n", 9) == 0) {
        endline = strstr(buf, "\n[end]\n");
        /* Exactly one [end] section: a second occurrence means payload
         * content is masquerading as structure (only possible if the
         * writer's marker hygiene was bypassed) - reject. */
        endline2 = endline ? strstr(endline + 1, "\n[end]\n") : NULL;
        hashfield = strstr(buf, "\npayload-sha256: ");
        if (endline && !endline2 && hashfield && hashfield > endline) {
            /* payload = everything before the newline that precedes [end] */
            size_t payload_len = (size_t)(endline - buf);
            hexstart = hashfield + strlen("\npayload-sha256: ");
            /* The hash field must be exactly 64 lowercase hex digits
             * followed by newline-or-EOF: a prefix-only comparison
             * would accept trailing garbage after a valid hash. */
            int hexok = 1;
            for (int i = 0; i < 64; i++) {
                char c = hexstart[i];
                if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) {
                    hexok = 0;
                    break;
                }
            }
            if (hexok && (hexstart[64] == '\n' || hexstart[64] == '\0')) {
                ds_sha256_hex(buf, payload_len, hex);
                if (strncmp(hexstart, hex, 64) == 0)
                    rc = 0;
            }
        }
    }
    free(buf);
    return rc;
}
