#ifndef DISPLAY_H
#define DISPLAY_H

#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#ifdef _WIN32
#include <io.h>
#define DISPLAY_ISATTY(fd) 0
#else
#include <unistd.h>
#include <stdlib.h>
#define DISPLAY_ISATTY(fd) isatty(fd)
#endif

/*
 * display.h — deterministic-stdout + live-TTY display split.
 *
 * RULE (record integrity): stdout must stay byte-deterministic across
 * runs, compilers (-march aside), and wall time. Therefore wall time,
 * progress bars, and heartbeats may ONLY go to stderr, and ONLY when
 * stderr is a TTY and CARBON_QUIET is unset. Record runs redirect both
 * streams to files (not TTYs) → silent, deterministic, s01-clean.
 * All functions are best-effort and never fail the program.
 */

/* FULL-AUDIT C31: re-read per call (isatty + getenv are negligible vs MD).
 * The old cache-once froze late CARBON_QUIET changes and hid threaded state. */
static int __attribute__((unused)) display_live(void) {
    return DISPLAY_ISATTY(STDERR_FILENO) && (getenv("CARBON_QUIET") == NULL);
}

static void __attribute__((unused)) progress(const char *fmt, ...) {
    va_list ap;
    if (!display_live()) return;
    va_start(ap, fmt);
    fputs("  [..] ", stderr);
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    va_end(ap);
}

/* Heartbeat fraction helper: returns 1 when iter hits a 10% mark. */
static int __attribute__((unused)) progress_mark(int iter, int max_iter) {
    int mark = max_iter / 10;
    if (mark < 1) mark = 1;
    return (iter % mark) == 0 || iter + 1 >= max_iter;
}

#endif /* DISPLAY_H */
