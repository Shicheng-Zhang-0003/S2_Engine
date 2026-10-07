#ifndef TUI_SCREEN_H
#define TUI_SCREEN_H

#include <stddef.h>
#include "types.h"
#include "tui_view.h"

/*
 * tui_screen.h — fullscreen ncurses viewport over the live world.
 *
 * Line mode stays the authority: every keystroke there is POSIX shell
 * text on stdin/stdout, piped and byte-stable, and the batch record
 * never sees a keypress. `screen` is the good-TUI surface for humans
 * on a real Linux terminal: alternate screen, terminfo colors, resize,
 * and the same viewport cells, side panel, status line and ps/vi
 * keybinds as the ANSI grid — both backends consume view_compute(),
 * so the two pictures cannot disagree about what is where.
 *
 * Availability is two-gated: compile-time S2_HAVE_CURSES (ncursesw,
 * optional — the build works without it) and runtime TTY + TERM.
 * Piped, dumb, missing-lib and init-failure all degrade to a typed
 * message with status 0, never a crash and never stdout DATA.
 */

/* 1 when this binary can open a fullscreen (compiled + TTY + TERM ok). */
int screen_available(void);

/* Open/close the alternate screen (termios-free; curses owns the TTY).
 * screen_init returns 0 on success. Safe to call shutdown unmatched. */
int screen_init(void);
void screen_shutdown(void);

/* Current fullscreen dimensions (or 80x24 when closed). */
void screen_size(int *cols, int *rows);

/* The shared keybind footer, so monitor and screen print the same keys. */
void screen_help(char *b1, size_t n1, char *b2, size_t n2);

/* Blit one computed frame plus status/help footer lines. */
void screen_draw(const ViewCells *vc, const char *status,
                 const char *help1, const char *help2, int sel);

/* Key codes: printable passthrough, 1001-1004 arrows (LIVE_KEY_*),
 * -1 timeout, -2 unavailable, -3 resize (redraw). */
int screen_getch(int timeout_ms);

/* One-line prompt: leaves curses, asks on the terminal, resumes. */
int screen_prompt(const char *prompt, char *buf, size_t cap);

#endif /* TUI_SCREEN_H */
