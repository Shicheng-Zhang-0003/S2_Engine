#ifndef TUI_VIEW_H
#define TUI_VIEW_H

#include <stdio.h>
#include "types.h"

/* tui_view.h — orthographic ASCII display grid for the live sim.
 *
 * Two backends share one pure cell model, so the piped ANSI picture
 * and the fullscreen curses picture can never disagree about what is
 * where: view_compute() projects atoms/bonds/restraints into cells,
 * view_render_to() emits the classic ANSI frame to a stream (the
 * record-safe path used by `render`, `ps` snapshots and pipes), and
 * the optional ncurses screen in tui_screen.c consumes the same cells
 * for its viewport pane. No graphics.
 */

typedef struct {
    double yaw_deg;    /* rotation about z (top-down spin) */
    double pitch_deg;  /* rotation about x (tilt) */
    double zoom;       /* scale multiplier */
    Vec3 center;       /* look-at point (A) */
    int has_center;    /* 0 = auto-fit bbox center */
    double slice;      /* half-thickness roasted z-slab (A), 0 = off */
    int auto_render;   /* redraw after commands when stdout is a TTY */
} ViewCam;

void view_cam_reset(ViewCam *cam);

/* Pure cell model: vw x vh grid plus the six side-panel lines.
 * Caller frees with view_cells_free(). Returns NULL on OOM (renders
 * as "empty"). Column/row fitting, glyphs, depth shading and the
 * restraint-x marks are identical for every backend. */
typedef struct {
    int vw, vh;        /* grid dimensions */
    char *ch;          /* vw*vh glyphs, row-major */
    int *fg;           /* vw*vh ANSI base colors 30-37 */
    int *br;           /* vw*vh bright flag (near = bright) */
    char panel[6][96]; /* side-panel lines (NUL-terminated, may be "") */
    char title[64];    /* frame title ("s2 grid" or "empty") */
} ViewCells;

ViewCells *view_compute(const Simulation *sim, const ViewCam *cam,
                        int cols, int rows);

/* Bounding-box center of finite atoms (the auto-fit anchor). Returns 0
 * with *out untouched when there is nothing finite to frame. */
int view_bbox_center(const Simulation *sim, Vec3 *out);

/* Bounded scalar for fixed-width status displays (see tui_view.c). */
void view_fmt_scalar(char *b, size_t n, double v, const char *fmt,
                     const char *unit);
void view_cells_free(ViewCells *c);

/* Render sim into the current terminal (or fixed 80x24 when piped).
 * force_ansi: emit ANSI clears even when not a TTY (for `render`). */
void view_render(const Simulation *sim, const ViewCam *cam, int force_ansi);

/* Same render, written to an explicit stream. The live simulator modes
 * draw to stderr so stdout stays DATA (the stream contract). */
void view_render_to(FILE *fp, const Simulation *sim, const ViewCam *cam,
                     int force_ansi);

#endif /* TUI_VIEW_H */
