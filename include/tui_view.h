#ifndef TUI_VIEW_H
#define TUI_VIEW_H

#include "types.h"

/* tui_view.h — Minecraft-style ASCII display grid for the live sim.
 * ANSI only (no ncurses): orthographic projection of atoms+bonds onto a
 * text grid with depth shading, plus a side panel. No graphics.
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

/* Render sim into the current terminal (or fixed 80x24 when piped).
 * force_ansi: emit ANSI clears even when not a TTY (for `render`). */
void view_render(const Simulation *sim, const ViewCam *cam, int force_ansi);

#endif /* TUI_VIEW_H */
