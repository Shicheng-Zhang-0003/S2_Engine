#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/ioctl.h>

#include "../include/tui_view.h"
#include "../include/periodic_table.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void view_cam_reset(ViewCam *cam) {
    if (!cam) return;
    cam->yaw_deg = 25.0;
    cam->pitch_deg = 18.0;
    cam->zoom = 1.0;
    cam->center = vec3_zero();
    cam->has_center = 0;
    cam->slice = 0.0;
    cam->auto_render = 1;
}

/* element glyph (case distinguishes N/n) + ANSI base color 30-37 */
static char glyph_for_z(int z) {
    switch (z) {
        case 1: return 'h';
        case 6: return 'c';
        case 7: return 'n';
        case 8: return 'o';
        case 11: return 'N';
        case 15: return 'p';
        case 19: return 'K';
        default: return '*';
    }
}

static int color_for_z(int z) {
    switch (z) {
        case 1: return 37; /* H white */
        case 6: return 32; /* C green */
        case 7: return 34; /* N blue */
        case 8: return 31; /* O red */
        case 11: return 33; /* Na yellow */
        case 15: return 35; /* P magenta */
        case 19: return 35; /* K magenta */
        default: return 36; /* other cyan */
    }
}

typedef struct { double x, y, z; int idx; } Proj;

/* Bounded scalar for fixed-width displays: raw %f of 1e100 prints a
 * hundred digits. fmt holds one %f/%e-style verb (e.g. "%.1f"). */
void view_fmt_scalar(char *b, size_t n, double v, const char *fmt,
                     const char *unit) {
    char f[16];
    if (!b || n == 0) return;
    if (!isfinite(v)) { snprintf(b, n, "---%s", unit ? unit : ""); return; }
    if (fabs(v) >= 1e6) snprintf(f, sizeof f, "%%.3e");
    else snprintf(f, sizeof f, "%s", fmt);
    {
        char nb[64];
        snprintf(nb, sizeof nb, f, v);
        snprintf(b, n, "%s%s", nb, unit ? unit : "");
    }
}

static void term_size(int *cols, int *rows) {
    struct winsize ws;
    if (cols) *cols = 80;
    if (rows) *rows = 24;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0) {
        if (ws.ws_col >= 40) *cols = ws.ws_col;
        if (ws.ws_row >= 15) *rows = ws.ws_row;
    }
}

ViewCells *view_compute(const Simulation *sim, const ViewCam *cam,
                        int cols, int rows) {
    ViewCells *vc;
    int panel_w = 30;
    int vw, vh;
    double yaw, pitch, cy, sy, cp, sp;
    int n, i, b, r, x, y;
    Proj *pr = NULL;
    int *order = NULL;
    int *skip = NULL;
    Vec3 ctr;
    ViewCam c;
    double zmin = 1e30, zmax = -1e30, fitsc = 1.0;

    if (cols < 60) panel_w = 20;
    vw = cols - panel_w - 3;
    if (vw < 20) vw = 20;
    vh = rows - 5;
    if (vh < 8) vh = 8;

    vc = (ViewCells *)calloc(1, sizeof(ViewCells));
    if (!vc) return NULL;
    vc->vw = vw;
    vc->vh = vh;
    vc->ch = (char *)malloc((size_t)(vw * vh));
    vc->fg = (int *)malloc(sizeof(int) * (size_t)(vw * vh));
    vc->br = (int *)malloc(sizeof(int) * (size_t)(vw * vh));
    if (!vc->ch || !vc->fg || !vc->br) {
        view_cells_free(vc);
        return NULL;
    }
    memset(vc->ch, ' ', (size_t)(vw * vh));
    for (i = 0; i < vw * vh; i++) { vc->fg[i] = 37; vc->br[i] = 0; }
    for (i = 0; i < 6; i++) vc->panel[i][0] = '\0';
    snprintf(vc->title, sizeof vc->title, "%s", sim ? "s2 grid" : "empty");

    if (cam) c = *cam;
    else view_cam_reset(&c);
    yaw = c.yaw_deg * M_PI / 180.0;
    pitch = c.pitch_deg * M_PI / 180.0;
    cy = cos(yaw); sy = sin(yaw);
    cp = cos(pitch); sp = sin(pitch);

    n = (sim && sim->atoms) ? sim->num_atoms : 0;
    if (n < 0) n = 0;
    if (n > 0) {
        pr = (Proj *)malloc(sizeof(Proj) * (size_t)n);
        order = (int *)malloc(sizeof(int) * (size_t)n);
        skip = (int *)calloc((size_t)n, sizeof(int));
        if (!pr || !order || !skip) {
    free(pr); free(order);
    /* blown count before freeing skip: panel reports excluded atoms */
    {
        int blown = 0;
        if (skip) for (int q = 0; q < n; q++) blown += skip[q] ? 1 : 0;
        if (blown < 0) blown = 0;
        vc->panel[5][0] = '\0';
        if (blown > 0)
            snprintf(vc->panel[5], sizeof vc->panel[5], "view %s%s +%d blown",
                     c.has_center ? "centered" : "fit",
                     c.slice > 0 ? " slice" : "", blown);
        else
            snprintf(vc->panel[5], sizeof vc->panel[5], "view %s%s",
                     c.has_center ? "centered" : "fit",
                     c.slice > 0 ? " slice" : "");
    }
    free(skip);
            return vc;
        }
    }
    ctr = c.center;
    if (!c.has_center && pr && sim && n > 0) {
        int seeded = 0;
        Vec3 mn = vec3_zero(), mx = vec3_zero();
        for (int i = 0; i < n; i++) {
            Vec3 p = sim->atoms[i].position;
            /* non-finite or escaped atoms never seed or widen the frame:
             * one NaN position otherwise poisons the center and flings
             * every glyph (including the neuron readout) to a corner,
             * and one parsec-distant atom collapses the fit to a line.
             * 1e6 A is four orders past any box; beyond it is blown. */
            if (!isfinite(p.x + p.y + p.z) ||
                fabs(p.x) > 1e6 || fabs(p.y) > 1e6 || fabs(p.z) > 1e6) continue;
            if (!seeded) { mn = mx = p; seeded = 1; continue; }
            if (p.x < mn.x) mn.x = p.x;
            if (p.x > mx.x) mx.x = p.x;
            if (p.y < mn.y) mn.y = p.y;
            if (p.y > mx.y) mx.y = p.y;
            if (p.z < mn.z) mn.z = p.z;
            if (p.z > mx.z) mx.z = p.z;
        }
        if (seeded) ctr = vec3_scale(vec3_add(mn, mx), 0.5);
    }
    if (pr && sim) {
        for (i = 0; i < n; i++) {
            Vec3 d = vec3_sub(sim->atoms[i].position, ctr);
            double x1, y1, y2, z2;
            if (!isfinite(d.x + d.y + d.z) ||
                fabs(d.x) > 1e6 || fabs(d.y) > 1e6 || fabs(d.z) > 1e6) {
                /* blown-up atom: never plotted, never widens the frame */
                skip[i] = 1;
                pr[i].x = 0; pr[i].y = 0; pr[i].z = 0; pr[i].idx = i;
                order[i] = i;
                continue;
            }
            x1 = d.x * cy - d.y * sy;
            y1 = d.x * sy + d.y * cy;
            y2 = y1 * cp - d.z * sp;
            z2 = y1 * sp + d.z * cp;
            pr[i].x = x1; pr[i].y = y2; pr[i].z = z2; pr[i].idx = i;
            order[i] = i;
            if (c.slice > 0.0 && fabs(z2) > c.slice) continue;
            if (z2 < zmin) zmin = z2;
            if (z2 > zmax) zmax = z2;
        }
        /* fit scale */
        {
            double bx = 1e-9, by = 1e-9;
            for (i = 0; i < n; i++) {
                if (skip[i]) continue;
                if (c.slice > 0.0 && fabs(pr[i].z) > c.slice) continue;
                if (fabs(pr[i].x) > bx) bx = fabs(pr[i].x);
                if (fabs(pr[i].y) > by) by = fabs(pr[i].y);
            }
            {
                double span = bx * 2.0 + 2.0, spanY = by * 2.0 + 2.0;
                double sc = vw / span;
                double scY = (vh * 2.0) / spanY; /* cells ~2x tall */
                if (scY < sc) sc = scY;
                sc *= c.zoom;
                if (!(sc > 1e-9) || !isfinite(sc)) sc = 1.0;
                fitsc = sc;
            }
        }
        for (i = 0; i < n; i++) { pr[i].x *= fitsc; pr[i].y *= fitsc; }
        /* far (small z) first */
        for (i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pr[order[j]].z < pr[order[i]].z) { int t = order[i]; order[i] = order[j]; order[j] = t; }
    }
    if (!(zmax >= zmin)) { zmin = 0; zmax = 1; }

    /* bonds first (dim dots) */
    if (sim && sim->bonds && pr) {
        for (b = 0; b < sim->num_bonds; b++) {
            int ia = sim->bonds[b].atom_a, ib = sim->bonds[b].atom_b;
            int x0, y0, x1, y1, dx, sx, dy, sy2, err, xx, yy, guard, e2;
            if (ia < 0 || ia >= n || ib < 0 || ib >= n) continue;
            if (skip[ia] || skip[ib]) continue;
            if (c.slice > 0.0 && (fabs(pr[ia].z) > c.slice || fabs(pr[ib].z) > c.slice)) continue;
            x0 = (int)lround(vw / 2.0 + pr[ia].x);
            y0 = (int)lround(vh / 2.0 - pr[ia].y / 2.0);
            x1 = (int)lround(vw / 2.0 + pr[ib].x);
            y1 = (int)lround(vh / 2.0 - pr[ib].y / 2.0);
            dx = abs(x1 - x0); sx = x0 < x1 ? 1 : -1;
            dy = -abs(y1 - y0); sy2 = y0 < y1 ? 1 : -1;
            err = dx + dy; xx = x0; yy = y0; guard = vw + vh + 4;
            while (guard-- > 0) {
                if (xx >= 0 && xx < vw && yy >= 0 && yy < vh && vc->ch[yy * vw + xx] == ' ') {
                    vc->ch[yy * vw + xx] = '.';
                    vc->fg[yy * vw + xx] = 37; vc->br[yy * vw + xx] = 0;
                }
                if (xx == x1 && yy == y1) break;
                e2 = 2 * err;
                if (e2 >= dy) { err += dy; xx += sx; }
                if (e2 <= dx) { err += dx; yy += sy2; }
            }
        }
    }
    /* atoms far->near */
    if (pr && sim) {
        for (int k = 0; k < n; k++) {
            int idx = order[k];
            double t;
            if (skip[idx]) continue;
            if (c.slice > 0.0 && fabs(pr[idx].z) > c.slice) continue;
            x = (int)lround(vw / 2.0 + pr[idx].x);
            y = (int)lround(vh / 2.0 - pr[idx].y / 2.0);
            if (x < 0 || x >= vw || y < 0 || y >= vh) continue;
            t = (zmax > zmin) ? (pr[idx].z - zmin) / (zmax - zmin) : 1.0;
            vc->ch[y * vw + x] = glyph_for_z(sim->atoms[idx].Z);
            vc->fg[y * vw + x] = color_for_z(sim->atoms[idx].Z);
            vc->br[y * vw + x] = (t > 0.66) ? 1 : 0; /* near = bright */
        }
    }
    /* restraint anchors as x (same fit scale) */
    /* FULL-AUDIT C33: check the arrays alongside the count (corrupt/partial
     * sims must not segfault the renderer). */
    if (sim && pr && sim->restraint_anchor && sim->restraint_k) {
        for (r = 0; r < sim->num_restraints; r++) {
            Vec3 d = vec3_sub(sim->restraint_anchor[r], ctr);
            double x1, y1, y2;
            if (!isfinite(d.x + d.y + d.z)) continue;
            x1 = d.x * cy - d.y * sy;
            y1 = d.x * sy + d.y * cy;
            y2 = y1 * cp - d.z * sp;
            x = (int)lround(vw / 2.0 + x1 * fitsc);
            y = (int)lround(vh / 2.0 - (y2 * fitsc) / 2.0);
            if (x >= 0 && x < vw && y >= 0 && y < vh && vc->ch[y * vw + x] == ' ') {
                vc->ch[y * vw + x] = 'x'; vc->fg[y * vw + x] = 33; vc->br[y * vw + x] = 0;
            }
        }
    }
    free(pr); free(order); free(skip);

    /* side panel (six lines, same text as the classic frame).
     * Scalars are bounded: a blown-up T=1e100 once printed a hundred
     * digits and buried the frame. Non-finite reads "---", absurd
     * magnitudes read scientific — the layout never overflows. */
    snprintf(vc->panel[0], sizeof vc->panel[0], "N=%d", sim ? sim->num_atoms : 0);
    if (sim) {
        char tb[32], eb[32];
        view_fmt_scalar(tb, sizeof tb, sim->temperature, "%.1f", "K");
        view_fmt_scalar(eb, sizeof eb, sim->potential_energy, "%.3f", "");
        snprintf(vc->panel[1], sizeof vc->panel[1], "step=%llu",
                 (unsigned long long)sim->step);
        snprintf(vc->panel[2], sizeof vc->panel[2], "T=%s E=%s",
                 tb, eb);
    }
    snprintf(vc->panel[3], sizeof vc->panel[3], "yaw=%.0f pit=%.0f z=%.2f",
             c.yaw_deg, c.pitch_deg, c.zoom);
    snprintf(vc->panel[4], sizeof vc->panel[4], "%s",
             "h c n o N p K =H C N O Na P K");
    return vc;
}

void view_cells_free(ViewCells *c) {
    if (!c) return;
    free(c->ch); free(c->fg); free(c->br);
    free(c);
}

int view_bbox_center(const Simulation *sim, Vec3 *out) {
    int seeded = 0;
    Vec3 mn = vec3_zero(), mx = vec3_zero();
    int i, n;
    if (!sim || !sim->atoms || !out) return 0;
    n = sim->num_atoms;
    for (i = 0; i < n; i++) {
        Vec3 p = sim->atoms[i].position;
        if (!isfinite(p.x + p.y + p.z) ||
            fabs(p.x) > 1e6 || fabs(p.y) > 1e6 || fabs(p.z) > 1e6) continue;
        if (!seeded) { mn = mx = p; seeded = 1; continue; }
        if (p.x < mn.x) mn.x = p.x;
        if (p.x > mx.x) mx.x = p.x;
        if (p.y < mn.y) mn.y = p.y;
        if (p.y > mx.y) mx.y = p.y;
        if (p.z < mn.z) mn.z = p.z;
        if (p.z > mx.z) mx.z = p.z;
    }
    if (!seeded) return 0;
    *out = vec3_scale(vec3_add(mn, mx), 0.5);
    return 1;
}

void view_render(const Simulation *sim, const ViewCam *cam, int force_ansi) {
    view_render_to(stdout, sim, cam, force_ansi);
}

void view_render_to(FILE *fp, const Simulation *sim, const ViewCam *cam,
                     int force_ansi) {
    int cols, rows, x, y;
    ViewCells *vc;
    int use_ansi;
    if (!fp) return;
    term_size(&cols, &rows);
    vc = view_compute(sim, cam, cols, rows);
    if (!vc) return;
    use_ansi = force_ansi || isatty(STDOUT_FILENO);

    if (use_ansi) fprintf(fp, "\x1b[H\x1b[2J");
    fputc('+', fp);
    for (x = 0; x < vc->vw; x++) fputc('-', fp);
    fprintf(fp, "+ %s\n", vc->title);
    for (y = 0; y < vc->vh; y++) {
        int curfg = -1, curbr = -1;
        fputc('|', fp);
        for (x = 0; x < vc->vw; x++) {
            int i = y * vc->vw + x;
            if (use_ansi && (vc->fg[i] != curfg || vc->br[i] != curbr)) {
                if (vc->br[i]) fprintf(fp, "\x1b[1;%dm", vc->fg[i]);
                else fprintf(fp, "\x1b[0;%dm", vc->fg[i]);
                curfg = vc->fg[i]; curbr = vc->br[i];
            }
            fputc(vc->ch[i], fp);
        }
        if (use_ansi) fprintf(fp, "\x1b[0m");
        /* side panel */
        if (y == 0) fprintf(fp, "| %s", vc->panel[0]);
        else if (y == 1 && vc->panel[1][0]) fprintf(fp, "| %s", vc->panel[1]);
        else if (y == 2 && vc->panel[2][0]) fprintf(fp, "| %s", vc->panel[2]);
        else if (y == 3) fprintf(fp, "| %s", vc->panel[3]);
        else if (y == 4) fprintf(fp, "| %s", vc->panel[4]);
        else if (y == 5) fprintf(fp, "| %s", vc->panel[5]);
        fputc('\n', fp);
    }
    fputc('+', fp);
    for (x = 0; x < vc->vw; x++) fputc('-', fp);
    fprintf(fp, "+\n");
    fflush(fp);

    view_cells_free(vc);
}
