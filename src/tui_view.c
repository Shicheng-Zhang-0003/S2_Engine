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

static void term_size(int *cols, int *rows) {
    struct winsize ws;
    if (cols) *cols = 80;
    if (rows) *rows = 24;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0) {
        if (ws.ws_col >= 40) *cols = ws.ws_col;
        if (ws.ws_row >= 15) *rows = ws.ws_row;
    }
}

void view_render(const Simulation *sim, const ViewCam *cam, int force_ansi) {
    view_render_to(stdout, sim, cam, force_ansi);
}

void view_render_to(FILE *fp, const Simulation *sim, const ViewCam *cam, int force_ansi) {
    if (!fp) return;
    int cols, rows;
    term_size(&cols, &rows);
    int use_ansi = force_ansi || isatty(STDOUT_FILENO);
    ViewCam c;
    if (cam) c = *cam;
    else view_cam_reset(&c);

    int panel_w = 30;
    if (cols < 60) panel_w = 20;
    int vw = cols - panel_w - 3;
    if (vw < 20) vw = 20;
    int vh = rows - 5;
    if (vh < 8) vh = 8;

    double yaw = c.yaw_deg * M_PI / 180.0;
    double pitch = c.pitch_deg * M_PI / 180.0;
    double cy = cos(yaw), sy = sin(yaw);
    double cp = cos(pitch), sp = sin(pitch);

    int n = (sim && sim->atoms) ? sim->num_atoms : 0;
    if (n < 0) n = 0;
    Proj *pr = NULL;
    int *order = NULL;
    if (n > 0) {
        pr = (Proj *)malloc(sizeof(Proj) * (size_t)n);
        order = (int *)malloc(sizeof(int) * (size_t)n);
    }
    Vec3 ctr = c.center;
    if (!c.has_center && pr) {
        Vec3 mn = sim->atoms[0].position, mx = mn;
        for (int i = 1; i < n; i++) {
            Vec3 p = sim->atoms[i].position;
            if (p.x < mn.x) mn.x = p.x;
            if (p.x > mx.x) mx.x = p.x;
            if (p.y < mn.y) mn.y = p.y;
            if (p.y > mx.y) mx.y = p.y;
            if (p.z < mn.z) mn.z = p.z;
            if (p.z > mx.z) mx.z = p.z;
        }
        ctr = vec3_scale(vec3_add(mn, mx), 0.5);
    }
    double zmin = 1e30, zmax = -1e30;
    double fitsc = 1.0;
    if (pr) {
        for (int i = 0; i < n; i++) {
            Vec3 d = vec3_sub(sim->atoms[i].position, ctr);
            double x1 = d.x * cy - d.y * sy;
            double y1 = d.x * sy + d.y * cy;
            double y2 = y1 * cp - d.z * sp;
            double z2 = y1 * sp + d.z * cp;
            pr[i].x = x1; pr[i].y = y2; pr[i].z = z2; pr[i].idx = i;
            order[i] = i;
            if (c.slice > 0.0 && fabs(z2) > c.slice) continue;
            if (z2 < zmin) zmin = z2;
            if (z2 > zmax) zmax = z2;
        }
        /* fit scale */
        double bx = 1e-9, by = 1e-9;
        for (int i = 0; i < n; i++) {
            if (c.slice > 0.0 && fabs(pr[i].z) > c.slice) continue;
            if (fabs(pr[i].x) > bx) bx = fabs(pr[i].x);
            if (fabs(pr[i].y) > by) by = fabs(pr[i].y);
        }
        double span = bx * 2.0 + 2.0, spanY = by * 2.0 + 2.0;
        double sc = vw / span;
        double scY = (vh * 2.0) / spanY; /* cells ~2x tall */
        if (scY < sc) sc = scY;
        sc *= c.zoom;
        if (!(sc > 1e-9) || !isfinite(sc)) sc = 1.0;
        fitsc = sc;
        for (int i = 0; i < n; i++) { pr[i].x *= sc; pr[i].y *= sc; }
        /* far (small z) first */
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (pr[order[j]].z < pr[order[i]].z) { int t = order[i]; order[i] = order[j]; order[j] = t; }
    }
    if (!(zmax >= zmin)) { zmin = 0; zmax = 1; }

    /* cell buffers */
    char *ch = (char *)malloc((size_t)(vw * vh));
    int *fg = (int *)malloc(sizeof(int) * (size_t)(vw * vh));
    int *br = (int *)malloc(sizeof(int) * (size_t)(vw * vh));
    if (!ch || !fg || !br) { free(pr); free(order); free(ch); free(fg); free(br); return; }
    memset(ch, ' ', (size_t)(vw * vh));
    for (int i = 0; i < vw * vh; i++) { fg[i] = 37; br[i] = 0; }

    /* bonds first (dim dots) */
    if (sim && sim->bonds && pr) {
        for (int b = 0; b < sim->num_bonds; b++) {
            int ia = sim->bonds[b].atom_a, ib = sim->bonds[b].atom_b;
            if (ia < 0 || ia >= n || ib < 0 || ib >= n) continue;
            if (c.slice > 0.0 && (fabs(pr[ia].z) > c.slice || fabs(pr[ib].z) > c.slice)) continue;
            int x0 = (int)lround(vw / 2.0 + pr[ia].x);
            int y0 = (int)lround(vh / 2.0 - pr[ia].y / 2.0);
            int x1 = (int)lround(vw / 2.0 + pr[ib].x);
            int y1 = (int)lround(vh / 2.0 - pr[ib].y / 2.0);
            int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
            int dy = -abs(y1 - y0), sy2 = y0 < y1 ? 1 : -1;
            int err = dx + dy, x = x0, y = y0, guard = vw + vh + 4;
            while (guard-- > 0) {
                if (x >= 0 && x < vw && y >= 0 && y < vh && ch[y * vw + x] == ' ') {
                    ch[y * vw + x] = '.';
                    fg[y * vw + x] = 37; br[y * vw + x] = 0;
                }
                if (x == x1 && y == y1) break;
                int e2 = 2 * err;
                if (e2 >= dy) { err += dy; x += sx; }
                if (e2 <= dx) { err += dx; y += sy2; }
            }
        }
    }
    /* atoms far->near */
    if (pr) {
        for (int k = 0; k < n; k++) {
            int i = order[k];
            if (c.slice > 0.0 && fabs(pr[i].z) > c.slice) continue;
            int x = (int)lround(vw / 2.0 + pr[i].x);
            int y = (int)lround(vh / 2.0 - pr[i].y / 2.0);
            if (x < 0 || x >= vw || y < 0 || y >= vh) continue;
            double t = (zmax > zmin) ? (pr[i].z - zmin) / (zmax - zmin) : 1.0;
            ch[y * vw + x] = glyph_for_z(sim->atoms[i].Z);
            fg[y * vw + x] = color_for_z(sim->atoms[i].Z);
            br[y * vw + x] = (t > 0.66) ? 1 : 0; /* near = bright */
        }
    }
    /* restraint anchors as x (same fit scale) */
    if (sim && pr) {
        for (int r = 0; r < sim->num_restraints; r++) {
            Vec3 d = vec3_sub(sim->restraint_anchor[r], ctr);
            double x1 = d.x * cy - d.y * sy;
            double y1 = d.x * sy + d.y * cy;
            double y2 = y1 * cp - d.z * sp;
            int x = (int)lround(vw / 2.0 + x1 * fitsc);
            int y = (int)lround(vh / 2.0 - (y2 * fitsc) / 2.0);
            if (x >= 0 && x < vw && y >= 0 && y < vh && ch[y * vw + x] == ' ') {
                ch[y * vw + x] = 'x'; fg[y * vw + x] = 33; br[y * vw + x] = 0;
            }
        }
    }

    if (use_ansi) fprintf(fp, "\x1b[H\x1b[2J");
    fputc('+', fp);
    for (int x = 0; x < vw; x++) fputc('-', fp);
    fprintf(fp, "+ %s\n", (sim ? "s2 grid" : "empty"));
    for (int y = 0; y < vh; y++) {
        fputc('|', fp);
        int curfg = -1, curbr = -1;
        for (int x = 0; x < vw; x++) {
            int i = y * vw + x;
            if (use_ansi && (fg[i] != curfg || br[i] != curbr)) {
                if (br[i]) fprintf(fp, "\x1b[1;%dm", fg[i]);
                else fprintf(fp, "\x1b[0;%dm", fg[i]);
                curfg = fg[i]; curbr = br[i];
            }
            fputc(ch[i], fp);
        }
        if (use_ansi) fprintf(fp, "\x1b[0m");
        /* side panel */
        if (y == 0) fprintf(fp, "| N=%d", sim ? sim->num_atoms : 0);
        else if (y == 1 && sim) fprintf(fp, "| step=%llu", (unsigned long long)sim->step);
        else if (y == 2 && sim) fprintf(fp, "| T=%.1fK E=%.3f", sim->temperature, sim->potential_energy);
        else if (y == 3) fprintf(fp, "| yaw=%.0f pit=%.0f z=%.2f", c.yaw_deg, c.pitch_deg, c.zoom);
        else if (y == 4) fprintf(fp, "| h c n o N p K =H C N O Na P K");
        else if (y == 5) fprintf(fp, "| view %s%s", c.has_center ? "centered" : "fit", c.slice > 0 ? " slice" : "");
        fputc('\n', fp);
    }
    fputc('+', fp);
    for (int x = 0; x < vw; x++) fputc('-', fp);
    fprintf(fp, "+\n");
    fflush(fp);

    free(pr); free(order); free(ch); free(fg); free(br);
}
