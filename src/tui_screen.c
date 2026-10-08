#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "../include/tui_screen.h"

#ifdef S2_HAVE_CURSES
#include <locale.h>
#include <curses.h>

/* ANSI 30-37 base colors map onto the eight curses colors 1:1. */
static int pair_for(int fg) {
    int c = fg - 30;
    if (c < 0 || c > 7) c = 7;
    return c + 1;
}

static int screen_inited = 0;

int screen_available(void) {
    const char *term, *env;
    if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO)) return 0;
    env = getenv("S2TUI_SCREEN");
    if (env && (!strcmp(env, "0") || !strcmp(env, "off"))) return 0;
    term = getenv("TERM");
    if (!term || !term[0]) return 0;
    if (!strcmp(term, "dumb") || !strcmp(term, "unknown")) return 0;
    return 1;
}

int screen_init(void) {
    if (screen_inited) return 0;
    if (!screen_available()) return 1;
    setlocale(LC_ALL, "");
    if (initscr() == NULL) return 1;
    /* FULL-AUDIT C33: curses setup returns are cosmetic-only (screen never
     * feeds the record); explicitly ignored so -Werror stays quiet. */
    (void)cbreak();
    (void)noecho();
    (void)keypad(stdscr, TRUE);
    (void)curs_set(0);
    set_escdelay(25);
    if (has_colors()) {
        start_color();
        use_default_colors();
        for (int i = 0; i < 8; i++)
            init_pair((short)(i + 1), (short)i, (short)-1);
    }
    screen_inited = 1;
    return 0;
}

void screen_shutdown(void) {
    if (!screen_inited) return;
    (void)curs_set(1);
    (void)endwin();
    screen_inited = 0;
}

void screen_size(int *cols, int *rows) {
    int r = 24, c = 80;
    if (screen_inited) getmaxyx(stdscr, r, c);
    if (cols) *cols = c;
    if (rows) *rows = r;
}

/* One status/help line pair along the bottom; kept in one place so the
 * monitor and any future screen surface print the same keybinds. */
void screen_help(char *b1, size_t n1, char *b2, size_t n2) {
    snprintf(b1, n1, "%s",
        "space run/pause  s step  +/- speed  ijkl rotate  arrows strafe  z/Z zoom  Tab sel  c center");
    snprintf(b2, n2, "%s",
        "w water N Na+ K K+ C Cl- u base a ala g gly d sugar x del f freeze p clone r repl H/L heat m min n neuron e cat  q quit");
}

void screen_draw(const ViewCells *vc, const char *status,
                 const char *help1, const char *help2, int sel) {
    int rows, cols, y, x;
    (void)sel;
    if (!screen_inited || !vc) return;
    getmaxyx(stdscr, rows, cols);
    erase();
    /* viewport grid, top-left */
    for (y = 0; y < vc->vh && y < rows - 3; y++) {
        for (x = 0; x < vc->vw && x < cols - 32; x++) {
            int i = y * vc->vw + x;
            int pair = pair_for(vc->fg[i]);
            attr_t at = COLOR_PAIR(pair);
            if (vc->br[i]) at |= A_BOLD;
            attron(at);
            mvaddch(y, x, (chtype)vc->ch[i]);
            attroff(at);
        }
    }
    /* side panel to the right of the grid */
    {
        int px = vc->vw + 2;
        if (px < cols - 2) {
            mvprintw(0, px, "%s", vc->title);
            for (int p = 0; p < 6 && 1 + p < rows - 2; p++)
                mvprintw(1 + p, px, "%s", vc->panel[p]);
        }
    }
    /* status + help along the bottom */
    if (rows - 3 >= 0) {
        attron(A_REVERSE);
        mvprintw(rows - 3, 0, "%-*.*s", cols, cols, status ? status : "");
        attroff(A_REVERSE);
    }
    if (rows - 2 >= 0) mvprintw(rows - 2, 0, "%-*.*s", cols, cols, help1 ? help1 : "");
    if (rows - 1 >= 0) mvprintw(rows - 1, 0, "%-*.*s", cols, cols, help2 ? help2 : "");
    refresh();
}

int screen_getch(int timeout_ms) {
    int k;
    if (!screen_inited) return -2;
    timeout(timeout_ms);
    k = getch();
    if (k == ERR) return -1;
    switch (k) {
        case KEY_UP: return 1001;    /* LIVE_KEY_UP */
        case KEY_DOWN: return 1002;  /* LIVE_KEY_DOWN */
        case KEY_RIGHT: return 1003; /* LIVE_KEY_RIGHT */
        case KEY_LEFT: return 1004;  /* LIVE_KEY_LEFT */
        case KEY_RESIZE: return -3;  /* just redraw */
        default: return k;
    }
}

/* Prompt for one line: leave curses, ask on the real terminal, resume. */
int screen_prompt(const char *prompt, char *buf, size_t cap) {
    int r = 1;
    if (!screen_inited || !buf || cap == 0) return 1;
    def_prog_mode();
    endwin();
    fprintf(stderr, "%s", prompt ? prompt : "");
    fflush(stderr);
    if (fgets(buf, (int)cap, stdin) != NULL) {
        buf[strcspn(buf, "\n")] = '\0';
        r = 0;
    }
    reset_prog_mode();
    refresh();
    return r;
}

#else /* !S2_HAVE_CURSES */

int screen_available(void) { return 0; }
int screen_init(void) { return 1; }
void screen_shutdown(void) {}
void screen_size(int *cols, int *rows) {
    if (cols) *cols = 80;
    if (rows) *rows = 24;
}
void screen_help(char *b1, size_t n1, char *b2, size_t n2) {
    if (b1 && n1) b1[0] = '\0';
    if (b2 && n2) b2[0] = '\0';
}
void screen_draw(const ViewCells *vc, const char *status,
                 const char *help1, const char *help2, int sel) {
    (void)vc; (void)status; (void)help1; (void)help2; (void)sel;
}
int screen_getch(int timeout_ms) { (void)timeout_ms; return -2; }
int screen_prompt(const char *prompt, char *buf, size_t cap) {
    (void)prompt; (void)buf; (void)cap; return 1;
}

#endif /* S2_HAVE_CURSES */
