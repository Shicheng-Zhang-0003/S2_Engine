#ifndef TUI_H
#define TUI_H

/*
 * tui.h — OpenWorm-in-the-TUI: text-mode live simulator control.
 *
 * No graphics. A line-oriented REPL over ONE live Simulation plus one
 * HH neuron. Batch mode (`./carbonsim`, byte-deterministic record) is
 * untouched: this lives in a separate `s2tui` binary so the record path
 * can never leak wall-time or keystroke timing into stdout.
 *
 * Each `test demo <id>` mirrors a batch demo's physics at reduced cost
 * (WHAM sampling shortened, knock landscape coarsened) and prints a
 * PASS/FAIL verdict, so every demo doubles as a live system test.
 */

int tui_main(void);

#endif /* TUI_H */
