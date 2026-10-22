# The interactive terminal — `s2tui`

*S2 Engine v9R4 · a POSIX shell over one live simulation.*

<a id="interactive-terminal-—-s2tui"></a>
### Interactive terminal — `s2tui`
A separate binary so the record path can never leak wall time or keystrokes
into stdout. It is a **fully POSIX command surface** over one live
`Simulation` plus one HH neuron: every command name is a POSIX command and
the shell grammar is POSIX (pipes, redirection, `&&`/`||`, `$( )`, globs,
`fc -l`, `tput`, `sh -c '...'` to reach the host), while operands are S2
words — species, atom indices, scopes. Non-POSIX names are refused with
status 127.

The centrepiece is the **pure simulator mode**: `dd if=petri of=world`
builds a 48 Å periodic dish at 310 K holding ~780 atoms of interacting
matter (waters, Na/K/Cl ions, organic monomers) plus an HH neuron, and two
keybind surfaces edit that same world:

* `ps` — live control monitor (space run/pause, `s` step, `+/-` speed,
  `ijkl` rotate, arrows strafe across the environment, and
  create/augment/adapt keys: `w` water, `N` Na⁺, `K` K⁺,
  `C` Cl⁻, `u` base, `a` alanine, `g` glycine, `d` sugar, `H`/`L`
  heat/cool, `f` freeze, `r` replace, `p` clone, `x` delete, `m`
  minimise, `n` neuron, `e` catalysis). Piped, it prints one snapshot.
* `vi world` — modal editor (hjkl pan, `w`/`b` select, `i` insert,
  `r` replace, `x` delete, `p` clone, `u` undo, `/` search, `:` command
  mode running any shell command against the live world).
* `screen` — fullscreen twin of `ps` on the alternate screen (ncursesw
  build, real terminal required): same world, same keybinds, terminfo
  colors, resize handling, status footer. Both surfaces consume the
  same viewport cells, so the two pictures cannot disagree. Piped,
  dumb-`TERM`, missing-lib or failed init degrades to a typed line
  with status 0; `S2TUI_SCREEN=0` forces line mode. Esc quits like `q`.

The POSIX mapping for the rest of the world: `touch` creates matter,
`ln` links, `ln -s` restrains, `unlink` clears, `rm atom <i>` deletes,
`fsck` detects topology, `sleep <n>` advances world time, `kill -SIG`
freezes / stimulates / heats, `nice` minimises, `df`/`du` report state,
`cp <species> x.mol` exports a template and `dd if=x.mol` instantiates it,
`make x.rxn` fires a reaction rule, `sync` saves, `man` documents (30 topics incl. test/help/tput/fc/ps-keys/vi-keys; `ps`/`vi` live keybinds also in this section). Demos
remain live system tests (`test all` = 15/15). TRUST: `!`/`sh -c`/`source` run with YOUR
privileges (host exec + filesystem write) — never source untrusted files;
set `S2_NO_HOST=1` to disable the host escape. `save`/`load` (`S2SAVE1`) is LOSSY:
dihedrals, explicit angle overrides, box/PBC, step/time, QM flags and RNG state are NOT
restored (a loud WARNING prints on every load); `S2SAVE2` with full state is future work.
The orthographic ASCII grid
has depth shading, bonds, restraint anchors and a z-slab cutaway. Stream
contract: stdout is DATA, stderr is diagnostics.

---

---

*This is part of the S2 Engine documentation. The entry point is*
*[`readme.md`](../readme.md); the full document map is there.*
