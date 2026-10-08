CC      = gcc
# FULL-AUDIT C24/C25: -ffp-contract=off pins FP codegen (the -march ban alone
# does not stop FMA fusion under -O3 on FMA hosts; -O3-vs-O0 parity was luck),
# and -Werror makes warnings fail the build instead of advisory text.
CFLAGS_BASE = -O3 -g -Wall -Wextra -Werror -ffp-contract=off -std=c11 -Iinclude
# [s41] datastream build identity: ds_open() records build-compiler and
# build-flags in every datastream [header] (DATASTREAM_SPEC.md). The
# strings are injected here so the binary knows how it was built.
# Builds that override CFLAGS on the command line (ASan, debug) fall
# back to "unknown" via the #ifndef defaults in datastream.c.
# Full compiler version (not major-only) for record precision.
CC_VERSION = $(shell $(CC) -dumpfullversion 2>/dev/null || $(CC) -dumpversion)
CFLAGS = $(CFLAGS_BASE) -DDS_COMPILER_ID='"$(CC) $(CC_VERSION)"' -DDS_CFLAGS_ID='"$(CFLAGS_BASE)"'
# PERIODIC_TABLE ground_config initializer flood, which s03-s05
# eliminated at the source (fully explicit {{{0}}, 0, 0}). Keeping
# it would only hide a future regression of the same class.
# -Wno-missing-braces removed by s06: dead suppression. The initializer flood it masked was eliminated at the source by s03-s05 (fully explicit initializers); keeping it would only hide a future regression.
# [s38] -march=native disabled for record portability. The record
# build must reproduce across machines; -march=native enables
# machine-specific instructions (e.g. FMA) that change FP
# contraction and thus low-order result digits. Re-enable only
# for performance runs, never for the record.
# CFLAGS += -march=native

# Uncomment for debug build with address sanitiser:
# CFLAGS = -g -O0 -Wall -std=c11 -Iinclude -fsanitize=address,undefined

# Dependency generation flags
DEPFLAGS = -MMD -MP -MF $(OBJ_DIR)/$*.d

SRC_DIR = src
OBJ_DIR = build
BIN     = carbonsim
TUI_BIN = s2tui

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
# tui.c carries its own main: keep it out of the batch record binary so
# ./carbonsim stays byte-deterministic. s2tui links engine objects minus main.
# tui_screen.o (the optional fullscreen) likewise never enters carbonsim
# or the test binaries: the record path stays libm-only by construction.
BIN_OBJS = $(filter-out $(OBJ_DIR)/tui.o $(OBJ_DIR)/tui_view.o $(OBJ_DIR)/tui_screen.o,$(OBJS))
ENG_OBJS = $(filter-out $(OBJ_DIR)/main.o $(OBJ_DIR)/tui.o $(OBJ_DIR)/tui_view.o $(OBJ_DIR)/tui_screen.o,$(OBJS))
TUI_OBJS = $(ENG_OBJS) $(OBJ_DIR)/tui.o $(OBJ_DIR)/tui_view.o $(OBJ_DIR)/tui_screen.o
DEPS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.d, $(SRCS))

.PHONY: all clean run selftest selftest-forces selftest-fire selftest-regression \
        selftest-external selftest-loop test deps audit audit-revert tui

# One command builds every executable: the batch record binary, the live
# TUI, and all six test binaries — `make` / `make -j$(nproc)` is enough,
# no separate `make tui` / `make selftest-*` runs needed. The selftest-*
# targets below remain as aliases (verify_scripts.sh calls them).
TEST_BINS = $(OBJ_DIR)/test_datastream $(OBJ_DIR)/test_forces \
            $(OBJ_DIR)/test_fire $(OBJ_DIR)/test_regression \
            $(OBJ_DIR)/test_external $(OBJ_DIR)/test_loop
# FULL-AUDIT O16: seal checker (calls ds_verify_file; not a gate itself).
VERIFY_BINS = $(OBJ_DIR)/verify_cvmds

all: $(OBJ_DIR) $(BIN) $(TUI_BIN) $(TEST_BINS) $(VERIFY_BINS)

$(OBJ_DIR)/verify_cvmds: tests/verify_cvmds.c $(OBJ_DIR)/datastream.o
	$(CC) $(CFLAGS) -o $@ tests/verify_cvmds.c $(OBJ_DIR)/datastream.o

# Live terminal (see readme `s2tui` section): separate binary with its own
# main, linked against engine objects minus main. Untracked tool, not record.
# Fullscreen `screen` is ncursesw-optional: detected via ncursesw6-config,
# linked into s2tui only. Without it the TUI still builds and every line
# command works; only the alternate-screen surface reports unavailable.
# carbonsim and all test binaries stay libm-only either way.
CURSES_LIBS := $(shell ncursesw6-config --libs 2>/dev/null)
CURSES_CFLAGS := $(shell ncursesw6-config --cflags 2>/dev/null)
tui: $(OBJ_DIR) $(TUI_BIN)

$(OBJ_DIR)/tui_screen.o: $(SRC_DIR)/tui_screen.c
ifneq ($(CURSES_LIBS),)
	$(CC) $(CFLAGS) $(CURSES_CFLAGS) -DS2_HAVE_CURSES=1 -MMD -MP -MF $(OBJ_DIR)/tui_screen.d -c -o $@ $<
else
	$(CC) $(CFLAGS) -MMD -MP -MF $(OBJ_DIR)/tui_screen.d -c -o $@ $<
endif

$(TUI_BIN): $(TUI_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm $(CURSES_LIBS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN): $(BIN_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

# Include auto-generated dependencies
-include $(DEPS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

# FULL-AUDIT O8: `make run` used to bypass ./run discipline (no staging, no
# digest, no gates). Now it execs ./run so there is exactly one record path.
run: all
	./run

# [s41] datastream selftest: tests/test_datastream.c links against
# src/datastream.o only — no dependency on the rest of the engine.
TEST_DIR = tests

$(OBJ_DIR)/test_datastream: $(OBJ_DIR) $(OBJ_DIR)/test_datastream.o $(OBJ_DIR)/datastream.o
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_datastream.o $(OBJ_DIR)/datastream.o -lm

selftest: $(OBJ_DIR)/test_datastream

# Test compilation uses TEST_DEPFLAGS for correct dependency paths
TEST_DEPFLAGS = -MMD -MP -MF $(OBJ_DIR)/$(*F).d

$(OBJ_DIR)/test_datastream.o: $(TEST_DIR)/test_datastream.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# forces selftest (audit P2): analytic dihedral vs FD oracle. Links the
# engine objects (minus main) since forces_dihedral lives in forces.o.
$(OBJ_DIR)/test_forces: $(OBJ_DIR) $(OBJ_DIR)/test_forces.o $(ENG_OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_forces.o $(ENG_OBJS) -lm

selftest-forces: $(OBJ_DIR)/test_forces

$(OBJ_DIR)/test_forces.o: $(TEST_DIR)/test_forces.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# fire selftest: FIRE vs steepest-descent minima agreement.
$(OBJ_DIR)/test_fire: $(OBJ_DIR) $(OBJ_DIR)/test_fire.o $(ENG_OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_fire.o $(ENG_OBJS) -lm

selftest-fire: $(OBJ_DIR)/test_fire

$(OBJ_DIR)/test_fire.o: $(TEST_DIR)/test_fire.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# audit regression suite: one check per defect found and fixed in the
# v9R4 audit. Independent oracles (finite differences, quadrature, NIST
# SHA-256 vectors, published reference values) so it can actually fail.
selftest-regression: $(OBJ_DIR)/test_regression

$(OBJ_DIR)/test_regression: $(OBJ_DIR) $(OBJ_DIR)/test_regression.o $(ENG_OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_regression.o $(ENG_OBJS) -lm

$(OBJ_DIR)/test_regression.o: $(TEST_DIR)/test_regression.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# external validation suite (audit fix E1): checks the engine against values
# obtained from OUTSIDE this repository - NIST CODATA, the original AMBER
# ff99 parm99.dat, FIPS 180-4 SHA-256 vectors, Griffiths' closed forms, and
# the deposited PDB 1K4C entry. test_regression.c re-derives its own
# oracles, which is self-referential: an AMBER R* compared as a sigma and a
# FIPS digest transcribed with one wrong character were both harness errors
# during this audit, and in both cases the engine was right.
selftest-external: $(OBJ_DIR)/test_external

$(OBJ_DIR)/test_external: $(OBJ_DIR) $(OBJ_DIR)/test_external.o $(ENG_OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_external.o $(ENG_OBJS) -lm

$(OBJ_DIR)/test_external.o: $(TEST_DIR)/test_external.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# loop-closure suite: bio/QC/QM consistency from SAME primaries.
# Nernst reversals, 1000:1 selectivity scale, Marcus/Shannon/LJ/Coulomb
# identities, fail-closed NaN, explicit abstraction inventory. Not part of
# the batch record; proves the three tracks share one physics.
selftest-loop: $(OBJ_DIR)/test_loop

$(OBJ_DIR)/test_loop: $(OBJ_DIR) $(OBJ_DIR)/test_loop.o $(OBJ_DIR)/tui_view.o $(ENG_OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJ_DIR)/test_loop.o $(OBJ_DIR)/tui_view.o $(ENG_OBJS) -lm

$(OBJ_DIR)/test_loop.o: $(TEST_DIR)/test_loop.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# Run every gate. `make test` is what CI should invoke.
# AUDIT FIX M5: the regression line used to be
#     @./build/test_regression | tail -3
# A pipeline reports the exit status of its LAST command, so make saw tail's 0
# and reported success on a red suite. Confirmed against the pre-fix solver:
# test_regression reported 4 FAIL lines and that recipe still exited 0.
# The log is now written to a file so the binary's own exit status reaches
# make unmasked, and the FAIL lines are grepped as well so the gate fails even
# if a future change to the suite decouples its exit code from its verdicts.
test: selftest selftest-forces selftest-fire selftest-regression selftest-external selftest-loop
	@./$(OBJ_DIR)/test_datastream   > $(OBJ_DIR)/datastream.log 2>&1; \
	  dc=$$?; if [ $$dc -ne 0 ]; then echo "  datastream   FAILED (rc=$$dc)"; cat $(OBJ_DIR)/datastream.log; exit 1; else echo "  datastream   OK"; fi
	@./$(OBJ_DIR)/test_forces      > $(OBJ_DIR)/forces.log 2>&1; \
	  fc=$$?; if [ $$fc -ne 0 ]; then echo "  forces       FAILED (rc=$$fc)"; cat $(OBJ_DIR)/forces.log; exit 1; else echo "  forces       OK"; fi
	@./$(OBJ_DIR)/test_fire        > $(OBJ_DIR)/fire.log 2>&1; \
	  fic=$$?; if [ $$fic -ne 0 ]; then echo "  fire         FAILED (rc=$$fic)"; cat $(OBJ_DIR)/fire.log; exit 1; else echo "  fire         OK"; fi
	@./$(OBJ_DIR)/test_loop        > $(OBJ_DIR)/loop.log 2>&1; \
	  lc=$$?; tail -3 $(OBJ_DIR)/loop.log; \
	  if [ $$lc -ne 0 ] || grep -q '^  FAIL' $(OBJ_DIR)/loop.log; then \
	    echo "  loop         FAILED (rc=$$lc)"; grep '^  FAIL' $(OBJ_DIR)/loop.log | head -20; \
	    exit 1; \
	  else \
	    echo "  loop         OK ($$(grep -c '^  PASS' $(OBJ_DIR)/loop.log) consistency checks)"; \
	  fi
	@./$(OBJ_DIR)/test_external    > $(OBJ_DIR)/external.log 2>&1; \
	  ex=$$?; tail -3 $(OBJ_DIR)/external.log; \
	  if [ $$ex -ne 0 ] || grep -q '^  FAIL' $(OBJ_DIR)/external.log; then \
	    echo "  external     FAILED (rc=$$ex)"; grep '^  FAIL' $(OBJ_DIR)/external.log | head -20; \
	    exit 1; \
	  else \
	    echo "  external     OK ($$(grep -c '^  PASS' $(OBJ_DIR)/external.log) vs outside references)"; \
	  fi
	@./$(OBJ_DIR)/test_regression  > $(OBJ_DIR)/regression.log; \
	  rc=$$?; tail -3 $(OBJ_DIR)/regression.log; \
	  if [ $$rc -ne 0 ] || grep -q '^  FAIL' $(OBJ_DIR)/regression.log; then \
	    echo "  regression  FAILED (rc=$$rc)"; grep '^  FAIL' $(OBJ_DIR)/regression.log | head -20; \
	    exit 1; \
	  else echo "  regression  OK"; fi

# The v9R4 audit's own working material: independent mathematical oracles,
# C harnesses linked against the engine objects, external-source checks, and
# -- with audit-revert -- the same regression suite linked against the
# PRE-FIX engine to prove it still catches the defects it claims to.
# Not part of `make test`: this is evidence, not a product gate.
audit:
	@./audit/run_audit.sh

audit-revert:
	@./audit/run_audit.sh --revert

clean:
	rm -rf $(OBJ_DIR) $(BIN) $(TUI_BIN) audit/build

# Print dependency info for debugging
deps:
	@echo "Sources: $(SRCS)"
	@echo "Objects: $(OBJS)"
	@echo "Deps: $(DEPS)"
