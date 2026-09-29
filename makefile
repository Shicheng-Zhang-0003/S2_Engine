CC      = gcc
CFLAGS_BASE = -O3 -g -Wall -Wextra -std=c11 -Iinclude
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

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
DEPS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.d, $(SRCS))

.PHONY: all clean run selftest selftest-forces selftest-fire deps

all: $(OBJ_DIR) $(BIN)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

# Include auto-generated dependencies
-include $(DEPS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

run: all
	./$(BIN)

# [s41] datastream selftest: tests/test_datastream.c links against
# src/datastream.o only — no dependency on the rest of the engine.
TEST_DIR = tests

selftest: $(OBJ_DIR) $(OBJ_DIR)/test_datastream.o $(OBJ_DIR)/datastream.o
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/test_datastream $(OBJ_DIR)/test_datastream.o $(OBJ_DIR)/datastream.o -lm

# Test compilation uses TEST_DEPFLAGS for correct dependency paths
TEST_DEPFLAGS = -MMD -MP -MF $(OBJ_DIR)/$(*F).d

$(OBJ_DIR)/test_datastream.o: $(TEST_DIR)/test_datastream.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# forces selftest (audit P2): analytic dihedral vs FD oracle. Links the
# engine objects (minus main) since forces_dihedral lives in forces.o.
selftest-forces: all $(OBJ_DIR)/test_forces.o
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/test_forces $(OBJ_DIR)/test_forces.o $(filter-out $(OBJ_DIR)/main.o,$(OBJS)) -lm

$(OBJ_DIR)/test_forces.o: $(TEST_DIR)/test_forces.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

# fire selftest: FIRE vs steepest-descent minima agreement.
selftest-fire: all $(OBJ_DIR)/test_fire.o
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/test_fire $(OBJ_DIR)/test_fire.o $(filter-out $(OBJ_DIR)/main.o,$(OBJS)) -lm

$(OBJ_DIR)/test_fire.o: $(TEST_DIR)/test_fire.c
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(BIN)

# Print dependency info for debugging
deps:
	@echo "Sources: $(SRCS)"
	@echo "Objects: $(OBJS)"
	@echo "Deps: $(DEPS)"
