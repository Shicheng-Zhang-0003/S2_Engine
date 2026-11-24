# side_fixes/ — a staging area for candidate engine fixes.
#
# These are NOT part of the record build and NOT part of `make test`.
# Each fix is a self-contained .c/.h pair linked against the engine's own
# objects, so the suite exercises the real engine code rather than a copy.
#
# Two policies are inherited from the root makefile on purpose:
#
#   1. -Werror. A warning in a side fix fails the build here exactly as it
#      does in the engine. The previous Makefile silently dropped -Werror,
#      so this directory could accumulate warnings the record build would
#      reject.
#   2. The engine objects are built by the root makefile, not assumed to
#      be present. The old Makefile listed ../build/*.o as bare
#      prerequisites with no rule, so `make -C side_fixes` on a freshly
#      cloned tree failed with "No rule to make target". Now `engine-objects`
#      delegates to the root makefile, which is incremental and a no-op
#      when the objects are already current.
#
# -I../include lives in CPPFLAGS rather than CFLAGS so that a caller who
# sets CFLAGS cannot silently drop the include path (the old
# `CFLAGS ?= ... -I../include` lost -I the moment CFLAGS was set in the
# environment, producing confusing "header not found" errors).

CC ?= gcc
ENGINE_DIR = ..

CPPFLAGS += -I$(ENGINE_DIR)/include
CFLAGS    = -O3 -g -Wall -Wextra -Werror -ffp-contract=off -std=c11
LDLIBS   += -lm

ENGINE_OBJS = $(addprefix $(ENGINE_DIR)/build/,\
	sim.o forces.o periodic_table.o quantum.o qm.o qm_eht_overlap.o \
	qm_eht_scf.o integrator.o kcsa_filter.o qm_fluctuation.o \
	aminoacids.o nucleobases.o datastream.o qm_born.o loop.o neuron.o)

SIDE_OBJS = side_fix_qeq_screened.o side_fix_pauli_exact.o \
	side_fix_rhf_calibrated.o side_fix_kcsa_stabilized.o \
	side_fix_cell_list_pbc.o test_side_fixes.o

TARGET = test_side_fixes

.PHONY: all test clean engine-objects

all: $(TARGET)

# Delegate to the root makefile so there is exactly one definition of how
# an engine object is built. Phony, so it is consulted on every run; the
# root makefile's own dependency tracking makes that cheap.
engine-objects:
	$(MAKE) -C $(ENGINE_DIR)

# The trailing empty recipe is what lets these resolve: without it make
# refuses a prerequisite that no rule in THIS makefile can create. The
# file is actually produced by `engine-objects`, which runs first.
$(ENGINE_OBJS): engine-objects ;

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(TARGET): $(SIDE_OBJS) $(ENGINE_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SIDE_OBJS) $(ENGINE_OBJS) $(LDLIBS)

# Exits non-zero if any check fails, so this is usable as a gate.
test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(SIDE_OBJS) $(TARGET)
