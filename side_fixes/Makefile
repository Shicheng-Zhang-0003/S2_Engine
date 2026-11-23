CC ?= gcc
CFLAGS ?= -O3 -g -Wall -Wextra -std=c11 -I../include -ffp-contract=off
LDFLAGS ?= -lm

ENGINE_OBJS = ../build/sim.o ../build/forces.o ../build/periodic_table.o \
              ../build/quantum.o ../build/qm.o ../build/qm_eht_overlap.o \
              ../build/qm_eht_scf.o ../build/integrator.o ../build/kcsa_filter.o \
              ../build/qm_fluctuation.o ../build/aminoacids.o ../build/nucleobases.o \
              ../build/datastream.o ../build/qm_born.o ../build/loop.o ../build/neuron.o

SIDE_OBJS = side_fix_qeq_screened.o side_fix_pauli_exact.o \
            side_fix_rhf_calibrated.o side_fix_kcsa_stabilized.o \
            side_fix_cell_list_pbc.o test_side_fixes.o

TARGET = test_side_fixes

.PHONY: all test clean

all: $(TARGET)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(SIDE_OBJS)
	$(CC) $(CFLAGS) -o $@ $(SIDE_OBJS) $(ENGINE_OBJS) $(LDFLAGS)

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(SIDE_OBJS) $(TARGET)
