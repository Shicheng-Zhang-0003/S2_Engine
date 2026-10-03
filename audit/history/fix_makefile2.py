#!/usr/bin/env python3
"""
Fix the makefile dependency paths.
"""

with open('TREE/makefile', 'r') as f:
    content = f.read()

# Fix the dependency generation for test files
old = '''DEPFLAGS = -MMD -MP -MF $(OBJ_DIR)/$*.d

SRC_DIR = src
OBJ_DIR = build
BIN     = carbonsim

SRCS = $(wildcard $(SRC_DIR)/%.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
DEPS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.d, $(SRCS))'''

new = '''DEPFLAGS = -MMD -MP -MF $(OBJ_DIR)/$(*F).d

SRC_DIR = src
OBJ_DIR = build
BIN     = carbonsim

SRCS = $(wildcard $(SRC_DIR)/%.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))
DEPS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.d, $(SRCS))

# Test sources and objects
TEST_DIR = tests
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_OBJS = $(patsubst $(TEST_DIR)/%.c, $(OBJ_DIR)/%.o, $(TEST_SRCS))
TEST_DEPS = $(patsubst $(TEST_DIR)/%.c, $(OBJ_DIR)/%.d, $(TEST_SRCS))

DEPS += $(TEST_DEPS)'''

content = content.replace(old, new)

# Fix the test compilation rules to use correct depflags
old_test = '''$(OBJ_DIR)/test_datastream.o: $(TEST_DIR)/test_datastream.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

# forces selftest (audit P2): analytic dihedral vs FD oracle. Links the
# engine objects (minus main) since forces_dihedral lives in forces.o.
selftest-forces: all $(OBJ_DIR)/test_forces.o
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/test_forces $(OBJ_DIR)/test_forces.o $(filter-out $(OBJ_DIR)/main.o,$(OBJS)) -lm

$(OBJ_DIR)/test_forces.o: $(TEST_DIR)/test_forces.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<

# fire selftest: FIRE vs steepest-descent minima agreement.
selftest-fire: all $(OBJ_DIR)/test_fire.o
	$(CC) $(CFLAGS) -o $(OBJ_DIR)/test_fire $(OBJ_DIR)/test_fire.o $(filter-out $(OBJ_DIR)/main.o,$(OBJS)) -lm

$(OBJ_DIR)/test_fire.o: $(TEST_DIR)/test_fire.c
	$(CC) $(CFLAGS) $(DEPFLAGS) -c -o $@ $<'''

new_test = '''# Test compilation uses TEST_DEPFLAGS for correct dependency paths
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
	$(CC) $(CFLAGS) $(TEST_DEPFLAGS) -c -o $@ $<'''

content = content.replace(old, new)
content = content.replace(old_test, new_test)

with open('TREE/makefile', 'w') as f:
    f.write(content)

print("Makefile fixed!")