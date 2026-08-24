CC      = gcc
CFLAGS  = -O3 -g -Wall -Wextra -std=c11 -Iinclude
# PERIODIC_TABLE ground_config initializer flood, which s03-s05
# eliminated at the source (fully explicit {{{0}}, 0, 0}). Keeping
# it would only hide a future regression of the same class.
# -Wno-missing-braces removed by s06: dead suppression. The initializer flood it masked was eliminated at the source by s03-s05 (fully explicit initializers); keeping it would only hide a future regression.
CFLAGS += -march=native

# Uncomment for debug build with address sanitiser:
# CFLAGS = -g -O0 -Wall -std=c11 -Iinclude -fsanitize=address,undefined

SRC_DIR = src
OBJ_DIR = build
BIN     = carbonsim

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

.PHONY: all clean run

all: $(OBJ_DIR) $(BIN)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: all
	./$(BIN)

clean:
	rm -rf $(OBJ_DIR) $(BIN)
