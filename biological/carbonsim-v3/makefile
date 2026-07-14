CC      = gcc
CFLAGS  = -O3 -march=native -Wall -Wextra -std=c11 -Iinclude
CFLAGS += -Wno-missing-braces -march=native -ffast-math

# Uncomment for debug build with address sanitiser:
# CFLAGS = -g -O0 -Wall -std=c11 -Iinclude -fsanitize=address,undefined -lm

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
	$(CC) -O3 -o $@ $^ -lm

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

run: all
	./$(BIN)

clean:
	rm -rf $(OBJ_DIR) $(BIN)
