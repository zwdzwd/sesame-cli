# sesamec -- POSIX make. No external build system; zlib is the only dependency.
#
#   make            build sesamec
#   make asan       build with ASan/UBSan
#   make test       run the golden tests (requires Rscript for the R oracle)
#   make clean

CC      ?= cc
CFLAGS  ?= -O2 -g
CFLAGS  += -std=c11 -Wall -Wextra -Wpedantic -Wshadow -Wconversion \
           -Wstrict-prototypes -Iinclude
LDLIBS  += -lz

SRC     := src/idat.c
CLI_SRC := cli/main.c
OBJ     := $(SRC:.c=.o)
CLI_OBJ := $(CLI_SRC:.c=.o)
BIN     := sesamec

.PHONY: all asan test fuzz fuzz-replay clean

all: $(BIN)

$(BIN): $(OBJ) $(CLI_OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.c include/sesame.h
	$(CC) $(CFLAGS) -c -o $@ $<

asan: clean
	$(MAKE) CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LDFLAGS="-fsanitize=address,undefined"

test: $(BIN)
	@tests/run_golden.sh

# Real fuzzing. Needs a clang with libFuzzer (Linux CI); Apple clang has none.
fuzz:
	$(CC) -std=c11 -g -O1 -Iinclude -fsanitize=fuzzer,address,undefined \
	    fuzz/fuzz_idat.c src/idat.c -lz -o fuzz_idat

# Corpus replayer under ASan/UBSan. Works everywhere, incl. Apple clang.
fuzz-replay:
	$(CC) -std=c11 -g -O1 -Iinclude -DSESAME_FUZZ_STANDALONE \
	    -fsanitize=address,undefined \
	    fuzz/fuzz_idat.c src/idat.c -lz -o fuzz_replay

clean:
	rm -f $(OBJ) $(CLI_OBJ) $(BIN)
	rm -rf $(BIN).dSYM
