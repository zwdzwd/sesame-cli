# sesame -- POSIX make. No external build system; zlib is the only dependency.
#
#   make            build sesame
#   make asan       build with ASan/UBSan
#   make test       run the golden tests (requires Rscript for the R oracle)
#   make clean

CC      ?= cc
CFLAGS  ?= -O2 -g
CFLAGS  += -DSESAME_HAVE_CURL -std=c11 -Wall -Wextra -Wpedantic -Wshadow -Wconversion \
           -Wstrict-prototypes -Iinclude
LDLIBS  += -lz -lm -lcurl

# --- YAME: linked directly for reading .cm masks (both AGPL). Built from the
#     pinned submodule so the static lib never goes stale. ---
YAME_DIR := YAME
YAME_LIB := $(YAME_DIR)/libyame.a
HTSLIB   := $(YAME_DIR)/htslib/libhts.a
YAME_INC := -I$(YAME_DIR)/src -I$(YAME_DIR)/htslib

SRC     := src/util.c src/sha256.c src/numerics.c src/idat.c src/index.c src/sigdf.c src/prep.c src/mask.c src/cache.c
CLI_SRC := cli/main.c
OBJ     := $(SRC:.c=.o)
CLI_OBJ := $(CLI_SRC:.c=.o)
BIN     := sesame

.PHONY: all asan test test-idat test-betas test-prep test-qmask index yame-lib fuzz fuzz-replay clean

all: $(BIN)

# Incrementally (re)build libyame.a from the checked-out submodule.
yame-lib:
	$(MAKE) -C $(YAME_DIR) lib
$(YAME_LIB) $(HTSLIB): yame-lib

$(BIN): $(OBJ) $(CLI_OBJ) $(YAME_LIB) $(HTSLIB)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(OBJ) $(CLI_OBJ) $(YAME_LIB) $(HTSLIB) -lpthread $(LDLIBS)

# mask.c includes YAME headers (C99/GNU style); relax pedantic/conversion for it
# only, and give it the YAME include path.
src/mask.o: src/mask.c include/sesame.h src/internal.h | $(YAME_LIB)
	$(CC) -O2 -g -std=gnu11 -Wall -Iinclude $(YAME_INC) -c -o $@ $<

%.o: %.c include/sesame.h src/internal.h src/registry.h
	$(CC) $(CFLAGS) -c -o $@ $<

asan: clean
	$(MAKE) CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LDFLAGS="-fsanitize=address,undefined"

test: test-idat test-betas test-prep test-qmask

test-idat: $(BIN)
	@tests/run_golden.sh

test-betas: $(BIN)
	@tests/run_betas.sh

test-prep: $(BIN)
	@tests/run_prep.sh

test-qmask: $(BIN)
	@tests/run_qmask.sh

# Export ordering tables from sesameData (bootstrap; needs Rscript + sesame).
PLATFORMS := HM450 EPIC EPICv2 MSA
index:
	@for p in $(PLATFORMS); do \
	    Rscript tools/export_ordering.R $$p testdata/$$p.ordering.tsv.gz; \
	done

# Real fuzzing. Needs a clang with libFuzzer (Linux CI); Apple clang has none.
fuzz:
	$(CC) -std=c11 -g -O1 -Iinclude -fsanitize=fuzzer,address,undefined \
	    fuzz/fuzz_idat.c src/idat.c src/util.c -lz -o fuzz_idat

# Corpus replayer under ASan/UBSan. Works everywhere, incl. Apple clang.
fuzz-replay:
	$(CC) -std=c11 -g -O1 -Iinclude -DSESAME_FUZZ_STANDALONE \
	    -fsanitize=address,undefined \
	    fuzz/fuzz_idat.c src/idat.c src/util.c -lz -o fuzz_replay

clean:
	rm -f $(OBJ) $(CLI_OBJ) $(BIN) fuzz_idat fuzz_replay
	rm -rf $(BIN).dSYM fuzz_replay.dSYM fuzz_idat.dSYM
