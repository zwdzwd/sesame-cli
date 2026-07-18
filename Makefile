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
# EXTRA_CFLAGS is appended to every compile (incl. mask.o) and EXTRA_LDFLAGS to
# every link, so `asan` can add sanitizers without clobbering the base flags a
# command-line CFLAGS= override would drop.
CFLAGS  += $(EXTRA_CFLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)
LDLIBS  += -lz -lm -lcurl

# --- YAME: linked directly for reading .cm masks (both AGPL). Built from the
#     pinned submodule so the static lib never goes stale. ---
YAME_DIR := YAME
YAME_LIB := $(YAME_DIR)/libyame.a
HTSLIB   := $(YAME_DIR)/htslib/libhts.a
YAME_INC := -I$(YAME_DIR)/src -I$(YAME_DIR)/htslib

SRC     := src/util.c src/sha256.c src/numerics.c src/idat.c src/index.c src/sigdf.c src/prep.c src/qc.c src/dml.c src/cnv.c src/cbs.c src/mask.c src/cgwrite.c src/attach.c src/cache.c
CLI_SRC := cli/main.c
OBJ     := $(SRC:.c=.o)
CLI_OBJ := $(CLI_SRC:.c=.o)
BIN     := sesame

.PHONY: all asan test test-idat test-betas test-prep test-qmask test-poobah test-noob test-batch test-qc test-dml test-cg test-attach test-cnv test-cbs index cnv-normals yame-lib fuzz fuzz-replay clean

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
	$(CC) -O2 -g -std=gnu11 -Wall -Iinclude $(YAME_INC) $(EXTRA_CFLAGS) -c -o $@ $<

# cgwrite.c also includes YAME headers; same relaxed rule as mask.o.
src/cgwrite.o: src/cgwrite.c include/sesame.h src/internal.h | $(YAME_LIB)
	$(CC) -O2 -g -std=gnu11 -Wall -Iinclude $(YAME_INC) $(EXTRA_CFLAGS) -c -o $@ $<

# attach.c also includes YAME headers; same relaxed rule as cgwrite.o.
src/attach.o: src/attach.c include/sesame.h src/internal.h | $(YAME_LIB)
	$(CC) -O2 -g -std=gnu11 -Wall -Iinclude $(YAME_INC) $(EXTRA_CFLAGS) -c -o $@ $<

%.o: %.c include/sesame.h src/internal.h src/registry.h
	$(CC) $(CFLAGS) -c -o $@ $<

asan: clean
	$(MAKE) EXTRA_CFLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
	        EXTRA_LDFLAGS="-fsanitize=address,undefined"

test: test-idat test-betas test-prep test-qmask test-poobah test-noob test-batch test-qc test-dml test-cg test-attach test-cnv test-cbs

test-idat: $(BIN)
	@tests/run_golden.sh

test-betas: $(BIN) pipeline_dump
	@tests/run_betas.sh

test-prep: $(BIN) pipeline_dump
	@tests/run_prep.sh

test-qmask: $(BIN) pipeline_dump
	@tests/run_qmask.sh

test-poobah: $(BIN) pipeline_dump
	@tests/run_poobah.sh

# The B unit harness links only numerics.o (no other deps).
normexp_test: tests/normexp_test.c src/numerics.o include/sesame.h src/internal.h
	$(CC) $(CFLAGS) -Isrc -o $@ tests/normexp_test.c src/numerics.o -lm

# CBS unit harness: reads a bin-signal vector, prints segments. numerics.o + cbs.o.
cbs_test: tests/cbs_test.c src/cbs.o src/numerics.o include/sesame.h src/internal.h
	$(CC) $(CFLAGS) -Isrc -o $@ tests/cbs_test.c src/cbs.o src/numerics.o -lm

# Validation harness: full-precision library pipeline (links libsesame + YAME),
# so differential tests keep their bit-identical/ULP gates that the .cg product
# (float32) cannot carry.
pipeline_dump: tests/pipeline_dump.c $(OBJ) $(YAME_LIB) $(HTSLIB) include/sesame.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ tests/pipeline_dump.c $(OBJ) \
	    $(YAME_LIB) $(HTSLIB) -lpthread $(LDLIBS)

# M/U TSV -> format-3 .cg, for packing the CNV normal reference (export_cnvnormals.R).
mu2cg: tools/mu2cg.c $(OBJ) $(YAME_LIB) $(HTSLIB) include/sesame.h
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ tools/mu2cg.c $(OBJ) \
	    $(YAME_LIB) $(HTSLIB) -lpthread $(LDLIBS)

test-noob: $(BIN) normexp_test pipeline_dump
	@tests/run_noob.sh

test-batch: $(BIN)
	@tests/run_batch.sh

test-qc: $(BIN)
	@tests/run_qc.sh

test-dml: $(BIN)
	@tests/run_dml.sh

test-cg: $(BIN) pipeline_dump
	@tests/run_cg.sh

test-attach: $(BIN) mu2cg
	@tests/run_attach.sh

test-cnv: $(BIN)
	@tests/run_cnv.sh

test-cbs: cbs_test
	@tests/run_cbs.sh

# Export ordering tables from sesameData (bootstrap; needs Rscript + sesame).
PLATFORMS := HM450 EPIC EPICv2 MSA
index:
	@for p in $(PLATFORMS); do \
	    Rscript tools/export_ordering.R $$p testdata/$$p.ordering.tsv.gz; \
	done

# Export the per-platform CNV normal references (format-3 .cg) from sesameData
# into the store. Genome-level genomeInfo is NOT here -- it lives in the separate
# zhou-lab/genomes repo (generate it with tools/export_genomeinfo.R) and is
# pulled by `sesame fetch genome <build>`. Needs Rscript.
cnv-normals: mu2cg
	@for p in EPIC EPICv2; do \
	    mkdir -p data/$$p; \
	    Rscript tools/export_cnvnormals.R $$p testdata/$$p.ordering.tsv.gz \
	        data/$$p/$$p.cnvnormals.mu.tsv.gz && \
	    ./mu2cg data/$$p/$$p.cnvnormals.mu.tsv.gz data/$$p/$$p.cnvnormals.cg && \
	    rm -f data/$$p/$$p.cnvnormals.mu.tsv.gz; \
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
	rm -f $(OBJ) $(CLI_OBJ) $(BIN) normexp_test cbs_test pipeline_dump mu2cg fuzz_idat fuzz_replay
	rm -rf $(BIN).dSYM normexp_test.dSYM fuzz_replay.dSYM fuzz_idat.dSYM
