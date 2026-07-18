#!/bin/bash
set -euo pipefail

# The vendored YAME/htslib static libs (and their .o) are committed for the
# author's platform; drop them so they rebuild for THIS target arch. sesame's
# own `make clean` doesn't touch the submodule, so clean it by hand too.
make clean >/dev/null 2>&1 || true
rm -f YAME/libyame.a YAME/htslib/libhts.a
find YAME -name '*.o' -delete 2>/dev/null || true

# conda supplies CC plus the sysroot / $PREFIX include+lib flags in CPPFLAGS /
# CFLAGS / LDFLAGS. The vendored YAME and htslib Makefiles set rigid CFLAGS and
# ignore CPPFLAGS, so bake conda's *compile* flags into CC -- every level of
# sub-make (sesame -> YAME -> htslib) uses CC verbatim, and command-line CC
# propagates to all sub-makes. LDFLAGS reaches the final link on its own:
# sesame's Makefile does `LDFLAGS += $(EXTRA_LDFLAGS)`, which appends to the
# LDFLAGS make already imports from the environment (leave EXTRA_LDFLAGS unset so
# it is not added twice).
make -j"${CPU_COUNT:-1}" CC="${CC} ${CPPFLAGS} ${CFLAGS}"

make install PREFIX="${PREFIX}"
