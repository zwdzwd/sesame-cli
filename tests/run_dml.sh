#!/bin/sh
# Level (DML): per-probe differential methylation vs R's DML. No data-lineage
# caveat -- the betas matrix is the input -- so the OLS estimates, t/F p-values,
# effect sizes, and BH adjustment must match R to ~1e-8. The R oracle generates a
# small dataset, runs `sesame dml`, and compares (see compare_dml.R).
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

if Rscript --vanilla "$here/compare_dml.R" "$bin"; then
    echo "passed 1, failed 0"
else
    echo "passed 0, failed 1"; exit 1
fi
