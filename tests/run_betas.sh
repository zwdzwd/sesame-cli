#!/bin/sh
# Golden: betas out of C must be bit-identical to R with prep="".
#
# Needs the exported ordering tables (make index) and Rscript with sesame.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesamec"
idats=${SESAMEC_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
[ -d "$idats" ] || { echo "SKIP: $idats not found"; exit 0; }

fail=0
pass=0

# platform <space> prefix-relative-to-$idats
while IFS=' ' read -r plat rel; do
    [ -n "$plat" ] || continue
    idx="$root/data/$plat.ordering.tsv.gz"
    pfx="$idats/$rel"
    if [ ! -f "$idx" ]; then
        echo "SKIP $plat: no $idx (run: make index)"
        continue
    fi
    if ! "$bin" betas --index "$idx" --f64 "$pfx" > "$work/c.f64" 2>"$work/c.err"; then
        echo "FAIL $plat: sesamec errored"; sed 's/^/    /' "$work/c.err" | head -3
        fail=$((fail+1)); continue
    fi
    if Rscript --vanilla "$here/compare_betas.R" "$plat" "$pfx" "$work/c.f64" \
         2>"$work/r.err"; then
        pass=$((pass+1))
    else
        sed 's/^/    /' "$work/r.err" | head -3
        fail=$((fail+1))
    fi
done <<'EOF'
HM450 HM450/3999492009_R01C01
EPIC EPIC/GSM2995280_201868590258_R01C01
EPICv2 EPICv2/206909630040_R03C01
EPICv2 EPICv2/206909630042_R08C01
EPICv2 EPICv2/MaleBlood/206909630014_R05C01
MSA MSA/207760740030_R01C03
EOF

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
