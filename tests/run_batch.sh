#!/bin/sh
# Batch: `preprocess` over many samples in one process must be deterministic
# across thread counts and degrade gracefully on a bad sample. Pure
# self-consistency; three distinct EPICv2 samples with prep="" (no mask store).
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
idx="$root/testdata/EPICv2.ordering.tsv.gz"
yame="$root/YAME/yame"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$yame" ] || make -C "$root/YAME" >/dev/null 2>&1 || true
P1="$idats/EPICv2/206909630040_R03C01"
P2="$idats/EPICv2/206909630042_R08C01"
P3="$idats/EPICv2/MaleBlood/206909630014_R05C01"
[ -f "$idx" ] || { echo "SKIP batch: no $idx"; exit 0; }
for p in "$P1" "$P2" "$P3"; do
    if [ ! -f "$p"_Grn.idat ] && [ ! -f "$p"_Grn.idat.gz ]; then
        echo "SKIP batch: no IDAT $p"; exit 0; fi
done

PASS=0; FAIL=0
check() { if [ "$2" -eq 0 ]; then echo "ok   $1"; PASS=$((PASS+1))
          else echo "FAIL $1"; FAIL=$((FAIL+1)); fi; }

pp() { "$bin" preprocess --index "$idx" --prep "" --output beta "$@" 2>/dev/null; }

# (1) determinism: --threads 1 and --threads 4 give a byte-identical beta.cg
# (each sample's cdata is deterministic and written in argv order)
mkdir "$work/t1" "$work/t4"
pp --threads 1 --out "$work/t1" "$P1" "$P2" "$P3"
pp --threads 4 --out "$work/t4" "$P1" "$P2" "$P3"
cmp -s "$work/t1/beta.cg" "$work/t4/beta.cg" && rc=0 || rc=1
check "--threads 1 == --threads 4 (byte-identical .cg)" "$rc"

# (2) values agree with a single-sample run, per sample, via yame named selection
mkdir "$work/s2"
pp --out "$work/s2" "$P2"
"$yame" unpack "$work/s2/beta.cg" 2>/dev/null > "$work/single2.txt" || true
"$yame" unpack "$work/t1/beta.cg" 206909630042_R08C01 2>/dev/null > "$work/batch2.txt" || true
if [ -x "$yame" ] && [ -s "$work/batch2.txt" ] && [ -s "$work/single2.txt" ]; then
    cmp -s "$work/single2.txt" "$work/batch2.txt" && rc=0 || rc=1
    check "batch sample values == single-sample (yame named)" "$rc"
fi

# (3) a bad sample in the middle -> exit 1, but all three blocks written (bad=NA)
pp --threads 2 --out "$work/bad" "$P1" "$idats/NOPE_BADSAMPLE" "$P3" && rc=0 || rc=$?
[ "${rc:-0}" -eq 1 ] && rc=0 || rc=1
check "bad sample -> exit 1" "$rc"
[ "$(wc -l < "$work/bad/beta.cg.idx")" -eq 3 ] && rc=0 || rc=1
check "bad batch still writes 3 blocks" "$rc"

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
