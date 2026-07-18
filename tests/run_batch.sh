#!/bin/sh
# Batch (R4): many prefixes in one process must equal N single-sample runs, be
# deterministic across thread counts, and degrade gracefully on a bad sample.
# Pure self-consistency -- no R oracle -- so it is fast.
#
# Uses three distinct EPICv2 samples with the exported ordering (prep=""), so it
# needs no mask store. Skips cleanly if the ordering or IDATs are absent.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
idx="$root/testdata/EPICv2.ordering.tsv.gz"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

P1="$idats/EPICv2/206909630040_R03C01"
P2="$idats/EPICv2/206909630042_R08C01"
P3="$idats/EPICv2/MaleBlood/206909630014_R05C01"

[ -f "$idx" ] || { echo "SKIP batch: no $idx"; exit 0; }
for p in "$P1" "$P2" "$P3"; do
    if [ ! -f "$p"_Grn.idat ] && [ ! -f "$p"_Grn.idat.gz ]; then
        echo "SKIP batch: no IDAT $p"; exit 0; fi
done

PASS=0; FAIL=0
check() { # msg  cond(0=ok)
    if [ "$2" -eq 0 ]; then echo "ok   $1"; PASS=$((PASS+1))
    else echo "FAIL $1"; FAIL=$((FAIL+1)); fi
}

# single-sample runs, concatenated = the expected sample-major matrix
for i in 1 2 3; do
    eval "p=\$P$i"
    "$bin" betas --index "$idx" --f64 "$p" 2>/dev/null > "$work/s$i.f64"
done
cat "$work/s1.f64" "$work/s2.f64" "$work/s3.f64" > "$work/expected.f64"

# (1) batch --f64 == concatenation of single runs (column identity + ordering)
"$bin" betas --index "$idx" --threads 4 --f64 "$P1" "$P2" "$P3" 2>/dev/null > "$work/b4.f64"
cmp -s "$work/expected.f64" "$work/b4.f64"; check "batch == N single runs (byte-identical)" $?

# (2) determinism across thread counts
"$bin" betas --index "$idx" --threads 1 --f64 "$P1" "$P2" "$P3" 2>/dev/null > "$work/b1.f64"
cmp -s "$work/b1.f64" "$work/b4.f64"; check "--threads 1 == --threads 4" $?

# (3) text matrix: correct header + dimensions
"$bin" betas --index "$idx" "$P1" "$P2" "$P3" 2>/dev/null > "$work/m.txt"
n=$(( $(gzip -dc "$idx" | wc -l) - 1 ))
python3 - "$work/m.txt" "$n" <<'PY'
import sys
lines=open(sys.argv[1]).read().splitlines(); n=int(sys.argv[2])
hdr=lines[0].split('\t')
ok = (hdr==['Probe_ID','206909630040_R03C01','206909630042_R08C01','206909630014_R05C01']
      and len(lines)==n+1 and all(len(l.split('\t'))==4 for l in lines[1:]))
sys.exit(0 if ok else 1)
PY
check "text matrix header + dimensions" $?

# (4) a bad prefix in the middle -> NA column, neighbors intact, exit 1
set +e
"$bin" betas --index "$idx" --threads 2 --f64 "$P1" "$idats/NOPE_BADSAMPLE" "$P3" \
    2>/dev/null > "$work/bad.f64"
rc=$?
set -e
[ "$rc" -eq 1 ]; check "bad sample -> exit 1" $?
python3 - "$work/bad.f64" "$work/s1.f64" "$work/s3.f64" <<'PY'
import sys, math
b=open(sys.argv[1],'rb').read(); s1=open(sys.argv[2],'rb').read(); s3=open(sys.argv[3],'rb').read()
n=len(s1)//8; col=lambda k: b[k*len(s1):(k+1)*len(s1)]
import struct
mid=struct.unpack('<%dd'%n, col(1))
ok = (col(0)==s1 and col(2)==s3 and all(math.isnan(x) for x in mid))
sys.exit(0 if ok else 1)
PY
check "bad column all-NA; neighbors byte-identical to single runs" $?

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
