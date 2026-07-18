#!/bin/sh
# .cg output: `sesame intensity --cg` must write a YAME format-4 .cg that the
# yame toolchain reads back to the same values. Builds the yame binary from the
# submodule for the round-trip; skips cleanly if that or the store is missing.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
yame="$root/YAME/yame"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
plat=MSA; rel=MSA/207760740030_R01C03; pfx="$idats/$rel"
[ -d "$store/$plat" ] || { echo "SKIP cg: no store $store/$plat"; exit 0; }
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP cg: no IDAT $pfx"; exit 0; fi
[ -x "$yame" ] || make -C "$root/YAME" >/dev/null 2>&1 || true
[ -x "$yame" ] || { echo "SKIP cg: could not build $yame"; exit 0; }

name=$(basename "$pfx")
# raw-signal reference (double), and the two .cg via preprocess (raw signal)
SESAME_INDEX_DIR="$store" "$dump" --prep "" --what total "$pfx" 2>/dev/null > "$work/tot.tsv"
SESAME_INDEX_DIR="$store" "$bin" preprocess --prep "" --raw-signal \
    --output intensity,total_intensity --out "$work" "$pfx" 2>/dev/null

"$yame" unpack -f -1 "$work/intensity.cg"       2>/dev/null > "$work/mu.txt"    # M<TAB>U (fmt3)
"$yame" unpack       "$work/total_intensity.cg" 2>/dev/null > "$work/f4.txt"    # float total (fmt4)

python3 - "$work/tot.tsv" "$work/mu.txt" "$work/f4.txt" "$work/intensity.cg.idx" "$name" <<'PY'
import sys, math
tot = [l.rstrip("\n").split("\t")[1] for l in open(sys.argv[1])]
mu  = [l.split("\t") for l in open(sys.argv[2]).read().splitlines()]
f4  = [l.strip() for l in open(sys.argv[3])]
idx = open(sys.argv[4]).read().split("\t")[0]
name = sys.argv[5]
if not (len(tot) == len(mu) == len(f4)):
    print(f"FAIL: nrow tot={len(tot)} mu={len(mu)} f4={len(f4)}"); sys.exit(1)
if idx != name: print(f"FAIL: idx name '{idx}' != '{name}'"); sys.exit(1)
bad3 = bad4 = 0
for t, (m, u), v4 in zip(tot, mu, f4):
    tn = math.nan if t == "NA" else float(t)
    s3 = float(m) + float(u)                 # format 3: M+U == total (NA -> 0,0 -> 0)
    v  = math.nan if v4 in ("NA","-1","-1.0","nan") else float(v4)
    if math.isnan(tn):
        if s3 != 0: bad3 += 1
        if not math.isnan(v): bad4 += 1
    else:
        if abs(s3 - tn) > 0.5: bad3 += 1
        if math.isnan(v) or abs(v - tn) > 0.5: bad4 += 1
if bad3 or bad4:
    print(f"FAIL: format3 M+U mismatches={bad3}, format4 mismatches={bad4}"); sys.exit(1)
print(f"ok   MSA/{name} .cg round-trip: {len(mu)} probes; fmt3 M+U and fmt4 total match")
PY
rc=$?

echo
echo "passed $([ $rc -eq 0 ] && echo 1 || echo 0), failed $([ $rc -eq 0 ] && echo 0 || echo 1)"
[ $rc -eq 0 ]
