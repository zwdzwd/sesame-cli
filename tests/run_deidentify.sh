#!/bin/sh
# deidentify: the SNP (rs) probe bead means are zeroed (default) or reversibly
# scrambled (--randomize), every other bead is byte-identical, and -r restores a
# scrambled file exactly. Needs a store with the EPICv2 ordering, an EPICv2 IDAT,
# and python3; skips cleanly otherwise.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "SKIP deidentify: no python3"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
pfx="$idats/$plat/206909630042_R08C01"
[ -f "$ord" ] || { echo "SKIP deidentify: missing $ord (fetch $plat)"; exit 0; }
g=""
[ -f "$pfx"_Grn.idat ] && g="$pfx"_Grn.idat
[ -f "$pfx"_Grn.idat.gz ] && g="$pfx"_Grn.idat.gz
[ -n "$g" ] || { echo "SKIP deidentify: no IDAT $pfx"; exit 0; }

"$bin" deidentify                       --index "$ord" "$g"           "$work/deid.idat" 2>/dev/null
"$bin" deidentify --randomize --seed 999 --index "$ord" "$g"          "$work/rnd.idat"  2>/dev/null
"$bin" deidentify -r          --seed 999 --index "$ord" "$work/rnd.idat" "$work/reid.idat" 2>/dev/null

python3 - "$bin" "$ord" "$g" "$work/deid.idat" "$work/rnd.idat" "$work/reid.idat" <<'PY'
import sys, gzip, subprocess, collections
b, ordp, origp, deidp, rndp, reidp = sys.argv[1:7]
snp = set()
for l in gzip.open(ordp, 'rt').read().splitlines()[1:]:
    f = l.split('\t')
    if f[0].startswith('rs'):
        for a in (f[1], f[2]):
            if a != 'NA': snp.add(int(a))
def dump(p):
    out = subprocess.check_output([b, "idat-dump", "--tsv", p]).decode()
    return {int(t[0]): int(t[1]) for t in (x.split() for x in out.splitlines())
            if len(t) == 4 and t[0].isdigit()}
o, d, r, ri = dump(origp), dump(deidp), dump(rndp), dump(reidp)
si = [a for a in snp if a in o]
zero  = bool(si) and all(d[a] == 0 for a in si) and sum(o[a] for a in si) > 0
intact = all(d.get(a) == o[a] for a in o if a not in snp)
scram = (collections.Counter(o[a] for a in si) == collections.Counter(r[a] for a in si)
         and not all(o[a] == r[a] for a in si))
rtrip = all(ri.get(a) == o[a] for a in o)
print("deidentify: %d SNP beads; zero=%s non-SNP-intact=%s scramble=%s roundtrip=%s"
      % (len(si), zero, intact, scram, rtrip))
if not (zero and intact and scram and rtrip):
    print("FAIL"); sys.exit(1)
print("ok   deidentify: SNP beads zeroed/scrambled, non-SNP intact, -r round-trips exactly")
PY
rc=$?

echo
echo "passed $([ $rc -eq 0 ] && echo 1 || echo 0), failed $([ $rc -eq 0 ] && echo 0 || echo 1)"
[ $rc -eq 0 ]
