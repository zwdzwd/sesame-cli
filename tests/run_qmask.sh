#!/bin/sh
# Level-3 (Q): sesame's qualityMask must mask EXACTLY the union of the recommended
# YAME mask tracks, applied to the ordering.
#
# This is a self-consistency gate, not a differential test against R: the
# published .cm is a newer mask lineage than sesameData's KYCG object (see
# NUMERICS.md, "mask lineage"), so it cannot be bit-identical to R's qualityMask.
# What it CAN pin exactly is that sesame reads and unions the .cm correctly and
# aligns it to the ordering.
#
# Needs: the `yame` binary on PATH, and a fetched store with the platform's .cm
# (e.g. data/MSA/). Skips cleanly if either is missing.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
dump="$root/pipeline_dump"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v yame >/dev/null 2>&1 || { echo "SKIP Q: yame not on PATH"; exit 0; }

# platform  sample-prefix  <recommended track names...>
run_one() {
    plat=$1; rel=$2; shift 2
    cm=$(ls "$store/$plat"/*.cm 2>/dev/null | head -1 || true)
    pfx="$idats/$rel"
    if [ -z "$cm" ]; then echo "SKIP $plat Q: no .cm in $store/$plat"; return; fi
    if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
        echo "SKIP $plat Q: no IDAT $pfx"; return; fi

    # expected: ordering Probe_IDs where the yame union of the given tracks is set
    printf '%s\n' "$@" > "$work/names.txt"
    gzcat "$store/$plat/$plat.ordering.tsv.gz" | tail -n +2 | cut -f1 > "$work/ids.txt"
    yame unpack -l "$work/names.txt" "$cm" > "$work/tab.txt" 2>/dev/null
    python3 - "$work/ids.txt" "$work/tab.txt" "$work/expected.txt" <<'PY'
import sys
ids=[l.rstrip("\n") for l in open(sys.argv[1])]
exp=set()
for i,row in enumerate(open(sys.argv[2])):
    if i>=len(ids): break
    if any(c!="0" for c in row.rstrip("\n").split("\t")): exp.add(ids[i])
open(sys.argv[3],"w").write("\n".join(sorted(exp))+"\n")
PY

    # actual: probes sesame masks (NA) under prep=Q
    SESAME_INDEX_DIR="$store" "$dump" --prep Q --what beta "$pfx" 2>/dev/null \
      | python3 -c "import sys; print('\n'.join(sorted(l.split('\t')[0] for l in sys.stdin if l.rstrip('\n').split('\t')[1]=='NA')))" \
      > "$work/actual.txt"

    if cmp -s "$work/expected.txt" "$work/actual.txt"; then
        echo "ok   $plat Q  $(wc -l < "$work/expected.txt" | tr -d ' ') probes = yame union of $# tracks"
        PASS=$((PASS+1))
    else
        echo "FAIL $plat Q: masked set != yame union"
        echo "  expected $(wc -l < "$work/expected.txt"|tr -d ' '), got $(wc -l < "$work/actual.txt"|tr -d ' ')"
        FAIL=$((FAIL+1))
    fi
}

PASS=0; FAIL=0
run_one MSA MSA/207760740030_R01C03 \
    M_1baseSwitchSNPcommon_5pt M_2extBase_SNPcommon_5pt M_mapping M_nonuniq M_SNPcommon_5pt

echo
echo "passed $PASS, failed $FAIL"
[ "$FAIL" -eq 0 ]
