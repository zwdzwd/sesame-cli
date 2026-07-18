#!/bin/sh
# cbs: validate the deterministic circular binary segmentation. Two checks:
#  (1) exactness on a clean synthetic step signal (no RNG, must be exact);
#  (2) agreement with DNAcopy::segment on REAL K562 copy-number bin signals --
#      C must recover >= 90% of R's segment breakpoints (within 1 bin). DNAcopy's
#      significance is permutation-based (and cnSegmentation never seeds it), so
#      exact match is impossible; recall of R's breakpoints is the meaningful gate.
# Needs Rscript + sesameData for (2); (1) always runs.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
cbs="$root/cbs_test"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$cbs" ] || { echo "FAIL: $cbs not built"; exit 1; }
fail=0

# (1) clean step: 60 x 0, 40 x -1.5, 80 x 0  ->  exactly three segments at 60,100,180
python3 -c "
for _ in range(60): print(0.0)
for _ in range(40): print(-1.5)
for _ in range(80): print(0.0)" | "$cbs" > "$work/step.out"
ends=$(cut -f2 "$work/step.out" | tr '\n' ',')
means=$(cut -f4 "$work/step.out" | tr '\n' ',')
if [ "$ends" != "60,100,180," ] || [ "$means" != "0,-1.5,0," ]; then
    echo "FAIL: clean step -> ends=$ends means=$means (want 60,100,180 / 0,-1.5,0)"; fail=1
else
    echo "ok   clean step signal: exact 3 segments (60,100,180)"
fi

# (2) real K562 vs DNAcopy
if command -v Rscript >/dev/null 2>&1 && Rscript -e 'quit(status=!requireNamespace("DNAcopy",quietly=TRUE) || !requireNamespace("sesameData",quietly=TRUE))' 2>/dev/null; then
    if Rscript "$here/compare_cbs.R" "$work/cbs" >/dev/null 2>&1; then
        python3 - "$cbs" "$work/cbs.sig" "$work/cbs.seg" <<'PY'
import sys, subprocess, collections
cbs, sigf, segf = sys.argv[1], sys.argv[2], sys.argv[3]
sig = collections.OrderedDict()
for l in open(sigf):
    c, v = l.rstrip("\n").split("\t"); sig.setdefault(c, []).append(v)
rseg = collections.defaultdict(list)
for l in open(segf):
    c, e = l.rstrip("\n").split("\t"); rseg[c].append(int(e))
R_int = C_int = recall = 0
for c, vals in sig.items():
    out = subprocess.run([cbs], input="\n".join(vals)+"\n", capture_output=True, text=True).stdout
    cb = [int(l.split("\t")[1]) for l in out.splitlines()][:-1]   # interior breakpoints
    rb = rseg[c][:-1]
    R_int += len(rb); C_int += len(cb)
    for b in rb:
        if any(abs(b-cc) <= 1 for cc in cb): recall += 1
rc = recall / R_int if R_int else 1.0
print(f"cbs vs DNAcopy (K562): R breakpoints={R_int} C breakpoints={C_int} recall={rc*100:.1f}%")
if rc < 0.90:
    print(f"FAIL: recall {rc*100:.1f}% < 90%"); sys.exit(1)
if C_int > 2 * R_int:
    print(f"FAIL: C over-segments ({C_int} > 2x {R_int})"); sys.exit(1)
print(f"ok   cbs: recovers {rc*100:.1f}% of DNAcopy's K562 breakpoints ({C_int} vs {R_int} segments)")
PY
        [ $? -eq 0 ] || fail=1
    else
        echo "SKIP cbs K562: compare_cbs.R failed (sesameData/EPICv2.8.SigDF unavailable)"
    fi
else
    echo "SKIP cbs K562: no Rscript with DNAcopy + sesameData"
fi

echo
echo "passed $([ $fail -eq 0 ] && echo 1 || echo 0), failed $([ $fail -eq 0 ] && echo 0 || echo 1)"
[ $fail -eq 0 ]
