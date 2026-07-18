#!/bin/sh
# attach-probe: prepend the ordering's Probe_ID to a positional file's rows.
# Fully self-contained -- builds a tiny ordering, a format-3 .cg (via mu2cg), and
# a text table, then checks the labeled output and the lineage-mismatch guard. No
# IDATs, no store, no R oracle.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
mu2cg="$root/mu2cg"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ]   || { echo "FAIL: $bin not built"; exit 1; }
[ -x "$mu2cg" ] || { echo "SKIP attach: $mu2cg not built"; exit 0; }

# --- fixtures: 5 probes, a 5-col ordering (mask inline) and a 4-col one -------
printf 'Probe_ID\tM\tU\tcol\tmask\n' > "$work/ord5.tsv"
printf 'Probe_ID\tM\tU\tcol\n'       > "$work/ord4.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\tNA\t%d\t2\t0\n' "$i" "$((1000+i))" >> "$work/ord5.tsv"
    printf 'cg%07d_BC21\tNA\t%d\t2\n'    "$i" "$((1000+i))" >> "$work/ord4.tsv"
    i=$((i+1))
done
gzip -f "$work/ord5.tsv" "$work/ord4.tsv"

# a format-3 .cg from a tiny M/U table (one sample "S1")
printf 'Probe_ID\tS1_M\tS1_U\n' > "$work/mu.tsv"
i=1
while [ $i -le 5 ]; do
    printf 'cg%07d_BC21\t%d\t%d\n' "$i" "$((100*i))" "$((10*i))" >> "$work/mu.tsv"
    i=$((i+1))
done
"$mu2cg" "$work/mu.tsv" "$work/s1.cg" >/dev/null 2>&1

# a text table (header + 5 data rows), coord-style, no Probe_ID column
printf 'chrm\tpos\n' > "$work/coord.tsv"
i=1
while [ $i -le 5 ]; do printf 'chr1\t%d\n' "$((5000+i))" >> "$work/coord.tsv"; i=$((i+1)); done

fail=0

# 1) text: Probe_ID prepended, header kept, order preserved
"$bin" attach-probe --index "$work/ord5.tsv.gz" "$work/coord.tsv" > "$work/coord.out" 2>/dev/null
head1=$(sed -n '1p' "$work/coord.out")
[ "$head1" = "$(printf 'Probe_ID\tchrm\tpos')" ] || { echo "FAIL: text header '$head1'"; fail=1; }
[ "$(sed -n '2p' "$work/coord.out")" = "$(printf 'cg0000001_BC21\tchr1\t5001')" ] || { echo "FAIL: text row1"; fail=1; }
[ "$(wc -l < "$work/coord.out")" -eq 6 ] || { echo "FAIL: text nrow"; fail=1; }

# 2) text against a 4-column ordering (mask moved to .cm) parses the same
"$bin" attach-probe --index "$work/ord4.tsv.gz" "$work/coord.tsv" > "$work/coord4.out" 2>/dev/null
cmp -s "$work/coord.out" "$work/coord4.out" || { echo "FAIL: 4-col vs 5-col ordering differ"; fail=1; }

# 3) fmt3 .cg default -> M<TAB>U columns, header "<name>_M <name>_U"
"$bin" attach-probe --index "$work/ord5.tsv.gz" "$work/s1.cg" > "$work/cg.out" 2>/dev/null
[ "$(sed -n '1p' "$work/cg.out")" = "$(printf 'Probe_ID\tS1_M\tS1_U')" ] || { echo "FAIL: cg header"; fail=1; }
[ "$(sed -n '2p' "$work/cg.out")" = "$(printf 'cg0000001_BC21\t100\t10')" ] || { echo "FAIL: cg M/U row1"; fail=1; }

# 4) fmt3 .cg --beta -> beta = M/(M+U); probe 1 = 100/110
"$bin" attach-probe --beta --index "$work/ord5.tsv.gz" "$work/s1.cg" > "$work/beta.out" 2>/dev/null
b=$(sed -n '2p' "$work/beta.out" | cut -f2)
case "$b" in 0.9090*) : ;; *) echo "FAIL: cg beta '$b' != 0.9090..."; fail=1 ;; esac

# 5) lineage mismatch: a 4-probe ordering vs the 5-row file must fail, no stdout
printf 'Probe_ID\tM\tU\tcol\n' > "$work/ord_short.tsv"
i=1; while [ $i -le 4 ]; do printf 'cgX%06d_BC21\tNA\t1\t2\n' "$i" >> "$work/ord_short.tsv"; i=$((i+1)); done
gzip -f "$work/ord_short.tsv"
if "$bin" attach-probe --index "$work/ord_short.tsv.gz" "$work/coord.tsv" > "$work/bad.out" 2>/dev/null; then
    echo "FAIL: mismatch did not error"; fail=1
fi
[ -s "$work/bad.out" ] && { echo "FAIL: mismatch emitted output before erroring"; fail=1; } || true

echo
if [ $fail -eq 0 ]; then
    echo "ok   attach-probe: text + 4/5-col ordering + fmt3 M/U + beta + mismatch guard"
    echo "passed 1, failed 0"
else
    echo "passed 0, failed 1"
fi
[ $fail -eq 0 ]
