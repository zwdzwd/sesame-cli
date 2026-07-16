#!/bin/sh
# Golden level-1: the C IDAT reader must agree with R's readIDAT() bit-for-bit
# on IlluminaID / Mean / SD / NBeads. No tolerance -- these are integers.
#
# Requires Rscript with sesame + sesameData installed (the R oracle).
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesamec"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }

extdata=$(Rscript -e 'cat(system.file("extdata","",package="sesameData"))' 2>/dev/null)
[ -n "$extdata" ] || { echo "SKIP: sesameData not installed"; exit 0; }

fail=0
pass=0

for f in "$extdata"/*.idat; do
    name=$(basename "$f")

    # R oracle: dump addr/mean/sd/nbeads in file order, no header.
    Rscript --vanilla -e '
        suppressMessages(library(sesame))
        f <- commandArgs(trailingOnly=TRUE)[1]
        r <- suppressWarnings(sesame:::readIDAT(f))
        q <- r$Quants
        cat(sprintf("%s\t%d\t%d\t%d",
            rownames(q), q[,"Mean"], q[,"SD"], q[,"NBeads"]), sep="\n")
    ' "$f" > "$work/r.tsv" 2>"$work/r.err" || {
        echo "FAIL $name: R oracle errored"; sed 's/^/    /' "$work/r.err"; fail=$((fail+1)); continue; }

    "$bin" idat-dump --tsv "$f" > "$work/c.tsv" || {
        echo "FAIL $name: sesamec errored"; fail=$((fail+1)); continue; }

    if cmp -s "$work/r.tsv" "$work/c.tsv"; then
        n=$(wc -l < "$work/c.tsv" | tr -d ' ')
        echo "ok   $name ($n records, bit-identical)"
        pass=$((pass+1))
    else
        echo "FAIL $name: C and R differ"
        echo "  first differing lines (R | C):"
        diff "$work/r.tsv" "$work/c.tsv" | head -6 | sed 's/^/    /'
        fail=$((fail+1))
    fi
done

# gz round-trip: the same file gzipped must read identically.
src=$(ls "$extdata"/*_Grn.idat 2>/dev/null | head -1)
if [ -n "${src:-}" ]; then
    gzip -c "$src" > "$work/t.idat.gz"
    "$bin" idat-dump --tsv "$src"          > "$work/plain.tsv"
    "$bin" idat-dump --tsv "$work/t.idat.gz" > "$work/gz.tsv"
    if cmp -s "$work/plain.tsv" "$work/gz.tsv"; then
        echo "ok   gz round-trip (plain == gz)"
        pass=$((pass+1))
    else
        echo "FAIL gz round-trip"
        fail=$((fail+1))
    fi
fi

echo
echo "passed $pass, failed $fail"
[ "$fail" -eq 0 ]
