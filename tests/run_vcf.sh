#!/bin/sh
# vcf: validate SNP genotyping (sesame formatVCF) against R on real EPICv2 data.
# GT (the call) and PVF (variant fraction) must match R exactly; GS (the quality
# score) matches except in the deep tail, where R's dbinom is more precise than
# the log-pmf here. Needs an EPICv2 ordering, the SNP annotation, an EPICv2 IDAT,
# and Rscript with sesame; skips cleanly otherwise.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(dirname "$here")
bin="$root/sesame"
store=${SESAME_INDEX_DIR:-$root/data}
idats=${SESAME_TEST_IDATS:-$HOME/repo/InfiniumTestIDATs}
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

[ -x "$bin" ] || { echo "FAIL: $bin not built"; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "SKIP vcf: no Rscript"; exit 0; }

plat=EPICv2
ord="$store/$plat/$plat.ordering.tsv.gz"
[ -f "$ord" ] || ord="$root/testdata/$plat.ordering.tsv.gz"
pfx="$idats/$plat/206909630042_R08C01"

# locate the SNP annotation (not yet in the fetch set): store, $SESAME_SNP_DIR, or a sibling InfiniumAnnotation checkout
snp=""
for c in "$store/$plat/$plat.hg38.snp.tsv.gz" \
         "${SESAME_SNP_DIR:-}/$plat.hg38.snp.tsv.gz" \
         "$HOME/repo/InfiniumAnnotation/$plat/$plat.hg38.snp.tsv.gz"; do
    [ -f "$c" ] && { snp="$c"; break; }
done
[ -n "$snp" ] || { echo "SKIP vcf: no $plat.hg38.snp.tsv.gz (set SESAME_SNP_DIR)"; exit 0; }
[ -f "$ord" ] || { echo "SKIP vcf: no ordering $ord"; exit 0; }
if [ ! -f "$pfx"_Grn.idat ] && [ ! -f "$pfx"_Grn.idat.gz ]; then
    echo "SKIP vcf: no IDAT $pfx"; exit 0; fi

# C: full VCF, parse INFO -> Probe_ID GT GS PVF
"$bin" vcf "$pfx" --index "$ord" --snp "$snp" 2>/dev/null | grep -v '^#' | \
  python3 -c "
import sys
print('Probe_ID\tGT\tGS\tPVF')
for l in sys.stdin:
    p=l.rstrip().split('\t'); d=dict(kv.split('=',1) for kv in p[7].split(';') if '=' in kv)
    print(d['Probe_ID'],d['GT'],p[5],d['PVF'],sep='\t')" > "$work/c.tsv"

# R oracle
Rscript "$here/compare_vcf.R" "$pfx" $plat "$snp" "$work/r.tsv" >/dev/null 2>&1

python3 - "$work/c.tsv" "$work/r.tsv" <<'PY'
import sys
def load(fn):
    d={}
    for l in open(fn):
        p=l.rstrip('\n').split('\t')
        if p[0]=='Probe_ID': continue
        d[p[0]]=(p[1],p[2],p[3])
    return d
C=load(sys.argv[1]); R=load(sys.argv[2])
common=set(C)&set(R)
if not common: print("FAIL: no shared probes"); sys.exit(1)
gt_mis=pvf_max=0.0; gt_mis=0; gs_off=0
for k in common:
    cg,cs,cp=C[k]; rg,rs,rp=R[k]
    if cg!=rg: gt_mis+=1
    try:
        if abs(int(cs)-int(rs))>1: gs_off+=1
    except ValueError: pass
    try: pvf_max=max(pvf_max, abs(float(cp)-float(rp)))
    except ValueError: pass
n=len(common)
print(f"vcf vs R formatVCF: n={n} GT-mismatch={gt_mis} PVF-maxdiff={pvf_max:.1e} GS-off(>1)={gs_off}")
ok = (gt_mis==0 and pvf_max < 1e-6 and gs_off <= n*0.001)
if not ok:
    print(f"FAIL: GT={gt_mis} PVF={pvf_max:.1e} GS-off={gs_off} (>{int(n*0.001)})"); sys.exit(1)
print(f"ok   vcf: {n} probes, GT + PVF exact vs R, GS off on {gs_off} (deep-tail dbinom)")
PY
rc=$?

echo
echo "passed $([ $rc -eq 0 ] && echo 1 || echo 0), failed $([ $rc -eq 0 ] && echo 0 || echo 1)"
[ $rc -eq 0 ]
