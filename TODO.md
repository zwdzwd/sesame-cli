# TODO

## QC plots (sesameQC_plot*)

Implement QC visualizations like the ones in the sesame vignette's Quality Control
section: <https://zhou-lab.github.io/sesame-docs/dev/supplemental.html#Quality_Control>

R exposes these as `sesameQC_plot*`; the CLI should emit plot-ready TSVs and feed
them to **cinderplot** (same pattern as `region`/`cnv`), not draw them itself:

- `sesameQC_plotBar` — a sample's QC metrics as a labeled bar panel (from `qc.tsv`).
- `sesameQC_plotBetaByDesign` — beta density split by Infinium design (I vs II).
- `sesameQC_plotIntensVsBetas` — per-probe total intensity vs beta (detection view).
- `sesameQC_plotRedGrnQQ` — Red-vs-Green intensity Q-Q (dye-bias diagnostic).
- `sesameQC_plotHeatSNPs` — rs-probe genotype heatmap across samples (sample ID/mixup).

Most of the underlying numbers are already computed (`qc.tsv`, `intensity.cg`,
betas, the SNP genotypes from `vcf`); this is mostly a data-shaping + cinderplot
recipe task. `plotBar` needs only `qc.tsv`; the others need the per-probe signal.
