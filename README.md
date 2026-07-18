# sesame-cli

![license](https://img.shields.io/badge/license-AGPL--3.0--or--later-blue)
![language](https://img.shields.io/badge/C-C11-00599C)
![arrays](https://img.shields.io/badge/arrays-EPIC%20%7C%20EPICv2%20%7C%20HM450%20%7C%20MSA-brightgreen)
![vs R](https://img.shields.io/badge/betas%20vs%20R-bit--identical-success)

A standalone C implementation of [sesame](https://github.com/zwdzwd/sesame)'s
basic Infinium DNA methylation preprocessing: **IDAT → betas**, with no R and no
Bioconductor. Builds one binary, `sesame`.

> **Status.** The IDAT reader and `betas` with no preprocessing are
> **bit-identical to R**, and the full **`QCDPB`** pipeline is implemented: `Q`
> (self-consistent vs R by mask lineage), `C` (channel-identical), `D` (within
> ~2 ULP), `P` (detection-p exact; residual is 0.05-boundary/lineage), and `B`
> (noob — `normExpSignal`/`huber` proven to a few ULP in isolation; residual is
> lineage only). Multi-sample **batch** runs are parallel and deterministic. See
> [Fidelity and validation](#fidelity-and-validation) and `NUMERICS.md`.

## Contents

- [Why sesame-cli](#why-sesame-cli)
- [Scope](#scope)
- [Installation](#installation)
- [Quickstart](#quickstart)
- [Commands](#commands)
- [The QCDPB pipeline](#the-qcdpb-pipeline)
- [Batch mode](#batch-mode)
- [Output formats](#output-formats)
- [Data files and the store](#data-files-and-the-store)
- [Configuration](#configuration)
- [Exit status](#exit-status)
- [Fidelity and validation](#fidelity-and-validation)
- [Development](#development)
- [License](#license)

## Why sesame-cli

`openSesame` needs R + Bioconductor + `sesameData` (ExperimentHub) +
`preprocessCore`/`MASS`/`BiocParallel`. That blocks:

- **Pipeline deployment** — a mostly-static binary for Nextflow/Snakemake/cloud.
- **Scale** — 10k–100k+ samples. The math is already compiled; the glue is not
  (data.frame copies per prep step, `match()` on ~937k *strings*).
- **Release cadence** — Bioconductor ships only release+devel, so patches never
  reach users on older setups.
- **Non-R users** — Python users currently reach for `methylprep` instead.

The R package is **not** modified by this project. It stays the reference
implementation and the validation oracle, permanently. See `NUMERICS.md` for
every intentional difference.

## Scope

**In:** `readIDATpair` → `prepSesame(sdf, "QCDPB")` → `getBetas`, the `sesameQC`
panel, and per-probe differential methylation (`DML`), for EPIC, EPICv2, HM450,
MSA.
**Out:** DMR (region calling — needs a genomic-coordinate annotation not yet
hosted), KYCG, CNV, visualization, and all inference (species/strain/age/sex/
ethnicity) — those stay in R.

## Installation

Dependencies: a C compiler + `make`, **zlib**, **libcurl**, and the **YAME**
submodule (bundled htslib) — linked in to read `.cm` mask files. Both sesame-cli
and YAME are AGPL-3.0, so linking is clean.

```sh
git clone --recurse-submodules <repo>          # or: git submodule update --init --recursive
cd sesame-cli
make                # builds libyame.a from the submodule, then ./sesame
```

`make` produces a single binary, `./sesame`. Put it on your `PATH` (e.g.
`cp sesame ~/.local/bin/`). At **runtime** it needs only zlib and libcurl
(YAME/htslib are static); the `yame` tool does not need to be installed.

Other build targets are covered under [Development](#development).

## Quickstart

```sh
# 1. Fetch a platform's ordering table + mask into the local store (one time).
sesame fetch MSA

# 2. Raw betas (no preprocessing) — platform auto-detected from the IDAT.
sesame betas path/to/sample > sample.betas.tsv      # reads sample_{Grn,Red}.idat[.gz]

# 3. The full default pipeline, same as openSesame(prep="QCDPB").
sesame betas --prep QCDPB path/to/sample > sample.betas.tsv

# 4. A whole cohort in one parallel process -> a Probe_ID x sample matrix.
#    A prefix is a path stem; derive them from the Grn files if needed.
sesame betas --prep QCDPB $(ls idats/*_Grn.idat.gz | sed 's/_Grn.idat.gz//') > cohort.betas.tsv

# 5. Per-sample QC metrics (detection success rate + the sesameQC panel).
sesame qc path/to/sample                            # one sample -> one TSV row
```

Output row `Probe_ID<TAB>beta`, `NA` for missing/masked probes. A `<prefix>`
resolves `<prefix>_Grn.idat[.gz]` and `<prefix>_Red.idat[.gz]`.

## Commands

```
sesame betas      [options] <prefix> [<prefix> ...]   # IDAT -> betas
sesame intensity  [options] <prefix> [<prefix> ...]   # IDAT -> total intensity (M+U); --cg for YAME
sesame qc         [options] <prefix> [<prefix> ...]   # IDAT -> QC metrics (TSV)
sesame dml        --betas <matrix> (--formula .. --meta .. | --design ..)  # diff. methylation
sesame fetch      [<platform>] [--force]              # download data into the store
sesame index-info                                     # show store + pinned tag + platforms
sesame idat-dump  [--head N] [--tsv] <file.idat[.gz]> # inspect a raw IDAT
sesame version
```

### `sesame betas`

Compute beta values from an IDAT pair. Equivalent to
`openSesame(prefix, prep=CODE, func=getBetas)`. One prefix prints a two-column
table; multiple prefixes print a matrix (see [Batch mode](#batch-mode)).

| flag | meaning |
|---|---|
| `--prep CODE` | preprocessing steps in order, e.g. `QCDPB` (default: none → raw betas) |
| `--index <path>` | use this ordering table directly, bypassing the store |
| `--platform P` | force the platform (`EPIC`/`EPICv2`/`HM450`/`MSA`) instead of bead-count detection |
| `--min-beads N` | mask probes with any bead count `< N` (default: off, matching R's `NULL`) |
| `--no-mask` | ignore the mask column; emit every beta |
| `--threads N`, `-t N` | worker threads for a batch (default: online CPUs) |
| `--f64` | write raw little-endian float64 instead of text (see [Output formats](#output-formats)) |
| `--dump-col` | emit `Probe_ID<TAB>col` (G/R/2) instead of betas — single prefix, for testing `C` |

`Q`, `P`, and `B` need the platform's `.cm` mask in the store; run
`sesame fetch <platform>` first.

### `sesame intensity`

Total signal intensity (M+U) per probe — R's `totalIntensities` — as a matrix,
same shape as `betas` and batch-parallel. It's the **CNV signal input** (and the
tool that generates the copy-number normal reference from normal IDATs).

`--cg <out.cg>` writes a **YAME `.cg`** (+ `<out.cg>.idx` of sample names) instead
of the TSV, so the output feeds straight into the `yame` toolchain. The default
is **format 3** (M and U counts) — from which yame derives both the beta
(`MU2beta`) and the coverage/total (`MU2cov`), so it subsumes the total-intensity
number; M/U round to integers (exact for raw IDAT signal). `--f4` instead writes
**format 4** (the total intensity as an exact float per probe).

```sh
sesame intensity sample                       # Probe_ID <TAB> M+U
sesame intensity --cg cohort.cg s1 s2 s3      # YAME format-3 (M/U), one block per sample
sesame intensity --cg cohort.cg --f4 sample   # YAME format-4 (total intensity float)
```

### `sesame qc`

Per-sample QC metrics — the `sesameQC_calcStats` panel — as a TSV: one row per
sample, one column per metric, headline being **detection success rate**
(`frac_dt`, the fraction of probes with pOOBAH detection p ≤ 0.05). Computed from
the raw signal; the detection and beta groups run pOOBAH internally, so the
platform's `.cm` mask must be in the store. Batch-parallel like `betas`; a failed
sample becomes an all-`NA` row. Flags: `--index`, `--platform`, `--min-beads`,
`--threads`.

```sh
sesame qc --threads 8 $(ls idats/*_Grn.idat.gz | sed 's/_Grn.idat.gz//') > cohort.qc.tsv
```

The panel covers detection, probe counts, signal intensity (in/out-of-band),
Infinium-I channel switches, dye bias (`RGdistort`), and the beta distribution —
mirroring R's `sesameQC_calcStats`. See `NUMERICS.md` for the two noted behaviors
(`num_dtna` under the D2 fix; the beta group computed after `D→B→P`).

### `sesame dml`

Per-probe **differential methylation**: for each probe, an OLS of its betas
across samples on a design, giving per-coefficient estimates and t-tests, a
holdout F-test per categorical variable, effect sizes, and BH-adjusted p-values —
sesame's `DML` / `summaryExtractTest`, as a TSV (one row per probe). The betas
matrix is `sesame betas` batch output; `--meta` is a TSV whose first column
matches the sample names.

```sh
# betas matrix from a cohort, then test each CpG against a group (+ covariate)
sesame betas --prep QCDPB $(ls idats/*_Grn.idat.gz|sed 's/_Grn.idat.gz//') > betas.tsv
sesame dml --betas betas.tsv --meta samples.tsv --formula '~ group + age' > dml.tsv
```

`--formula` takes **main effects** only (a `~` and metadata column names joined by
`+`): categorical columns are auto-dummied to match R's `model.matrix` (treatment
contrasts, alphabetical levels), continuous columns enter as-is, with an
intercept. For interactions/splines, build the design elsewhere and pass it with
`--design <numeric.tsv>` (first column = sample id, remaining columns = terms).
Runs in parallel (`--threads`). Because DML consumes the betas matrix, it has no
data-lineage caveat — it matches R's `lm` to ~1e-9 (see `NUMERICS.md`). Region
calling (DMR) is not yet included: it needs per-probe genomic coordinates.

### `sesame fetch`

Download a platform's data (ordering table + `.cm` mask + `SHA256SUMS`) into the
[store](#data-files-and-the-store), verifying every file against a digest
compiled into the build. With no platform, fetches all published platforms.
`--force` re-downloads even if a file is already present and matches. **This is
the only command that touches the network**, and it never prompts.

### `sesame index-info`

Print the store location, which annotation tag this build pins, and which
platforms are present locally versus still need fetching.

### `sesame idat-dump`

Inspect a single raw `.idat`/`.idat.gz`. Default prints a summary header; `--tsv`
emits `addr<TAB>mean<TAB>sd<TAB>nbeads` with no header; `--head N` limits rows.

## The QCDPB pipeline

`--prep` applies steps in the **order given**; `QCDPB` is the `openSesame`
default. Each is an independent, self-contained port of the R function
(referenced in `NUMERICS.md`).

| code | step | what it does |
|---|---|---|
| `Q` | `qualityMask` | OR the platform's recommended probe mask (from the `.cm`) into the SigDF mask |
| `C` | `inferInfiniumIChannel` | reassign each Infinium-I probe to its brighter color channel |
| `D` | `dyeBiasNL` | nonlinear dye-bias correction — pull the red and green Inf-I distributions to a common curve |
| `P` | `pOOBAH` | detection-p masking from the out-of-band + negative-control background; mask probes with `p > 0.05` |
| `B` | `noob` | normal-exponential background subtraction |

Pass any subset in any order (e.g. `--prep CD` for just channel + dye bias, or
`--prep ""`/omit for raw betas). `Q`, `P`, and `B` read the `.cm` mask, so the
platform must be present in the store.

## Batch mode

Passing **multiple prefixes** runs them in one process. The index and masks are
parsed **once** and shared read-only; the independent per-sample work runs across
a pthread pool (`--threads`, default = online CPUs).

- **All prefixes must be the same platform** (the row space of the matrix).
- Output is a `Probe_ID` matrix, one column per sample named by basename — or a
  sample-major stream with `--f64`.
- **Deterministic:** a sample's column is its position on the command line, so
  output order never depends on scheduling, and `--threads 1` is byte-identical
  to `--threads N`.
- **Robust:** a sample that fails (missing/corrupt IDAT, wrong platform) becomes
  an all-`NA` column with a warning and sets [exit status](#exit-status) 1 — one
  bad IDAT never aborts the run.
- **Scales past RAM:** the result store is an unlinked temp file that the OS pages
  to disk, so a matrix larger than memory still works. For whole-cohort scale,
  shard the sample list across invocations (pipelines already do this).

```sh
sesame betas --prep QCDPB --threads 8 run/S1 run/S2 run/S3 ... > cohort.betas.tsv
```

On three EPICv2 samples this is ~3× faster than three separate invocations (index
parsed once, work parallelized); the gap widens with cohort size and core count.

## Output formats

- **Single prefix, text** (default): `Probe_ID<TAB>beta` per line, `NA` for
  missing/masked probes.
- **Multiple prefixes, text**: a header `Probe_ID<TAB><name1><TAB>…` (names are
  prefix basenames) followed by one row per probe.
- **`--f64`**: raw little-endian `double` (NA as `NaN`), no IDs — one block of
  `nprobes` values per sample, in command-line order. Probe order is the ordering
  table; sample order is argv.

Use `--f64` for lossless comparison. Do **not** diff betas via text: R's parser
does not correctly round 17-digit decimals, which manufactures a phantom ~1e-16
disagreement (see `NUMERICS.md`).

## Data files and the store

Ordering tables and masks live in the annotation repo
[`zhou-lab/InfiniumAnnotation`](https://github.com/zhou-lab/InfiniumAnnotation),
one subfolder per platform, versioned by **git tag**:

```
https://github.com/zhou-lab/InfiniumAnnotation/raw/<tag>/<platform>/<file>
                                                   ^^^^^ the version    e.g. v1/MSA/MSA.ordering.tsv.gz
```

Each platform folder holds its ordering table, its `.cm` mask (+ `.idx`), and a
`SHA256SUMS`. This build pins tag **v1**. *(Only MSA is published at v1 so far;
the others error cleanly until published, and still auto-detect by bead count.)*

`sesame fetch <platform>` downloads that platform's whole folder into the local
**store**, mirroring the remote exactly:

```
<store>/MSA/SHA256SUMS            byte-identical copy of the remote
<store>/MSA/MSA.ordering.tsv.gz
<store>/MSA/MSA.hg38.mask.cm(.idx)
```

so `cd <store>/MSA && shasum -a 256 -c SHA256SUMS` verifies the store by hand.
Fetch first pulls `SHA256SUMS`, verifies it against a digest compiled into this
build (a hard trust anchor), then verifies every file against it; a file already
present with the right digest is skipped.

**Store location** (`sesame fetch` writes here) — resolved in this order:

| | |
|---|---|
| `$SESAME_INDEX_DIR` | explicit |
| `<dir of the binary>/data` | a checkout — found from **any** working directory |
| `$XDG_CACHE_HOME/sesame`, `~/Library/Caches/sesame`, `~/.cache/sesame` | fallback |

The binary-relative default is why a source checkout finds `data/` from any cwd,
while an installed binary (no `data/` beside it) falls through to the XDG store.

**Lookup order for an ordering table**: `--index <path>` >
`<store>/<platform>/<platform>.ordering.tsv.gz` > `./<platform>.ordering.tsv.gz`.

**sesame never downloads implicitly and never prompts** — a prompt would hang
forever, silently, in a Nextflow job or a Docker build. `sesame fetch` is the
only path that touches the network, so behaviour is identical with or without a
TTY. If an index is missing, you get an error naming the exact command to run.

## Configuration

sesame-cli is configured by flags and a few environment variables; there is no
config file.

| variable | effect |
|---|---|
| `SESAME_INDEX_DIR` | the store — where `fetch` writes and where indices/masks are looked up |
| `XDG_CACHE_HOME` | fallback store root when `SESAME_INDEX_DIR` is unset and no `data/` sits beside the binary |
| `TMPDIR` | where a batch's temp-file-backed result matrix is created (falls back to `/tmp`) |
| `NO_COLOR` | disables ANSI color in `index-info` (also auto-off when stdout is not a TTY) |
| `SESAME_TEST_IDATS` | test IDAT directory for `make test` (default `~/repo/InfiniumTestIDATs`) |

## Exit status

| code | meaning |
|---|---|
| `0` | success |
| `1` | a runtime error, or (in a batch) at least one sample failed — the matrix still prints, with `NA` columns for the failures |
| `2` | usage error (bad flags/arguments) |

## Fidelity and validation

R is the oracle, permanently. Betas are compared as raw float64, never text. The
golden ladder (`make test`):

| level | gate | status |
|---|---|---|
| 1. IDAT reader | `IlluminaID`/`Mean`/`SD`/`NBeads` **bit-identical**, no tolerance | ✅ 33/33 files, ~20.4M records, 9 platforms, plain + gz |
| 3. Per-step `C` | channel calls identical | ✅ 5/5 samples |
| 3. Per-step `D` | max relative Δ ≤ 1e-12 (D8 clean-room qnorm) | ✅ ~2 ULP |
| 3. Per-step `Q` | masks exactly the recommended-track union of the `.cm` | ✅ self-consistent (0.982 Jaccard vs R — mask lineage, `NUMERICS.md`) |
| 3. Per-step `P` | every C-vs-R disagreement is 0.05-boundary or D2 | ✅ R-only 0; all 12 flips at the cutoff (`NUMERICS.md`) |
| 3. Per-step `B` | `normExpSignal`/`huber` vs R on identical inputs; betas on raw-identical probes | ✅ arithmetic ≤ few ULP; betas median 6.5e-6, max 1.6e-3 (lineage) |
| 4. Betas `prep=""` | **bit-identical** | ✅ 6/6 samples, 4.4M betas |
| 5. Batch | each column == its single-sample run; `--threads 1 == N`; bad sample → NA column | ✅ byte-identical; ThreadSanitizer-clean |
| 6. QC panel | every `sesameQC_calcStats` metric within lineage scale | ✅ 65/65 metrics; worst 1.75e-2 (small count, `NUMERICS.md`) |
| 7. DML | vs R `DML`/`summaryExtractTest` (no lineage) | ✅ Est/Pval/FPval/Eff/BH match to ~5e-10 |
| 8. .cg output | `intensity --cg` round-trips through `yame unpack` | ✅ format 3 (M/U) + format 4; names + values match |

"Lineage" means the published ordering/mask is a newer data version than the
installed `sesameData`, so a handful of probes differ for reasons that predate any
arithmetic. Every such case is quantified and explained in `NUMERICS.md`, which
carries one row per intentional difference from R (the `D1`–`D8` register).

The IDAT parser eats untrusted files from GEO — `nFields`, per-field `byteOffset`,
and `nSNPsRead` are all attacker-controlled, and the field table is a 10-byte
*unpadded* record. It is fuzzed and run under ASan/UBSan.

## Development

```sh
make                # build ./sesame
make test           # the full golden ladder vs the R oracle
make asan           # rebuild sesame's own objects under ASan + UBSan
make fuzz-replay    # replay the IDAT-parser corpus under ASan/UBSan
make index          # regenerate testdata/*.ordering.tsv.gz from sesameData (needs R)
make clean
```

`make test` needs `Rscript` with the sesame R package installed (the differential
oracle) and test IDATs at `$SESAME_TEST_IDATS` (default
`~/repo/InfiniumTestIDATs`); the `Q` sub-test also wants the `yame` binary on
`PATH` for its oracle and skips cleanly without it. Individual levels:
`make test-idat test-betas test-prep test-qmask test-poobah test-noob test-batch`.

Sanitizers compose via additive flags, so you can target anything, e.g. a
ThreadSanitizer build of the batch path:

```sh
make clean
make EXTRA_CFLAGS="-fsanitize=thread -g -O1" EXTRA_LDFLAGS="-fsanitize=thread"
```

Layout: `cli/` (the `sesame` command), `src/` (the `libsesame` core — IDAT
reader, index, SigDF, prep steps, numerics, fetch/cache), `include/sesame.h` (the
public API), `tests/` (shell drivers + R comparison oracles), `YAME/` (submodule).
`NUMERICS.md` documents every intentional divergence from R and is required
reading before changing a prep step.

## License

**AGPL-3.0-or-later.** See `LICENSE`. This matches the wider Zhou Lab toolchain
(YAME is also AGPL-3.0), which is what lets sesame-cli link those tools directly —
e.g. linking YAME (a git submodule) to read `.cm` mask files in-process — without
a licence conflict.

sesame (the R package) is MIT and is a *separate program*; sesame-cli linking or
invoking it, or vice versa, is unaffected by this choice.

The quantile-normalization step (`dyeBiasNL`) is a clean-room reimplementation
rather than a use of `preprocessCore` (`LGPL >= 2`). Under AGPL we *could* now
vendor `qnorm.c` (LGPL is upward-compatible with GPL/AGPL) for bit-exactness, but
the clean-room already agrees to ~2 ULP and avoids the extra dependency, so it
stays. See `NUMERICS.md` (D8).
