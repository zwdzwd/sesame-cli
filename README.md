# sesame-cli

A standalone C implementation of [sesame](https://github.com/zwdzwd/sesame)'s
basic Infinium DNA methylation preprocessing: **IDAT → betas**, with no R and no
Bioconductor. Builds one binary, `sesame`.

> **Status.** The IDAT reader and `betas` with no preprocessing are
> **bit-identical to R**, and the full **`QCDPB`** pipeline is implemented: `Q`
> (self-consistent vs R by mask lineage), `C` (channel-identical), `D` (within
> ~2 ULP), `P` (detection-p exact; residual is 0.05-boundary/lineage), and `B`
> (noob — `normExpSignal`/`huber` proven to a few ULP in isolation; residual is
> lineage only). See the validation table and `NUMERICS.md`.

## Why

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

**In:** `readIDATpair` → `prepSesame(sdf, "QCDPB")` → `getBetas`, for EPIC,
EPICv2, HM450, MSA.
**Out:** DML, KYCG, CNV, visualization, and all inference (species/strain/age/
sex/ethnicity) — those stay in R.

## Build

Dependencies: a C compiler + `make`, **zlib**, **libcurl**, and the **YAME**
submodule (bundled htslib) — linked in to read `.cm` mask files. Both sesame-cli
and YAME are AGPL-3.0, so linking is clean.

```sh
git clone --recurse-submodules <repo>          # or: git submodule update --init --recursive
make                # builds libyame.a from the submodule, then ./sesame
make test           # golden tests vs the R oracle (needs Rscript + the sesame R package;
                    # the Q test also needs the `yame` binary on PATH)
                    # uses $SESAME_TEST_IDATS (default ~/repo/InfiniumTestIDATs)
make asan           # sesame's own objects under ASan/UBSan
make fuzz-replay    # IDAT-parser corpus replayer under ASan/UBSan
```

At **runtime** the binary needs only zlib and libcurl (YAME/htslib are static);
`yame` does not need to be installed.

## Data files

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

**Store location** (`sesame fetch` writes here) — one variable:

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

## Use

```sh
sesame idat-dump sample_Grn.idat            # summary of an IDAT
sesame idat-dump --tsv sample_Grn.idat.gz   # addr<TAB>mean<TAB>sd<TAB>nbeads

sesame fetch MSA                            # download MSA's ordering + mask into the store
sesame fetch                                # all published platforms
sesame index-info                           # store location, pinned tag, which platforms present

# Betas. Platform is auto-detected from the IDAT bead count.
sesame betas <prefix>                       # prep="" — equivalent to openSesame(prep="")
sesame betas --prep QCDPB <prefix>          # full pipeline: Q, C, D, P, B in order
sesame betas --index <ordering.tsv.gz> <prefix>   # bypass the store
sesame betas --f64 <prefix> > betas.f64     # raw float64 (lossless)
```

`<prefix>` resolves `<prefix>_Grn.idat[.gz]` and `<prefix>_Red.idat[.gz]`. `Q`,
`P`, and `B` require the platform's `.cm` mask in the store
(`sesame fetch <platform>`).

Use `--f64` (raw little-endian float64) for lossless output. Do **not** compare
betas via text: R's parser does not correctly round 17-digit decimals, which
manufactures a phantom ~1e-16 disagreement (see `NUMERICS.md`).

## Validation

R is the oracle, permanently. The golden ladder (`make test`):

| level | gate | status |
|---|---|---|
| 1. IDAT reader | `IlluminaID`/`Mean`/`SD`/`NBeads` **bit-identical**, no tolerance | ✅ 33/33 files, ~20.4M records, 9 platforms, plain + gz |
| 3. Per-step `C` | channel calls identical | ✅ 5/5 samples |
| 3. Per-step `D` | max relative Δ ≤ 1e-12 (D8 clean-room qnorm) | ✅ ~2 ULP |
| 3. Per-step `Q` | masks exactly the recommended-track union of the `.cm` | ✅ self-consistent (0.982 Jaccard vs R — mask lineage, `NUMERICS.md`) |
| 3. Per-step `P` | every C-vs-R disagreement is 0.05-boundary or D2 | ✅ R-only 0; all 12 flips at the cutoff (`NUMERICS.md`) |
| 3. Per-step `B` | `normExpSignal`/`huber` vs R on identical inputs; betas on raw-identical probes | ✅ arithmetic ≤ few ULP; betas median 6.5e-6, max 1.6e-3 (lineage) |
| 4. Betas `prep=""` | **bit-identical** | ✅ 6/6 samples, 4.4M betas |

The IDAT parser eats untrusted files from GEO — `nFields`, per-field
`byteOffset`, and `nSNPsRead` are all attacker-controlled, and the field table is
a 10-byte *unpadded* record. It is fuzzed and run under ASan/UBSan.

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
