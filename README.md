# sesame-cli

`sesame` — a standalone C implementation of Infinium preprocessing.

A standalone C implementation of [sesame](https://github.com/zwdzwd/sesame)'s
basic Infinium DNA methylation preprocessing: **IDAT → betas**, with no R and no
Bioconductor. Builds one binary, `sesame`.

> Status: **P0 done, P1 partial.** `sesame betas` produces beta values that are
> **bit-identical to R** with `prep=""`. Preprocessing (QCDPB) is NOT yet
> implemented -- that is P2-P4.

## Why

`openSesame` needs R + Bioconductor + `sesameData` (ExperimentHub) +
`preprocessCore`/`MASS`/`BiocParallel`. That blocks:

- **Pipeline deployment** — a static binary for Nextflow/Snakemake/cloud.
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

Only dependency is zlib.

```sh
make            # build ./sesame
make test       # golden tests vs the R oracle (needs Rscript + sesame)
                # uses $SESAME_TEST_IDATS (default ~/repo/InfiniumTestIDATs)
                # for real full-size arrays, in addition to sesameData extdata
make asan       # build with ASan/UBSan
make fuzz       # libFuzzer target (Linux/clang; Apple clang has no libFuzzer)
make fuzz-replay  # corpus replayer under ASan/UBSan (works everywhere)
```

## Index files

Ordering tables are published as **per-platform release assets**:

```
https://github.com/zwdzwd/sesame-cli/releases/download/EPICv2.ordering.v1/EPICv2.ordering.tsv.gz
                                                       ^^^^^^^^^^^^^^^ tag carries the version
```

The filename is stable across versions — the tag is the version, so updating one
platform means one new tag and one registry row; the others are untouched and are
not re-downloaded. Assets are immutable: new data always means a new tag, never a
re-upload under an existing one.

Retrieval is `(tag, file)`, and the cache mirrors the URL so pinned versions coexist:

```
~/.cache/sesame/EPICv2.ordering.v1/EPICv2.ordering.tsv.gz
~/.cache/sesame/EPICv2.ordering.v2/EPICv2.ordering.tsv.gz
```

Resolution order: `--index <path>` > `$SESAME_INDEX_DIR` > `./` > `./data` > cache
(`$SESAME_CACHE` | `$XDG_CACHE_HOME/sesame` | `~/Library/Caches/sesame` | `~/.cache/sesame`).

**sesame never downloads implicitly and never prompts** — a prompt would hang
forever, silently, in a Nextflow job or a Docker build. `sesame fetch` is the only
path that touches the network, so behaviour is identical with or without a TTY.
Downloads are verified against a pinned sha256; a mismatch is fatal.

## Use

```sh
sesame idat-dump sample_Grn.idat            # summary
sesame idat-dump --tsv sample_Grn.idat.gz   # addr<TAB>mean<TAB>sd<TAB>nbeads

sesame index-info                           # cache dir + which platforms are present
sesame fetch EPICv2                         # download the pinned index
sesame fetch --tag EPICv2.ordering.v2 EPICv2.ordering.tsv.gz   # explicit (tag, file)

# Betas (no preprocessing yet -- equivalent to openSesame(prefix, prep=""))
sesame betas <prefix>                       # platform auto-detected from bead count
sesame betas --index data/EPICv2.ordering.tsv.gz <prefix>
sesame betas --f64 <prefix> > betas.f64     # lossless
```

Both plain `.idat` and gzipped `.idat.gz` are read through the same path.

Use `--f64` (raw little-endian float64) for lossless output. Do not compare
betas via text: R's parser does not correctly round 17-digit decimals, which
manufactures a phantom ~1e-16 disagreement. See `NUMERICS.md`.

## Validation

R is the oracle forever, not for one release. The golden ladder:

| level | gate |
|---|---|
| 1. IDAT reader | `IlluminaID`/`Mean`/`SD`/`NBeads` **bit-identical**. No tolerance. ✅ passing — 33/33 files, ~20.4M bead records, 9 platforms (EPIC, EPIC+, EPICv2, HM27, HM450, Mammal40, MM285, MSA), plain and gzipped |
| 2. Index | exact set + order equality vs `sesameAnno_buildAddressFile()` / `getMask()` / `backgroundMask()` — ✅ ordering table exported from R and verified via the beta gate below |
| 3. Per-step `Q C D P B` | applied independently against a fixed SigDF |
| 4. End-to-end betas (`prep=""`) | ✅ passing — **bit-identical**, 6/6 samples, 4.4M betas across HM450/EPIC/EPICv2/MSA |
| 4b. End-to-end betas (`QCDPB`) | blocked on P2-P4 |

The IDAT parser eats untrusted files from GEO — `nFields`, per-field
`byteOffset`, and `nSNPsRead` are all attacker-controlled, and the field table
is a 10-byte *unpadded* record. It is fuzzed and run under ASan/UBSan.

## License

MIT, matching sesame. Note `preprocessCore` is LGPL-2 — its `qnorm.c` must not
be copied; the quantile-normalization replacement is clean-room.
