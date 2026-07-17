# NUMERICS.md — divergence register

sesame is an independent implementation of sesame's basic preprocessing. It is
**not** bit-exact with the R package by mandate: where the R code has a defect or
a numerically fragile construction, sesame fixes it and records the difference
here.

The R package (`zwdzwd/sesame`) remains the reference oracle and is never
modified by this project. Every difference below must be *justified*, not merely
reported — for numerical reformulations that means a proof or an
arbitrary-precision oracle showing sesame is closer to truth.

| id | site | R behavior | sesame behavior | status |
|----|------|-----------|------------------|--------|
| D1 | `R/sesame.R:293-298` `readIDAT1` | Merges Grn/Red by **positional `cbind`**, rownames taken from Grn. No join on IlluminaID. The assumption *holds for well-formed files* (verified: Grn/Red address vectors are identical, sorted and unique on every test array). It fails silently and array-wide only when mismatched IDATs are paired. | Verifies the address vectors match and errors otherwise. Confirmed: pairing an EPICv2 Grn with an HM450 Red is rejected, where R would silently emit garbage. | **done (P1)** |
| D2 | `R/detection.R:163-165` `pOOBAH` | `pmax(..., na.rm=TRUE)` but `pmin(...)` **without** `na.rm` → an NA in one channel yields `p=NA`, and `NA > 0.05` is NA, so the probe is **not** masked. Sibling `ELBAR` at `R/detection.R:60` does `pvals[is.na(pvals)] <- 1.0`, evidence this is an oversight. | Both channels NA → `p=1` (mask). One channel NA → use the other. | planned (P4) |
| D3 | `R/readIDAT.R:22` vs `:178` | Header parsed **twice**: `readIDAT()` reads `version` as int32 to dispatch, then `readIDAT_nonenc()` reopens the file and re-reads it as int64. Works only because the value is small and little-endian. | Parsed once, as int64. | **done (P0)** |
| D4 | `R/readIDAT.R:44-46` `readLong` | `readBin(what="integer", size=8)` — base R does not support 8-byte integers; this works by accident. | Real `int64_t`. | **done (P0)** |
| D5 | `R/background.R:143-167` `normExpSignal` | Computes `exp(dnorm(...,log=TRUE) - pnorm(...,lower.tail=FALSE,log.p=TRUE))`; a difference of logs, which loses precision in the left tail where the `1e-6` floor fires. | Evaluates the inverse Mills ratio `λ(t)=φ(t)/Φ(t)` directly (`erfc` for `t > -5`, Laplace continued fraction below). Mathematically identical, better conditioned, never differences logs. | planned (P4) |
| D6 | `MASS::huber` via `R/background.R:132,135` | **Errors** when `MAD == 0`. | Falls back to `s = max(mad, IQR/1.349, 1e-8)` and sets a status bit — a single degenerate array must not abort a 100k-sample run. | planned (P3) |
| D7 | `R/dye_bias.R:122`, `R/background.R:98`, `R/detection.R:155`, `R/mask.R:173` | Escape hatches are **silent**: `dyeBiasNL`→`maskIG` when RGdistort is NA or >10; `noob` returns unchanged if background <100; `pOOBAH` substitutes a `1:1000` prior; `backgroundMask` returns NULL for unknown platforms (silently changing both `noob` and `pOOBAH` background composition). | Each sets a bit in a per-sample status word exposed by the API, so pipelines can count fallbacks instead of discovering them in the betas. Done for `dyeBiasNL` (`SESAME_STAT_DYEBIAS_FAILED`); rest planned. | partial |
| D8 | `preprocessCore::normalize.quantiles.use.target` via `R/dye_bias.R:135,150` | preprocessCore's own compiled C (`qnorm.c`). | Clean-room reimplementation — see below. Algebraically identical; **not bit-identical**. | **done (D)** |

### D8 — the one difference I could not eliminate, and why

`normalize.quantiles.use.target` is the only compiled third-party dependency in
the QCDPB path. preprocessCore is `LGPL (>= 2)`; sesame-cli is AGPL-3.0-or-later.

**Since sesame-cli went AGPL (2026-07-17), vendoring `qnorm.c` is now legally
open:** LGPL is upward-compatible with GPL/AGPL (a recipient may elect LGPL-3,
which is GPL-3 plus permissions and thus absorbable into an AGPL-3 work), and
AGPL already requires shipping complete source, so LGPL's relink obligation is
moot. So this entry could be *erased* by vendoring preprocessCore's `qnorm.c`
for bit-exactness.

**It is deliberately kept clean-room anyway**, for engineering reasons rather
than licensing ones:

- the clean-room already agrees with R to ~2 ULP (measured below) — biologically
  and numerically meaningless;
- vendoring adds an LGPL component and pulls in preprocessCore's `include/` +
  stubs, which are `LinkingTo:` plumbing for R packages, not a clean library a
  standalone C program links against — i.e. real extraction work;
- the clean-room is ~40 lines with no dependency.

If exact reproduction of R's last bit ever matters more than that, vendoring is
now the sanctioned route. Until then the clean-room stands. It was characterized
purely by black-box probing:

| probe | observed |
|---|---|
| `n == m`, no ties | `sorted(target)[rank(x)]` |
| `n != m` | `quantile(target, (rank-1)/(n-1), type=7)` |
| ties in `x` | quantile at the **average** rank |
| NA in target | dropped before ranking |
| NA in `x` | preserved; ranking uses non-NA only |

The tie rule was pinned with a deliberately non-linear target, which separates
"quantile at the average rank" (actual: `100`) from "mean of the quantiles at
each tied rank" (`106.67`). The former wins.

**The algebra is settled; the last bit is not.** I could not find an arithmetic
form that reproduces `qnorm.c` bit-for-bit. Both `(1-h)*lo + h*hi` and
`lo + h*(hi-lo)`, across three index formulations, top out at ~92% bitwise
agreement — so preprocessCore evaluates it in some third order I would have to
read the source to learn, which is precisely what the licence forbids.

Measured consequence on end-to-end betas (`prep="CD"`, NA patterns identical
everywhere):

| platform | betas differing | max abs | max rel |
|---|---|---|---|
| HM450 | 1.05% | 5.6e-16 | 2.5e-15 |
| EPIC | 1.47% | 1.0e-15 | 4.0e-15 |
| EPICv2 | 1.82% | 7.8e-16 | 3.0e-15 |
| MSA | 4.39% | 3.3e-16 | 1.1e-15 |

~2 ULP on a value in [0,1]. The acceptance gate for `D` is max relative
difference ≤ 1e-12, which this clears by three orders of magnitude — and the
project's own stated gate was 1e-8. This is accepted as last-bit arithmetic
ordering, not an algorithmic difference: it is not attributable to any decision
in this implementation, and neither result is "more correct" than the other.

## Observations about the R implementation (not divergences)

**Ordering files are sorted by Probe_ID under R's *locale* collation, not byte
order.** All four platforms are in `order()` sequence under `en_US.UTF-8`.
HM450/EPIC/EPICv2 happen to coincide with byte order; **MSA does not**. It is the
only platform carrying `rs` IDs where one number prefixes another, e.g.

```
locale : rs12051_TC21 < rs12051548_TC21     <- what the files use
radix  : rs12051548_TC21 < rs12051_TC21
```

ICU ignores the `_` and compares `rs12051TC21` vs `rs12051548TC21` (`5` < `T`);
byte order compares `_` (0x5F) vs `5` (0x35) and flips it. Exactly 11 such pairs
in MSA.

This does not affect sesame-cli: Probe_IDs are consumed in file order and the
only binary search is over integer addresses. Two hazards it *does* create:

1. **Rebuild reproducibility.** `order()` is locale-dependent, so regenerating
   the `.address` object under `LC_COLLATE=C` (Docker, CI, minimal images) would
   reorder MSA, change its sha256, and change the output row order — which is
   load-bearing (`R/sesame.R:504`). `tools/export_ordering.R` preserves the
   object's existing order and is therefore safe; anyone rebuilding *upstream*
   must pin the locale.
2. **Never binary-search Probe_IDs with `strcmp`.** It would silently fail on
   those 11 MSA probes. Use a hash if ID lookup is ever needed.

**The ordering table's `mask` column is dead data.** `.address$ordering` carries
a `mask` column built by `sesameAnno_buildAddressFile` from
`create_default_mask()$ref_issue` (`R/sesameAnno.R:170`). Nothing ever reads it:
`chipAddressToSignal` explicitly sets `mask=FALSE` (`R/sesame.R:475,498`) and a
grep over `R/` finds no consumer anywhere in the pipeline. sesame parses the
column but deliberately ignores it, matching R. Seeding the SigDF mask from it
would wrongly mask ~257 (HM450) / ~2263 (EPICv2) / ~3120 (MSA) probes that R
keeps. Design masking is `qualityMask`'s job (the `Q` step), which uses a
different, curated vocabulary. *This is worth a look from the sesame side: the
column is either vestigial or an intent that was never wired up.*

**R's text parser is not correctly rounded — do not compare betas via text.**
`as.numeric("0.96236179722418258")` returns a double one ULP away from the value
that produced that string, even though C's own `strtod` round-trips it exactly.
R uses its own `R_strtod` rather than a correctly-rounded parser. Comparing
through `read.table` manufactured a ~2.2e-16 "divergence" across ~19% of probes
that did not exist in the computation: with raw float64 the betas are
**bit-identical**. The golden harness therefore compares binary (`--f64`), never
text.

## Deliberately *not* divergent

- **Row order.** The SigDF row order follows the manifest ordering table
  (`R/sesame.R:504`). This is load-bearing and is preserved exactly.
- **`approx(ties=mean)`** in `dyeBiasNL` (`R/dye_bias.R:143,158`). The tie
  collapse is load-bearing — intensities are integer-valued with massive
  duplication and tied `x` map to different `y`. Replicated exactly.
- **`quantile` type 7** at `R/channel_inference.R:32` and inside `IQR` at
  `R/background.R:103-104`. Replicated exactly.
- **`getBetas` clamps** `pmax(M,1)/pmax(M+U,2)` — note sesame uses **no** `+100`
  Illumina offset. Replicated exactly.
