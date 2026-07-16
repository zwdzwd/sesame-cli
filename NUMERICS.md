# NUMERICS.md — divergence register

sesamec is an independent implementation of sesame's basic preprocessing. It is
**not** bit-exact with the R package by mandate: where the R code has a defect or
a numerically fragile construction, sesamec fixes it and records the difference
here.

The R package (`zwdzwd/sesame`) remains the reference oracle and is never
modified by this project. Every difference below must be *justified*, not merely
reported — for numerical reformulations that means a proof or an
arbitrary-precision oracle showing sesamec is closer to truth.

| id | site | R behavior | sesamec behavior | status |
|----|------|-----------|------------------|--------|
| D1 | `R/sesame.R:293-298` `readIDAT1` | Merges Grn/Red by **positional `cbind`**, rownames taken from Grn. No join on IlluminaID. If the two files' address sets differ or are ordered differently, R silently produces wrong betas array-wide. | Join on `addr`. Mismatch is an error under strict mode, otherwise NA + a status bit. | planned (P1) |
| D2 | `R/detection.R:163-165` `pOOBAH` | `pmax(..., na.rm=TRUE)` but `pmin(...)` **without** `na.rm` → an NA in one channel yields `p=NA`, and `NA > 0.05` is NA, so the probe is **not** masked. Sibling `ELBAR` at `R/detection.R:60` does `pvals[is.na(pvals)] <- 1.0`, evidence this is an oversight. | Both channels NA → `p=1` (mask). One channel NA → use the other. | planned (P4) |
| D3 | `R/readIDAT.R:22` vs `:178` | Header parsed **twice**: `readIDAT()` reads `version` as int32 to dispatch, then `readIDAT_nonenc()` reopens the file and re-reads it as int64. Works only because the value is small and little-endian. | Parsed once, as int64. | **done (P0)** |
| D4 | `R/readIDAT.R:44-46` `readLong` | `readBin(what="integer", size=8)` — base R does not support 8-byte integers; this works by accident. | Real `int64_t`. | **done (P0)** |
| D5 | `R/background.R:143-167` `normExpSignal` | Computes `exp(dnorm(...,log=TRUE) - pnorm(...,lower.tail=FALSE,log.p=TRUE))`; a difference of logs, which loses precision in the left tail where the `1e-6` floor fires. | Evaluates the inverse Mills ratio `λ(t)=φ(t)/Φ(t)` directly (`erfc` for `t > -5`, Laplace continued fraction below). Mathematically identical, better conditioned, never differences logs. | planned (P4) |
| D6 | `MASS::huber` via `R/background.R:132,135` | **Errors** when `MAD == 0`. | Falls back to `s = max(mad, IQR/1.349, 1e-8)` and sets a status bit — a single degenerate array must not abort a 100k-sample run. | planned (P3) |
| D7 | `R/dye_bias.R:122`, `R/background.R:98`, `R/detection.R:155`, `R/mask.R:173` | Escape hatches are **silent**: `dyeBiasNL`→`maskIG` when RGdistort is NA or >10; `noob` returns unchanged if background <100; `pOOBAH` substitutes a `1:1000` prior; `backgroundMask` returns NULL for unknown platforms (silently changing both `noob` and `pOOBAH` background composition). | Each sets a bit in a per-sample status word exposed by the API, so pipelines can count fallbacks instead of discovering them in the betas. | planned (P2–P4) |

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
