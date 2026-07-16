/* prep.c -- preprocessing steps (the prepSesame codes).
 *
 * Implemented so far: C (inferInfiniumIChannel).
 * Outstanding: Q, D, P, B.
 *
 * SPDX-License-Identifier: MIT
 */
#include "sesame.h"
#include "internal.h"

#include <math.h>
#include <stdlib.h>

/* ---------------------------------------------------------------- C ---
 *
 * inferInfiniumIChannel (R/channel_inference.R:20-55).
 *
 * For every Infinium-I probe, take the larger of the two alleles in each
 * channel and assign the probe to whichever channel is brighter. The
 * "background" is then the pooled *out-of-band* signal under that new
 * assignment, and any probe whose brighter channel does not clear the 95th
 * percentile of it is deemed failed.
 *
 * Four R semantics that a naive port gets wrong, all verified against the
 * source:
 *
 *   1. `pmax(MR, UR)` has NO na.rm (:26-27) -- an NA allele makes the max NA.
 *      C's fmax() would instead ignore the NaN and return the other operand.
 *      sesame__pmax2 replicates R.
 *   2. `ifelse(red_max > grn_max, "R", "G")` (:28) yields **NA**, not "G", when
 *      the comparison is NA. Those probes are neither R nor G at this stage.
 *   3. `sdf1[new_col == "R",]` (:30-31) indexes with a logical containing NAs,
 *      which yields NA-*filled rows*; `quantile(..., na.rm=TRUE)` then drops
 *      them. Net effect: NA probes contribute nothing to the background pool.
 *      Not an accident to be tidied -- it is the behaviour.
 *   4. `red_max > grn_max` is strict, so exact ties go to **G**.
 *
 * With the defaults (switch_failed=0, mask_failed=0) failed probes revert to
 * the manifest channel and are not masked.
 */
int sesame_prep_infer_channel(sesame_sigdf_t *s, int switch_failed,
                              int mask_failed, sesame_err_t *err)
{
    int32_t n, i, n_inf1 = 0, n_pool = 0, n_ok;
    double *rmax = NULL, *gmax = NULL, *pool = NULL, *buf = NULL;
    uint8_t *newc = NULL;
    double bg_max;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;

    rmax = (double  *)malloc((size_t)n * sizeof(double));
    gmax = (double  *)malloc((size_t)n * sizeof(double));
    newc = (uint8_t *)malloc((size_t)n);
    /* pool holds at most 2 values per Inf-I probe */
    pool = (double  *)malloc((size_t)n * 2 * sizeof(double));
    buf  = (double  *)malloc((size_t)n * 2 * sizeof(double));
    if (!rmax || !gmax || !newc || !pool || !buf) {
        rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom in C step");
        goto out;
    }

    /* Pass 1: per-probe maxima and the tentative assignment.
     * newc is SESAME_COL_II for "NA" here -- Inf-I probes can never legitimately
     * be II, so it is an unambiguous sentinel for R's NA. */
    for (i = 0; i < n; i++) {
        if (s->col[i] == SESAME_COL_II) { newc[i] = SESAME_COL_II; continue; }
        n_inf1++;
        rmax[i] = sesame__pmax2(s->MR[i], s->UR[i]);
        gmax[i] = sesame__pmax2(s->MG[i], s->UG[i]);
        if (isnan(rmax[i]) || isnan(gmax[i]))
            newc[i] = SESAME_COL_II;               /* R's NA (semantic 2) */
        else
            newc[i] = (rmax[i] > gmax[i]) ? SESAME_COL_R  /* strict: ties -> G */
                                          : SESAME_COL_G;
    }
    if (n_inf1 == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no Infinium-I probes");
        goto out;
    }

    /* Pass 2: background = out-of-band signal under the NEW assignment.
     * Probes whose newc is the NA sentinel contribute nothing (semantic 3). */
    for (i = 0; i < n; i++) {
        if (newc[i] == SESAME_COL_R) {          /* inferred red -> green is oob */
            pool[n_pool++] = s->MG[i];
            pool[n_pool++] = s->UG[i];
        } else if (newc[i] == SESAME_COL_G && s->col[i] != SESAME_COL_II) {
            pool[n_pool++] = s->MR[i];          /* inferred green -> red is oob */
            pool[n_pool++] = s->UR[i];
        }
    }

    n_ok = sesame__drop_na(pool, n_pool, buf);
    if (n_ok == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
                          "no out-of-band signal for channel inference");
        goto out;
    }
    bg_max = sesame__quantile7(buf, n_ok, 0.95);

    /* Pass 3: failed probes revert to the manifest channel (default). */
    for (i = 0; i < n; i++) {
        int failed;
        if (s->col[i] == SESAME_COL_II) continue;

        failed = isnan(rmax[i]) || isnan(gmax[i]) ||
                 (sesame__pmax2(rmax[i], gmax[i]) < bg_max);

        if (newc[i] == SESAME_COL_II) {
            /* R's NA. With switch_failed=TRUE R would leave col=NA -- a
             * degenerate state our encoding cannot represent (and which would
             * silently turn an Inf-I probe into Inf-II). Revert instead. Moot
             * at the default switch_failed=FALSE, where `failed` is already
             * TRUE for these and reverts them anyway. */
            newc[i] = s->col[i];
        } else if (failed && !switch_failed) {
            newc[i] = s->col[i];
        }
        if (failed && mask_failed) s->mask[i] = 1;
        s->col[i] = newc[i];
    }

out:
    free(rmax); free(gmax); free(newc); free(pool); free(buf);
    return rc;
}
