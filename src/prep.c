/* prep.c -- preprocessing steps (the prepSesame codes).
 *
 * Implemented so far: C (inferInfiniumIChannel).
 * Outstanding: Q, D, P, B.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---------------------------------------------------------------- Q ---
 *
 * qualityMask (R/mask.R:286-297): OR a precomputed recommended-mask vector into
 * the SigDF mask. The vector is loaded once per platform via sesame_quality_mask
 * (which shells out to yame); here we just apply it. It is aligned to the
 * ordering, so mask[i] corresponds to probe i. */
int sesame_prep_quality_mask(sesame_sigdf_t *s, const uint8_t *qmask,
                             int32_t qn, sesame_err_t *err)
{
    int32_t i;
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s || !qmask) return sesame__fail(err, SESAME_ERR_IO, "null argument");
    if (qn != s->n)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "mask length %d != probe count %d (ordering/mask mismatch)", qn, s->n);
    for (i = 0; i < s->n; i++) if (qmask[i]) s->mask[i] = 1;
    return SESAME_OK;
}

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

/* ---------------------------------------------------------------- D ---
 *
 * dyeBiasNL (R/dye_bias.R:118-167).
 *
 * Pull the Red and Green Infinium-I distributions to their common midpoint:
 * quantile-normalize Red onto Green, average that with Red-as-is, and use the
 * resulting IR1 -> IRmid map as a lookup applied to every Red value. Green is
 * handled symmetrically. Both maps are built from the ORIGINAL signals before
 * either is applied (R builds both closures first, R/dye_bias.R:164-165).
 */

/* A monotone lookup table: knots (x, y), linearly interpolated. */
typedef struct { double *x, *y; int32_t n; } knots_t;

/* Build knots from a sorted src and its aligned midpoint values, collapsing
 * tied src values. This is R's approx(ties=mean) (R/dye_bias.R:143): approx
 * needs unique x, so tied x collapse to one knot carrying the mean of their y.
 * Here the mean is a formality -- qnorm_use_target gives every element of a
 * tied run the same value (average rank), so y is constant across the run.
 * Verified on real data: 0 of 11259 tied runs on HM450 had non-constant y. */
static void knots_build(knots_t *k, const double *src, const double *mid,
                        int32_t n)
{
    int32_t a = 0;
    k->n = 0;
    while (a < n) {
        int32_t b = a, i;
        double sum = 0.0;
        while (b + 1 < n && src[b + 1] == src[a]) b++;
        for (i = a; i <= b; i++) sum += mid[i];
        k->x[k->n] = src[a];
        k->y[k->n] = sum / (double)(b - a + 1);
        k->n++;
        a = b + 1;
    }
}

/* Linear interpolation on the knots. Callers guarantee k->x[0] <= v <=
 * k->x[n-1], which is R's approx(rule=1) guarded by the explicit in-support
 * test, so no extrapolation is possible here. */
static double knots_interp(const knots_t *k, double v)
{
    int32_t lo = 0, hi = k->n - 1;

    if (v <= k->x[0])          return k->y[0];
    if (v >= k->x[k->n - 1])   return k->y[k->n - 1];

    while (hi - lo > 1) {                      /* largest lo with x[lo] <= v */
        int32_t mid = lo + ((hi - lo) >> 1);
        if (k->x[mid] <= v) lo = mid; else hi = mid;
    }
    if (k->x[lo] == v) return k->y[lo];
    return k->y[lo] + (v - k->x[lo]) / (k->x[lo + 1] - k->x[lo]) *
                      (k->y[lo + 1] - k->y[lo]);
}

/* R's fitfun (R/dye_bias.R:139-147): interpolate inside the support, shift
 * additively above it, scale multiplicatively below it. NA passes through. */
static double fitfun(const knots_t *k, double v,
                     double lo_in, double hi_in, double lo_mid, double hi_mid)
{
    if (isnan(v))   return v;
    if (v > hi_in)  return v - hi_in + hi_mid;          /* oversupp  */
    if (v < lo_in)  return lo_mid / lo_in * v;          /* undersupp */
    return knots_interp(k, v);                          /* insupp    */
}

/* Total intensity per probe, as totalIntensities(sdf, mask=FALSE) via signalMU:
 * col G -> MG+UG, col R -> MR+UR, col II -> UG+UR (R/sesame.R:111-115,
 * R/SigDFMethods.R:69-81). */
static double total_intensity(const sesame_sigdf_t *s, int32_t i)
{
    switch (s->col[i]) {
    case SESAME_COL_G: return s->MG[i] + s->UG[i];
    case SESAME_COL_R: return s->MR[i] + s->UR[i];
    default:           return s->UG[i] + s->UR[i];
    }
}

/* RGdistort = (topR/topG) / (medR/medG), where med* is the median total
 * intensity of that channel's Inf-I probes and top* the median of its largest
 * 20 (R/QC.R:359-366). NA or > 10 means the green channel has failed, and
 * dyeBiasNL gives up and masks all Inf-I green probes. */
static double rg_distort(const sesame_sigdf_t *s, double *buf)
{
    double medR, medG, topR, topG;
    int32_t i, nR = 0, nG = 0;

    for (i = 0; i < s->n; i++)
        if (s->col[i] == SESAME_COL_R) {
            double v = total_intensity(s, i);
            if (!isnan(v)) buf[nR++] = v;
        }
    sesame__sort(buf, nR);
    if (nR == 0) return NAN;
    medR = sesame__median_sorted(buf, nR);
    topR = sesame__median_sorted(buf + (nR > 20 ? nR - 20 : 0),
                                 nR > 20 ? 20 : nR);

    for (i = 0; i < s->n; i++)
        if (s->col[i] == SESAME_COL_G) {
            double v = total_intensity(s, i);
            if (!isnan(v)) buf[nG++] = v;
        }
    sesame__sort(buf, nG);
    if (nG == 0) return NAN;
    medG = sesame__median_sorted(buf, nG);
    topG = sesame__median_sorted(buf + (nG > 20 ? nG - 20 : 0),
                                 nG > 20 ? 20 : nG);

    return (topR / topG) / (medR / medG);
}

int sesame_prep_dye_bias_nl(sesame_sigdf_t *s, sesame_err_t *err)
{
    int32_t n, i, nIG = 0, nIR = 0, nIG1, nIR1;
    double *IG0 = NULL, *IR0 = NULL, *IG1 = NULL, *IR1 = NULL;
    double *IG2 = NULL, *IR2 = NULL, *IGmid = NULL, *IRmid = NULL, *tbuf = NULL;
    knots_t kR = {0}, kG = {0};
    double maxIG, minIG, maxIR, minIR, maxIGmid, minIGmid, maxIRmid, minIRmid;
    double dist;
    int rc = SESAME_OK;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!s) return sesame__fail(err, SESAME_ERR_IO, "null sigdf");
    n = s->n;

    tbuf = (double *)malloc((size_t)n * sizeof(double));
    if (!tbuf) return sesame__fail(err, SESAME_ERR_NOMEM, "oom in D step");

    /* Guard: green channel failed completely -> mask Inf-I green and stop
     * (R/dye_bias.R:122-123, maskIG at :94-97). R does this silently; we set a
     * status bit so a pipeline can count it. */
    dist = rg_distort(s, tbuf);
    if (isnan(dist) || dist > 10.0) {
        for (i = 0; i < n; i++)
            if (s->col[i] == SESAME_COL_G) s->mask[i] = 1;
        s->status |= SESAME_STAT_DYEBIAS_FAILED;
        free(tbuf);
        return SESAME_OK;
    }

    /* Pools. mask=TRUE is the default and means "include masked probes"
     * (R/dye_bias.R:126) -- so this is every Inf-I probe, masked or not. */
    IG0 = (double *)malloc((size_t)n * 2 * sizeof(double));
    IR0 = (double *)malloc((size_t)n * 2 * sizeof(double));
    if (!IG0 || !IR0) { rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out; }
    for (i = 0; i < n; i++) {
        if (s->col[i] == SESAME_COL_G) {
            IG0[nIG++] = s->MG[i]; IG0[nIG++] = s->UG[i];
        } else if (s->col[i] == SESAME_COL_R) {
            IR0[nIR++] = s->MR[i]; IR0[nIR++] = s->UR[i];
        }
    }
    if (nIG == 0 || nIR == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
                          "no Infinium-I %s probes for dye bias",
                          nIG == 0 ? "green" : "red");
        goto out;
    }

    /* min/max over the raw pools, na.rm=TRUE (R/dye_bias.R:130-131). */
    maxIG = -INFINITY; minIG = INFINITY; maxIR = -INFINITY; minIR = INFINITY;
    for (i = 0; i < nIG; i++) if (!isnan(IG0[i])) {
        if (IG0[i] > maxIG) maxIG = IG0[i];
        if (IG0[i] < minIG) minIG = IG0[i];
    }
    for (i = 0; i < nIR; i++) if (!isnan(IR0[i])) {
        if (IR0[i] > maxIR) maxIR = IR0[i];
        if (IR0[i] < minIR) minIR = IR0[i];
    }
    if (maxIG <= 0.0 || maxIR <= 0.0) {   /* R/dye_bias.R:133 */
        rc = SESAME_OK;
        goto out;
    }

    IG1 = (double *)malloc((size_t)nIG * sizeof(double));
    IR1 = (double *)malloc((size_t)nIR * sizeof(double));
    IG2 = (double *)malloc((size_t)nIG * sizeof(double));
    IR2 = (double *)malloc((size_t)nIR * sizeof(double));
    IGmid = (double *)malloc((size_t)nIG * sizeof(double));
    IRmid = (double *)malloc((size_t)nIR * sizeof(double));
    kR.x = (double *)malloc((size_t)nIR * sizeof(double));
    kR.y = (double *)malloc((size_t)nIR * sizeof(double));
    kG.x = (double *)malloc((size_t)nIG * sizeof(double));
    kG.y = (double *)malloc((size_t)nIG * sizeof(double));
    if (!IG1 || !IR1 || !IG2 || !IR2 || !IGmid || !IRmid ||
        !kR.x || !kR.y || !kG.x || !kG.y) {
        rc = sesame__fail(err, SESAME_ERR_NOMEM, "oom"); goto out;
    }

    /* sort() drops NA (na.last=NA), so IR1 can be shorter than IR0. */
    nIR1 = sesame__drop_na(IR0, nIR, IR1); sesame__sort(IR1, nIR1);
    nIG1 = sesame__drop_na(IG0, nIG, IG1); sesame__sort(IG1, nIG1);
    if (nIR1 == 0 || nIG1 == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no non-NA Inf-I signal");
        goto out;
    }

    /* Red onto Green, and Green onto Red. R re-sorts the result; the output is
     * already non-decreasing, so that sort is a no-op. tgt is sorted in place,
     * so pass a scratch copy of the raw pool. */
    memcpy(tbuf, IG0, (size_t)nIG * sizeof(double));
    {
        int32_t mt = sesame__drop_na(tbuf, nIG, tbuf);   /* NA target dropped */
        sesame__qnorm_use_target(IR1, nIR1, tbuf, mt, IR2);
    }
    for (i = 0; i < nIR1; i++) IRmid[i] = (IR1[i] + IR2[i]) / 2.0;

    memcpy(tbuf, IR0, (size_t)nIR * sizeof(double));
    {
        int32_t mt = sesame__drop_na(tbuf, nIR, tbuf);
        sesame__qnorm_use_target(IG1, nIG1, tbuf, mt, IG2);
    }
    for (i = 0; i < nIG1; i++) IGmid[i] = (IG1[i] + IG2[i]) / 2.0;

    minIRmid = INFINITY; maxIRmid = -INFINITY;
    for (i = 0; i < nIR1; i++) {
        if (IRmid[i] < minIRmid) minIRmid = IRmid[i];
        if (IRmid[i] > maxIRmid) maxIRmid = IRmid[i];
    }
    minIGmid = INFINITY; maxIGmid = -INFINITY;
    for (i = 0; i < nIG1; i++) {
        if (IGmid[i] < minIGmid) minIGmid = IGmid[i];
        if (IGmid[i] > maxIGmid) maxIGmid = IGmid[i];
    }

    knots_build(&kR, IR1, IRmid, nIR1);
    knots_build(&kG, IG1, IGmid, nIG1);

    /* Both maps are built from the original signals before either is applied. */
    for (i = 0; i < n; i++) {
        s->MR[i] = fitfun(&kR, s->MR[i], minIR, maxIR, minIRmid, maxIRmid);
        s->UR[i] = fitfun(&kR, s->UR[i], minIR, maxIR, minIRmid, maxIRmid);
        s->MG[i] = fitfun(&kG, s->MG[i], minIG, maxIG, minIGmid, maxIGmid);
        s->UG[i] = fitfun(&kG, s->UG[i], minIG, maxIG, minIGmid, maxIGmid);
    }

out:
    free(IG0); free(IR0); free(IG1); free(IR1); free(IG2); free(IR2);
    free(IGmid); free(IRmid); free(tbuf);
    free(kR.x); free(kR.y); free(kG.x); free(kG.y);
    return rc;
}
