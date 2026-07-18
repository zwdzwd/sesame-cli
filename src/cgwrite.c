/* cgwrite.c -- write per-probe values as a YAME format-4 .cg file (one float32
 * per probe, NA encoded as -1.0, RLE over NA runs), plus a companion .idx of
 * "<sample>\t<bgzf-offset>" lines. Links YAME directly (both AGPL), same as the
 * mask reader.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include "cfile.h"    /* cdata_write1, BGZF */
#include "cdata.h"    /* cdata_t, cdata_compress, bgzf_open2/tell/close */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Write the companion "<path>.idx" of "<name>\t<bgzf-offset>" lines. */
static void write_idx(const char *path, char *const *names,
                      const int64_t *offs, int32_t nsamp)
{
    size_t z = strlen(path) + 5;
    char *idxpath = (char *)malloc(z);
    FILE *ix;
    int32_t j;
    if (!idxpath) return;
    snprintf(idxpath, z, "%s.idx", path);
    if ((ix = fopen(idxpath, "w"))) {
        for (j = 0; j < nsamp; j++)
            fprintf(ix, "%s\t%lld\n", names[j], (long long)offs[j]);
        fclose(ix);
    }
    free(idxpath);
}

/* Format 3 (M/U counts), the default. Rounds M,U to integers -- exact for raw
 * IDAT signal (uint16 means), sub-integer rounding for preprocessed floats
 * (noise). NA (either allele NaN) is stored as (0,0), which yame reads as
 * beta=NA / cov=0. yame derives beta (MU2beta) and total intensity (MU2cov).
 * matM/matU are sample-major: sample j, probe i at mat[j*nprobe+i]. */
int sesame_write_cg_mu(const char *path, const double *matM, const double *matU,
                       int32_t nprobe, int32_t nsamp, char *const *names,
                       sesame_err_t *err)
{
    BGZF *fp;
    int64_t *offs;
    int32_t j, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!(fp = bgzf_open2(path, "w")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s for writing", path);
    if (!(offs = (int64_t *)malloc((size_t)nsamp * sizeof(int64_t)))) {
        bgzf_close(fp); return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    for (j = 0; j < nsamp; j++) {
        cdata_t c;
        memset(&c, 0, sizeof c);
        c.fmt = '3'; c.compressed = 0; c.n = (uint64_t)nprobe; c.unit = 8;
        c.s = (uint8_t *)calloc((size_t)nprobe, 8);
        if (!c.s) { free(offs); bgzf_close(fp);
                    return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
        for (i = 0; i < nprobe; i++) {
            double m = matM[(size_t)j*(size_t)nprobe + (size_t)i];
            double u = matU[(size_t)j*(size_t)nprobe + (size_t)i];
            if (isnan(m) || isnan(u)) f3_set_mu(&c, (uint64_t)i, 0, 0);   /* NA */
            else f3_set_mu(&c, (uint64_t)i, (uint64_t)llround(m), (uint64_t)llround(u));
        }
        offs[j] = bgzf_tell(fp);
        cdata_compress(&c);
        cdata_write1(fp, &c);
        free(c.s);
    }
    bgzf_close(fp);
    write_idx(path, names, offs, nsamp);
    free(offs);
    return SESAME_OK;
}

/* Format 4 (one float32 per probe, NA=-1.0), on request -- exact values (e.g.
 * betas or total intensity) with no rounding. */
int sesame_write_cg(const char *path, const double *mat, int32_t nprobe,
                    int32_t nsamp, char *const *names, sesame_err_t *err)
{
    BGZF *fp;
    int64_t *offs;
    int32_t j, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!(fp = bgzf_open2(path, "w")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s for writing", path);
    if (!(offs = (int64_t *)malloc((size_t)nsamp * sizeof(int64_t)))) {
        bgzf_close(fp); return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    }
    for (j = 0; j < nsamp; j++) {
        cdata_t c;
        float *s = (float *)malloc((size_t)nprobe * sizeof(float));
        if (!s) { free(offs); bgzf_close(fp);
                  return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
        for (i = 0; i < nprobe; i++) {
            double v = mat[(size_t)j * (size_t)nprobe + (size_t)i];
            s[i] = isnan(v) ? -1.0f : (float)v;     /* NA -> negative (bit31 set) */
        }
        memset(&c, 0, sizeof c);
        c.fmt = '4'; c.compressed = 0; c.n = (uint64_t)nprobe;
        c.unit = sizeof(float); c.s = (uint8_t *)s;
        offs[j] = bgzf_tell(fp);
        cdata_compress(&c);                         /* frees s, sets c.s compressed */
        cdata_write1(fp, &c);
        free(c.s);
    }
    bgzf_close(fp);
    write_idx(path, names, offs, nsamp);
    free(offs);
    return SESAME_OK;
}
