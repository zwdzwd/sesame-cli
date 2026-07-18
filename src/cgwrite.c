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

/* mat is sample-major: sample j, probe i at mat[j*nprobe + i]. Writes one
 * format-4 cdata block per sample to <path>, and <path>.idx with each sample's
 * name and the bgzf virtual offset of its block (for yame's named access). */
int sesame_write_cg(const char *path, const double *mat, int32_t nprobe,
                    int32_t nsamp, char *const *names, sesame_err_t *err)
{
    BGZF *fp;
    int64_t *offs;
    int32_t j, i;
    char *idxpath;

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
        offs[j] = bgzf_tell(fp);                    /* block start, for the index */
        cdata_compress(&c);                         /* frees s, sets c.s compressed */
        cdata_write1(fp, &c);
        free(c.s);
    }
    bgzf_close(fp);

    idxpath = (char *)malloc(strlen(path) + 5);
    if (idxpath) {
        FILE *ix;
        snprintf(idxpath, strlen(path) + 5, "%s.idx", path);
        if ((ix = fopen(idxpath, "w"))) {
            for (j = 0; j < nsamp; j++)
                fprintf(ix, "%s\t%lld\n", names[j], (long long)offs[j]);
            fclose(ix);
        }
        free(idxpath);
    }
    free(offs);
    return SESAME_OK;
}
