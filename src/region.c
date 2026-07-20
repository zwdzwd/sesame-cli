/* region.c -- extract a genomic region's betas from a multi-sample beta .cg as
 * long-form (tidy) TSV: one row per (probe in region) x (sample), columns
 *   chrom  beg  end  Probe_ID  beta  sample_name
 * The beta .cg, the ordering, and the per-probe coordinate table are all
 * positional in the same ordering, so probe i lines up across all three. This is
 * the plot-ready feed for a locus / region view (e.g. cinderplot).
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
#include <zlib.h>

/* one logical gz line into buf, newline stripped. 1 on a line, 0 EOF, -1 OOM. */
static int gzline(gzFile f, char **buf, size_t *cap)
{
    size_t len;
    if (!gzgets(f, *buf, (int)*cap)) return 0;
    len = strlen(*buf);
    while (len + 1 == *cap && (*buf)[len-1] != '\n') {
        char *n = (char *)realloc(*buf, *cap * 2);
        if (!n) return -1;
        *buf = n; *cap *= 2;
        if (!gzgets(f, *buf + len, (int)(*cap - len))) break;
        len += strlen(*buf + len);
    }
    (*buf)[strcspn(*buf, "\r\n")] = '\0';
    return 1;
}

/* Load the coordinate table (CpG_chrm CpG_beg ...), positional in the ordering:
 * chrom[i] strdup'd ("" if unmapped), pos[i] the CpG position (or -1). */
static int load_coords(const char *path, int32_t np, char ***chrom_out,
                       long **pos_out, sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf, *tab, *tab2; size_t cap = 1 << 16;
    char **chrom; long *pos;
    int32_t row = 0, r;

    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    chrom = (char **)malloc((size_t)np * sizeof(char *));
    pos = (long *)malloc((size_t)np * sizeof(long));
    if (!chrom || !pos) { free(buf); free(chrom); free(pos); gzclose(f); return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }

    gzline(f, &buf, &cap);                                   /* header */
    while ((r = gzline(f, &buf, &cap)) == 1) {
        if (row >= np) { row++; continue; }
        tab = strchr(buf, '\t'); if (tab) *tab = '\0';
        if (buf[0] == '\0' || !strcmp(buf, "*") || !strcmp(buf, "NA")) {
            chrom[row] = strdup(""); pos[row] = -1;
        } else {
            chrom[row] = strdup(buf);
            tab2 = tab ? strchr(tab + 1, '\t') : NULL;
            if (tab2) *tab2 = '\0';
            pos[row] = tab ? strtol(tab + 1, NULL, 10) : -1;
        }
        row++;
    }
    free(buf); gzclose(f);
    if (r < 0 || row != np) {
        for (int32_t i = 0; i < row && i < np; i++) free(chrom[i]);
        free(chrom); free(pos);
        if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d data rows, ordering has %d -- lineage mismatch", path, row, np);
    }
    *chrom_out = chrom; *pos_out = pos;
    return SESAME_OK;
}

int sesame_region_extract(const char *beta_cg, const sesame_index_t *ix,
                          const char *coords_path, const char *chrom,
                          long beg, long end, FILE *out, sesame_err_t *err)
{
    double *mat = NULL;
    char **names = NULL, **cchrom = NULL;
    long *cpos = NULL;
    int32_t np = 0, ns = 0, i, s, nprobe_ix, kept = 0;
    int rc = SESAME_ERR_IO;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }

    if (sesame_read_cg(beta_cg, &mat, &np, &ns, &names, err) != SESAME_OK) return err ? err->code : SESAME_ERR_IO;
    nprobe_ix = sesame_index_nprobes(ix);
    if (np != nprobe_ix) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d probes, ordering has %d", beta_cg, np, nprobe_ix);
        goto done;
    }
    if (load_coords(coords_path, np, &cchrom, &cpos, err) != SESAME_OK) goto done;

    fputs("chrom\tbeg\tend\tProbe_ID\tbeta\tsample_name\n", out);
    for (i = 0; i < np; i++) {
        if (cpos[i] < 0 || strcmp(cchrom[i], chrom) != 0) continue;
        if (cpos[i] < beg || cpos[i] > end) continue;
        kept++;
        for (s = 0; s < ns; s++) {
            double b = mat[(size_t)s * (size_t)np + (size_t)i];
            fprintf(out, "%s\t%ld\t%ld\t%s\t", cchrom[i], cpos[i], cpos[i] + 2,
                    sesame_index_probe_id(ix, i));
            if (isnan(b)) fputs("NA", out); else fprintf(out, "%.4f", b);
            fprintf(out, "\t%s\n", names[s]);
        }
    }
    if (kept == 0)
        fprintf(stderr, "sesame: no probes in %s:%ld-%ld\n", chrom, beg, end);
    rc = SESAME_OK;

done:
    free(mat);
    if (names) { for (s = 0; s < ns; s++) free(names[s]); free(names); }
    if (cchrom) { for (i = 0; i < np; i++) free(cchrom[i]); free(cchrom); }
    free(cpos);
    return rc;
}
