/* attach.c -- attach the ordering's Probe_IDs to a positional data file.
 *
 * A YAME .cg/.cm/.cx stores one value per probe in ordering order, with NO probe
 * id inside the container (row names live in the ordering, not the data). The
 * per-probe genomic coordinate tables (<platform>.hg38.coord.tsv.gz) are the
 * same: positional, no Probe_ID column. This turns either into a labeled TSV by
 * pairing row i with the ordering's i-th Probe_ID -- so the output is directly
 * greppable / joinable. The lineage must match (same platform + tag that
 * produced the file); a row-count mismatch is a hard error, not silent
 * misalignment.
 *
 * The YAME per-format rendering mirrors `yame unpack`'s print_cdata1 (fmt0 mask
 * bit, fmt3 M/U or beta, fmt4 float, fmt5 ternary, fmt1/2 raw), except floats
 * print at full precision rather than 3 decimals since the .cg already stores
 * float32. Links YAME directly (both AGPL), same as cgwrite.c.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include "cfile.h"    /* open_cfile, read_cdata1 */
#include "cdata.h"    /* cdata_t, decompress, free_cdata, f3_get_mu, f2_get_string */
#include "index.h"    /* loadSampleNamesFromIndex, cleanSampleNames2, snames_t */

#include <ctype.h>
#include <inttypes.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

static int has_suffix(const char *s, const char *suf)
{
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

/* A YAME container, by extension. Everything else is treated as text. */
static int is_yame(const char *path)
{
    return has_suffix(path, ".cg") || has_suffix(path, ".cm") ||
           has_suffix(path, ".cx") || has_suffix(path, ".cr");
}

/* Read one logical (possibly long) line into the grown buffer, stripping the
 * trailing newline. Returns 1 on a line, 0 at EOF, -1 on OOM. */
static int read_gzline(gzFile f, char **buf, size_t *cap)
{
    size_t len = 0;
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

/* Count data rows (total lines minus the header, unless no_header). */
static int count_text_rows(const char *path, int no_header, int32_t *nrow,
                           sesame_err_t *err)
{
    gzFile f = gzopen(path, "rb");
    char *buf; size_t cap = 1 << 16;
    long lines = 0; int r;
    if (!f) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }
    while ((r = read_gzline(f, &buf, &cap)) == 1) lines++;
    free(buf); gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    *nrow = (int32_t)(no_header ? lines : (lines > 0 ? lines - 1 : 0));
    return SESAME_OK;
}

/* ------------------------------------------------------------- text ------
 *
 * gzopen reads plain and gzipped text alike. The file's own header is kept and
 * prefixed with "Probe_ID"; each data line gets its positional Probe_ID. Row
 * count is checked up front so a lineage mismatch fails before any output. */
static int attach_text(const char *path, const sesame_index_t *ix,
                       const sesame_attach_opt_t *opt, FILE *out,
                       sesame_err_t *err)
{
    gzFile f;
    char *buf = NULL;
    size_t cap = 1 << 16;
    int32_t nid = sesame_index_nprobes(ix), nrow = 0, row = 0;
    int first = 1, r;

    if (count_text_rows(path, opt->no_header, &nrow, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;
    if (nrow != nid)
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %d data rows, ordering has %d -- lineage mismatch",
            path, nrow, nid);

    if (!(f = gzopen(path, "rb")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    if (!(buf = (char *)malloc(cap))) { gzclose(f);
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom"); }

    while ((r = read_gzline(f, &buf, &cap)) == 1) {
        if (first && !opt->no_header) {
            fprintf(out, "Probe_ID\t%s\n", buf);
            first = 0;
            continue;
        }
        first = 0;
        fprintf(out, "%s\t%s\n", sesame_index_probe_id(ix, row), buf);
        row++;
    }
    free(buf);
    gzclose(f);
    if (r < 0) return sesame__fail(err, SESAME_ERR_NOMEM, "oom");
    return SESAME_OK;
}

/* ------------------------------------------------------------- YAME ------ */

/* One value of decompressed record d at row i. Mirrors yame unpack. */
static void render1(FILE *out, cdata_t *d, uint64_t i,
                    const sesame_attach_opt_t *opt)
{
    switch (d->fmt) {
    case '0':                              /* mask bit */
        fputc(((d->s[i>>3] >> (i&0x7)) & 0x1) + '0', out);
        break;
    case '1':                              /* raw byte */
        fputc(d->s[i], out);
        break;
    case '2':                              /* state label */
        fputs(f2_get_string(d, i), out);
        break;
    case '3': {                            /* M/U counts */
        uint64_t mu = f3_get_mu(d, i);
        uint64_t M = mu >> 32, U = (mu << 32) >> 32;
        if (opt->beta) {
            if (M == 0 && U == 0) fputs("NA", out);
            else fprintf(out, "%.6g", (double)M / (double)(M + U));
        } else {
            fprintf(out, "%" PRIu64 "\t%" PRIu64, M, U);
        }
        break;
    }
    case '4': {                            /* float32, negative = NA */
        float v = ((float *)d->s)[i];
        if (v < 0) fputs("NA", out); else fprintf(out, "%.6g", (double)v);
        break;
    }
    case '5':                              /* ternary, 2 = NA */
        if (d->s[i] == 2) fputs("NA", out); else fputc(d->s[i] + '0', out);
        break;
    default:
        fputs("NA", out);
        break;
    }
}

/* Header cell(s) for a sample under format fmt (fmt3 M/U is two columns). */
static void render_head(FILE *out, char fmt, const char *name,
                        const sesame_attach_opt_t *opt)
{
    if (fmt == '3' && !opt->beta) fprintf(out, "%s_M\t%s_U", name, name);
    else                          fputs(name, out);
}

static int attach_yame(const char *path, const sesame_index_t *ix,
                       const sesame_attach_opt_t *opt, FILE *out,
                       sesame_err_t *err)
{
    cfile_t cf;
    snames_t sn;
    cdata_t *recs = NULL;
    char pathbuf[4096];
    int32_t nid = sesame_index_nprobes(ix), nrec = 0, cap = 0, j;
    int64_t np = -1, i;
    int rc = SESAME_OK;

    snprintf(pathbuf, sizeof pathbuf, "%s", path);
    cf = open_cfile(pathbuf);
    if (!cf.fh) return sesame__fail(err, SESAME_ERR_IO, "cannot open %s", path);
    sn = loadSampleNamesFromIndex(pathbuf);

    for (;;) {
        cdata_t c = read_cdata1(&cf), d;
        if (c.n == 0) break;                            /* EOF */
        d = decompress(c);
        free_cdata(&c);
        if (np < 0) np = (int64_t)d.n;
        else if ((int64_t)d.n != np) {
            free_cdata(&d);
            rc = sesame__fail(err, SESAME_ERR_FORMAT,
                              "inconsistent record length in %s", path);
            goto done;
        }
        if (nrec >= cap) {
            cap = cap ? cap * 2 : 8;
            recs = (cdata_t *)realloc(recs, (size_t)cap * sizeof(cdata_t));
        }
        recs[nrec++] = d;
        if (!opt->all) break;                           /* first record only */
    }

    if (nrec == 0) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT, "no records in %s", path);
        goto done;
    }
    if (np != nid) {
        rc = sesame__fail(err, SESAME_ERR_FORMAT,
            "%s has %" PRId64 " probes, ordering has %d -- lineage mismatch",
            path, np, nid);
        goto done;
    }

    if (!opt->no_header) {
        fputs("Probe_ID", out);
        for (j = 0; j < nrec; j++) {
            const char *nm = (j < sn.n) ? sn.s[j] : "";
            fputc('\t', out);
            render_head(out, recs[j].fmt, (nm && *nm) ? nm : "V", opt);
        }
        fputc('\n', out);
    }
    for (i = 0; i < np; i++) {
        fputs(sesame_index_probe_id(ix, (int32_t)i), out);
        for (j = 0; j < nrec; j++) {
            fputc('\t', out);
            render1(out, &recs[j], (uint64_t)i, opt);
        }
        fputc('\n', out);
    }

done:
    for (j = 0; j < nrec; j++) free_cdata(&recs[j]);
    free(recs);
    cleanSampleNames2(sn);
    bgzf_close(cf.fh);
    return rc;
}

int sesame_attach_probe(const char *path, const sesame_index_t *ix,
                        const sesame_attach_opt_t *opt, FILE *out,
                        sesame_err_t *err)
{
    static const sesame_attach_opt_t deflt = { 0, 0, 0 };
    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!opt) opt = &deflt;
    return is_yame(path) ? attach_yame(path, ix, opt, out, err)
                         : attach_text(path, ix, opt, out, err);
}
