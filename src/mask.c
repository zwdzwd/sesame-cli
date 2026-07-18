/* mask.c -- quality masking (the Q step), reading YAME .cm mask files.
 *
 * The recommended mask sets live in a YAME-packed .cm (a family of format-0
 * bit-tracks, one bit per ordering row). We link YAME directly (both projects
 * are AGPL-3.0, same author) and read the .cm in-process: no `yame` binary at
 * runtime, no fork/exec.
 *
 * A probe is quality-masked iff it is set in ANY recommended track -- the union
 * of recommendedMaskNames(platform) (R/mask.R:197-220).
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

/* YAME (linked): cdata/cfile I/O, the sample index, the fmt0 bit accessor, and
 * prepare_mask (normalizes a mask block to fmt0). */
#include "cfile.h"    /* open_cfile, read_cdata1, cdata_t, cfile_t, bgzf_close */
#include "index.h"    /* loadSampleNamesFromIndex, snames_t, cleanSampleNames2 */
#include "summary.h"  /* prepare_mask */
#include "cdata.h"    /* FMT0_IN_SET, free_cdata */

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* recommendedMaskNames() (R/mask.R:197-220); NULL-terminated. */
static const char *const REC_EPIC[] = {
    "mapping", "channel_switch", "snp5_GMAF1p", "extension", "sub30_copy", NULL };
static const char *const REC_MSA[] = {   /* also EPICv2 */
    "M_1baseSwitchSNPcommon_5pt", "M_2extBase_SNPcommon_5pt",
    "M_mapping", "M_nonuniq", "M_SNPcommon_5pt", NULL };

static const char *const *recommended_names(const char *platform)
{
    if (!platform) return NULL;
    if (strcmp(platform, "EPIC") == 0 || strcmp(platform, "HM450") == 0)
        return REC_EPIC;
    if (strcmp(platform, "MSA") == 0 || strcmp(platform, "EPICv2") == 0)
        return REC_MSA;
    return NULL;
}

static int is_recommended(const char *name, const char *const *rec)
{
    for (; *rec; rec++) if (strcmp(name, *rec) == 0) return 1;
    return 0;
}

/* Locate the platform's .cm mask in the store: <store>/<platform>/ *.cm
 * (ending in ".cm", not ".cm.idx"). 0 on success. */
static int find_cm(const char *platform, char *out, size_t n)
{
    char store[4096], pdir[4096];
    DIR *d;
    struct dirent *de;
    int found = 0;

    sesame_store_dir(store, sizeof store);
    snprintf(pdir, sizeof pdir, "%s/%s", store, platform);
    if (!(d = opendir(pdir))) return -1;
    while ((de = readdir(d))) {
        size_t L = strlen(de->d_name);
        if (L > 3 && strcmp(de->d_name + L - 3, ".cm") == 0) {
            snprintf(out, n, "%s/%s", pdir, de->d_name);
            found = 1;
            break;
        }
    }
    closedir(d);
    return found ? 0 : -1;
}

/* Load the recommended quality mask for a platform as a 0/1 vector aligned to
 * the ordering (mask[i] == 1 -> probe i is masked). *out is malloc'd (caller
 * frees), *out_n is the probe count. */
int sesame_quality_mask(const char *platform, uint8_t **out, int32_t *out_n,
                        sesame_err_t *err)
{
    const char *const *names = recommended_names(platform);
    char cm[4096];
    cfile_t cf;
    snames_t sn;
    uint8_t *mask = NULL;
    int32_t n = 0;
    int k = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!names)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "no recommended mask set for platform '%s'", platform);
    if (find_cm(platform, cm, sizeof cm) != 0)
        return sesame__fail(err, SESAME_ERR_IO,
            "no .cm mask for %s in the store -- run: sesame fetch %s",
            platform, platform);

    /* All 26 tracks come back in index order; keep the union of the recommended
     * ones. The file is small, so reading every block is cheap. */
    cf = open_cfile(cm);
    sn = loadSampleNamesFromIndex(cm);
    for (;;) {
        cdata_t c = read_cdata1(&cf);
        const char *name;
        if (c.n == 0) break;                 /* EOF (methscope pattern) */
        prepare_mask(&c);                    /* normalize to fmt0 bit form */
        if (!mask) {
            n = (int32_t)c.n;                /* bit count == probe count */
            mask = (uint8_t *)calloc((size_t)n, 1);
        }
        name = (k < sn.n) ? sn.s[k] : "";
        if (mask && is_recommended(name, names)) {
            int32_t i;
            for (i = 0; i < n; i++) if (FMT0_IN_SET(c, i)) mask[i] = 1;
        }
        free_cdata(&c);
        k++;
    }
    bgzf_close(cf.fh);
    cleanSampleNames2(sn);

    if (!mask)
        return sesame__fail(err, SESAME_ERR_FORMAT, "empty .cm mask for %s", platform);

    *out = mask;
    *out_n = n;
    return SESAME_OK;
}
