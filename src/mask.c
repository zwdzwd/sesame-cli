/* mask.c -- quality masking (the Q step), reading YAME .cm mask files.
 *
 * The recommended mask sets live in a YAME-packed .cm (a family of format-0
 * bit-tracks, one bit per ordering row). YAME is AGPL and its format code stays
 * out of this build: we shell out to the `yame` binary and read its text output.
 * Process separation, not linking -- so YAME's licence does not reach here and
 * its format stays its own concern.
 *
 * `yame unpack -l <names> <file.cm>` emits a tab-delimited table, one row per
 * probe (in ordering order) and one 0/1 column per requested track. A probe is
 * quality-masked iff it is set in ANY recommended track -- the union.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

/* recommendedMaskNames() (R/mask.R:197-220). The union of these is what
 * qualityMask applies by default. Terminated by NULL. */
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
    return NULL;   /* unknown platform: caller reports it */
}

/* Locate the platform's .cm mask in the store: <store>/<platform>/ *.cm
 * (the file ending in ".cm", not ".cm.idx"). 0 on success. */
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
 * the ordering (mask[i] == 1 -> probe i is masked). *out is malloc'd; caller
 * frees. *out_n is the probe count. Shells out to `yame`. */
int sesame_quality_mask(const char *platform, uint8_t **out, int32_t *out_n,
                        sesame_err_t *err)
{
    const char *const *names = recommended_names(platform);
    char cm[4096], listfile[] = "/tmp/sesame_recmask_XXXXXX";
    char **argv;
    int nnames = 0, i, pipefd[2], fd, status;
    pid_t pid;
    FILE *lf, *rd;
    uint8_t *mask = NULL;
    int32_t cap = 1 << 16, n = 0;
    char line[512];

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!names)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "no recommended mask set for platform '%s'", platform);
    if (find_cm(platform, cm, sizeof cm) != 0)
        return sesame__fail(err, SESAME_ERR_IO,
            "no .cm mask for %s in the store -- run: sesame fetch %s",
            platform, platform);
    while (names[nnames]) nnames++;

    /* Track names go in a list file (yame -l), never on the command line -- no
     * shell, no injection surface. */
    if ((fd = mkstemp(listfile)) < 0)
        return sesame__fail(err, SESAME_ERR_IO, "cannot create temp file");
    lf = fdopen(fd, "w");
    for (i = 0; i < nnames; i++) fprintf(lf, "%s\n", names[i]);
    fclose(lf);

    /* argv: yame unpack -l <listfile> <cm> */
    argv = (char **)calloc(6, sizeof(char *));
    argv[0] = (char *)"yame"; argv[1] = (char *)"unpack";
    argv[2] = (char *)"-l";   argv[3] = listfile;
    argv[4] = cm;             argv[5] = NULL;

    if (pipe(pipefd) != 0) {
        free(argv); remove(listfile);
        return sesame__fail(err, SESAME_ERR_IO, "pipe failed");
    }
    pid = fork();
    if (pid < 0) {
        close(pipefd[0]); close(pipefd[1]); free(argv); remove(listfile);
        return sesame__fail(err, SESAME_ERR_IO, "fork failed");
    }
    if (pid == 0) {                          /* child: yame -> pipe */
        dup2(pipefd[1], STDOUT_FILENO);
        close(pipefd[0]); close(pipefd[1]);
        execvp("yame", argv);
        _exit(127);                          /* execvp failed (yame not found) */
    }
    close(pipefd[1]);
    free(argv);

    mask = (uint8_t *)malloc((size_t)cap);
    rd = fdopen(pipefd[0], "r");
    while (rd && fgets(line, sizeof line, rd)) {
        int masked = 0;
        const char *p = line;
        /* union across the tab-separated 0/1 columns */
        for (;;) {
            while (*p == ' ' || *p == '\t') p++;
            if (*p == '\0' || *p == '\n') break;
            if (*p != '0') masked = 1;
            while (*p && *p != '\t' && *p != '\n') p++;
        }
        if (n >= cap) {
            uint8_t *nm = (uint8_t *)realloc(mask, (size_t)(cap *= 2));
            if (!nm) { free(mask); mask = NULL; break; }
            mask = nm;
        }
        if (mask) mask[n++] = (uint8_t)masked;
    }
    if (rd) fclose(rd); else close(pipefd[0]);
    waitpid(pid, &status, 0);
    remove(listfile);

    if (!mask)
        return sesame__fail(err, SESAME_ERR_NOMEM, "oom reading mask");
    if (WIFEXITED(status) && WEXITSTATUS(status) == 127) {
        free(mask);
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "the 'yame' tool is required to read .cm masks but was not found on "
            "PATH\n  install it (e.g. conda install -c bioconda yame)");
    }
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
        free(mask);
        return sesame__fail(err, SESAME_ERR_IO, "yame unpack failed for %s", cm);
    }

    *out = mask;
    *out_n = n;
    return SESAME_OK;
}
