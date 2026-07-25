/* cache.c -- index discovery and fetching.
 *
 * The store is SHARED with every other tool in the suite (kycg, methscope-cli)
 * and lives in YAME, which all of them already link: see YAME/src/assets.h. It
 * mirrors the remote exactly, one version, keyed by upstream repo and platform
 * so a file one tool downloads is found by the next:
 *
 *   <store>/InfiniumAnnotation/MSA/SHA256SUMS   byte-identical copy of the remote
 *   <store>/InfiniumAnnotation/MSA/MSA.ordering.tsv.gz
 *   <store>/InfiniumAnnotation/MSA/MSA.hg38.mask.cm
 *   <store>/InfiniumAnnotation/MSA/KYCG/...     kycg's knowledgebase sets
 *   <store>/genomes/hg38/seqinfo.tsv.gz
 *
 * Before this, sesame kept its own cache under ~/.cache/sesame and kycg kept
 * one under ~/.cache/kycg, and the two of them downloaded the identical
 * InfiniumAnnotation files separately.
 *
 * Mirroring means the store is self-describing and hand-verifiable per platform
 * (cd <store>/InfiniumAnnotation/MSA && shasum -a 256 -c SHA256SUMS), with no
 * local bookkeeping -- and that same stored SHA256SUMS is how a directory says
 * which upstream tag filled it, so two tools pinned to different tags cannot
 * silently overwrite each other. See yame_assets_pin_check().
 *
 * No tag subdirectories, no `current` symlink, no local version management.
 * Most users only ever want one version; anyone who genuinely needs two can
 * point SESAME_INDEX_DIR at two directories, which is simpler than anything we
 * could build for them.
 *
 * De-duplication falls out of that: fetch pulls the platform's SHA256SUMS first
 * (a few hundred bytes), compares each digest against what is already on disk,
 * and downloads only the files that differ. When the build's pinned tag bumps,
 * unchanged files are not re-downloaded; upstream, git is content-addressed, so
 * they cost nothing in the new tag either.
 *
 * The trust anchor is a per-platform sha256(SHA256SUMS) compiled into this build
 * (see registry.h). Fetch verifies the downloaded SHA256SUMS against it, then
 * verifies every file against a digest from that SHA256SUMS -- a hard chain.
 *
 * Store location (fetch writes here):
 *
 *   1. $SESAME_INDEX_DIR                explicit, still sesame's own
 *   2. <dir of the binary>/data         a checkout: found from any cwd
 *   3. the shared root, yame_assets_root(): $YAME_DATA_HOME, else
 *      ${XDG_DATA_HOME:-~/.local/share}/yame
 *
 * Step 3 replaced ~/.cache/sesame. The data tier rather than a cache tier is
 * deliberate -- these are multi-GB references that should survive somebody
 * clearing ~/.cache. An existing ~/.cache/sesame is detected and named in a
 * one-line notice, never read from and never moved.
 *
 * Resolution order (explicit always wins; same rule for library, CLI and any
 * future binding):
 *
 *   1. --index <path>                                 (caller-supplied)
 *   2. <store>/InfiniumAnnotation/<platform>/<platform>.ordering.tsv.gz
 *   3. ./<platform>.ordering.tsv.gz
 *
 * sesame NEVER prompts and NEVER downloads implicitly. If the index is missing,
 * the caller gets an error naming the exact command to run. An interactive
 * prompt would hang forever in a Nextflow job or a Docker build -- silently, on
 * question one, across a whole run. `sesame fetch` is the only download path, so
 * behaviour is identical whether or not a TTY is attached.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "internal.h"
#include "registry.h"
#include "assets.h"                  /* YAME: the shared store and downloader */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

/* The store layout, keyed by upstream repo so every tool agrees on it. */
#define SESAME_IA_SUB      "InfiniumAnnotation"
#define SESAME_GENOME_SUB  "genomes"

static int is_file(const char *p)
{
    return p && *p && yame_assets_is_file(p);
}

const sesame_reg_t *sesame__reg_for_platform(const char *platform)
{
    if (!platform) return NULL;
    for (const sesame_reg_t *r = SESAME_REGISTRY; r->platform; r++)
        if (strcmp(r->platform, platform) == 0) return r;
    return NULL;
}

const char *sesame_default_tag(void) { return SESAME_DEFAULT_TAG; }

const char *sesame_platform_from_beads(int32_t beads)
{
    for (const sesame_reg_t *r = SESAME_REGISTRY; r->platform; r++)
        if (r->beads && r->beads == beads) return r->platform;
    return NULL;
}

static int is_dir(const char *p)
{
    return p && *p && yame_assets_is_dir(p);
}

/* Absolute directory containing the running executable, symlinks resolved.
 * Returns 0 on success.
 *
 * This is what makes the store follow the *binary* rather than the working
 * directory: a checkout keeps data/ beside ./sesame, so it is found from
 * anywhere, whereas a CWD-relative ./data evaporates the moment a pipeline
 * cd's into a task scratch dir. */
static int exe_dir(char *out, size_t n)
{
    char raw[4096], real[4096];
    char *slash;

#ifdef __APPLE__
    uint32_t sz = (uint32_t)sizeof raw;
    if (_NSGetExecutablePath(raw, &sz) != 0) return -1;
#elif defined(__linux__)
    ssize_t k = readlink("/proc/self/exe", raw, sizeof raw - 1);
    if (k <= 0) return -1;
    raw[k] = '\0';
#else
    return -1;                       /* unknown platform: fall back to XDG */
#endif
    if (!realpath(raw, real)) return -1;   /* resolve any symlink to the binary */
    slash = strrchr(real, '/');
    if (!slash) return -1;
    *slash = '\0';
    snprintf(out, n, "%s", real);
    return 0;
}

/* <dir of the binary>/data, if it exists. 0 on success. Exposed so the CLI can
 * explain which default is in play. */
int sesame__exe_data_dir(char *out, size_t n)
{
    char exe[4096], cand[4096];
    if (exe_dir(exe, sizeof exe) != 0) return -1;
    snprintf(cand, sizeof cand, "%s/data", exe);
    if (!is_dir(cand)) return -1;
    snprintf(out, n, "%s", cand);
    return 0;
}

/* Where fetch writes and where the active links live. */
const char *sesame_store_dir(char *out, size_t n)
{
    const char *e;
    char exe[4096], cand[4096];

    if ((e = getenv("SESAME_INDEX_DIR")) && *e) {
        snprintf(out, n, "%s", e);
    } else if (exe_dir(exe, sizeof exe) == 0 &&
               (snprintf(cand, sizeof cand, "%s/data", exe), is_dir(cand))) {
        /* A checkout: data/ sits next to the binary. Found from any cwd. An
         * installed binary has no data/ beside it, so this simply does not
         * match and we fall through to the shared store. */
        snprintf(out, n, "%s", cand);
    } else {
        /* The store every tool in the suite shares. Passing NULL for the
         * per-tool variable because SESAME_INDEX_DIR was already handled above
         * -- it stays first in precedence, as it always was. */
        yame_assets_root(NULL, NULL, out, n);
    }
    return out;
}


int sesame_index_locate(const char *platform, char *out, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];

    if (!platform) return -1;

    (void)reg;
    sesame_store_dir(dir, sizeof dir);
    /* Mirrors the remote under the shared key:
     * <store>/InfiniumAnnotation/<platform>/<platform>.ordering.tsv.gz */
    snprintf(out, n, "%s/%s/%s/%s.ordering.tsv.gz",
             dir, SESAME_IA_SUB, platform, platform);
    if (is_file(out)) return 0;

    snprintf(out, n, "./%s.ordering.tsv.gz", platform);   /* cwd convenience */
    if (is_file(out)) return 0;

    out[0] = '\0';
    return -1;
}

int sesame_genome_locate(const char *genome, const char *file, char *out, size_t n)
{
    char dir[4096];
    if (!genome || !file) return -1;
    sesame_store_dir(dir, sizeof dir);
    /* Mirrors the remote: <store>/genomes/<genome>/<file>. The directory is
     * the upstream repo name (zhou-lab/genomes), which is the rule the whole
     * shared layout follows; it was singular "genome/" before the move. */
    snprintf(out, n, "%s/%s/%s/%s", dir, SESAME_GENOME_SUB, genome, file);
    if (is_file(out)) return 0;

    snprintf(out, n, "./%s", file);   /* cwd convenience */
    if (is_file(out)) return 0;

    out[0] = '\0';
    return -1;
}

/* Fills msg with the "here is how to fix it" text used when no index is found.
 * Deliberately verbose: this is the error a first-time user will hit. */
void sesame_index_missing_help(const char *platform, char *msg, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no index found for platform %s\n"
        "  searched:\n"
        "    %s/" SESAME_IA_SUB "/%s/%s\n"
        "    ./%s.ordering.tsv.gz\n"
        "  fix, any of:\n"
        "    sesame fetch                 download the pinned tag (%s)\n"
        "    sesame betas --index <path> ...\n"
        "    export SESAME_INDEX_DIR=<dir>",
        platform,
        dir, platform, reg ? reg->ordering : "<platform>.ordering.tsv.gz",
        platform,
        SESAME_DEFAULT_TAG);
}

/* The searched path in that message must match sesame_index_locate(). */

/* ---------------------------------------------------------------- fetching
 *
 * A thin layer over YAME's shared downloader (YAME/src/assets.h): this file
 * supplies the trust anchors out of registry.h and the progress rendering,
 * and yame_assets_fetch_subtree() does the transfer, the verification, and
 * the pin check that stops one tool overwriting a directory another tool
 * populated from a different upstream tag.
 *
 * What used to live here -- a curl handle, a SHA256SUMS parser, a
 * download-verify-rename dance -- was the third copy of the same code in the
 * suite. Two things it did wrong are fixed by the move: the temp file was a
 * fixed "<dest>.part" opened without O_EXCL, so two concurrent fetches into one
 * store wrote the same file (kycg had already hit this and moved to a per-pid
 * name), and manifest entry names were joined into a path with no traversal
 * check, so a "../.." entry would have written outside the store.
 *
 * There is also no longer a compile-time split on libcurl: YAME reports at run
 * time whether it has it, which is one code path instead of two.
 */

/* Human-readable byte size: 913 -> "913 B", 156784 -> "153 KB", etc. */
static void human_size(double b, char *out, size_t n)
{
    static const char *u[] = { "B", "KB", "MB", "GB" };
    int i = 0;
    while (b >= 1024.0 && i < 3) { b /= 1024.0; i++; }
    if (i == 0) snprintf(out, n, "%.0f %s", b, u[i]);
    else        snprintf(out, n, "%.1f %s", b, u[i]);
}

/* Live download progress, rendered in place on a TTY only. */
typedef struct {
    char label[128];
    int  last_pct;
    int  tty;
    int  ndl;                 /* files actually transferred */
} progress_t;

static void on_begin(void *ud, const char *name, uint64_t total)
{
    progress_t *pr = (progress_t *)ud;
    (void)total;
    snprintf(pr->label, sizeof pr->label, "%s", name ? name : "");
    pr->last_pct = -1;
}

static void on_progress(void *ud, uint64_t now, uint64_t total)
{
    progress_t *pr = (progress_t *)ud;
    int pct, fill, i;
    char bar[3 * 24 + 1], cur[24], tot[24];
    const int W = 24;

    if (!pr->tty || total == 0) return;
    pct = (int)((now * 100) / total);
    if (pct == pr->last_pct) return;          /* redraw only on a percent change */
    pr->last_pct = pct;

    fill = pct * W / 100;
    bar[0] = '\0';
    for (i = 0; i < W; i++) strcat(bar, i < fill ? "\xe2\x96\x88"   /* U+2588 */
                                                 : "\xe2\x96\x91"); /* U+2591 */
    human_size((double)now, cur, sizeof cur);
    human_size((double)total, tot, sizeof tot);
    fprintf(stderr, "\r  %-24s %s %3d%%  %s / %s\033[K",
            pr->label, bar, pct, cur, tot);
    fflush(stderr);
}

static void on_done(void *ud, const char *name, uint64_t bytes, int ok)
{
    progress_t *pr = (progress_t *)ud;
    char sz[24];

    if (pr->tty) fprintf(stderr, "\r\033[K");    /* clear the in-place bar */
    if (!ok) return;
    pr->ndl++;
    human_size((double)bytes, sz, sizeof sz);
    fprintf(stderr, "  \xe2\x9c\x93 %-24s %s\n", name ? name : "", sz);
}

/* Fetch a SHA256SUMS-anchored subtree into <store_sub>/, mirroring the remote.
 * `noun` labels the summary line. Files already present and correct print
 * nothing, so a fetch with nothing to do stays quiet. */
static int fetch_subtree(const char *base, const char *tag,
                         const char *remote_sub, const char *store_sub,
                         const char *anchor_sha, const char *noun, int force,
                         sesame_err_t *err)
{
    yame_fetch_opt_t opt;
    progress_t pr;
    char *yerr = NULL;
    int rc;

    memset(&pr, 0, sizeof pr);
    pr.last_pct = -1;
    pr.tty = isatty(STDERR_FILENO);

    memset(&opt, 0, sizeof opt);
    opt.force = force;
    opt.on_begin = on_begin;
    opt.on_progress = on_progress;
    opt.on_done = on_done;
    opt.ud = &pr;

    rc = yame_assets_fetch_subtree(base, tag, remote_sub, store_sub,
                                   anchor_sha, &opt, &yerr);
    if (rc != 0) {
        int code = sesame__fail(err, SESAME_ERR_IO, "%s: %s", noun,
                                yerr ? yerr : "fetch failed");
        free(yerr);
        return code;
    }
    free(yerr);

    if (pr.ndl == 0)
        fprintf(stderr, "sesame: %s up to date (%s)\n", noun, tag);
    else
        fprintf(stderr, "sesame: %s ready (%s, %d file%s downloaded)\n",
                noun, tag, pr.ndl, pr.ndl == 1 ? "" : "s");
    return SESAME_OK;
}

/* Fetch one platform at the pinned tag into
 * <store>/InfiniumAnnotation/<platform>/. out_path receives the ordering
 * table's path. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char store[4096], pdir[4096];
    int rc;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (out_path && out_n) out_path[0] = '\0';
    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);
    if (!reg->sums_sha256)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "%s is not published at tag %s yet", platform, SESAME_DEFAULT_TAG);

    sesame_store_dir(store, sizeof store);
    if (!yame_assets_root_writable(store))
        return sesame__fail(err, SESAME_ERR_IO,
            "store %s is not writable; set SESAME_INDEX_DIR or YAME_DATA_HOME "
            "to a directory you own", store);

    snprintf(pdir, sizeof pdir, "%s/%s/%s", store, SESAME_IA_SUB, platform);
    rc = fetch_subtree(SESAME_BASE_URL, SESAME_DEFAULT_TAG, platform, pdir,
                       reg->sums_sha256, platform, force, err);
    if (rc != SESAME_OK) return rc;

    if (out_path && out_n) {
        snprintf(out_path, out_n, "%s/%s", pdir, reg->ordering);
        if (!is_file(out_path)) out_path[0] = '\0';
    }
    return SESAME_OK;
}

static const sesame_genome_reg_t *sesame__genome_reg_for(const char *genome)
{
    if (!genome) return NULL;
    for (const sesame_genome_reg_t *g = SESAME_GENOME_REGISTRY; g->genome; g++)
        if (strcmp(g->genome, genome) == 0) return g;
    return NULL;
}

/* Fetch one genome's genome-level annotation into <store>/genomes/<genome>/. */
int sesame_fetch_genome(const char *genome, int force, sesame_err_t *err)
{
    const sesame_genome_reg_t *g = sesame__genome_reg_for(genome);
    char store[4096], gdir[4096];

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!g)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown genome '%s' (known: hg38, mm10, mm39)", genome);
    if (!g->sums_sha256)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "genome %s is not published at tag %s yet", genome, SESAME_GENOME_TAG);

    sesame_store_dir(store, sizeof store);
    if (!yame_assets_root_writable(store))
        return sesame__fail(err, SESAME_ERR_IO,
            "store %s is not writable; set SESAME_INDEX_DIR or YAME_DATA_HOME "
            "to a directory you own", store);

    snprintf(gdir, sizeof gdir, "%s/%s/%s", store, SESAME_GENOME_SUB, genome);
    return fetch_subtree(SESAME_GENOME_BASE_URL, SESAME_GENOME_TAG, genome, gdir,
                         g->sums_sha256, genome, force, err);
}

/* Fetch every platform that is published at the pinned tag. */
int sesame_fetch_all(int force, sesame_err_t *err)
{
    const sesame_reg_t *r;
    char path[4096];
    int any = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    for (r = SESAME_REGISTRY; r->platform; r++) {
        if (!r->sums_sha256) continue;    /* not published at the pinned tag */
        if (sesame_fetch_index(r->platform, force, path, sizeof path, err)
                != SESAME_OK)
            return err ? err->code : SESAME_ERR_IO;
        any = 1;
    }
    if (!any)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "no platforms are published at tag %s", SESAME_DEFAULT_TAG);
    return SESAME_OK;
}
