/* cache.c -- index discovery and fetching.
 *
 * The store mirrors the remote exactly, one version, one subfolder per platform:
 *
 *   <store>/MSA/SHA256SUMS         byte-identical copy of the remote
 *   <store>/MSA/MSA.ordering.tsv.gz
 *   <store>/MSA/MSA.hg38.mask.cm
 *   <store>/MSA/KYCG/...           room to grow, no cross-platform conflicts
 *
 * Mirroring means the store is self-describing and hand-verifiable per platform
 * (cd <store>/MSA && shasum -a 256 -c SHA256SUMS), with no local bookkeeping.
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
 * Store location (fetch writes here). ONE sesame variable, deliberately:
 *
 *   1. $SESAME_INDEX_DIR                explicit
 *   2. <dir of the binary>/data         a checkout: found from any cwd
 *   3. $XDG_CACHE_HOME/sesame           standard fallback (not ours)
 *   4. ~/Library/Caches/sesame (macOS) | ~/.cache/sesame
 *
 * Resolution order (explicit always wins; same rule for library, CLI and any
 * future binding):
 *
 *   1. --index <path>                                 (caller-supplied)
 *   2. <store>/<platform>/<platform>.ordering.tsv.gz
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef SESAME_HAVE_CURL
#include <curl/curl.h>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

static int is_file(const char *p)
{
    struct stat st;
    return p && *p && stat(p, &st) == 0 && S_ISREG(st.st_mode);
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
    struct stat st;
    return p && *p && stat(p, &st) == 0 && S_ISDIR(st.st_mode);
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
         * match and we fall through to the XDG store. */
        snprintf(out, n, "%s", cand);
    } else if ((e = getenv("XDG_CACHE_HOME")) && *e) {
        snprintf(out, n, "%s/sesame", e);
    } else if ((e = getenv("HOME")) && *e) {
#ifdef __APPLE__
        snprintf(out, n, "%s/Library/Caches/sesame", e);
#else
        snprintf(out, n, "%s/.cache/sesame", e);
#endif
    } else {
        snprintf(out, n, ".");
    }
    return out;
}

/* mkdir -p */
static int mkdirs(const char *path)
{
    char tmp[4096];
    size_t len;
    snprintf(tmp, sizeof tmp, "%s", path);
    len = strlen(tmp);
    if (len && tmp[len-1] == '/') tmp[len-1] = '\0';
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0) { /* EEXIST is fine */ }
            *p = '/';
        }
    }
    return (mkdir(tmp, 0755) == 0 || 1) ? 0 : -1;
}


int sesame_index_locate(const char *platform, char *out, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];

    if (!platform) return -1;

    (void)reg;
    sesame_store_dir(dir, sizeof dir);
    /* The store mirrors the remote: <store>/<platform>/<platform>.ordering.tsv.gz */
    snprintf(out, n, "%s/%s/%s.ordering.tsv.gz", dir, platform, platform);
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
    /* Mirrors the remote: <store>/genome/<genome>/<file> */
    snprintf(out, n, "%s/genome/%s/%s", dir, genome, file);
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
        "    %s/%s/%s\n"
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

#ifdef SESAME_HAVE_CURL
static size_t wr(void *p, size_t sz, size_t nm, void *ud)
{
    return fwrite(p, sz, nm, (FILE *)ud);
}

/* --- SHA256SUMS: the per-file digests for a tag ---------------------------
 *
 * Lines are the coreutils/shasum format: "<64 hex><2 spaces><filename>".
 * Fetching this one small file first is what lets verification and
 * de-duplication work for ANY tag rather than only the tag compiled in. */
#define SUMS_MAX 64
typedef struct { char file[192]; char sha[65]; } sums_t;

static int sums_load(const char *path, sums_t *out, int max)
{
    FILE *f = fopen(path, "rb");
    char line[512];
    int n = 0;

    if (!f) return -1;
    while (n < max && fgets(line, sizeof line, f)) {
        char *p = line, *name;
        size_t k;
        while (*p == ' ' || *p == '\t') p++;
        if (strlen(p) < 66) continue;
        for (k = 0; k < 64; k++)
            if (!isxdigit((unsigned char)p[k])) break;
        if (k != 64) continue;                      /* not a digest line */
        name = p + 64;
        while (*name == ' ' || *name == '*' || *name == '\t') name++;
        name[strcspn(name, "\r\n")] = '\0';
        if (!*name) continue;
        memcpy(out[n].sha, p, 64);
        out[n].sha[64] = '\0';
        snprintf(out[n].file, sizeof out[n].file, "%s", name);
        n++;
    }
    fclose(f);
    return n;
}


/* Fetch <file> at <tag> into <store>/<file>.
 *
 * If it is already on disk with the right digest, nothing is downloaded. That
 * one check is the whole de-duplication story. want_sha may be NULL only for the
 * SHA256SUMS of a tag this build does not pin. */
/* Download url into dest, verifying want_sha (may be NULL for the anchor-less
 * case). If dest already matches want_sha and !force, skips the network. */
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
typedef struct { const char *label; int last_pct; int tty; } progress_t;

static int on_xfer(void *p, curl_off_t dltotal, curl_off_t dlnow,
                   curl_off_t ultotal, curl_off_t ulnow)
{
    progress_t *pr = (progress_t *)p;
    int pct, fill, i;
    char bar[3 * 24 + 1], cur[24], tot[24];
    const int W = 24;

    (void)ultotal; (void)ulnow;
    if (!pr->tty || dltotal <= 0) return 0;
    pct = (int)((dlnow * 100) / dltotal);
    if (pct == pr->last_pct) return 0;       /* redraw only on a percent change */
    pr->last_pct = pct;

    fill = pct * W / 100;
    bar[0] = '\0';
    for (i = 0; i < W; i++) strcat(bar, i < fill ? "\xe2\x96\x88"   /* U+2588 */
                                                 : "\xe2\x96\x91"); /* U+2591 */
    human_size((double)dlnow, cur, sizeof cur);
    human_size((double)dltotal, tot, sizeof tot);
    fprintf(stderr, "\r  %-24s %s %3d%%  %s / %s\033[K",
            pr->label, bar, pct, cur, tot);
    fflush(stderr);
    return 0;
}

/* label == NULL: fetch silently (used for the tiny SHA256SUMS). Otherwise show a
 * progress bar on a TTY while downloading, then a one-line "got it" for files
 * actually downloaded -- files already present and correct print nothing, so
 * a fetch that has nothing to do stays quiet. *downloaded (if non-NULL) reports
 * whether a transfer happened. */
static int download_verify(const char *url, const char *want_sha,
                           const char *dest, const char *label,
                           int force, int *downloaded, sesame_err_t *err)
{
    char tmp[4096], got[65], sz[24];
    CURL *cu;
    FILE *f;
    CURLcode rc;
    long code = 0;
    progress_t pr;
    struct stat st;

    if (downloaded) *downloaded = 0;

    if (!force && want_sha && is_file(dest) &&
        sesame__sha256_file(dest, got) == 0 && strcmp(got, want_sha) == 0)
        return SESAME_OK;                 /* already correct; nothing to do */

    snprintf(tmp, sizeof tmp, "%s.part", dest);
    if (!(f = fopen(tmp, "wb")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot write %s", tmp);

    pr.label = label; pr.last_pct = -1;
    pr.tty = (label != NULL) && isatty(STDERR_FILENO);

    curl_global_init(CURL_GLOBAL_DEFAULT);
    if (!(cu = curl_easy_init())) {
        fclose(f); remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "curl init failed");
    }
    curl_easy_setopt(cu, CURLOPT_URL, url);
    curl_easy_setopt(cu, CURLOPT_WRITEFUNCTION, wr);
    curl_easy_setopt(cu, CURLOPT_WRITEDATA, f);
    curl_easy_setopt(cu, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(cu, CURLOPT_FAILONERROR, 1L);
    curl_easy_setopt(cu, CURLOPT_NOPROGRESS, pr.tty ? 0L : 1L);
    curl_easy_setopt(cu, CURLOPT_XFERINFOFUNCTION, on_xfer);
    curl_easy_setopt(cu, CURLOPT_XFERINFODATA, &pr);
    curl_easy_setopt(cu, CURLOPT_USERAGENT, "sesame-cli");
    rc = curl_easy_perform(cu);
    curl_easy_getinfo(cu, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(cu);
    fclose(f);
    if (pr.tty) fprintf(stderr, "\r\033[K");   /* clear the in-place bar */

    if (rc != CURLE_OK) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "download failed (%s): %s",
                            curl_easy_strerror(rc), url);
    }
    if (want_sha) {
        if (sesame__sha256_file(tmp, got) != 0) {
            remove(tmp);
            return sesame__fail(err, SESAME_ERR_IO, "cannot hash %s", tmp);
        }
        if (strcmp(got, want_sha) != 0) {
            remove(tmp);
            return sesame__fail(err, SESAME_ERR_FORMAT,
                "sha256 mismatch for %s\n  expected %s\n  got      %s\n"
                "  refusing to install a file that does not match its pinned digest",
                dest, want_sha, got);
        }
    }
    if (rename(tmp, dest) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "cannot install %s", dest);
    }
    if (downloaded) *downloaded = 1;
    if (label) {
        if (stat(dest, &st) == 0) { human_size((double)st.st_size, sz, sizeof sz); }
        else sz[0] = '\0';
        fprintf(stderr, "  \xe2\x9c\x93 %-24s %s\n", label, sz);
    }
    return SESAME_OK;
}

/* mkdir -p the parent directory of a file path. */
static void mkdir_parent(const char *path)
{
    char tmp[4096];
    char *slash;
    snprintf(tmp, sizeof tmp, "%s", path);
    slash = strrchr(tmp, '/');
    if (slash) { *slash = '\0'; mkdirs(tmp); }
}

/* Fetch a SHA256SUMS-anchored subtree from <base>/<tag>/<remote_sub>/ into
 * <store_sub>/, mirroring the remote exactly.
 *
 * Repo layout: <base>/<tag>/<remote_sub>/{SHA256SUMS, <files...>}, where files
 * may themselves be nested (e.g. KYCG/foo.cm). We pull SHA256SUMS first, verify
 * it against the digest compiled into this build (a hard trust anchor) and KEEP
 * it -- so <store_sub>/SHA256SUMS is a byte-identical copy of the remote and
 * `shasum -a 256 -c SHA256SUMS` verifies the subtree by hand. Then every file it
 * lists is fetched, a matching one skipped. `noun` labels the summary line. If
 * match_file is non-NULL, the on-disk path of the file whose name equals it is
 * written to out_match. Mirroring means no naming conflicts, room to grow
 * subfolders, and no local merge logic. Shared by the platform index and the
 * genome-annotation fetch. */
static int fetch_subtree(const char *base, const char *tag,
                         const char *remote_sub, const char *store_sub,
                         const char *anchor_sha, const char *noun, int force,
                         const char *match_file, char *out_match, size_t out_n,
                         sesame_err_t *err)
{
    char url[4096], sums_path[4096];
    sums_t sums[SUMS_MAX];
    int nsums, i, ndl = 0, dl;

    mkdirs(store_sub);
    if (out_match && out_n) out_match[0] = '\0';

    /* SHA256SUMS: verified against the anchor, and kept in place. */
    snprintf(url, sizeof url, "%s/%s/%s/%s", base, tag, remote_sub, SESAME_SUMS_FILE);
    snprintf(sums_path, sizeof sums_path, "%s/%s", store_sub, SESAME_SUMS_FILE);
    if (download_verify(url, anchor_sha, sums_path, NULL, 1, NULL, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;

    nsums = sums_load(sums_path, sums, SUMS_MAX);
    if (nsums <= 0)
        return sesame__fail(err, SESAME_ERR_FORMAT,
                            "empty or unreadable SHA256SUMS for %s", noun);

    for (i = 0; i < nsums; i++) {
        char dest[4096];
        snprintf(dest, sizeof dest, "%s/%s", store_sub, sums[i].file);
        snprintf(url,  sizeof url,  "%s/%s/%s/%s", base, tag, remote_sub, sums[i].file);
        mkdir_parent(dest);                /* nested files, e.g. KYCG/foo.cm */
        if (download_verify(url, sums[i].sha, dest, sums[i].file, force, &dl, err)
                != SESAME_OK)
            return err ? err->code : SESAME_ERR_IO;
        ndl += dl;
        if (match_file && out_match && strcmp(sums[i].file, match_file) == 0)
            snprintf(out_match, out_n, "%s", dest);
    }

    /* One clear summary. When nothing was downloaded, say so plainly instead of
     * repeating a "(cached)" line for every file. */
    if (ndl == 0)
        fprintf(stderr, "sesame: %s up to date (%s, %d files)\n", noun, tag, nsums);
    else
        fprintf(stderr, "sesame: %s ready (%s, %d file%s downloaded)\n",
                noun, tag, ndl, ndl == 1 ? "" : "s");
    return SESAME_OK;
}

/* Fetch one platform at the pinned tag into <store>/<platform>/. out_path
 * receives the ordering table's path. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char store[4096], pdir[4096];

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);
    if (!reg->sums_sha256)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "%s is not published at tag %s yet", platform, SESAME_DEFAULT_TAG);

    sesame_store_dir(store, sizeof store);
    snprintf(pdir, sizeof pdir, "%s/%s", store, platform);
    return fetch_subtree(SESAME_BASE_URL, SESAME_DEFAULT_TAG, platform, pdir,
                         reg->sums_sha256, platform, force,
                         reg->ordering, out_path, out_n, err);
}

static const sesame_genome_reg_t *sesame__genome_reg_for(const char *genome)
{
    if (!genome) return NULL;
    for (const sesame_genome_reg_t *g = SESAME_GENOME_REGISTRY; g->genome; g++)
        if (strcmp(g->genome, genome) == 0) return g;
    return NULL;
}

/* Fetch one genome's genome-level annotation into <store>/genome/<genome>/. */
int sesame_fetch_genome(const char *genome, int force, sesame_err_t *err)
{
    const sesame_genome_reg_t *g = sesame__genome_reg_for(genome);
    char store[4096], gdir[4096];

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!g)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown genome '%s' (known: hg38)", genome);
    if (!g->sums_sha256)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "genome %s is not published at tag %s yet", genome, SESAME_GENOME_TAG);

    sesame_store_dir(store, sizeof store);
    snprintf(gdir, sizeof gdir, "%s/genome/%s", store, genome);
    return fetch_subtree(SESAME_GENOME_BASE_URL, SESAME_GENOME_TAG, genome, gdir,
                         g->sums_sha256, genome, force, NULL, NULL, 0, err);
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
#else
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    (void)force; (void)out_path; (void)out_n;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support; download %s/%s/%s/ manually",
        SESAME_BASE_URL, SESAME_DEFAULT_TAG, reg ? reg->platform : "<platform>");
}

int sesame_fetch_all(int force, sesame_err_t *err)
{
    (void)force;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support");
}

int sesame_fetch_genome(const char *genome, int force, sesame_err_t *err)
{
    (void)force;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support; download %s/%s/%s/ manually",
        SESAME_GENOME_BASE_URL, SESAME_GENOME_TAG, genome ? genome : "<genome>");
}
#endif
