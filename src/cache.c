/* cache.c -- index discovery and fetching.
 *
 * The store holds exactly ONE version of the annotation, flat:
 *
 *   <store>/EPIC.ordering.tsv.gz
 *   <store>/EPICv2.ordering.tsv.gz
 * (SHA256SUMS is fetched per platform, verified, parsed, and discarded -- it is
 * not kept in the store. The tag is fixed by this build.)
 *
 * No tag subdirectories, no `current` symlink, no local version management.
 * Most users only ever want one version; anyone who genuinely needs two can
 * point SESAME_INDEX_DIR at two directories, which is simpler than anything we
 * could build for them.
 *
 * De-duplication falls out of that: fetch pulls the tag's SHA256SUMS first (a
 * few hundred bytes), compares each digest against what is already on disk, and
 * downloads only the files that differ. Moving v1 -> v2 when only EPICv2 changed
 * therefore transfers only EPICv2. Upstream, git is content-addressed, so the
 * unchanged files cost nothing in the new tag either.
 *
 * SHA256SUMS for the pinned tag is verified against a digest compiled into this
 * build -- a hard trust anchor. For any other tag it is trusted on first use,
 * which still gives per-file integrity within that tag.
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
 *   1. --index <path>                        (caller-supplied)
 *   2. <store>/<platform>.ordering.tsv.gz
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
    snprintf(out, n, "%s/%s.ordering.tsv.gz", dir, platform);
    if (is_file(out)) return 0;

    snprintf(out, n, "./%s.ordering.tsv.gz", platform);
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
        "    %s/%s\n"
        "    ./%s.ordering.tsv.gz\n"
        "  fix, any of:\n"
        "    sesame fetch                 download the pinned tag (%s)\n"
        "    sesame betas --index <path> ...\n"
        "    export SESAME_INDEX_DIR=<dir>",
        platform,
        dir, reg ? reg->ordering : "<platform>.ordering.tsv.gz",
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
static int download_verify(const char *url, const char *want_sha,
                           const char *dest, int force, sesame_err_t *err)
{
    char tmp[4096], got[65];
    CURL *cu;
    FILE *f;
    CURLcode rc;
    long code = 0;

    if (!force && want_sha && is_file(dest) &&
        sesame__sha256_file(dest, got) == 0 && strcmp(got, want_sha) == 0)
        return SESAME_OK;                 /* already correct; no download */

    snprintf(tmp, sizeof tmp, "%s.part", dest);
    if (!(f = fopen(tmp, "wb")))
        return sesame__fail(err, SESAME_ERR_IO, "cannot write %s", tmp);

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
    curl_easy_setopt(cu, CURLOPT_NOPROGRESS, 0L);
    curl_easy_setopt(cu, CURLOPT_USERAGENT, "sesame-cli");
    rc = curl_easy_perform(cu);
    curl_easy_getinfo(cu, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(cu);
    fclose(f);

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
    return SESAME_OK;
}

/* Fetch one platform at the pinned tag into the flat store.
 *
 * Repo layout is <base>/<tag>/<platform>/{SHA256SUMS, <files...>}. We pull the
 * platform's SHA256SUMS first and verify it against the digest compiled into
 * this build -- a hard trust anchor -- then fetch every file it lists (the
 * ordering table plus, for later steps, the .cm mask). A file already on disk
 * with the right digest is skipped, which is the whole de-duplication story. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    const char *tag = SESAME_DEFAULT_TAG;
    char dir[4096], url[4096], sums_path[4096], ordering[4096] = "";
    sums_t sums[SUMS_MAX];
    int nsums, i;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);
    if (!reg->sums_sha256)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "%s is not published at tag %s yet", platform, tag);

    sesame_store_dir(dir, sizeof dir);
    mkdirs(dir);

    snprintf(url, sizeof url, "%s/%s/%s/%s",
             SESAME_BASE_URL, tag, platform, SESAME_SUMS_FILE);
    snprintf(sums_path, sizeof sums_path, "%s/.%s.SHA256SUMS.tmp", dir, platform);
    if (download_verify(url, reg->sums_sha256, sums_path, 1, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;

    nsums = sums_load(sums_path, sums, SUMS_MAX);
    remove(sums_path);                    /* transient; not part of the store */
    if (nsums <= 0)
        return sesame__fail(err, SESAME_ERR_FORMAT,
                            "empty or unreadable SHA256SUMS for %s", platform);

    for (i = 0; i < nsums; i++) {
        char dest[4096];
        snprintf(dest, sizeof dest, "%s/%s", dir, sums[i].file);
        snprintf(url,  sizeof url,  "%s/%s/%s/%s",
                 SESAME_BASE_URL, tag, platform, sums[i].file);
        if (download_verify(url, sums[i].sha, dest, force, err) != SESAME_OK)
            return err ? err->code : SESAME_ERR_IO;
        if (strcmp(sums[i].file, reg->ordering) == 0)
            snprintf(ordering, sizeof ordering, "%s", dest);
    }
    if (out_path) snprintf(out_path, out_n, "%s", ordering);
    return SESAME_OK;
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
#endif
