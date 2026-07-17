/* cache.c -- index discovery and fetching.
 *
 * The store holds exactly ONE version of the annotation, flat:
 *
 *   <store>/EPIC.ordering.tsv.gz
 *   <store>/EPICv2.ordering.tsv.gz
 *   <store>/SHA256SUMS     per-file digests, with the tag in a leading
 *                          "# <tag>" comment -- so there is no separate TAG
 *                          file that could drift out of sync with them
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

/* Which tag the store holds, read from the leading "# <tag>" comment of
 * <store>/SHA256SUMS. 0 on success.
 *
 * The tag lives in that comment rather than a third column because `shasum -c`
 * parses everything after the two spaces as the filename -- a third column
 * would break hand-verification with standard tools. A comment keeps the file
 * both self-describing and checkable, and means there is no separate TAG file
 * to drift out of sync with the digests it belongs to. */
int sesame_index_active(const char *store, char *tag_out, size_t n)
{
    char path[4096], buf[256];
    FILE *f;
    char *p;

    snprintf(path, sizeof path, "%s/%s", store, SESAME_SUMS_FILE);
    if (!(f = fopen(path, "rb"))) return -1;
    if (!fgets(buf, sizeof buf, f)) { fclose(f); return -1; }
    fclose(f);
    buf[strcspn(buf, "\r\n")] = '\0';
    p = buf;
    while (*p == ' ' || *p == '\t') p++;
    if (*p != '#') return -1;
    p++;
    while (*p == ' ' || *p == '\t') p++;
    if (!*p) return -1;
    snprintf(tag_out, n, "%s", p);
    return 0;
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
        dir, reg ? reg->file : "<platform>.ordering.tsv.gz",
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

static const char *sums_lookup(const sums_t *s, int n, const char *file)
{
    int i;
    for (i = 0; i < n; i++)
        if (strcmp(s[i].file, file) == 0) return s[i].sha;
    return NULL;
}

/* Fetch <file> at <tag> into <store>/<file>.
 *
 * If it is already on disk with the right digest, nothing is downloaded. That
 * one check is the whole de-duplication story. want_sha may be NULL only for the
 * SHA256SUMS of a tag this build does not pin. */
int sesame_fetch(const char *tag, const char *file, const char *want_sha,
                 int force, char *out_path, size_t out_n, sesame_err_t *err)
{
    char dir[4096], dest[4096], tmp[4096], url[4096], got[65];
    CURL *cu = NULL;
    FILE *f = NULL;
    CURLcode rc;
    long code = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!tag || !file)
        return sesame__fail(err, SESAME_ERR_IO, "fetch needs a tag and a file");

    sesame_store_dir(dir, sizeof dir);
    mkdirs(dir);
    snprintf(dest, sizeof dest, "%s/%s", dir, file);
    snprintf(tmp,  sizeof tmp,  "%s/.%s.part", dir, file);
    snprintf(url,  sizeof url,  "%s/%s/%s", SESAME_BASE_URL, tag, file);

    /* Already correct -> the remote copy is byte-identical; skip the network. */
    if (!force && want_sha && is_file(dest) &&
        sesame__sha256_file(dest, got) == 0 && strcmp(got, want_sha) == 0) {
        snprintf(out_path, out_n, "%s", dest);
        return SESAME_OK;
    }

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

    /* Verify before publishing into the cache. A wrong index silently changes
     * betas, so a digest mismatch is fatal, never a warning. */
    if (want_sha) {
        if (sesame__sha256_file(tmp, got) != 0) {
            remove(tmp);
            return sesame__fail(err, SESAME_ERR_IO, "cannot hash %s", tmp);
        }
        if (strcmp(got, want_sha) != 0) {
            remove(tmp);
            return sesame__fail(err, SESAME_ERR_FORMAT,
                "sha256 mismatch for %s\n  expected %s\n  got      %s\n"
                "  refusing to install an index that does not match the pinned digest",
                file, want_sha, got);
        }
    }
    if (rename(tmp, dest) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "cannot install %s", dest);
    }
    snprintf(out_path, out_n, "%s", dest);
    return SESAME_OK;
}

/* Convenience: one platform at the pinned tag, digest via SHA256SUMS. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char sums_path[4096];
    sums_t sums[SUMS_MAX];
    const char *want;
    int nsums;

    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);
    if (sesame_fetch(SESAME_DEFAULT_TAG, SESAME_SUMS_FILE, SESAME_SUMS_SHA256,
                     force, sums_path, sizeof sums_path, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;
    nsums = sums_load(sums_path, sums, SUMS_MAX);
    want = (nsums > 0) ? sums_lookup(sums, nsums, reg->file) : NULL;
    if (!want)
        return sesame__fail(err, SESAME_ERR_FORMAT,
                            "no digest for %s in %s", reg->file, sums_path);
    return sesame_fetch(SESAME_DEFAULT_TAG, reg->file, want, force,
                        out_path, out_n, err);
}

/* Fetch every platform at `tag` and make it current.
 *
 * SHA256SUMS is fetched first: for the pinned tag it is verified against the
 * digest compiled into this build (a hard trust anchor); for any other tag it
 * is trusted on first use, which still gives per-file integrity within that tag
 * and still enables de-duplication. Files whose content already exists under
 * another tag are hardlinked rather than downloaded. */
int sesame_fetch_tag(const char *tag, int force, sesame_err_t *err)
{
    char dir[4096], path[4096], sums_path[4096];
    sums_t sums[SUMS_MAX];
    const sesame_reg_t *r;
    int pinned, nsums;

    if (!tag || !*tag) tag = SESAME_DEFAULT_TAG;
    pinned = (strcmp(tag, SESAME_DEFAULT_TAG) == 0);
    sesame_store_dir(dir, sizeof dir);

    if (sesame_fetch(tag, SESAME_SUMS_FILE,
                     pinned ? SESAME_SUMS_SHA256 : NULL,
                     force, sums_path, sizeof sums_path, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;

    nsums = sums_load(sums_path, sums, SUMS_MAX);
    if (nsums <= 0)
        return sesame__fail(err, SESAME_ERR_FORMAT,
                            "%s has no usable digest lines", sums_path);

    for (r = SESAME_REGISTRY; r->platform; r++) {
        const char *want = sums_lookup(sums, nsums, r->file);
        if (!want)
            return sesame__fail(err, SESAME_ERR_FORMAT,
                "%s lists no digest for %s -- tag %s does not carry this platform",
                SESAME_SUMS_FILE, r->file, tag);
        if (sesame_fetch(tag, r->file, want, force,
                         path, sizeof path, err) != SESAME_OK)
            return err ? err->code : SESAME_ERR_IO;
    }
    return SESAME_OK;   /* the tag is recorded in SHA256SUMS itself */
}
#else
int sesame_fetch(const char *tag, const char *file, const char *want_sha,
                 int force, char *out_path, size_t out_n, sesame_err_t *err)
{
    (void)want_sha; (void)force; (void)out_path; (void)out_n;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support; download %s/%s/%s manually",
        SESAME_BASE_URL, tag ? tag : "<tag>", file ? file : "<file>");
}

int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    (void)force; (void)out_path; (void)out_n;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support; download %s/%s/%s manually",
        SESAME_BASE_URL, reg ? reg->tag : "<tag>",
        reg ? reg->file : "<platform>.ordering.tsv.gz");
}
#endif
