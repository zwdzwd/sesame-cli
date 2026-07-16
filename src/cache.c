/* cache.c -- index discovery, platform detection, and fetching.
 *
 * Resolution order (explicit always wins, and it is the same order for the
 * library, the CLI and any future binding -- one cache, one rule):
 *
 *   1. --index <path>                              (caller-supplied)
 *   2. $SESAMEC_INDEX_DIR/<platform>.ordering.tsv.gz
 *   3. ./<platform>.ordering.tsv.gz
 *   4. <cache>/<platform>.ordering.tsv.gz
 *        cache = $SESAMEC_CACHE
 *              | $XDG_CACHE_HOME/sesamec
 *              | ~/Library/Caches/sesamec   (macOS)
 *              | ~/.cache/sesamec
 *
 * sesamec NEVER prompts and NEVER downloads implicitly. If the index is
 * missing, the caller gets an error naming the exact command to run. An
 * interactive prompt would hang forever in a Nextflow job or a Docker build --
 * silently, on question one, across a whole run. `sesamec fetch` is the only
 * download path, so behaviour is identical whether or not a TTY is attached.
 *
 * SPDX-License-Identifier: MIT
 */
#include "sesame.h"
#include "internal.h"
#include "registry.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef SESAMEC_HAVE_CURL
#include <curl/curl.h>
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

const char *sesame_platform_from_beads(int32_t beads)
{
    for (const sesame_reg_t *r = SESAME_REGISTRY; r->platform; r++)
        if (r->beads && r->beads == beads) return r->platform;
    return NULL;
}

const char *sesame_cache_dir(char *out, size_t n)
{
    const char *e;
    if ((e = getenv("SESAMEC_CACHE")) && *e) {
        snprintf(out, n, "%s", e);
    } else if ((e = getenv("XDG_CACHE_HOME")) && *e) {
        snprintf(out, n, "%s/sesamec", e);
    } else if ((e = getenv("HOME")) && *e) {
#ifdef __APPLE__
        snprintf(out, n, "%s/Library/Caches/sesamec", e);
#else
        snprintf(out, n, "%s/.cache/sesamec", e);
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
    const char *e;
    char dir[4096];

    if (!platform) return -1;

    if ((e = getenv("SESAMEC_INDEX_DIR")) && *e) {
        snprintf(out, n, "%s/%s.ordering.tsv.gz", e, platform);
        if (is_file(out)) return 0;
    }
    snprintf(out, n, "./%s.ordering.tsv.gz", platform);
    if (is_file(out)) return 0;

    /* data/ next to the cwd -- convenient in a source checkout. */
    snprintf(out, n, "./data/%s.ordering.tsv.gz", platform);
    if (is_file(out)) return 0;

    sesame_cache_dir(dir, sizeof dir);
    snprintf(out, n, "%s/%s.ordering.tsv.gz", dir, platform);
    if (is_file(out)) return 0;

    out[0] = '\0';
    return -1;
}

/* Fills msg with the "here is how to fix it" text used when no index is found.
 * Deliberately verbose: this is the error a first-time user will hit. */
void sesame_index_missing_help(const char *platform, char *msg, size_t n)
{
    const char *idxdir = getenv("SESAMEC_INDEX_DIR");
    char dir[4096];
    sesame_cache_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no index found for platform %s\n"
        "  searched:\n"
        "    $SESAMEC_INDEX_DIR  %s\n"
        "    ./%s.ordering.tsv.gz\n"
        "    ./data/%s.ordering.tsv.gz\n"
        "    %s/%s.ordering.tsv.gz\n"
        "  fix, any of:\n"
        "    sesamec fetch %s            download it to the cache\n"
        "    sesamec betas --index <path> ...\n"
        "    export SESAMEC_INDEX_DIR=<dir>",
        platform,
        (idxdir && *idxdir) ? idxdir : "(unset)",
        platform, platform, dir, platform, platform);
}

#ifdef SESAMEC_HAVE_CURL
static size_t wr(void *p, size_t sz, size_t nm, void *ud)
{
    return fwrite(p, sz, nm, (FILE *)ud);
}

int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096], dest[4096], tmp[4096], url[4096], got[65];
    CURL *cu = NULL;
    FILE *f = NULL;
    CURLcode rc;
    long code = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);

    sesame_cache_dir(dir, sizeof dir);
    mkdirs(dir);
    snprintf(dest, sizeof dest, "%s/%s", dir, reg->file);
    snprintf(tmp,  sizeof tmp,  "%s/.%s.part", dir, reg->file);
    snprintf(url,  sizeof url,  "%s/%s", SESAME_INDEX_BASE_URL, reg->file);

    if (!force && is_file(dest)) {
        if (sesame__sha256_file(dest, got) == 0 &&
            strcmp(got, reg->sha256) == 0) {
            snprintf(out_path, out_n, "%s", dest);
            return SESAME_OK;   /* already have it, verified */
        }
        /* present but wrong: fall through and refetch */
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
    curl_easy_setopt(cu, CURLOPT_USERAGENT, "sesamec/" SESAME_INDEX_VERSION);
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
    if (sesame__sha256_file(tmp, got) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "cannot hash %s", tmp);
    }
    if (strcmp(got, reg->sha256) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_FORMAT,
            "sha256 mismatch for %s\n  expected %s\n  got      %s\n"
            "  refusing to install a index that does not match the pinned digest",
            reg->file, reg->sha256, got);
    }
    if (rename(tmp, dest) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "cannot install %s", dest);
    }
    snprintf(out_path, out_n, "%s", dest);
    return SESAME_OK;
}
#else
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    (void)force; (void)out_path; (void)out_n;
    return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
        "this build has no network support; download %s/%s manually",
        SESAME_INDEX_BASE_URL, reg ? reg->file : "<platform>.ordering.tsv.gz");
}
#endif
