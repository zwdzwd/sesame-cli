/* cache.c -- index discovery, platform detection, and fetching.
 *
 * The store keeps every version side by side and marks the active one with a
 * symlink, so versions coexist and switching is instant:
 *
 *   <store>/EPICv2.ordering.v1/EPICv2.ordering.tsv.gz
 *   <store>/EPICv2.ordering.v2/EPICv2.ordering.tsv.gz
 *   <store>/EPICv2.ordering.tsv.gz -> EPICv2.ordering.v2/EPICv2.ordering.tsv.gz
 *
 * The symlink is what the flat lookup finds, so "drop a file in a directory and
 * it is used" still works -- a plain file put there by hand is equally valid.
 * The link target is relative, so the whole store can be moved or mounted
 * elsewhere. `sesame use <platform> <tag>` just repoints it.
 *
 * Store location (fetch writes here). ONE sesame variable, not two -- the store
 * is a managed, versioned thing you switch with `sesame use`, so it is not a
 * "cache" and does not need a second name:
 *
 *   1. $SESAME_INDEX_DIR                       explicit
 *   2. <dir of the binary>/data                a checkout: found from any cwd
 *   3. $XDG_CACHE_HOME/sesame                  standard fallback (not ours)
 *   4. ~/Library/Caches/sesame (macOS) | ~/.cache/sesame
 *
 * Resolution order (explicit always wins; same rule for library, CLI and any
 * future binding):
 *
 *   1. --index <path>                       (caller-supplied)
 *   2. <store>/<platform>.ordering.tsv.gz   (the active symlink, or a hand-placed file)
 *   3. ./<platform>.ordering.tsv.gz
 *   4. ./data/<platform>.ordering.tsv.gz
 *   5. <store>/<tag>/<platform>.ordering.tsv.gz   (the registry-pinned version)
 *
 * sesame NEVER prompts and NEVER downloads implicitly. If the index is
 * missing, the caller gets an error naming the exact command to run. An
 * interactive prompt would hang forever in a Nextflow job or a Docker build --
 * silently, on question one, across a whole run. `sesame fetch` is the only
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

/* Point <store>/<file> at <tag>/<file>. The target is RELATIVE so the store can
 * be moved or mounted elsewhere. Repointing is atomic: symlink to a temp name,
 * then rename over -- a concurrent reader sees either the old or the new target,
 * never a missing one.
 *
 * Refuses to clobber a regular file: if someone dropped their own
 * <platform>.ordering.tsv.gz there by hand, that is a deliberate choice and we
 * do not silently replace it with a link. */
int sesame_index_activate(const char *store, const char *tag, const char *file,
                          sesame_err_t *err)
{
    char link[4096], tmp[4096], target[4096];
    struct stat st;

    if (!store || !tag || !file)
        return sesame__fail(err, SESAME_ERR_IO, "activate needs store/tag/file");

    snprintf(link,   sizeof link,   "%s/%s", store, file);
    snprintf(tmp,    sizeof tmp,    "%s/.%s.link.tmp", store, file);
    snprintf(target, sizeof target, "%s/%s", tag, file);   /* relative */

    if (lstat(link, &st) == 0 && !S_ISLNK(st.st_mode))
        return sesame__fail(err, SESAME_ERR_IO,
            "%s is a regular file, not a version link -- refusing to replace it "
            "(remove it to let sesame manage this platform)", link);

    remove(tmp);
    if (symlink(target, tmp) != 0)
        return sesame__fail(err, SESAME_ERR_IO, "cannot create link %s", tmp);
    if (rename(tmp, link) != 0) {
        remove(tmp);
        return sesame__fail(err, SESAME_ERR_IO, "cannot activate %s", link);
    }
    return SESAME_OK;
}

/* Which tag is currently active, i.e. what the symlink points at. Returns 0 and
 * fills tag_out on success. */
int sesame_index_active(const char *store, const char *file,
                        char *tag_out, size_t n)
{
    char link[4096], buf[4096];
    ssize_t k;
    char *slash;

    snprintf(link, sizeof link, "%s/%s", store, file);
    k = readlink(link, buf, sizeof buf - 1);
    if (k < 0) return -1;
    buf[k] = '\0';
    slash = strchr(buf, '/');
    if (!slash) return -1;
    *slash = '\0';
    snprintf(tag_out, n, "%s", buf);
    return 0;
}

int sesame_index_locate(const char *platform, char *out, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    char dir[4096];

    if (!platform) return -1;

    /* The active symlink (or a file placed by hand). is_file() stats through
     * the link, so a dangling one simply falls through to the next candidate. */
    sesame_store_dir(dir, sizeof dir);
    snprintf(out, n, "%s/%s.ordering.tsv.gz", dir, platform);
    if (is_file(out)) return 0;

    snprintf(out, n, "./%s.ordering.tsv.gz", platform);
    if (is_file(out)) return 0;

    /* data/ next to the cwd -- convenient in a source checkout. */
    snprintf(out, n, "./data/%s.ordering.tsv.gz", platform);
    if (is_file(out)) return 0;

    /* The registry-pinned version, even if no symlink points at it. */
    if (reg) {
        snprintf(out, n, "%s/%s/%s", dir, reg->tag, reg->file);
        if (is_file(out)) return 0;
    }

    out[0] = '\0';
    return -1;
}

/* Fills msg with the "here is how to fix it" text used when no index is found.
 * Deliberately verbose: this is the error a first-time user will hit. */
void sesame_index_missing_help(const char *platform, char *msg, size_t n)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    const char *idxdir = getenv("SESAME_INDEX_DIR");
    char dir[4096];
    sesame_store_dir(dir, sizeof dir);
    snprintf(msg, n,
        "no index found for platform %s\n"
        "  searched:\n"
        "    $SESAME_INDEX_DIR  %s\n"
        "    ./%s.ordering.tsv.gz\n"
        "    ./data/%s.ordering.tsv.gz\n"
        "    %s/%s/%s\n"
        "  fix, any of:\n"
        "    sesame fetch %s             download it to the cache\n"
        "    sesame betas --index <path> ...\n"
        "    export SESAME_INDEX_DIR=<dir>",
        platform,
        (idxdir && *idxdir) ? idxdir : "(unset)",
        platform, platform,
        dir, reg ? reg->tag : "<tag>", reg ? reg->file : "<file>",
        platform);
}

#ifdef SESAME_HAVE_CURL
static size_t wr(void *p, size_t sz, size_t nm, void *ud)
{
    return fwrite(p, sz, nm, (FILE *)ud);
}

/* The retrieval primitive: fetch (tag, file) into <cache>/<tag>/<file>.
 * want_sha may be NULL for an unpinned explicit fetch. */
int sesame_fetch(const char *tag, const char *file, const char *want_sha,
                 int force, char *out_path, size_t out_n, sesame_err_t *err)
{
    char dir[4096], sub[4096], dest[4096], tmp[4096], url[4096], got[65];
    CURL *cu = NULL;
    FILE *f = NULL;
    CURLcode rc;
    long code = 0;

    if (err) { err->code = SESAME_OK; err->msg[0] = '\0'; }
    if (!tag || !file)
        return sesame__fail(err, SESAME_ERR_IO, "fetch needs a tag and a file");

    sesame_store_dir(dir, sizeof dir);
    snprintf(sub,  sizeof sub,  "%s/%s", dir, tag);
    mkdirs(sub);
    snprintf(dest, sizeof dest, "%s/%s", sub, file);
    snprintf(tmp,  sizeof tmp,  "%s/.%s.part", sub, file);
    snprintf(url,  sizeof url,  "%s/%s/%s", SESAME_BASE_URL, tag, file);

    if (!force && is_file(dest)) {
        if (!want_sha || (sesame__sha256_file(dest, got) == 0 &&
                          strcmp(got, want_sha) == 0)) {
            sesame_index_activate(dir, tag, file, NULL);  /* make it current */
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
    /* Freshly fetched version becomes the active one. */
    if (sesame_index_activate(dir, tag, file, err) != SESAME_OK)
        return err ? err->code : SESAME_ERR_IO;
    snprintf(out_path, out_n, "%s", dest);
    return SESAME_OK;
}

/* Convenience over the primitive: look the platform's pinned (tag, file,
 * sha256) up in the registry and fetch that. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err)
{
    const sesame_reg_t *reg = sesame__reg_for_platform(platform);
    if (!reg)
        return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
            "unknown platform '%s' (known: EPIC, EPICv2, HM450, MSA)", platform);
    return sesame_fetch(reg->tag, reg->file, reg->sha256, force,
                        out_path, out_n, err);
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
