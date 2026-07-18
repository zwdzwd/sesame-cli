/* main.c -- sesame command line interface
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "../src/internal.h"

#include <dirent.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

static int usage(void)
{
    fputs(
      "sesame -- Infinium DNA methylation preprocessing\n"
      "\n"
      "usage:\n"
      "  sesame idat-dump [--head N] [--tsv] <file.idat|file.idat.gz>\n"
      "      Print IDAT contents. Default prints a summary header;\n"
      "      --tsv emits addr<TAB>mean<TAB>sd<TAB>nbeads with no header.\n"
      "      --head N limits to the first N records (default: all for --tsv,\n"
      "      5 for summary).\n"
      "\n"
      "  sesame betas [--index <ordering.tsv.gz>] [--platform P] [--prep CODE]\n"
      "                [--min-beads N] [--no-mask] [--threads N] [--f64]\n"
      "                <prefix> [<prefix> ...]\n"
      "      Compute beta values from <prefix>_Grn.idat[.gz] and\n"
      "      <prefix>_Red.idat[.gz]. Equivalent to\n"
      "      openSesame(prefix, prep=CODE, func=getBetas).\n"
      "      One prefix  -> Probe_ID<TAB>beta (NA for missing).\n"
      "      Many        -> a matrix: header Probe_ID + one column per sample\n"
      "                     (named by basename); the index and masks are parsed\n"
      "                     once and samples run in parallel. All prefixes must\n"
      "                     be the same platform. A sample that fails becomes an\n"
      "                     NA column (with a warning) and the exit code is 1.\n"
      "      --prep CODE    preprocessing steps to apply (QCDPB; e.g. \"QCDPB\")\n"
      "      --min-beads N  mask probes with any bead count < N (default: off)\n"
      "      --no-mask      ignore the mask column; emit every beta\n"
      "      --threads N    worker threads (default: online CPUs)\n"
      "      --dump-col     emit Probe_ID<TAB>col (G/R/2) instead of betas\n"
      "                     (single prefix; for differential-testing a prep step)\n"
      "      --f64          write raw little-endian float64 to stdout instead\n"
      "                     of text (NA as NaN), sample-major for a batch. For\n"
      "                     lossless comparison: R's text parser does not\n"
      "                     correctly round 17-digit decimals, losing up to 1 ULP.\n"
      "      If --index is omitted the platform is detected from the IDAT bead\n"
      "      count and the index looked up in $SESAME_INDEX_DIR, ., ./data,\n"
      "      then the cache. sesame never downloads implicitly.\n"
      "\n"
      "  sesame intensity [--index ..] [--platform P] [--min-beads N]\n"
      "            [--threads N] [--f64 | --cg <out.cg>] <prefix> [<prefix> ...]\n"
      "      Total signal intensity (M+U) per probe -- the CNV signal input, and\n"
      "      how the normal reference is generated. Same matrix output as betas\n"
      "      (one column per sample); no preprocessing. --cg writes a YAME\n"
      "      format-4 .cg (float per probe, + .idx of sample names) instead of\n"
      "      the TSV.\n"
      "\n"
      "  sesame qc [--index <ordering.tsv.gz>] [--platform P] [--min-beads N]\n"
      "            [--threads N] <prefix> [<prefix> ...]\n"
      "      Per-sample QC metrics (the sesameQC panel) as a TSV: one row per\n"
      "      sample, one column per metric, headline being detection success\n"
      "      rate (%% Detection Success). Computed from the raw signal; runs\n"
      "      pOOBAH internally, so it needs the platform's .cm mask in the store.\n"
      "      Batch-parallel like betas.\n"
      "\n"
      "  sesame dml --betas <matrix.tsv> (--formula '~a+b' --meta <s.tsv>\n"
      "             | --design <X.tsv>) [--threads N]\n"
      "      Per-probe differential methylation: OLS of each probe's betas on the\n"
      "      design, with per-coefficient t-tests, a holdout F-test per\n"
      "      categorical variable, effect sizes, and BH-adjusted p-values. The\n"
      "      betas matrix is `sesame betas` batch output (Probe_ID + one column\n"
      "      per sample); --meta is a TSV whose first column matches those sample\n"
      "      names. --formula takes main effects only (categorical auto-dummied);\n"
      "      use --design to pass a pre-built numeric design for interactions.\n"
      "\n"
      "  sesame fetch [<platform>] [--force]\n"
      "      Download <platform> (default: all published) at the pinned tag.\n"
      "      Verifies every file against its digest; a file already present and\n"
      "      matching is skipped. The ONLY path that touches the network.\n"
      "      Never prompts.\n"
      "\n"
      "  sesame index-info\n"
      "      Show the store location, which tag it holds, and each platform.\n"
      "\n"
      "  sesame version\n",
      stderr);
    return 2;
}

/* Resolve <prefix>_Grn.idat, trying .gz as R does (R/sesame.R:331-345). */
static int resolve_idat(const char *prefix, const char *chan,
                        char *out, size_t outsz)
{
    FILE *f;
    snprintf(out, outsz, "%s_%s.idat", prefix, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    snprintf(out, outsz, "%s_%s.idat.gz", prefix, chan);
    if ((f = fopen(out, "rb"))) { fclose(f); return 0; }
    return -1;
}

static int cmd_fetch(int argc, char **argv)
{
    const char *platform = NULL;
    int force = 0, i;
    char dir[4096], path[4096];
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--force") == 0) force = 1;
        else if (argv[i][0] == '-' && argv[i][1] != '\0') return usage();
        else platform = argv[i];
    }

    (void)dir;
    if (platform) {
        if (sesame_fetch_index(platform, force, path, sizeof path, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg);
            return 1;
        }
    } else if (sesame_fetch_all(force, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg);
        return 1;
    }
    return 0;
}

/* ANSI only when stdout is a terminal, and never when NO_COLOR is set
 * (no-color.org). Piping index-info into grep must stay plain. */
static int use_color(void)
{
    const char *nc = getenv("NO_COLOR");
    if (nc && *nc) return 0;
    return isatty(STDOUT_FILENO);
}
#define C(code) (use_color() ? (code) : "")
#define C_RESET C("\x1b[0m")
#define C_BOLD  C("\x1b[1m")
#define C_DIM   C("\x1b[2m")
#define C_GRN   C("\x1b[32m")
#define C_YEL   C("\x1b[33m")

static int cmd_index_info(int argc, char **argv)
{
    char dir[4096], path[4096], exe[4096];
    static const char *plats[] = { "EPIC", "EPICv2", "HM450", "MSA", NULL };
    const char *e = getenv("SESAME_INDEX_DIR");
    (void)argc; (void)argv;

    sesame_store_dir(dir, sizeof dir);

    printf("%sstore%s  %s%s%s\n", C_BOLD, C_RESET, C_BOLD, dir, C_RESET);
    if (e && *e)
        printf("       %sfrom SESAME_INDEX_DIR%s\n", C_DIM, C_RESET);
    else if (sesame__exe_data_dir(exe, sizeof exe) == 0)
        printf("       %sSESAME_INDEX_DIR unset - using data/ beside the binary%s\n",
               C_DIM, C_RESET);
    else
        printf("       %sSESAME_INDEX_DIR unset - using the XDG store%s\n",
               C_DIM, C_RESET);

    printf("%stag%s    %s%s%s   %spinned by this build%s\n",
           C_BOLD, C_RESET, C_GRN, sesame_default_tag(), C_RESET,
           C_DIM, C_RESET);

    printf("\n%s%-9s %s%s\n", C_BOLD, "PLATFORM", "RESOLVED", C_RESET);
    for (int i = 0; plats[i]; i++) {
        if (sesame_index_locate(plats[i], path, sizeof path) == 0)
            printf("%-9s %s\n", plats[i], path);
        else
            printf("%-9s %smissing%s  %s(sesame fetch)%s\n",
                   plats[i], C_YEL, C_RESET, C_DIM, C_RESET);
    }
    return 0;
}

/* --------------------------------------------------------------- betas ---
 *
 * Batch-capable: many <prefix> in one process share the opened index and the
 * loaded masks, and independent per-sample work runs across a pthread pool. */

/* Basename of a path, without modifying it (basename(3) may). */
static const char *path_basename(const char *p)
{
    const char *s = strrchr(p, '/');
    return s ? s + 1 : p;
}

/* Apply the prep code string to a SigDF; masks are preloaded by the caller. */
static int apply_prep(sesame_sigdf_t *s, const char *prep,
                      const uint8_t *qmask, int32_t qn,
                      const uint8_t *bgmask, int32_t bgn, sesame_err_t *err)
{
    for (const char *c = prep; *c; c++) {
        int rc;
        switch (*c) {
        case 'Q': rc = sesame_prep_quality_mask(s, qmask, qn, err); break;
        case 'C': rc = sesame_prep_infer_channel(s, 0, 0, err); break;
        case 'D': rc = sesame_prep_dye_bias_nl(s, err); break;
        case 'P': rc = sesame_prep_poobah(s, bgmask, bgn, 0.05, 1, err); break;
        case 'B': rc = sesame_prep_noob(s, bgmask, bgn, 1, 15.0, err); break;
        default:
            return sesame__fail(err, SESAME_ERR_UNSUPPORTED,
                "prep code '%c' not implemented (have: Q, C, D, P, B)", *c);
        }
        if (rc != SESAME_OK) return rc;
    }
    return SESAME_OK;
}

/* Read one prefix's IDAT pair and build its SigDF against the shared index. If
 * expect_plat is non-NULL the sample's detected platform must match it -- this
 * guards a mixed-platform batch from silently producing a misaligned column. */
static int build_sigdf_for(const char *prefix, const sesame_index_t *ix,
                           const char *expect_plat, int min_beads,
                           sesame_sigdf_t **out, sesame_err_t *err)
{
    char gpath[4096], rpath[4096];
    sesame_idat_t *g = NULL, *r = NULL;
    sesame_sigdf_t *s = NULL;
    int rc;

    if (resolve_idat(prefix, "Grn", gpath, sizeof gpath) ||
        resolve_idat(prefix, "Red", rpath, sizeof rpath))
        return sesame__fail(err, SESAME_ERR_IO, "IDAT pair not found for %s", prefix);

    if ((rc = sesame_idat_read(gpath, &g, err)) != SESAME_OK ||
        (rc = sesame_idat_read(rpath, &r, err)) != SESAME_OK)
        goto done;

    if (expect_plat) {
        const char *p = sesame_platform_from_beads(g->n);
        if (!p || strcmp(p, expect_plat) != 0) {
            rc = sesame__fail(err, SESAME_ERR_FORMAT,
                "%s is %s (%d beads), not the batch platform %s",
                prefix, p ? p : "unknown", g->n, expect_plat);
            goto done;
        }
    }

    if (!(s = sesame_sigdf_from_idats(g, r, ix, min_beads, err))) {
        rc = err ? err->code : SESAME_ERR_IO; goto done;
    }
    *out = s; s = NULL; rc = SESAME_OK;

done:
    sesame_idat_free(g);
    sesame_idat_free(r);
    sesame_sigdf_free(s);
    return rc;
}

/* Compute one sample's betas into out (n = index nprobes); *status_out gets the
 * SigDF status word. */
static int betas_one(const char *prefix, const sesame_index_t *ix,
                     const char *expect_plat, const char *prep,
                     const uint8_t *qmask, int32_t qn,
                     const uint8_t *bgmask, int32_t bgn,
                     int min_beads, int apply_mask,
                     double *out, uint32_t *status_out, sesame_err_t *err)
{
    sesame_sigdf_t *s = NULL;
    int rc = build_sigdf_for(prefix, ix, expect_plat, min_beads, &s, err);
    if (rc != SESAME_OK) return rc;
    if ((rc = apply_prep(s, prep, qmask, qn, bgmask, bgn, err)) == SESAME_OK)
        rc = sesame_get_betas(s, apply_mask, out, err);
    if (status_out) *status_out = s->status;
    sesame_sigdf_free(s);
    return rc;
}

/* A samples x probes result store, sample-major (column j at offset j*n). Backed
 * by an unlinked temp file so the OS can page a matrix larger than RAM to disk;
 * falls back to malloc. Either backing is just a flat double * to callers. */
typedef struct {
    double *data;
    size_t  bytes;
    int     fd;         /* -1 when malloc-backed */
} betas_matrix;

static int betas_matrix_init(betas_matrix *m, int32_t n, int32_t nsamp)
{
    size_t bytes = (size_t)n * (size_t)nsamp * sizeof(double);
    const char *tmp = getenv("TMPDIR");
    char path[4096];
    int fd;

    m->data = NULL; m->bytes = bytes; m->fd = -1;
    if (bytes == 0) return -1;

    snprintf(path, sizeof path, "%s/sesame-betas-XXXXXX",
             tmp && *tmp ? tmp : "/tmp");
    if ((fd = mkstemp(path)) >= 0) {
        unlink(path);                        /* auto-clean on close/crash */
        if (ftruncate(fd, (off_t)bytes) == 0) {
            void *p = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
            if (p != MAP_FAILED) { m->data = (double *)p; m->fd = fd; return 0; }
        }
        close(fd);
    }
    m->data = (double *)malloc(bytes);       /* fallback: anonymous RAM */
    return m->data ? 0 : -1;
}

static void betas_matrix_free(betas_matrix *m)
{
    if (!m->data) return;
    if (m->fd >= 0) { munmap(m->data, m->bytes); close(m->fd); }
    else free(m->data);
    m->data = NULL;
}

/* Shared state for the worker pool. Only `next`, `status_or`, and `nfail` are
 * mutable, all guarded by `lock`; everything else is read-only or a disjoint
 * per-column write into `mat`. */
typedef struct {
    char           **prefixes;
    int32_t          nsamp, n;
    const sesame_index_t *ix;
    const char      *plat, *prep;
    const uint8_t   *qmask, *bgmask;
    int32_t          qn, bgn;
    int              min_beads, apply_mask;
    double          *mat;
    pthread_mutex_t  lock;
    int32_t          next, nfail;
    uint32_t         status_or;
} batch_ctx;

static void *batch_worker(void *arg)
{
    batch_ctx *c = (batch_ctx *)arg;
    for (;;) {
        int32_t j, k;
        uint32_t st = 0;
        double *col;
        sesame_err_t e;
        int rc;

        pthread_mutex_lock(&c->lock);
        j = c->next < c->nsamp ? c->next++ : -1;
        pthread_mutex_unlock(&c->lock);
        if (j < 0) break;

        col = c->mat + (size_t)j * (size_t)c->n;
        rc = betas_one(c->prefixes[j], c->ix, c->plat, c->prep,
                       c->qmask, c->qn, c->bgmask, c->bgn,
                       c->min_beads, c->apply_mask, col, &st, &e);
        if (rc != SESAME_OK)
            for (k = 0; k < c->n; k++) col[k] = NAN;   /* NA column, keep shape */

        pthread_mutex_lock(&c->lock);
        if (rc != SESAME_OK) {
            c->nfail++;
            fprintf(stderr, "sesame: %s: %s\n",
                    path_basename(c->prefixes[j]), e.msg);
        } else {
            c->status_or |= st;
        }
        pthread_mutex_unlock(&c->lock);
    }
    return NULL;
}

static int cmd_betas(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *prep = "", *plat = NULL;
    int min_beads = 0, apply_mask = 1, f64 = 0, dump_col = 0, i, rc = 1;
    int nthreads = 0;                        /* 0 => auto (online CPUs) */
    char **prefixes = NULL;
    int32_t nsamp = 0, n = 0, qn = 0, bgn = 0;
    char resolved[4096];
    sesame_index_t *ix = NULL;
    uint8_t *qmask = NULL, *bgmask = NULL;
    betas_matrix mat = { NULL, 0, -1 };
    int have_mat = 0, nfail = 0;
    uint32_t status_or = 0;
    sesame_err_t e;

    prefixes = (char **)malloc((size_t)(argc + 1) * sizeof(char *));
    if (!prefixes) { fprintf(stderr, "sesame: out of memory\n"); return 1; }

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i + 1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i + 1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i + 1 < argc)
            min_beads = (int)strtol(argv[++i], NULL, 10);
        else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0)
                 && i + 1 < argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--no-mask") == 0) apply_mask = 0;
        else if (strcmp(argv[i], "--prep") == 0 && i + 1 < argc) prep = argv[++i];
        else if (strcmp(argv[i], "--dump-col") == 0) dump_col = 1;
        else if (strcmp(argv[i], "--f64") == 0) f64 = 1;
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            free(prefixes); return usage();
        } else prefixes[nsamp++] = argv[i];
    }
    if (nsamp == 0) { free(prefixes); return usage(); }
    if (dump_col && nsamp != 1) {
        fprintf(stderr, "sesame: --dump-col takes a single prefix\n");
        free(prefixes); return 1;
    }

    /* Determine the platform and open the index ONCE. Auto-detect from the first
     * sample's Grn bead count when neither --index nor --platform is given. */
    if (!idxpath) {
        if (!platform) {
            char gp[4096]; sesame_idat_t *g0 = NULL; int32_t beads;
            if (resolve_idat(prefixes[0], "Grn", gp, sizeof gp)) {
                fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefixes[0]);
                goto out;
            }
            if (sesame_idat_read(gp, &g0, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            beads = g0->n;
            platform = sesame_platform_from_beads(beads);
            sesame_idat_free(g0);
            if (!platform) {
                fprintf(stderr,
                    "sesame: cannot identify platform from %d beads.\n"
                    "  pass --platform <P> or --index <path>.\n", beads);
                goto out;
            }
            fprintf(stderr, "sesame: detected %s (%d beads)\n", platform, beads);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024];
            sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); goto out;
        }
        idxpath = resolved;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    n = sesame_index_nprobes(ix);
    plat = platform;   /* NULL only when --index is given without --platform */

    /* Masks: loaded once (they read the .cm), shared read-only across samples. */
    if (strchr(prep, 'Q') || strchr(prep, 'P') || strchr(prep, 'B')) {
        if (!plat) {
            fprintf(stderr, "sesame: --prep Q/P/B needs a known platform; "
                            "pass --platform with --index.\n");
            goto out;
        }
        if (strchr(prep, 'Q') &&
            sesame_quality_mask(plat, &qmask, &qn, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out;
        }
        if ((strchr(prep, 'P') || strchr(prep, 'B')) &&
            sesame_background_mask(plat, &bgmask, &bgn, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out;
        }
    }

    if (dump_col) {                          /* single-sample debugging aid */
        sesame_sigdf_t *s = NULL;
        static const char *nm[] = { "2", "G", "R" };
        if (build_sigdf_for(prefixes[0], ix, plat, min_beads, &s, &e) != SESAME_OK ||
            apply_prep(s, prep, qmask, qn, bgmask, bgn, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg);
            sesame_sigdf_free(s); goto out;
        }
        for (int32_t k = 0; k < s->n; k++)
            printf("%s\t%s\n", sesame_index_probe_id(ix, k), nm[s->col[k]]);
        sesame_sigdf_free(s);
        rc = 0; goto out;
    }

    if (betas_matrix_init(&mat, n, nsamp) != 0) {
        fprintf(stderr, "sesame: cannot allocate %d x %d beta matrix\n", nsamp, n);
        goto out;
    }
    have_mat = 1;

    /* Run the pool: T-1 spawned workers plus this thread, all draining a shared
     * counter. Output column = argv position, so results are deterministic
     * regardless of scheduling; --threads 1 == --threads N byte-for-byte. */
    {
        batch_ctx ctx;
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1);
        pthread_t *th = NULL;
        int spawned = 0, t;

        if (T > nsamp) T = nsamp;
        if (T < 1) T = 1;

        ctx.prefixes = prefixes; ctx.nsamp = nsamp; ctx.n = n;
        ctx.ix = ix; ctx.plat = plat; ctx.prep = prep;
        ctx.qmask = qmask; ctx.qn = qn; ctx.bgmask = bgmask; ctx.bgn = bgn;
        ctx.min_beads = min_beads; ctx.apply_mask = apply_mask;
        ctx.mat = mat.data; ctx.next = 0; ctx.nfail = 0; ctx.status_or = 0;
        pthread_mutex_init(&ctx.lock, NULL);

        if (T > 1) th = (pthread_t *)malloc((size_t)(T - 1) * sizeof(pthread_t));
        if (th)
            for (t = 0; t < T - 1; t++)
                if (pthread_create(&th[spawned], NULL, batch_worker, &ctx) == 0)
                    spawned++;
        batch_worker(&ctx);                  /* this thread participates */
        for (t = 0; t < spawned; t++) pthread_join(th[t], NULL);
        free(th);
        pthread_mutex_destroy(&ctx.lock);

        nfail = ctx.nfail; status_or = ctx.status_or;
    }

    if (f64) {
        size_t tot = (size_t)n * (size_t)nsamp;
        if (fwrite(mat.data, sizeof(double), tot, stdout) != tot) {
            fprintf(stderr, "sesame: short write\n"); goto out;
        }
    } else if (nsamp == 1) {                 /* unchanged single-sample format */
        for (int32_t k = 0; k < n; k++) {
            const char *id = sesame_index_probe_id(ix, k);
            if (isnan(mat.data[k])) printf("%s\tNA\n", id);
            else                    printf("%s\t%.17g\n", id, mat.data[k]);
        }
    } else {                                 /* matrix: Probe_ID + one col/sample */
        fputs("Probe_ID", stdout);
        for (int32_t j = 0; j < nsamp; j++)
            printf("\t%s", path_basename(prefixes[j]));
        putchar('\n');
        for (int32_t k = 0; k < n; k++) {
            fputs(sesame_index_probe_id(ix, k), stdout);
            for (int32_t j = 0; j < nsamp; j++) {
                double b = mat.data[(size_t)j * (size_t)n + (size_t)k];
                if (isnan(b)) fputs("\tNA", stdout);
                else          printf("\t%.17g", b);
            }
            putchar('\n');
        }
    }

    if (status_or & SESAME_STAT_ADDR_MISSING)
        fprintf(stderr, "sesame: note: some probes had addresses absent from an IDAT\n");
    if (status_or & SESAME_STAT_DYEBIAS_FAILED)
        fprintf(stderr, "sesame: note: dye-bias correction gave up on a sample (green channel failed)\n");
    if (status_or & SESAME_STAT_NOOB_SKIPPED)
        fprintf(stderr, "sesame: note: noob skipped on a sample (insufficient background)\n");
    if (status_or & SESAME_STAT_NOOB_MAD0)
        fprintf(stderr, "sesame: note: noob MAD==0 fallback fired on a sample\n");
    if (nfail)
        fprintf(stderr, "sesame: %d of %d samples failed (NA columns)\n", nfail, nsamp);
    rc = nfail ? 1 : 0;

out:
    if (have_mat) betas_matrix_free(&mat);
    free(qmask);
    free(bgmask);
    free(prefixes);
    sesame_index_close(ix);
    return rc;
}

/* ---------------------------------------------------------- intensity ---
 *
 * Total signal intensity (M+U) per probe -- the CNV signal input, and the tool
 * that generates the normal reference from normal IDATs. Batch-parallel with the
 * same matrix output as `betas`; no prep, no masks. */

typedef struct {
    char           **prefixes;
    int32_t          nsamp, n;
    const sesame_index_t *ix;
    const char      *plat;
    int              min_beads;
    double          *mat;
    pthread_mutex_t  lock;
    int32_t          next, nfail;
    uint32_t         status_or;
} intens_ctx;

static void *intens_worker(void *arg)
{
    intens_ctx *c = (intens_ctx *)arg;
    for (;;) {
        int32_t j, k;
        sesame_sigdf_t *s = NULL;
        sesame_err_t e;
        uint32_t st = 0;
        double *col;
        int rc;

        pthread_mutex_lock(&c->lock);
        j = c->next < c->nsamp ? c->next++ : -1;
        pthread_mutex_unlock(&c->lock);
        if (j < 0) break;

        col = c->mat + (size_t)j * (size_t)c->n;
        rc = build_sigdf_for(c->prefixes[j], c->ix, c->plat, c->min_beads, &s, &e);
        if (rc == SESAME_OK) { st = s->status; rc = sesame_total_intensities(s, col, &e); }
        sesame_sigdf_free(s);
        if (rc != SESAME_OK)
            for (k = 0; k < c->n; k++) col[k] = NAN;

        pthread_mutex_lock(&c->lock);
        if (rc != SESAME_OK) {
            c->nfail++;
            fprintf(stderr, "sesame: %s: %s\n", path_basename(c->prefixes[j]), e.msg);
        } else c->status_or |= st;
        pthread_mutex_unlock(&c->lock);
    }
    return NULL;
}

static int cmd_intensity(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *plat = NULL, *cgpath = NULL;
    int min_beads = 0, nthreads = 0, f64 = 0, i, rc = 1;
    char **prefixes = NULL;
    int32_t nsamp = 0, n = 0;
    char resolved[4096];
    sesame_index_t *ix = NULL;
    betas_matrix mat = { NULL, 0, -1 };
    int have_mat = 0, nfail = 0;
    sesame_err_t e;

    prefixes = (char **)malloc((size_t)(argc + 1) * sizeof(char *));
    if (!prefixes) { fprintf(stderr, "sesame: out of memory\n"); return 1; }
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i+1 < argc)
            min_beads = (int)strtol(argv[++i], NULL, 10);
        else if ((strcmp(argv[i], "--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--f64") == 0) f64 = 1;
        else if (strcmp(argv[i], "--cg") == 0 && i+1 < argc) cgpath = argv[++i];
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            free(prefixes); return usage();
        } else prefixes[nsamp++] = argv[i];
    }
    if (nsamp == 0) { free(prefixes); return usage(); }

    if (!idxpath) {
        if (!platform) {
            char gp[4096]; sesame_idat_t *g0 = NULL; int32_t beads;
            if (resolve_idat(prefixes[0], "Grn", gp, sizeof gp)) {
                fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefixes[0]); goto out; }
            if (sesame_idat_read(gp, &g0, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
            beads = g0->n; platform = sesame_platform_from_beads(beads); sesame_idat_free(g0);
            if (!platform) { fprintf(stderr,
                "sesame: cannot identify platform from %d beads.\n"
                "  pass --platform <P> or --index <path>.\n", beads); goto out; }
            fprintf(stderr, "sesame: detected %s (%d beads)\n", platform, beads);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); goto out;
        }
        idxpath = resolved;
    }
    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    n = sesame_index_nprobes(ix); plat = platform;

    if (betas_matrix_init(&mat, n, nsamp) != 0) {
        fprintf(stderr, "sesame: cannot allocate %d x %d matrix\n", nsamp, n); goto out; }
    have_mat = 1;

    {
        intens_ctx ctx;
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1);
        pthread_t *th = NULL; int spawned = 0, t;
        if (T > nsamp) T = nsamp;
        if (T < 1) T = 1;
        ctx.prefixes=prefixes; ctx.nsamp=nsamp; ctx.n=n; ctx.ix=ix; ctx.plat=plat;
        ctx.min_beads=min_beads; ctx.mat=mat.data; ctx.next=0; ctx.nfail=0; ctx.status_or=0;
        pthread_mutex_init(&ctx.lock, NULL);
        if (T > 1) th = (pthread_t *)malloc((size_t)(T-1)*sizeof(pthread_t));
        if (th) for (t=0;t<T-1;t++) if (pthread_create(&th[spawned],NULL,intens_worker,&ctx)==0) spawned++;
        intens_worker(&ctx);
        for (t=0;t<spawned;t++) pthread_join(th[t],NULL);
        free(th); pthread_mutex_destroy(&ctx.lock);
        nfail = ctx.nfail;
        if (ctx.status_or & SESAME_STAT_ADDR_MISSING)
            fprintf(stderr, "sesame: note: some probes had addresses absent from an IDAT\n");
    }

    if (cgpath) {                                /* YAME format-4 .cg (+ .idx) */
        char **names = (char **)malloc((size_t)nsamp * sizeof(char *));
        int32_t j;
        if (!names) { fprintf(stderr, "sesame: oom\n"); goto out; }
        for (j = 0; j < nsamp; j++) names[j] = (char *)path_basename(prefixes[j]);
        if (sesame_write_cg(cgpath, mat.data, n, nsamp, names, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); free(names); goto out; }
        free(names);
        fprintf(stderr, "sesame: wrote %s (%d sample%s, format 4) + %s.idx\n",
                cgpath, nsamp, nsamp==1?"":"s", cgpath);
    } else if (f64) {
        size_t tot = (size_t)n * (size_t)nsamp;
        if (fwrite(mat.data, sizeof(double), tot, stdout) != tot) {
            fprintf(stderr, "sesame: short write\n"); goto out; }
    } else if (nsamp == 1) {
        for (int32_t k = 0; k < n; k++) {
            const char *id = sesame_index_probe_id(ix, k);
            if (isnan(mat.data[k])) printf("%s\tNA\n", id);
            else                    printf("%s\t%.10g\n", id, mat.data[k]);
        }
    } else {
        fputs("Probe_ID", stdout);
        for (int32_t j = 0; j < nsamp; j++) printf("\t%s", path_basename(prefixes[j]));
        putchar('\n');
        for (int32_t k = 0; k < n; k++) {
            fputs(sesame_index_probe_id(ix, k), stdout);
            for (int32_t j = 0; j < nsamp; j++) {
                double b = mat.data[(size_t)j*(size_t)n + (size_t)k];
                if (isnan(b)) fputs("\tNA", stdout); else printf("\t%.10g", b);
            }
            putchar('\n');
        }
    }
    if (nfail) fprintf(stderr, "sesame: %d of %d samples failed (NA columns)\n", nfail, nsamp);
    rc = nfail ? 1 : 0;

out:
    if (have_mat) betas_matrix_free(&mat);
    free(prefixes);
    sesame_index_close(ix);
    return rc;
}

/* ------------------------------------------------------------------ qc ---
 *
 * The sesameQC panel per sample, as a TSV (one row per sample). Batch-parallel
 * like `betas`; the per-sample result is a small struct, so no matrix store. */

typedef struct {
    char           **prefixes;
    int32_t          nsamp;
    const sesame_index_t *ix;
    const char      *plat;
    const uint8_t   *bgmask;
    int32_t          bgn;
    int              min_beads;
    sesame_qc_t     *res;
    pthread_mutex_t  lock;
    int32_t          next, nfail;
} qc_ctx;

static void *qc_worker(void *arg)
{
    qc_ctx *c = (qc_ctx *)arg;
    for (;;) {
        int32_t j;
        sesame_sigdf_t *s = NULL;
        sesame_err_t e;
        int rc;

        pthread_mutex_lock(&c->lock);
        j = c->next < c->nsamp ? c->next++ : -1;
        pthread_mutex_unlock(&c->lock);
        if (j < 0) break;

        rc = build_sigdf_for(c->prefixes[j], c->ix, c->plat, c->min_beads, &s, &e);
        if (rc == SESAME_OK)
            rc = sesame_qc_calc(s, c->bgmask, c->bgn, &c->res[j], &e);
        sesame_sigdf_free(s);

        if (rc != SESAME_OK) {                 /* leave the row as pre-filled NaN */
            pthread_mutex_lock(&c->lock);
            c->nfail++;
            fprintf(stderr, "sesame: %s: %s\n", path_basename(c->prefixes[j]), e.msg);
            pthread_mutex_unlock(&c->lock);
        }
    }
    return NULL;
}

static int cmd_qc(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *plat = NULL;
    int min_beads = 0, nthreads = 0, i, rc = 1;
    char **prefixes = NULL;
    int32_t nsamp = 0, bgn = 0;
    char resolved[4096];
    sesame_index_t *ix = NULL;
    uint8_t *bgmask = NULL;
    sesame_qc_t *res = NULL;
    sesame_err_t e;

    prefixes = (char **)malloc((size_t)(argc + 1) * sizeof(char *));
    if (!prefixes) { fprintf(stderr, "sesame: out of memory\n"); return 1; }

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i + 1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i + 1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i + 1 < argc)
            min_beads = (int)strtol(argv[++i], NULL, 10);
        else if ((strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0)
                 && i + 1 < argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            free(prefixes); return usage();
        } else prefixes[nsamp++] = argv[i];
    }
    if (nsamp == 0) { free(prefixes); return usage(); }

    if (!idxpath) {
        if (!platform) {
            char gp[4096]; sesame_idat_t *g0 = NULL; int32_t beads;
            if (resolve_idat(prefixes[0], "Grn", gp, sizeof gp)) {
                fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefixes[0]);
                goto out;
            }
            if (sesame_idat_read(gp, &g0, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            beads = g0->n;
            platform = sesame_platform_from_beads(beads);
            sesame_idat_free(g0);
            if (!platform) {
                fprintf(stderr, "sesame: cannot identify platform from %d beads.\n"
                    "  pass --platform <P> or --index <path>.\n", beads);
                goto out;
            }
            fprintf(stderr, "sesame: detected %s (%d beads)\n", platform, beads);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024];
            sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); goto out;
        }
        idxpath = resolved;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    plat = platform;

    /* QC's detection and beta groups run pOOBAH internally, so the mask is
     * mandatory (unlike `betas` with prep=""). */
    if (!plat) {
        fprintf(stderr, "sesame: qc needs a known platform for the mask; "
                        "pass --platform with --index.\n");
        goto out;
    }
    if (sesame_background_mask(plat, &bgmask, &bgn, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    res = (sesame_qc_t *)malloc((size_t)nsamp * sizeof(sesame_qc_t));
    if (!res) { fprintf(stderr, "sesame: out of memory\n"); goto out; }
    /* pre-fill every field NaN so a failed sample prints an all-NA row */
    {
        size_t k, nd = sizeof(sesame_qc_t) / sizeof(double);
        double *p = (double *)res;
        for (k = 0; k < (size_t)nsamp * nd; k++) p[k] = NAN;
    }

    {
        qc_ctx ctx;
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1);
        pthread_t *th = NULL;
        int spawned = 0, t;

        if (T > nsamp) T = nsamp;
        if (T < 1) T = 1;

        ctx.prefixes = prefixes; ctx.nsamp = nsamp; ctx.ix = ix; ctx.plat = plat;
        ctx.bgmask = bgmask; ctx.bgn = bgn; ctx.min_beads = min_beads;
        ctx.res = res; ctx.next = 0; ctx.nfail = 0;
        pthread_mutex_init(&ctx.lock, NULL);

        if (T > 1) th = (pthread_t *)malloc((size_t)(T - 1) * sizeof(pthread_t));
        if (th)
            for (t = 0; t < T - 1; t++)
                if (pthread_create(&th[spawned], NULL, qc_worker, &ctx) == 0) spawned++;
        qc_worker(&ctx);
        for (t = 0; t < spawned; t++) pthread_join(th[t], NULL);
        free(th);
        pthread_mutex_destroy(&ctx.lock);

        printf("sample\t%s\n", sesame_qc_header());
        for (i = 0; i < nsamp; i++) {
            char row[8192];
            if (sesame_qc_format_row(&res[i], row, sizeof row) < 0)
                { fprintf(stderr, "sesame: qc row too long\n"); goto out; }
            printf("%s\t%s\n", path_basename(prefixes[i]), row);
        }
        if (ctx.nfail)
            fprintf(stderr, "sesame: %d of %d samples failed (NA rows)\n", ctx.nfail, nsamp);
        rc = ctx.nfail ? 1 : 0;
    }

out:
    free(res);
    free(bgmask);
    free(prefixes);
    sesame_index_close(ix);
    return rc;
}

/* ----------------------------------------------------------------- dml ---
 *
 * Per-probe differential methylation from a betas matrix + sample metadata. The
 * heavy per-probe OLS lives in src/dml.c; here we parse the inputs, build the
 * design (a simple main-effects formula, or an explicit matrix), fit in
 * parallel, BH-adjust, and emit a TSV. */

/* Read one \n-terminated line (newline stripped) into *buf (grown as needed).
 * Returns length, or -1 at EOF. Portable substitute for getline. */
static long read_line(FILE *f, char **buf, size_t *cap)
{
    size_t len = 0;
    if (!*buf) { *cap = 65536; *buf = (char *)malloc(*cap); if (!*buf) return -1; }
    for (;;) {
        if (len + 1 >= *cap) {
            char *nb = (char *)realloc(*buf, *cap * 2); if (!nb) return -1;
            *buf = nb; *cap *= 2;
        }
        if (!fgets(*buf + len, (int)(*cap - len), f)) break;
        len += strlen(*buf + len);
        if (len && (*buf)[len-1] == '\n') { (*buf)[--len] = '\0'; return (long)len; }
        if (feof(f)) break;
    }
    return len ? (long)len : -1;
}

/* Split s in place on tabs; store up to max field pointers in tok. Returns the
 * field count (may exceed max, but only max are stored). */
static int split_tabs(char *s, char **tok, int max)
{
    int n = 0;
    for (;;) {
        char *t = strchr(s, '\t');
        if (n < max) tok[n] = s;
        n++;
        if (!t) break;
        *t = '\0'; s = t + 1;
    }
    return n;
}

static double parse_num(const char *s)
{
    char *end;
    double v;
    if (!s || !*s) return NAN;
    if (strcmp(s, "NA") == 0 || strcmp(s, "NaN") == 0) return NAN;
    v = strtod(s, &end);
    return end == s ? NAN : v;
}

/* Load a betas matrix (Probe_ID + one column per sample). Fills probe ids, an
 * nprobe*nsamp row-major matrix, and the sample names. 0 on success. */
static int load_betas(const char *path, char ***ids_out, double **mat_out,
                      int32_t *nprobe_out, int32_t *nsamp_out, char ***samp_out)
{
    FILE *f = fopen(path, "rb");
    char *line = NULL, **tok = NULL;
    size_t cap = 0, cap_rows = 0;
    int32_t nsamp = 0, np = 0, i;
    char **ids = NULL, **samp = NULL;
    double *mat = NULL;
    int maxtok;

    if (!f) return -1;
    if (read_line(f, &line, &cap) < 0) { fclose(f); free(line); return -1; }
    maxtok = (int)(cap);                            /* generous upper bound */
    tok = (char **)malloc((size_t)maxtok * sizeof(char *));
    if (!tok) { fclose(f); free(line); return -1; }
    nsamp = (int32_t)split_tabs(line, tok, maxtok) - 1;   /* minus Probe_ID col */
    if (nsamp < 1) { fclose(f); free(line); free(tok); return -1; }
    samp = (char **)malloc((size_t)nsamp * sizeof(char *));
    for (i = 0; i < nsamp; i++) samp[i] = strdup(tok[i+1]);

    while (read_line(f, &line, &cap) >= 0) {
        int nf = split_tabs(line, tok, maxtok);
        if (nf < 1 || tok[0][0] == '\0') continue;
        if ((size_t)np >= cap_rows) {
            size_t nc = cap_rows ? cap_rows * 2 : 8192;
            char **ni = (char **)realloc(ids, nc * sizeof(char *));
            double *nm = (double *)realloc(mat, nc * (size_t)nsamp * sizeof(double));
            if (!ni || !nm) { free(ni?ni:ids); free(nm?nm:mat); fclose(f);
                              free(line); free(tok); free(samp); return -1; }
            ids = ni; mat = nm; cap_rows = nc;
        }
        ids[np] = strdup(tok[0]);
        for (i = 0; i < nsamp; i++)
            mat[(size_t)np*(size_t)nsamp + (size_t)i] = (i+1 < nf) ? parse_num(tok[i+1]) : NAN;
        np++;
    }
    fclose(f); free(line); free(tok);
    *ids_out = ids; *mat_out = mat; *nprobe_out = np; *nsamp_out = nsamp;
    *samp_out = samp;
    return 0;
}

/* Benjamini-Hochberg adjust p[0..n-1] (NaN preserved) into out; scratch holds n
 * pairs. */
typedef struct { double p; int32_t i; } bh_pair;
static int bh_cmp(const void *a, const void *b)
{
    double x = ((const bh_pair *)a)->p, y = ((const bh_pair *)b)->p;
    return x < y ? -1 : x > y ? 1 : 0;
}
static void bh_adjust(const double *p, int32_t n, double *out, bh_pair *pr)
{
    int32_t cnt = 0, r;
    double prev = INFINITY;
    for (r = 0; r < n; r++) { out[r] = NAN; if (!isnan(p[r])) { pr[cnt].p = p[r]; pr[cnt].i = r; cnt++; } }
    if (cnt == 0) return;
    qsort(pr, (size_t)cnt, sizeof(bh_pair), bh_cmp);
    for (r = cnt; r >= 1; r--) {                     /* from largest p down: cummin */
        double a = pr[r-1].p * (double)cnt / (double)r;
        if (a < prev) prev = a;
        out[pr[r-1].i] = prev > 1.0 ? 1.0 : prev;
    }
}

/* The built design plus the labels needed to write the header. */
typedef struct {
    double  *X;                     /* nsamp * p, row-major */
    int32_t  p, nvar;
    char   **coef;                  /* p coefficient names */
    char   **vname;                 /* nvar categorical-variable names */
    int32_t *var_off, *var_col;     /* grouping for F-tests / effect size */
} design_t;

static void design_free(design_t *d)
{
    int32_t j;
    if (d->coef) { for (j = 0; j < d->p; j++) free(d->coef[j]); free(d->coef); }
    if (d->vname) { for (j = 0; j < d->nvar; j++) free(d->vname[j]); free(d->vname); }
    free(d->X); free(d->var_off); free(d->var_col);
    memset(d, 0, sizeof *d);
}

static char *join2(const char *a, const char *b)
{
    size_t la = strlen(a), lb = strlen(b);
    char *s = (char *)malloc(la + lb + 1);
    if (s) { memcpy(s, a, la); memcpy(s + la, b, lb + 1); }
    return s;
}

static int cmp_str(const void *a, const void *b)
{ return strcmp(*(const char *const *)a, *(const char *const *)b); }

/* One static partition of probes for a fitting thread. */
typedef struct {
    const sesame_dml_design_t *d;
    const double *betas;
    int32_t lo, hi;
    double *est, *pval, *fpval, *eff;
    int32_t *nobs;
} dml_part_t;

static void *dml_worker(void *arg)
{
    dml_part_t *P = (dml_part_t *)arg;
    sesame_dml_work_t *w = sesame_dml_work_new(P->d->m, P->d->p);
    int32_t i, nv = P->d->nvar ? P->d->nvar : 1, p = P->d->p, m = P->d->m;
    if (!w) return NULL;
    for (i = P->lo; i < P->hi; i++)
        P->nobs[i] = sesame_dml_fit(P->d, w, P->betas + (size_t)i*(size_t)m,
            P->est + (size_t)i*(size_t)p, P->pval + (size_t)i*(size_t)p,
            P->fpval + (size_t)i*(size_t)nv, P->eff + (size_t)i*(size_t)nv);
    sesame_dml_work_free(w);
    return NULL;
}

static int cmd_dml(int argc, char **argv)
{
    const char *betapath = NULL, *metapath = NULL, *formula = NULL, *designpath = NULL;
    int nthreads = 0, i, j, rc = 1;
    char **ids = NULL, **samp = NULL;
    double *betas = NULL;
    int32_t nprobe = 0, nsamp = 0;
    design_t D; memset(&D, 0, sizeof D);
    double *est = NULL, *pval = NULL, *fpval = NULL, *eff = NULL, *padj = NULL, *fpadj = NULL;
    int32_t *nobs = NULL;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--betas") == 0 && i+1 < argc) betapath = argv[++i];
        else if (strcmp(argv[i], "--meta") == 0 && i+1 < argc) metapath = argv[++i];
        else if (strcmp(argv[i], "--formula") == 0 && i+1 < argc) formula = argv[++i];
        else if (strcmp(argv[i], "--design") == 0 && i+1 < argc) designpath = argv[++i];
        else if ((strcmp(argv[i], "--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else { fprintf(stderr, "sesame: unknown/incomplete option %s\n", argv[i]); return usage(); }
    }
    if (!betapath || (!formula && !designpath)) {
        fprintf(stderr, "sesame: dml needs --betas and either --formula (with --meta) or --design\n");
        return usage();
    }
    if (load_betas(betapath, &ids, &betas, &nprobe, &nsamp, &samp) != 0) {
        fprintf(stderr, "sesame: cannot read betas matrix %s\n", betapath); return 1;
    }

    /* ---- build the design (samples aligned to the betas columns) ---- */
    {
        FILE *mf = NULL; char *line = NULL, **tok = NULL; size_t cap = 0;
        char **mid = NULL, ***mcell = NULL, **mhdr = NULL;
        int32_t mrows = 0, mcols = 0, *smap = NULL, s;
        const char *srcpath = designpath ? designpath : metapath;
        int maxtok, ok = 1;

        if (!srcpath) { fprintf(stderr, "sesame: --formula needs --meta\n"); goto dclean; }
        mf = fopen(srcpath, "rb");
        if (!mf) { fprintf(stderr, "sesame: cannot open %s\n", srcpath); goto dclean; }
        if (read_line(mf, &line, &cap) < 0) { fprintf(stderr,"sesame: empty %s\n",srcpath); goto dclean; }
        maxtok = (int)cap;
        tok = (char **)malloc((size_t)maxtok * sizeof *tok);
        mcols = (int32_t)split_tabs(line, tok, maxtok);          /* col0 = sample id */
        mhdr = (char **)malloc((size_t)mcols * sizeof *mhdr);
        for (j = 0; j < mcols; j++) mhdr[j] = strdup(tok[j]);
        { size_t crow = 0;
          while (read_line(mf, &line, &cap) >= 0) {
            int nf = split_tabs(line, tok, maxtok);
            if (nf < 1 || tok[0][0]=='\0') continue;
            if ((size_t)mrows >= crow) {
                size_t nc = crow ? crow*2 : 256;
                mid = realloc(mid, nc*sizeof *mid); mcell = realloc(mcell, nc*sizeof *mcell);
                crow = nc;
            }
            mid[mrows] = strdup(tok[0]);
            mcell[mrows] = (char **)malloc((size_t)mcols * sizeof(char *));
            for (j = 0; j < mcols; j++) mcell[mrows][j] = strdup(j < nf ? tok[j] : "");
            mrows++;
          } }
        fclose(mf); mf = NULL;

        /* map each betas sample to its metadata row */
        smap = (int32_t *)malloc((size_t)nsamp * sizeof(int32_t));
        for (s = 0; s < nsamp; s++) {
            int32_t r, found = -1;
            for (r = 0; r < mrows; r++) if (strcmp(mid[r], samp[s]) == 0) { found = r; break; }
            if (found < 0) { fprintf(stderr, "sesame: sample %s not in %s\n", samp[s], srcpath); ok = 0; break; }
            smap[s] = found;
        }
        if (!ok) goto dclean;

        if (designpath) {                          /* explicit numeric design */
            D.p = mcols - 1; D.nvar = 0;
            D.X = (double *)malloc((size_t)nsamp * (size_t)D.p * sizeof(double));
            D.coef = (char **)malloc((size_t)D.p * sizeof(char *));
            for (j = 0; j < D.p; j++) D.coef[j] = strdup(mhdr[j+1]);
            for (s = 0; s < nsamp; s++)
                for (j = 0; j < D.p; j++)
                    D.X[(size_t)s*(size_t)D.p + (size_t)j] = parse_num(mcell[smap[s]][j+1]);
        } else {                                   /* main-effects formula */
            char fbuf[1024]; char *terms[64]; int nterm = 0, ok2 = 1;
            const char *fp = formula;
            while (*fp == ' ' || *fp == '~') fp++;
            if (strpbrk(fp, "*:^(|")) {
                fprintf(stderr, "sesame: --formula supports main effects only "
                    "(no * : ^ () ). Use --design for interactions/splines.\n");
                goto dclean;
            }
            snprintf(fbuf, sizeof fbuf, "%s", fp);
            { char *t = strtok(fbuf, "+");
              while (t && nterm < 64) {
                while (*t==' ') t++;
                { char *e = t + strlen(t); while (e>t && e[-1]==' ') *--e = '\0'; }
                if (*t && strcmp(t,"1")!=0) terms[nterm++] = t;
                t = strtok(NULL, "+");
              } }
            /* pass 1: classify each term, collect levels; count columns */
            {
                int *is_cat = (int *)calloc((size_t)nterm, sizeof(int));
                int *mcol   = (int *)malloc((size_t)nterm * sizeof(int));
                char ***levs = (char ***)calloc((size_t)nterm, sizeof(char **));
                int *nlev = (int *)calloc((size_t)nterm, sizeof(int));
                int32_t p = 1, nvar = 0, k;
                for (k = 0; k < nterm; k++) {
                    int c = -1; int cat = 0;
                    for (j = 1; j < mcols; j++) if (strcmp(mhdr[j], terms[k])==0) { c = j; break; }
                    if (c < 0) { fprintf(stderr,"sesame: no metadata column '%s'\n",terms[k]); ok2=0; break; }
                    mcol[k] = c;
                    for (s = 0; s < nsamp; s++) {
                        const char *v = mcell[smap[s]][c];
                        if (*v && isnan(parse_num(v))) { cat = 1; break; }
                    }
                    is_cat[k] = cat;
                    if (cat) {                     /* collect sorted unique levels */
                        char **lv = (char **)malloc((size_t)nsamp * sizeof(char *));
                        int nl = 0, u;
                        for (s = 0; s < nsamp; s++) {
                            const char *v = mcell[smap[s]][c]; int seen = 0;
                            for (u = 0; u < nl; u++) if (strcmp(lv[u],v)==0){seen=1;break;}
                            if (!seen) lv[nl++] = (char *)v;
                        }
                        qsort(lv, (size_t)nl, sizeof(char *), cmp_str);
                        levs[k] = lv; nlev[k] = nl;
                        p += nl - 1; nvar++;
                    } else p += 1;
                }
                if (!ok2) { free(is_cat);free(mcol);free(nlev);
                            for(k=0;k<nterm;k++) free(levs[k]); free(levs); goto dclean; }
                /* pass 2: fill X, coef names, variable groupings */
                D.p = p; D.nvar = nvar;
                D.X = (double *)calloc((size_t)nsamp * (size_t)p, sizeof(double));
                D.coef = (char **)malloc((size_t)p * sizeof(char *));
                D.vname = nvar ? (char **)malloc((size_t)nvar*sizeof(char*)) : NULL;
                D.var_off = (int32_t *)malloc((size_t)(nvar+1)*sizeof(int32_t));
                D.var_col = (int32_t *)malloc((size_t)(p)*sizeof(int32_t));
                { int32_t col = 0, vv = 0, vc = 0, lvl;
                  for (s = 0; s < nsamp; s++) D.X[(size_t)s*(size_t)p + 0] = 1.0;   /* intercept */
                  D.coef[0] = strdup("(Intercept)");
                  D.var_off[0] = 0;
                  col = 1;
                  for (k = 0; k < nterm; k++) {
                    int c = mcol[k];
                    if (!is_cat[k]) {
                        for (s = 0; s < nsamp; s++)
                            D.X[(size_t)s*(size_t)p + (size_t)col] = parse_num(mcell[smap[s]][c]);
                        D.coef[col] = strdup(terms[k]);
                        col++;
                    } else {
                        D.vname[vv] = strdup(terms[k]);
                        for (lvl = 1; lvl < nlev[k]; lvl++) {
                            for (s = 0; s < nsamp; s++) {
                                const char *v = mcell[smap[s]][c];
                                D.X[(size_t)s*(size_t)p + (size_t)col] = (strcmp(v, levs[k][lvl])==0) ? 1.0 : 0.0;
                            }
                            D.coef[col] = join2(terms[k], levs[k][lvl]);
                            D.var_col[vc++] = col;
                            col++;
                        }
                        vv++; D.var_off[vv] = vc;
                    }
                  }
                }
                free(is_cat); free(mcol); free(nlev);
                for (k = 0; k < nterm; k++) free(levs[k]);
                free(levs);
            }
            if (!ok2) goto dclean;
        }
    dclean:
        if (mf) fclose(mf);
        free(line); free(tok);
        if (mid) { for (j=0;j<mrows;j++) free(mid[j]); free(mid); }
        if (mcell) { int32_t r; for (r=0;r<mrows;r++){ if(mcell[r]){int q;for(q=0;q<mcols;q++)free(mcell[r][q]);free(mcell[r]);} } free(mcell); }
        if (mhdr) { for (j=0;j<mcols;j++) free(mhdr[j]); free(mhdr); }
        free(smap);
        if (!D.X) goto out;
    }

    /* ---- fit every probe in parallel ---- */
    est   = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    pval  = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    fpval = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    eff   = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    nobs  = (int32_t *)malloc((size_t)nprobe * sizeof(int32_t));
    if (!est || !pval || !fpval || !eff || !nobs) { fprintf(stderr,"sesame: oom\n"); goto out; }

    {
        sesame_dml_design_t dz = { nsamp, D.p, D.nvar, D.X, D.var_off, D.var_col };
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1);
        dml_part_t *parts; pthread_t *th; int t;
        if (T > nprobe) T = nprobe;
        if (T < 1) T = 1;
        parts = (dml_part_t *)malloc((size_t)T * sizeof *parts);
        th = (pthread_t *)malloc((size_t)T * sizeof *th);
        if (!parts || !th) { free(parts); free(th); fprintf(stderr,"sesame: oom\n"); goto out; }
        for (t = 0; t < T; t++) {
            parts[t].d = &dz; parts[t].betas = betas;
            parts[t].lo = (int32_t)((long)nprobe*t/T);
            parts[t].hi = (int32_t)((long)nprobe*(t+1)/T);
            parts[t].est=est; parts[t].pval=pval; parts[t].fpval=fpval;
            parts[t].eff=eff; parts[t].nobs=nobs;
        }
        for (t = 1; t < T; t++) pthread_create(&th[t], NULL, dml_worker, &parts[t]);
        dml_worker(&parts[0]);
        for (t = 1; t < T; t++) pthread_join(th[t], NULL);
        free(parts); free(th);
    }

    /* ---- BH-adjust each p-value column across probes ---- */
    padj  = (double *)malloc((size_t)nprobe * (size_t)D.p * sizeof(double));
    fpadj = (double *)malloc((size_t)nprobe * (size_t)(D.nvar?D.nvar:1) * sizeof(double));
    {
        bh_pair *pr = (bh_pair *)malloc((size_t)nprobe * sizeof(bh_pair));
        double *colb = (double *)malloc((size_t)nprobe * sizeof(double));
        double *cola = (double *)malloc((size_t)nprobe * sizeof(double));
        int32_t k;
        if (pr && colb && cola) {
            for (j = 0; j < D.p; j++) {
                for (k = 0; k < nprobe; k++) colb[k] = pval[(size_t)k*(size_t)D.p+(size_t)j];
                bh_adjust(colb, nprobe, cola, pr);
                for (k = 0; k < nprobe; k++) padj[(size_t)k*(size_t)D.p+(size_t)j] = cola[k];
            }
            for (j = 0; j < D.nvar; j++) {
                for (k = 0; k < nprobe; k++) colb[k] = fpval[(size_t)k*(size_t)D.nvar+(size_t)j];
                bh_adjust(colb, nprobe, cola, pr);
                for (k = 0; k < nprobe; k++) fpadj[(size_t)k*(size_t)D.nvar+(size_t)j] = cola[k];
            }
        }
        free(pr); free(colb); free(cola);
    }

    /* ---- emit TSV ---- */
    fputs("Probe_ID\tN", stdout);
    for (j = 0; j < D.p; j++) printf("\tEst_%s", D.coef[j]);
    for (j = 0; j < D.p; j++) printf("\tPval_%s", D.coef[j]);
    for (j = 0; j < D.nvar; j++) printf("\tFPval_%s", D.vname[j]);
    for (j = 0; j < D.nvar; j++) printf("\tEff_%s", D.vname[j]);
    for (j = 0; j < D.p; j++) printf("\tPadj_%s", D.coef[j]);
    for (j = 0; j < D.nvar; j++) printf("\tFPadj_%s", D.vname[j]);
    putchar('\n');
    for (i = 0; i < nprobe; i++) {
        printf("%s\t%d", ids[i], nobs[i]);
        for (j = 0; j < D.p; j++) { double x=est[(size_t)i*(size_t)D.p+(size_t)j];   if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.p; j++) { double x=pval[(size_t)i*(size_t)D.p+(size_t)j];  if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.nvar; j++){ double x=fpval[(size_t)i*(size_t)D.nvar+(size_t)j];if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        for (j = 0; j < D.nvar; j++){ double x=eff[(size_t)i*(size_t)D.nvar+(size_t)j]; if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        for (j = 0; j < D.p; j++) { double x=padj[(size_t)i*(size_t)D.p+(size_t)j];  if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x); }
        for (j = 0; j < D.nvar; j++){ double x=fpadj[(size_t)i*(size_t)D.nvar+(size_t)j];if(isnan(x))fputs("\tNA",stdout);else printf("\t%.10g",x);}
        putchar('\n');
    }
    rc = 0;

out:
    free(est); free(pval); free(fpval); free(eff); free(padj); free(fpadj); free(nobs);
    design_free(&D);
    if (ids) { for (i=0;i<nprobe;i++) free(ids[i]); free(ids); }
    if (samp) { for (i=0;i<nsamp;i++) free(samp[i]); free(samp); }
    free(betas);
    return rc;
}

static int cmd_idat_dump(int argc, char **argv)
{
    const char *path = NULL;
    long head = -1;
    int tsv = 0, i;
    sesame_idat_t *d = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--tsv") == 0) {
            tsv = 1;
        } else if (strcmp(argv[i], "--head") == 0 && i + 1 < argc) {
            head = strtol(argv[++i], NULL, 10);
        } else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            return usage();
        } else {
            path = argv[i];
        }
    }
    if (!path) return usage();

    if (sesame_idat_read(path, &d, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s: %s\n", path, e.msg);
        return 1;
    }

    if (tsv) {
        long lim = (head < 0 || head > d->n) ? d->n : head;
        for (long k = 0; k < lim; k++)
            printf("%u\t%u\t%u\t%u\n", d->addr[k], d->mean[k],
                   d->sd[k], d->nbeads[k]);
    } else {
        long lim = (head < 0) ? 5 : (head > d->n ? d->n : head);
        printf("file       %s\n", path);
        printf("version    %lld\n", (long long)d->version);
        printf("nFields    %d\n", d->n_fields);
        printf("nSNPsRead  %d\n", d->n);
        printf("\n%-12s %8s %8s %7s\n", "addr", "mean", "sd", "nbeads");
        for (long k = 0; k < lim; k++)
            printf("%-12u %8u %8u %7u\n", d->addr[k], d->mean[k],
                   d->sd[k], d->nbeads[k]);
        if (lim < d->n) printf("... (%d more)\n", d->n - (int)lim);
    }

    sesame_idat_free(d);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) return usage();
    if (strcmp(argv[1], "idat-dump") == 0)
        return cmd_idat_dump(argc - 2, argv + 2);
    if (strcmp(argv[1], "betas") == 0)
        return cmd_betas(argc - 2, argv + 2);
    if (strcmp(argv[1], "intensity") == 0)
        return cmd_intensity(argc - 2, argv + 2);
    if (strcmp(argv[1], "qc") == 0)
        return cmd_qc(argc - 2, argv + 2);
    if (strcmp(argv[1], "dml") == 0)
        return cmd_dml(argc - 2, argv + 2);
    if (strcmp(argv[1], "fetch") == 0)
        return cmd_fetch(argc - 2, argv + 2);
    if (strcmp(argv[1], "index-info") == 0)
        return cmd_index_info(argc - 2, argv + 2);
    if (strcmp(argv[1], "version") == 0) {
        puts("sesame 0.0.1-dev");
        return 0;
    }
    return usage();
}
