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
      "  sesame preprocess [--output LIST] [--prep QCDPB] [--raw-signal]\n"
      "            [--index <ordering.tsv.gz>] [--platform P] [--min-beads N]\n"
      "            [--out DIR] [--tmp DIR] [--threads N] <prefix> [<prefix> ...]\n"
      "      The pipeline: apply --prep (default QCDPB) per sample and write YAME\n"
      "      .cg outputs to DIR (default .), one indexed file per output over the\n"
      "      cohort. --output is a comma list from intensity,total_intensity,beta,\n"
      "      pval,qc (default beta,intensity,pval,qc): beta.cg (fmt4), intensity.cg\n"
      "      (fmt3 M/U -- yame derives beta and coverage), total_intensity.cg\n"
      "      (fmt4), pval.cg (fmt4), qc.tsv. Signal outputs reflect the prep;\n"
      "      --raw-signal takes intensity/total from the raw signal. Many prefixes\n"
      "      -> one indexed .cg per output; a failed sample -> NA + exit 1. If\n"
      "      --index is omitted the platform is detected from the bead count and\n"
      "      the index looked up in $SESAME_INDEX_DIR, ., ./data, then the cache.\n"
      "\n"
      "  sesame dml --betas <beta.cg|matrix.tsv> [--index <ordering.tsv.gz>]\n"
      "             (--formula '~a+b' --meta <s.tsv> | --design <X.tsv>) [--threads N]\n"
      "      Per-probe differential methylation: OLS of each probe's betas on the\n"
      "      design, with per-coefficient t-tests, a holdout F-test per categorical\n"
      "      variable, effect sizes, and BH-adjusted p-values. --betas is a\n"
      "      preprocess beta.cg (needs --index for probe IDs) or a Probe_ID matrix\n"
      "      TSV; --meta's first column matches the sample names. --formula takes\n"
      "      main effects (categorical auto-dummied); --design passes a numeric\n"
      "      design for interactions.\n"
      "\n"
      "  sesame cnv --target <total_intensity.cg> [--normals <cnvnormals.cg>]\n"
      "             [--platform P | --index <ordering>] [--genome hg38]\n"
      "             [--coords <coord.tsv.gz>] [--probes|--bins] [--out <file>]\n"
      "      Copy-number: regress the target's per-probe total intensity on a\n"
      "      panel of normals (OLS), then log2(target/fitted) per probe, binned\n"
      "      along the genome (median per bin, default --bins). normals/coords/\n"
      "      genome default to the fetched store for --platform. CBS segmentation\n"
      "      is not included. --target may hold several samples (one profile each).\n"
      "\n"
      "  sesame attach-probe [--index <ordering.tsv.gz> | --platform P]\n"
      "                      [--all] [--beta] [--no-header] <file>\n"
      "      Prepend the ordering's Probe_ID to a positional file's rows, as TSV\n"
      "      on stdout. <file> is a YAME .cg/.cm/.cx (fmt0 mask, fmt3 M/U or\n"
      "      --beta, fmt4 float) or a text .tsv[.gz] (e.g. a .hg38.coord.tsv.gz).\n"
      "      --all emits every sample column; the ordering is --index, else\n"
      "      --platform, else inferred from the filename prefix. Row count must\n"
      "      match the ordering (same platform + tag that made the file).\n"
      "\n"
      "  sesame idat-dump [--head N] [--tsv] <file.idat|file.idat.gz>\n"
      "      Print IDAT contents. Default prints a summary header; --tsv emits\n"
      "      addr<TAB>mean<TAB>sd<TAB>nbeads. --head N limits records.\n"
      "\n"
      "  sesame fetch [<platform>] [--force]\n"
      "  sesame fetch genome [<build>] [--force]\n"
      "      Download <platform> (default: all published), or a genome build's\n"
      "      genome-level annotation (default: hg38), at the pinned tag. Verifies\n"
      "      every file against its digest; a file already present and matching is\n"
      "      skipped. The ONLY path that touches the network. Never prompts.\n"
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
    const char *platform = NULL, *genome = NULL;
    int force = 0, want_genome = 0, i;
    char path[4096];
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--force") == 0) force = 1;
        else if (strcmp(argv[i], "genome") == 0) want_genome = 1;
        else if (argv[i][0] == '-' && argv[i][1] != '\0') return usage();
        else if (want_genome && !genome) genome = argv[i];
        else platform = argv[i];
    }

    if (want_genome) {                    /* sesame fetch genome [<build>] */
        if (sesame_fetch_genome(genome ? genome : "hg38", force, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg);
            return 1;
        }
    } else if (platform) {
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

/* ---------------------------------------------------------- preprocess ---
 *
 * The one pipeline command: apply --prep (default QCDPB) per sample and emit any
 * of {intensity, total_intensity, beta, pval} as YAME .cg (one indexed file per
 * output over the whole cohort) plus a qc.tsv. Signal-derived outputs reflect the
 * prep by default; --raw-signal takes intensity/total from the raw SigDF. */

typedef struct {
    char           **prefixes;
    int32_t          nsamp, n;
    const sesame_index_t *ix;
    const char      *plat, *prep;
    const uint8_t   *qmask, *bgmask;
    int32_t          qn, bgn;
    int              min_beads, raw_signal;
    int              want_beta, want_int, want_tot, want_pval, want_qc;
    double          *matBeta, *matM, *matU, *matTot, *matPval;   /* NULL if unwanted */
    sesame_qc_t     *qcres;
    pthread_mutex_t  lock;
    int32_t          next, nfail;
    uint32_t         status_or;
} pp_ctx;

static void fill_na(double *col, int32_t n) { int32_t k; if (col) for (k = 0; k < n; k++) col[k] = NAN; }

static void *pp_worker(void *arg)
{
    pp_ctx *c = (pp_ctx *)arg;
    int32_t n = c->n;
    for (;;) {
        int32_t j;
        sesame_sigdf_t *raw = NULL;
        sesame_err_t e;
        uint32_t st = 0;
        int rc;
        double *beta, *M, *U, *tot, *pval;

        pthread_mutex_lock(&c->lock);
        j = c->next < c->nsamp ? c->next++ : -1;
        pthread_mutex_unlock(&c->lock);
        if (j < 0) break;

        beta = c->matBeta ? c->matBeta + (size_t)j*(size_t)n : NULL;
        M    = c->matM    ? c->matM    + (size_t)j*(size_t)n : NULL;
        U    = c->matU    ? c->matU    + (size_t)j*(size_t)n : NULL;
        tot  = c->matTot  ? c->matTot  + (size_t)j*(size_t)n : NULL;
        pval = c->matPval ? c->matPval + (size_t)j*(size_t)n : NULL;

        rc = build_sigdf_for(c->prefixes[j], c->ix, c->plat, c->min_beads, &raw, &e);
        if (rc == SESAME_OK) st = raw->status;
        if (rc == SESAME_OK && c->want_qc)
            rc = sesame_qc_calc(raw, c->bgmask, c->bgn, &c->qcres[j], &e);
        if (rc == SESAME_OK)
            rc = sesame_pipeline(raw, c->prep, c->qmask, c->qn, c->bgmask, c->bgn,
                                 c->raw_signal, beta, M, U, tot, pval, NULL, &e);
        sesame_sigdf_free(raw);

        if (rc != SESAME_OK) {
            fill_na(beta, n); fill_na(M, n); fill_na(U, n); fill_na(tot, n); fill_na(pval, n);
        }
        pthread_mutex_lock(&c->lock);
        if (rc != SESAME_OK) {
            c->nfail++;
            fprintf(stderr, "sesame: %s: %s\n", path_basename(c->prefixes[j]), e.msg);
        } else c->status_or |= st;
        pthread_mutex_unlock(&c->lock);
    }
    return NULL;
}

static int pp_parse_output(const char *s, pp_ctx *c)
{
    char buf[256], *t;
    snprintf(buf, sizeof buf, "%s", s);
    c->want_beta = c->want_int = c->want_tot = c->want_pval = c->want_qc = 0;
    for (t = strtok(buf, ","); t; t = strtok(NULL, ",")) {
        if      (!strcmp(t, "beta"))            c->want_beta = 1;
        else if (!strcmp(t, "intensity"))       c->want_int = 1;
        else if (!strcmp(t, "total_intensity")) c->want_tot = 1;
        else if (!strcmp(t, "pval"))            c->want_pval = 1;
        else if (!strcmp(t, "qc"))              c->want_qc = 1;
        else { fprintf(stderr, "sesame: unknown --output '%s'\n", t); return -1; }
    }
    return 0;
}

static int cmd_preprocess(int argc, char **argv)
{
    const char *idxpath = NULL, *platform = NULL, *plat = NULL, *prep = "QCDPB";
    const char *outlist = "beta,intensity,pval,qc", *outdir = ".", *tmpdir = NULL;
    int min_beads = 0, nthreads = 0, raw_signal = 0, i, rc = 1;
    char **prefixes = NULL, **names = NULL;
    int32_t nsamp = 0, n = 0, qn = 0, bgn = 0;
    char resolved[4096], path[4096];
    sesame_index_t *ix = NULL;
    uint8_t *qmask = NULL, *bgmask = NULL;
    sesame_qc_t *qcres = NULL;
    betas_matrix mB={NULL,0,-1}, mM={NULL,0,-1}, mU={NULL,0,-1}, mT={NULL,0,-1}, mP={NULL,0,-1};
    pp_ctx ctx; memset(&ctx, 0, sizeof ctx);
    sesame_err_t e;

    prefixes = (char **)malloc((size_t)(argc + 1) * sizeof(char *));
    if (!prefixes) { fprintf(stderr, "sesame: out of memory\n"); return 1; }
    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--prep") == 0 && i+1 < argc) prep = argv[++i];
        else if (strcmp(argv[i], "--output") == 0 && i+1 < argc) outlist = argv[++i];
        else if (strcmp(argv[i], "--out") == 0 && i+1 < argc) outdir = argv[++i];
        else if (strcmp(argv[i], "--tmp") == 0 && i+1 < argc) tmpdir = argv[++i];
        else if (strcmp(argv[i], "--min-beads") == 0 && i+1 < argc) min_beads = (int)strtol(argv[++i],NULL,10);
        else if ((strcmp(argv[i],"--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc) nthreads = (int)strtol(argv[++i],NULL,10);
        else if (strcmp(argv[i], "--raw-signal") == 0) raw_signal = 1;
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]); free(prefixes); return usage();
        } else prefixes[nsamp++] = argv[i];
    }
    if (nsamp == 0) { free(prefixes); return usage(); }
    if (pp_parse_output(outlist, &ctx) != 0) { free(prefixes); return 1; }
    if (tmpdir) setenv("TMPDIR", tmpdir, 1);
    mkdir(outdir, 0777);                              /* ok if it already exists */

    if (!idxpath) {
        if (!platform) {
            char gp[4096]; sesame_idat_t *g0 = NULL; int32_t beads;
            if (resolve_idat(prefixes[0], "Grn", gp, sizeof gp)) {
                fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefixes[0]); goto out; }
            if (sesame_idat_read(gp, &g0, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
            beads = g0->n; platform = sesame_platform_from_beads(beads); sesame_idat_free(g0);
            if (!platform) { fprintf(stderr, "sesame: cannot identify platform from %d beads.\n"
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

    /* masks: Q needs qmask; P/B/pval/qc need bgmask */
    if (strchr(prep, 'Q')) {
        if (!plat) { fprintf(stderr, "sesame: Q needs a known platform; pass --platform with --index\n"); goto out; }
        if (sesame_quality_mask(plat, &qmask, &qn, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    }
    if (strchr(prep,'P') || strchr(prep,'B') || ctx.want_pval || ctx.want_qc) {
        if (!plat) { fprintf(stderr, "sesame: P/B/pval/qc need a known platform; pass --platform with --index\n"); goto out; }
        if (sesame_background_mask(plat, &bgmask, &bgn, &e) != SESAME_OK) { fprintf(stderr, "sesame: %s\n", e.msg); goto out; }
    }

    if ((ctx.want_beta && betas_matrix_init(&mB, n, nsamp)) ||
        (ctx.want_int  && (betas_matrix_init(&mM, n, nsamp) || betas_matrix_init(&mU, n, nsamp))) ||
        (ctx.want_tot  && betas_matrix_init(&mT, n, nsamp)) ||
        (ctx.want_pval && betas_matrix_init(&mP, n, nsamp))) {
        fprintf(stderr, "sesame: cannot allocate output matrices\n"); goto out;
    }
    if (ctx.want_qc) {
        size_t k, nd = sizeof(sesame_qc_t)/sizeof(double); double *q;
        qcres = (sesame_qc_t *)malloc((size_t)nsamp * sizeof(sesame_qc_t));
        if (!qcres) { fprintf(stderr, "sesame: oom\n"); goto out; }
        q = (double *)qcres; for (k = 0; k < (size_t)nsamp*nd; k++) q[k] = NAN;
    }

    ctx.prefixes=prefixes; ctx.nsamp=nsamp; ctx.n=n; ctx.ix=ix; ctx.plat=plat; ctx.prep=prep;
    ctx.qmask=qmask; ctx.qn=qn; ctx.bgmask=bgmask; ctx.bgn=bgn;
    ctx.min_beads=min_beads; ctx.raw_signal=raw_signal;
    ctx.matBeta=mB.data; ctx.matM=mM.data; ctx.matU=mU.data; ctx.matTot=mT.data; ctx.matPval=mP.data;
    ctx.qcres=qcres; ctx.next=0; ctx.nfail=0; ctx.status_or=0;
    pthread_mutex_init(&ctx.lock, NULL);
    {
        long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
        int T = nthreads > 0 ? nthreads : (ncpu > 0 ? (int)ncpu : 1), t, spawned = 0;
        pthread_t *th = NULL;
        if (T > nsamp) T = nsamp;
        if (T < 1) T = 1;
        if (T > 1) th = (pthread_t *)malloc((size_t)(T-1)*sizeof(pthread_t));
        if (th) for (t=0;t<T-1;t++) if (pthread_create(&th[spawned],NULL,pp_worker,&ctx)==0) spawned++;
        pp_worker(&ctx);
        for (t=0;t<spawned;t++) pthread_join(th[t],NULL);
        free(th);
    }
    pthread_mutex_destroy(&ctx.lock);

    names = (char **)malloc((size_t)nsamp * sizeof(char *));
    if (!names) { fprintf(stderr, "sesame: oom\n"); goto out; }
    for (i = 0; i < nsamp; i++) names[i] = (char *)path_basename(prefixes[i]);

    if (ctx.want_beta) { snprintf(path,sizeof path,"%s/beta.cg",outdir);
        if (sesame_write_cg(path, mB.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_int)  { snprintf(path,sizeof path,"%s/intensity.cg",outdir);
        if (sesame_write_cg_mu(path, mM.data, mU.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_tot)  { snprintf(path,sizeof path,"%s/total_intensity.cg",outdir);
        if (sesame_write_cg(path, mT.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_pval) { snprintf(path,sizeof path,"%s/pval.cg",outdir);
        if (sesame_write_cg(path, mP.data, n, nsamp, names, &e)) { fprintf(stderr,"sesame: %s\n",e.msg); goto out; } }
    if (ctx.want_qc) {
        FILE *qf; snprintf(path,sizeof path,"%s/qc.tsv",outdir);
        if (!(qf = fopen(path, "w"))) { fprintf(stderr, "sesame: cannot write %s\n", path); goto out; }
        fprintf(qf, "sample\t%s\n", sesame_qc_header());
        for (i = 0; i < nsamp; i++) {
            char row[8192];
            if (sesame_qc_format_row(&qcres[i], row, sizeof row) < 0) { fclose(qf); fprintf(stderr,"sesame: qc row too long\n"); goto out; }
            fprintf(qf, "%s\t%s\n", names[i], row);
        }
        fclose(qf);
    }
    fprintf(stderr, "sesame: preprocess wrote %s to %s/ (%d sample%s, prep=%s)\n",
            outlist, outdir, nsamp, nsamp==1?"":"s", *prep?prep:"\"\"");
    if (ctx.nfail) fprintf(stderr, "sesame: %d of %d samples failed (NA)\n", ctx.nfail, nsamp);
    rc = ctx.nfail ? 1 : 0;

out:
    betas_matrix_free(&mB); betas_matrix_free(&mM); betas_matrix_free(&mU);
    betas_matrix_free(&mT); betas_matrix_free(&mP);
    free(qcres); free(names); free(qmask); free(bgmask); free(prefixes);
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
    const char *idxpath = NULL;
    int nthreads = 0, i, j, rc = 1;
    char **ids = NULL, **samp = NULL;
    double *betas = NULL;
    int32_t nprobe = 0, nsamp = 0;
    design_t D; memset(&D, 0, sizeof D);
    double *est = NULL, *pval = NULL, *fpval = NULL, *eff = NULL, *padj = NULL, *fpadj = NULL;
    int32_t *nobs = NULL;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--betas") == 0 && i+1 < argc) betapath = argv[++i];
        else if (strcmp(argv[i], "--meta") == 0 && i+1 < argc) metapath = argv[++i];
        else if (strcmp(argv[i], "--formula") == 0 && i+1 < argc) formula = argv[++i];
        else if (strcmp(argv[i], "--design") == 0 && i+1 < argc) designpath = argv[++i];
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) idxpath = argv[++i];
        else if ((strcmp(argv[i], "--threads")==0||strcmp(argv[i],"-t")==0) && i+1<argc)
            nthreads = (int)strtol(argv[++i], NULL, 10);
        else { fprintf(stderr, "sesame: unknown/incomplete option %s\n", argv[i]); return usage(); }
    }
    if (!betapath || (!formula && !designpath)) {
        fprintf(stderr, "sesame: dml needs --betas and either --formula (with --meta) or --design\n");
        return usage();
    }
    { size_t bl = strlen(betapath);
      if (bl >= 3 && strcmp(betapath + bl - 3, ".cg") == 0) {   /* YAME beta.cg */
        double *smat = NULL; int32_t np = 0, ns = 0, s, p; sesame_index_t *ix;
        if (!idxpath) { fprintf(stderr, "sesame: dml on a .cg needs --index (for probe IDs)\n"); return 1; }
        if (sesame_read_cg(betapath, &smat, &np, &ns, &samp, &e)) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
        if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); free(smat); return 1; }
        if (sesame_index_nprobes(ix) != np) {
            fprintf(stderr, "sesame: %s has %d probes, index has %d\n", betapath, np, sesame_index_nprobes(ix));
            sesame_index_close(ix); free(smat); return 1; }
        nprobe = np; nsamp = ns;
        ids = (char **)malloc((size_t)np * sizeof(char *));
        betas = (double *)malloc((size_t)np * (size_t)ns * sizeof(double));
        for (p = 0; p < np; p++) ids[p] = strdup(sesame_index_probe_id(ix, p));
        for (s = 0; s < ns; s++) for (p = 0; p < np; p++)   /* sample-major -> probe-major */
            betas[(size_t)p*(size_t)ns + (size_t)s] = smat[(size_t)s*(size_t)np + (size_t)p];
        free(smat); sesame_index_close(ix);
      } else if (load_betas(betapath, &ids, &betas, &nprobe, &nsamp, &samp) != 0) {
        fprintf(stderr, "sesame: cannot read betas matrix %s\n", betapath); return 1;
      }
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
                { char *ep = t + strlen(t); while (ep>t && ep[-1]==' ') *--ep = '\0'; }
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

static int cmd_cnv(int argc, char **argv)
{
    const char *target = NULL, *normals = NULL, *idxpath = NULL, *platform = NULL;
    const char *coords = NULL, *genome = "hg38", *outpath = NULL;
    int tilewidth = 50000, min_probes = 20, probes = 0, i, rc = 1;
    char store[4096], resolved[4096], cbuf[4096], nbuf[4096], sbuf[4096], gbuf[4096];
    sesame_index_t *ix = NULL;
    sesame_cnv_t *res = NULL;
    int32_t nres = 0, s, k;
    FILE *out = stdout;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--target") == 0 && i+1 < argc) target = argv[++i];
        else if (strcmp(argv[i], "--normals") == 0 && i+1 < argc) normals = argv[++i];
        else if (strcmp(argv[i], "--index") == 0 && i+1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--coords") == 0 && i+1 < argc) coords = argv[++i];
        else if (strcmp(argv[i], "--genome") == 0 && i+1 < argc) genome = argv[++i];
        else if (strcmp(argv[i], "--tilewidth") == 0 && i+1 < argc) tilewidth = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--min-probes") == 0 && i+1 < argc) min_probes = (int)strtol(argv[++i], NULL, 10);
        else if (strcmp(argv[i], "--probes") == 0) probes = 1;
        else if (strcmp(argv[i], "--bins") == 0) probes = 0;
        else if (strcmp(argv[i], "--out") == 0 && i+1 < argc) outpath = argv[++i];
        else if (argv[i][0] == '-' && argv[i][1] != '\0') { fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage(); }
        else if (!target) target = argv[i];
        else { fprintf(stderr, "sesame: unexpected argument %s\n", argv[i]); return usage(); }
    }
    if (!target) { fprintf(stderr, "sesame: cnv needs --target <total_intensity.cg>\n"); return usage(); }
    if (tilewidth < 1 || min_probes < 1) { fprintf(stderr, "sesame: --tilewidth and --min-probes must be positive\n"); return 1; }

    /* Defaults for normals/coords/index come from the store, keyed on platform. */
    sesame_store_dir(store, sizeof store);
    if (!idxpath) {
        if (!platform) { fprintf(stderr, "sesame: cnv needs --platform (or --index)\n"); return 1; }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024]; sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help); return 1;
        }
        idxpath = resolved;
    }
    if (!normals) {
        if (!platform) { fprintf(stderr, "sesame: cnv needs --normals or --platform\n"); return 1; }
        snprintf(nbuf, sizeof nbuf, "%s/%s/%s.cnvnormals.cg", store, platform, platform);
        normals = nbuf;
    }
    if (!coords) {
        if (!platform) { fprintf(stderr, "sesame: cnv needs --coords or --platform\n"); return 1; }
        snprintf(cbuf, sizeof cbuf, "%s/%s/%s.%s.coord.tsv.gz", store, platform, platform, genome);
        coords = cbuf;
    }
    if (sesame_genome_locate(genome, "seqinfo.tsv.gz", sbuf, sizeof sbuf) != 0 ||
        sesame_genome_locate(genome, "gaps.tsv.gz", gbuf, sizeof gbuf) != 0) {
        fprintf(stderr, "sesame: no genome info for %s in the store\n"
            "  fix: sesame fetch genome %s\n", genome, genome);
        return 1;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) { fprintf(stderr, "sesame: %s\n", e.msg); return 1; }
    if (sesame_cnv_run(target, normals, ix, coords, sbuf, gbuf,
                       tilewidth, min_probes, &res, &nres, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    if (outpath && !(out = fopen(outpath, "w"))) {
        fprintf(stderr, "sesame: cannot write %s\n", outpath); goto out;
    }
    if (probes) {
        fputs("sample\tProbe_ID\tchrom\tpos\tlog2ratio\n", out);
        for (s = 0; s < nres; s++)
            for (k = 0; k < res[s].n_probe; k++)
                fprintf(out, "%s\t%s\t%s\t%d\t%.6g\n", res[s].sample, res[s].probe_id[k],
                        res[s].chrom[k], res[s].pos[k], res[s].log2ratio[k]);
    } else {
        fputs("sample\tchrom\tstart\tend\tnprobes\tlog2ratio\n", out);
        for (s = 0; s < nres; s++)
            for (k = 0; k < res[s].n_bin; k++)
                fprintf(out, "%s\t%s\t%d\t%d\t%d\t%.6g\n", res[s].sample, res[s].bin_chrom[k],
                        res[s].bin_start[k], res[s].bin_end[k], res[s].bin_nprobe[k], res[s].bin_log2ratio[k]);
    }
    for (s = 0; s < nres; s++)
        fprintf(stderr, "sesame: cnv %s -- %d probes on %d normals, %d bins\n",
                res[s].sample, res[s].n_used, res[s].n_normal, res[s].n_bin);
    rc = 0;
out:
    if (out && out != stdout) fclose(out);
    sesame_cnv_free_array(res, nres);
    sesame_index_close(ix);
    return rc;
}

/* Infer a platform from a filename like "MSA.hg38.coord.tsv.gz" -> "MSA". Returns
 * a static registry string, or NULL. */
static const char *platform_from_basename(const char *path)
{
    static const char *plats[] = { "EPICv2", "EPIC", "HM450", "MSA", NULL };
    const char *base = path_basename(path);
    for (int i = 0; plats[i]; i++) {
        size_t l = strlen(plats[i]);
        if (strncmp(base, plats[i], l) == 0 && base[l] == '.') return plats[i];
    }
    return NULL;
}

static int cmd_attach_probe(int argc, char **argv)
{
    const char *path = NULL, *idxpath = NULL, *platform = NULL;
    sesame_attach_opt_t opt; memset(&opt, 0, sizeof opt);
    char resolved[4096];
    sesame_index_t *ix = NULL;
    sesame_err_t e;
    int i, rc = 1;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i+1 < argc) idxpath = argv[++i];
        else if (strcmp(argv[i], "--platform") == 0 && i+1 < argc) platform = argv[++i];
        else if (strcmp(argv[i], "--all") == 0 || strcmp(argv[i], "-a") == 0) opt.all = 1;
        else if (strcmp(argv[i], "--beta") == 0) opt.beta = 1;
        else if (strcmp(argv[i], "--no-header") == 0) opt.no_header = 1;
        else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]); return usage();
        } else path = argv[i];
    }
    if (!path) return usage();

    /* Probe IDs come from the ordering: --index wins, else --platform, else
     * inferred from the filename prefix (e.g. MSA.hg38.coord.tsv.gz). */
    if (!idxpath) {
        if (!platform) platform = platform_from_basename(path);
        if (!platform) {
            fprintf(stderr, "sesame: cannot tell which ordering to use for %s\n"
                "  pass --index <ordering.tsv.gz> or --platform <EPIC|EPICv2|HM450|MSA>\n",
                path);
            return 1;
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            char help[1024];
            sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help);
            return 1;
        }
        idxpath = resolved;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); return 1;
    }
    if (sesame_attach_probe(path, ix, &opt, stdout, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    rc = 0;
out:
    sesame_index_close(ix);
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
    if (strcmp(argv[1], "preprocess") == 0)
        return cmd_preprocess(argc - 2, argv + 2);
    if (strcmp(argv[1], "dml") == 0)
        return cmd_dml(argc - 2, argv + 2);
    if (strcmp(argv[1], "cnv") == 0)
        return cmd_cnv(argc - 2, argv + 2);
    if (strcmp(argv[1], "attach-probe") == 0)
        return cmd_attach_probe(argc - 2, argv + 2);
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
