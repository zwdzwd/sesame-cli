/* main.c -- sesame command line interface
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "sesame.h"
#include "../src/internal.h"

#include <dirent.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
      "  sesame betas [--index <ordering.tsv.gz>] [--platform P]\n"
      "                [--min-beads N] [--no-mask] [--f64] <prefix>\n"
      "      Compute beta values from <prefix>_Grn.idat[.gz] and\n"
      "      <prefix>_Red.idat[.gz]. Emits Probe_ID<TAB>beta (NA for missing).\n"
      "      Equivalent to openSesame(prefix, prep=\"\", func=getBetas).\n"
      "      NOTE: no preprocessing yet -- QCDPB is not implemented.\n"
      "      --min-beads N  mask probes with any bead count < N (default: off)\n"
      "      --no-mask      ignore the mask column; emit every beta\n"
      "      --prep CODE    preprocessing steps to apply (Q, C, D, P so far)\n"
      "      --dump-col     emit Probe_ID<TAB>col (G/R/2) instead of betas,\n"
      "                     for differential-testing a prep step\n"
      "      --f64          write raw little-endian float64 to stdout instead\n"
      "                     of text (NA as NaN). For lossless comparison: R's\n"
      "                     text parser does not correctly round 17-digit\n"
      "                     decimals, so text round-trips lose up to 1 ULP.\n"
      "      If --index is omitted the platform is detected from the IDAT bead\n"
      "      count and the index looked up in $SESAME_INDEX_DIR, ., ./data,\n"
      "      then the cache. sesame never downloads implicitly.\n"
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

static int cmd_betas(int argc, char **argv)
{
    const char *idxpath = NULL, *prefix = NULL, *platform = NULL, *prep = "";
    int min_beads = 0, apply_mask = 1, f64 = 0, dump_col = 0, i, rc = 1;
    char gpath[4096], rpath[4096], resolved[4096];
    sesame_index_t *ix = NULL;
    sesame_idat_t *g = NULL, *r = NULL;
    sesame_sigdf_t *s = NULL;
    double *betas = NULL;
    uint8_t *qmask = NULL, *bgmask = NULL;
    int32_t qn = 0, bgn = 0;
    sesame_err_t e;

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "--index") == 0 && i + 1 < argc) {
            idxpath = argv[++i];
        } else if (strcmp(argv[i], "--platform") == 0 && i + 1 < argc) {
            platform = argv[++i];
        } else if (strcmp(argv[i], "--min-beads") == 0 && i + 1 < argc) {
            min_beads = (int)strtol(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--no-mask") == 0) {
            apply_mask = 0;
        } else if (strcmp(argv[i], "--prep") == 0 && i + 1 < argc) {
            prep = argv[++i];
        } else if (strcmp(argv[i], "--dump-col") == 0) {
            dump_col = 1;
        } else if (strcmp(argv[i], "--f64") == 0) {
            f64 = 1;
        } else if (argv[i][0] == '-' && argv[i][1] != '\0') {
            fprintf(stderr, "sesame: unknown option %s\n", argv[i]);
            return usage();
        } else {
            prefix = argv[i];
        }
    }
    if (!prefix) return usage();

    if (resolve_idat(prefix, "Grn", gpath, sizeof gpath)) {
        fprintf(stderr, "sesame: Grn IDAT does not exist for %s\n", prefix);
        return 1;
    }
    if (resolve_idat(prefix, "Red", rpath, sizeof rpath)) {
        fprintf(stderr, "sesame: Red IDAT does not exist for %s\n", prefix);
        return 1;
    }

    if (sesame_idat_read(gpath, &g, &e) != SESAME_OK ||
        sesame_idat_read(rpath, &r, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    /* Resolve the index: explicit path > explicit platform > bead-count
     * detection. Never downloads, never prompts -- see src/cache.c. */
    if (!idxpath) {
        char help[1024];
        if (!platform) {
            platform = sesame_platform_from_beads(g->n);
            if (!platform) {
                fprintf(stderr,
                    "sesame: cannot identify platform from %d beads.\n"
                    "  pass --platform <P> or --index <path>.\n", g->n);
                goto out;
            }
            fprintf(stderr, "sesame: detected %s (%d beads)\n", platform, g->n);
        }
        if (sesame_index_locate(platform, resolved, sizeof resolved) != 0) {
            sesame_index_missing_help(platform, help, sizeof help);
            fprintf(stderr, "sesame: %s\n", help);
            goto out;
        }
        idxpath = resolved;
    }

    if (!(ix = sesame_index_open(idxpath, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }
    if (!(s = sesame_sigdf_from_idats(g, r, ix, min_beads, &e))) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    /* Q needs the recommended mask for the platform; load it once (it shells
     * out to yame) and reuse across every step/sample. */
    if (strchr(prep, 'Q') || strchr(prep, 'P')) {
        const char *plat = platform ? platform
                                    : sesame_platform_from_beads(g->n);
        if (strchr(prep, 'Q') &&
            sesame_quality_mask(plat, &qmask, &qn, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out;
        }
        if (strchr(prep, 'P') &&
            sesame_background_mask(plat, &bgmask, &bgn, &e) != SESAME_OK) {
            fprintf(stderr, "sesame: %s\n", e.msg); goto out;
        }
    }

    /* Apply prep steps in the order given. */
    for (const char *c = prep; *c; c++) {
        switch (*c) {
        case 'Q':
            if (sesame_prep_quality_mask(s, qmask, qn, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            break;
        case 'C':
            if (sesame_prep_infer_channel(s, 0, 0, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            break;
        case 'D':
            if (sesame_prep_dye_bias_nl(s, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            break;
        case 'P':
            if (sesame_prep_poobah(s, bgmask, bgn, 0.05, 1, &e) != SESAME_OK) {
                fprintf(stderr, "sesame: %s\n", e.msg); goto out;
            }
            break;
        default:
            fprintf(stderr, "sesame: prep code '%c' not implemented "
                            "(have: Q, C, D, P)\n", *c);
            goto out;
        }
    }

    if (dump_col) {
        static const char *nm[] = { "2", "G", "R" };
        for (int32_t k = 0; k < s->n; k++)
            printf("%s\t%s\n", sesame_index_probe_id(ix, k), nm[s->col[k]]);
        rc = 0;
        goto out;
    }

    betas = (double *)malloc((size_t)s->n * sizeof(double));
    if (!betas) { fprintf(stderr, "sesame: out of memory\n"); goto out; }
    if (sesame_get_betas(s, apply_mask, betas, &e) != SESAME_OK) {
        fprintf(stderr, "sesame: %s\n", e.msg); goto out;
    }

    if (f64) {
        /* Raw float64 so comparisons are lossless. */
        if (fwrite(betas, sizeof(double), (size_t)s->n, stdout) != (size_t)s->n) {
            fprintf(stderr, "sesame: short write\n"); goto out;
        }
    } else {
        for (int32_t k = 0; k < s->n; k++) {
            const char *id = sesame_index_probe_id(ix, k);
            if (isnan(betas[k])) printf("%s\tNA\n", id);
            else                 printf("%s\t%.17g\n", id, betas[k]);
        }
    }

    if (s->status & SESAME_STAT_ADDR_MISSING)
        fprintf(stderr, "sesame: note: %d probes had addresses absent from the IDAT\n",
                s->n_addr_missing);
    rc = 0;

out:
    free(qmask);
    free(bgmask);
    free(betas);
    sesame_sigdf_free(s);
    sesame_idat_free(g);
    sesame_idat_free(r);
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
    if (strcmp(argv[1], "betas") == 0)
        return cmd_betas(argc - 2, argv + 2);
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
