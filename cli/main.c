/* main.c -- sesamec command line interface
 * SPDX-License-Identifier: MIT
 */
#include "sesame.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int usage(void)
{
    fputs(
      "sesamec -- Infinium DNA methylation preprocessing\n"
      "\n"
      "usage:\n"
      "  sesamec idat-dump [--head N] [--tsv] <file.idat|file.idat.gz>\n"
      "      Print IDAT contents. Default prints a summary header;\n"
      "      --tsv emits addr<TAB>mean<TAB>sd<TAB>nbeads with no header.\n"
      "      --head N limits to the first N records (default: all for --tsv,\n"
      "      5 for summary).\n"
      "\n"
      "  sesamec version\n",
      stderr);
    return 2;
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
            fprintf(stderr, "sesamec: unknown option %s\n", argv[i]);
            return usage();
        } else {
            path = argv[i];
        }
    }
    if (!path) return usage();

    if (sesame_idat_read(path, &d, &e) != SESAME_OK) {
        fprintf(stderr, "sesamec: %s: %s\n", path, e.msg);
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
    if (strcmp(argv[1], "version") == 0) {
        puts("sesamec 0.0.1-dev");
        return 0;
    }
    return usage();
}
