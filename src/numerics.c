/* numerics.c -- the handful of R numerical primitives the prep steps rely on.
 *
 * Each one replicates a specific R semantic exactly; the comments say which,
 * because "close enough" here shows up as a different beta.
 *
 * SPDX-License-Identifier: MIT
 */
#include "internal.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    if (x < y) return -1;
    if (x > y) return 1;
    return 0;
}

/* Copy the non-NaN values of x into out (caller-allocated, >= n). Returns the
 * count. Mirrors R's na.rm=TRUE / sort()'s na.last=NA: NAs are dropped *before*
 * n is taken, so they never affect the quantile position. */
int32_t sesame__drop_na(const double *x, int32_t n, double *out)
{
    int32_t m = 0, i;
    for (i = 0; i < n; i++)
        if (!isnan(x[i])) out[m++] = x[i];
    return m;
}

/* R's quantile(x, p, type=7) on already-NaN-filtered data. Sorts in place.
 *
 *   h  = (n-1) * p
 *   lo = floor(h)
 *   Q  = x[lo] + (h - lo) * (x[lo+1] - x[lo])
 *
 * Returns NaN for n == 0 (R errors; callers guard). For n == 1 or p == 1 the
 * interpolation degenerates to the single/last order statistic -- guarded so we
 * never index x[n]. */
double sesame__quantile7(double *x, int32_t n, double p)
{
    double h, frac;
    int32_t lo;

    if (n <= 0) return NAN;
    if (n == 1) return x[0];

    qsort(x, (size_t)n, sizeof(double), cmp_double);

    h  = (double)(n - 1) * p;
    lo = (int32_t)floor(h);
    if (lo >= n - 1) return x[n - 1];
    frac = h - (double)lo;
    if (frac == 0.0) return x[lo];
    return x[lo] + frac * (x[lo + 1] - x[lo]);
}

/* R's pmax(a, b) with NO na.rm: NA in either operand propagates.
 *
 * C's fmax() does the opposite -- fmax(NaN, x) returns x -- which would
 * silently disagree with R at exactly the probes that matter (missing
 * addresses). Hence the explicit test. */
double sesame__pmax2(double a, double b)
{
    if (isnan(a) || isnan(b)) return NAN;
    return a > b ? a : b;
}
