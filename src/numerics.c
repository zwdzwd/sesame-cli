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

/* R's quantile(x, p, type=7) on data already sorted ascending and NaN-free.
 *
 * Transcribed from stats:::quantile.default, which for type 7 does:
 *
 *     index <- 1 + max(n - 1, 0) * probs
 *     lo <- floor(index);  hi <- ceiling(index)
 *     qs <- x[lo]
 *     i  <- which(index > lo & x[hi] != qs)
 *     h  <- (index - lo)[i]
 *     qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
 *
 * The interpolation MUST be written `(1-h)*a + h*b`, not the algebraically
 * equivalent `a + h*(b-a)`: the two disagree in the last bit for ~24% of random
 * inputs, which showed up as ~1% of betas differing by 1-2 ULP. Likewise the
 * `x[hi] != qs` guard is load-bearing -- when the two order statistics are equal
 * R returns qs untouched rather than recombining them.
 *
 * Returns NaN for n == 0 (R errors; callers guard). */
double sesame__quantile7_sorted(const double *x, int32_t n, double p)
{
    double index, h, qs;
    int32_t lo, hi;

    if (n <= 0) return NAN;
    if (n == 1) return x[0];

    index = 1.0 + (double)(n - 1) * p;   /* R is 1-based */
    lo    = (int32_t)floor(index);
    hi    = (int32_t)ceil(index);
    if (lo < 1) lo = 1;
    if (hi > n) hi = n;

    qs = x[lo - 1];                       /* to 0-based */
    if ((double)lo < index && x[hi - 1] != qs) {
        h  = index - (double)lo;
        qs = (1.0 - h) * qs + h * x[hi - 1];
    }
    return qs;
}

/* Same, but sorts x in place first. */
double sesame__quantile7(double *x, int32_t n, double p)
{
    if (n <= 0) return NAN;
    qsort(x, (size_t)n, sizeof(double), cmp_double);
    return sesame__quantile7_sorted(x, n, p);
}

void sesame__sort(double *x, int32_t n)
{
    if (n > 1) qsort(x, (size_t)n, sizeof(double), cmp_double);
}

/* R's median() on sorted NaN-free data. Identical to quantile type 7 at p=0.5
 * (even n -> mean of the two middle order statistics), so it is not a separate
 * implementation. */
double sesame__median_sorted(const double *x, int32_t n)
{
    return sesame__quantile7_sorted(x, n, 0.5);
}

/* preprocessCore::normalize.quantiles.use.target(matrix(x), target).
 *
 * CLEAN-ROOM. preprocessCore is LGPL-2 and sesame is MIT, so qnorm.c must not
 * be read or copied. This was characterized purely by black-box probing:
 *
 *   n==m, no ties  -> sorted(target)[rank(x)]
 *   n!=m           -> quantile(target, (rank-1)/(n-1), type=7)
 *   ties in x      -> quantile at the AVERAGE rank  (verified against the
 *                     competing "mean of the quantiles at each tied rank" rule
 *                     with a non-linear target: actual 100, rule-B 106.67)
 *   NA in target   -> dropped before ranking
 *   NA in x        -> preserved; ranking uses non-NA only
 *
 * This implementation requires x already sorted ascending and NaN-free, which
 * is how dyeBiasNL calls it (R/dye_bias.R:135). Because ties take the average
 * rank, every element of a tied run maps to the identical output -- confirmed
 * on real data (0 of 11259 tied runs on HM450 had non-constant output). That is
 * also why the subsequent approx(ties=mean) never actually averages anything.
 *
 * tgt is sorted in place. out must hold n doubles. Output is non-decreasing,
 * which is why R's subsequent sort() of the result is a no-op. */
void sesame__qnorm_use_target(const double *x_sorted, int32_t n,
                              double *tgt, int32_t m, double *out)
{
    int32_t a = 0;

    sesame__sort(tgt, m);
    while (a < n) {
        int32_t b = a, i;
        double r, p, q;
        while (b + 1 < n && x_sorted[b + 1] == x_sorted[a]) b++;
        /* average of the 1-based ranks spanned by the tied run [a, b] */
        r = ((double)(a + 1) + (double)(b + 1)) / 2.0;
        p = (n == 1) ? 0.0 : (r - 1.0) / (double)(n - 1);
        q = sesame__quantile7_sorted(tgt, m, p);
        for (i = a; i <= b; i++) out[i] = q;
        a = b + 1;
    }
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
