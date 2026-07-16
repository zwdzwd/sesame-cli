/* internal.h -- shared across libsesame TUs. Not installed, not public.
 * SPDX-License-Identifier: MIT
 */
#ifndef SESAME_INTERNAL_H
#define SESAME_INTERNAL_H

#include "sesame.h"

/* Fills err (if non-NULL) and returns code, so callers can `return
 * sesame__fail(...)`. Never longjmps, never writes to stderr. */
int sesame__fail(sesame_err_t *err, int code, const char *fmt, ...);

/* SHA-256 hex digest of a file; out must hold 65 bytes. 0 on success. */
int sesame__sha256_file(const char *path, char out[65]);

/* --- numerics: each replicates a specific R semantic exactly --- */

/* Copy non-NaN values of x into out (>= n); returns the count. Mirrors R's
 * na.rm=TRUE / sort(na.last=NA): NAs drop *before* n is taken. */
int32_t sesame__drop_na(const double *x, int32_t n, double *out);

/* R's quantile(x, p, type=7). _sorted takes ascending NaN-free data; the other
 * sorts x in place first. */
double sesame__quantile7_sorted(const double *x, int32_t n, double p);
double sesame__quantile7(double *x, int32_t n, double p);

/* R's median() -- identical to quantile type 7 at p=0.5. */
double sesame__median_sorted(const double *x, int32_t n);

void sesame__sort(double *x, int32_t n);

/* preprocessCore::normalize.quantiles.use.target, CLEAN-ROOM (preprocessCore is
 * LGPL-2, sesame is MIT: qnorm.c must not be read or copied). Characterized by
 * black-box probing; see the comment in numerics.c. x_sorted must be ascending
 * and NaN-free; tgt is sorted in place; out holds n doubles. */
void sesame__qnorm_use_target(const double *x_sorted, int32_t n,
                              double *tgt, int32_t m, double *out);

/* R's pmax(a, b) with NO na.rm -- NA propagates. Note C's fmax() does the
 * opposite and would silently disagree. */
double sesame__pmax2(double a, double b);

/* Registry lookup (registry.h must be included by the TU that uses the type). */
struct sesame_reg_t;

/* Index accessors -- keep the on-disk layout private to index.c. */
const uint32_t *sesame__index_M(const sesame_index_t *ix);
const uint32_t *sesame__index_U(const sesame_index_t *ix);
const uint8_t  *sesame__index_col(const sesame_index_t *ix);
const uint8_t  *sesame__index_mask(const sesame_index_t *ix);

#endif /* SESAME_INTERNAL_H */
