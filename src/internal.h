/* internal.h -- shared across libsesamec TUs. Not installed, not public.
 * SPDX-License-Identifier: MIT
 */
#ifndef SESAME_INTERNAL_H
#define SESAME_INTERNAL_H

#include "sesame.h"

/* Fills err (if non-NULL) and returns code, so callers can `return
 * sesame__fail(...)`. Never longjmps, never writes to stderr. */
int sesame__fail(sesame_err_t *err, int code, const char *fmt, ...);

/* Index accessors -- keep the on-disk layout private to index.c. */
const uint32_t *sesame__index_M(const sesame_index_t *ix);
const uint32_t *sesame__index_U(const sesame_index_t *ix);
const uint8_t  *sesame__index_col(const sesame_index_t *ix);
const uint8_t  *sesame__index_mask(const sesame_index_t *ix);

#endif /* SESAME_INTERNAL_H */
