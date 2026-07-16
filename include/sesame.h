/* sesame.h -- public API for libsesamec
 *
 * A standalone C implementation of sesame's basic Infinium preprocessing.
 * See NUMERICS.md for documented divergences from the R implementation.
 *
 * SPDX-License-Identifier: MIT
 */
#ifndef SESAME_H
#define SESAME_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SESAME_OK              0
#define SESAME_ERR_IO          1
#define SESAME_ERR_FORMAT      2
#define SESAME_ERR_UNSUPPORTED 3
#define SESAME_ERR_NOMEM       4

/* Errors are returned by status code and described in a caller-owned struct.
 * The library never longjmps and never writes to stderr. */
typedef struct {
    int  code;
    char msg[256];
} sesame_err_t;

const char *sesame_strerror(int code);

/* ---------------------------------------------------------------- IDAT ---
 *
 * Illumina IDAT v3 (non-encrypted). Only the four numeric sections are read
 * (IlluminaID/Mean/SD/NBeads); string and RunInfo sections are deliberately
 * skipped because they are unreliable in gzipped and de-identified files.
 * This mirrors R/readIDAT.R:244-246.
 */
typedef struct {
    int32_t   n;        /* nSNPsRead */
    uint32_t *addr;     /* field 102, IlluminaID (bead address) */
    uint16_t *mean;     /* field 104 */
    uint16_t *sd;       /* field 103 */
    uint8_t  *nbeads;   /* field 107 */
    int64_t   version;
    int32_t   n_fields;
} sesame_idat_t;

/* Reads a .idat or .idat.gz. On success *out is a heap object owned by the
 * caller; free with sesame_idat_free. On failure returns non-zero, *out is
 * untouched, and err (if non-NULL) is filled in. */
int  sesame_idat_read(const char *path, sesame_idat_t **out, sesame_err_t *err);
void sesame_idat_free(sesame_idat_t *idat);

#ifdef __cplusplus
}
#endif

#endif /* SESAME_H */
