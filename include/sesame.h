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

/* --------------------------------------------------------------- index ---
 *
 * The per-platform ordering table: Probe_ID, the M/U bead addresses, the
 * Infinium-I colour channel, and the default design mask.
 *
 * Opaque by design -- the on-disk format is an implementation detail. Today it
 * parses a TSV; a mmap-able binary can be swapped in without touching callers.
 *
 * Row order is load-bearing: it defines the canonical SigDF row order and
 * mirrors R/sesame.R:504. Never sort it.
 */
typedef struct sesame_index_t sesame_index_t;

/* Infinium design type / in-band channel. Mirrors sesame's col factor
 * levels c("G","R","2"). */
#define SESAME_COL_II 0  /* Infinium-II  (no M address; UG=meth, UR=unmeth) */
#define SESAME_COL_G  1  /* Infinium-I, green in-band */
#define SESAME_COL_R  2  /* Infinium-I, red in-band   */

sesame_index_t *sesame_index_open(const char *path, sesame_err_t *err);
void            sesame_index_close(sesame_index_t *ix);
int32_t         sesame_index_nprobes(const sesame_index_t *ix);
const char     *sesame_index_probe_id(const sesame_index_t *ix, int32_t i);

/* --------------------------------------------------------------- sigdf ---
 *
 * Signals, one row per index probe, in index order. NaN encodes NA.
 */

/* Status bits. R signals these conditions silently or not at all; we surface
 * them so a pipeline can count fallbacks instead of finding them in the betas. */
#define SESAME_STAT_ADDR_MISSING  (1u << 0) /* manifest address absent from IDAT */
#define SESAME_STAT_PAIR_MISMATCH (1u << 1) /* Grn/Red address vectors disagree  */

typedef struct {
    const sesame_index_t *ix;
    int32_t   n;
    double   *MG, *MR, *UG, *UR;
    uint8_t  *col;    /* SESAME_COL_* */
    uint8_t  *mask;
    uint32_t  status; /* SESAME_STAT_* */
    int32_t   n_addr_missing;
} sesame_sigdf_t;

/* Assembles a SigDF from a Grn/Red IDAT pair against the index.
 *
 * Divergence D1: R merges the two IDATs by *position* with no join on
 * IlluminaID (R/sesame.R:293-298). That assumption holds for well-formed files
 * -- verified across every test array -- but fails silently and array-wide when
 * mismatched files are paired. We verify it and error instead.
 *
 * min_beads: 0 disables bead filtering (R's default is NULL/off). Otherwise a
 * probe is masked when any relevant bead count is missing or < min_beads. */
sesame_sigdf_t *sesame_sigdf_from_idats(const sesame_idat_t *grn,
                                        const sesame_idat_t *red,
                                        const sesame_index_t *ix,
                                        int min_beads,
                                        sesame_err_t *err);
void sesame_sigdf_free(sesame_sigdf_t *sdf);

/* Beta values, one per probe, in index order. NaN for NA and for masked probes
 * when apply_mask is non-zero. out must hold sesame_index_nprobes() doubles. */
int sesame_get_betas(const sesame_sigdf_t *sdf, int apply_mask,
                     double *out, sesame_err_t *err);

#ifdef __cplusplus
}
#endif

#endif /* SESAME_H */
