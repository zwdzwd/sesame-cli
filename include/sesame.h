/* sesame.h -- public API for libsesame
 *
 * A standalone C implementation of sesame's basic Infinium preprocessing.
 * See NUMERICS.md for documented divergences from the R implementation.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
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

/* Platform auto-detection from the IDAT bead count (nSNPsRead). Returns NULL
 * when the count is not one we have verified -- callers must then be told the
 * platform explicitly rather than have it guessed. Avoids depending on
 * sesameData's idatSignature. */
const char *sesame_platform_from_beads(int32_t beads);

/* The annotation tag this build pins, and whose digests it can verify. */
const char *sesame_default_tag(void);

/* The index store -- where fetch writes, and where the active-version links
 * live. One variable, since the store is managed (see sesame_index_activate),
 * not a cache:
 *   $SESAME_INDEX_DIR | <dir of the binary>/data | $XDG_CACHE_HOME/sesame
 *                     | ~/Library/Caches/sesame (macOS) | ~/.cache/sesame
 * The binary-relative default means a checkout is found from any cwd, while an
 * installed binary (no data/ beside it) falls through to the XDG store.
 * Returns out. */
const char *sesame_store_dir(char *out, size_t n);

/* Finds an existing index for platform. 0 and fills out on success, -1 if
 * absent. Never downloads, never prompts. */
int sesame_index_locate(const char *platform, char *out, size_t n);

/* Point <store>/current at <tag> (relative link, atomic rename), selecting that
 * annotation snapshot for every platform at once. */
int sesame_index_activate(const char *store, const char *tag, sesame_err_t *err);

/* The tag `current` points at. 0 on success. */

/* Fills msg with actionable "no index found, here is how to fix it" text. */
void sesame_index_missing_help(const char *platform, char *msg, size_t n);

/* Fetch one platform at the pinned tag into the store: pull its SHA256SUMS
 * (verified against a digest compiled into this build), then every file it
 * lists (ordering table + the .cm mask). Files already present with the right
 * digest are skipped. out_path receives the ordering table's path.
 *
 * This and sesame_fetch_all are the ONLY paths that touch the network -- nothing
 * downloads implicitly and nothing ever prompts, so behaviour is identical with
 * or without a TTY. */
int sesame_fetch_index(const char *platform, int force,
                       char *out_path, size_t out_n, sesame_err_t *err);

/* Fetch every platform published at the pinned tag. */
int sesame_fetch_all(int force, sesame_err_t *err);

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
#define SESAME_STAT_ADDR_MISSING   (1u << 0) /* manifest address absent from IDAT */
#define SESAME_STAT_PAIR_MISMATCH  (1u << 1) /* Grn/Red address vectors disagree  */
#define SESAME_STAT_DYEBIAS_FAILED (1u << 2) /* D gave up: green channel failed   */
#define SESAME_STAT_NOOB_SKIPPED   (1u << 3) /* B skipped: <100 background probes  */
#define SESAME_STAT_NOOB_MAD0      (1u << 4) /* B: huber MAD==0 fallback (D6)       */

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

/* --------------------------------------------------------------- prep ---
 *
 * prepSesame steps. Order matters (R/open.R:30-35). Only C is implemented.
 */

/* Q -- load the recommended quality mask for a platform as a 0/1 vector aligned
 * to the ordering (1 = masked). Shells out to the `yame` binary to read the .cm
 * mask in the store; *out is malloc'd (caller frees), *out_n is the probe count.
 * yame is invoked as a separate process, so its AGPL code never links here. */
int sesame_quality_mask(const char *platform, uint8_t **out, int32_t *out_n,
                        sesame_err_t *err);

/* Q -- apply a mask vector from sesame_quality_mask to a SigDF (mask |= qmask).
 * qn must equal the probe count. */
int sesame_prep_quality_mask(sesame_sigdf_t *sdf, const uint8_t *qmask,
                             int32_t qn, sesame_err_t *err);

/* P -- the background mask (backgroundMask names) as a 0/1 vector aligned to the
 * ordering (1 = excluded from background estimation). Reads the same .cm. */
int sesame_background_mask(const char *platform, uint8_t **out, int32_t *out_n,
                           sesame_err_t *err);

/* P -- pOOBAH (R/detection.R:141-170). Masks probes whose detection p-value
 * exceeds pval_threshold, from the ecdf of out-of-band + negative-control
 * background (probes in bgmask excluded). combine_neg adds negative controls.
 * Implements D2 (NA-channel probes are masked, not silently kept). */
int sesame_prep_poobah(sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                       double pval_threshold, int combine_neg, sesame_err_t *err);

/* C -- inferInfiniumIChannel (R/channel_inference.R:20-55). Reassigns each
 * Infinium-I probe to its brighter channel, judged against the 95th percentile
 * of the pooled out-of-band signal under the new assignment. Defaults
 * (switch_failed=0, mask_failed=0) revert failed probes to the manifest channel
 * and do not mask them. */
int sesame_prep_infer_channel(sesame_sigdf_t *sdf, int switch_failed,
                              int mask_failed, sesame_err_t *err);

/* D -- dyeBiasNL (R/dye_bias.R:118-167). Pulls the Red and Green Infinium-I
 * distributions to their common midpoint. If the green channel has failed
 * (RGdistort NA or > 10) it masks all Inf-I green probes, sets
 * SESAME_STAT_DYEBIAS_FAILED and returns -- R does this silently. */
int sesame_prep_dye_bias_nl(sesame_sigdf_t *sdf, sesame_err_t *err);

/* B -- noob (R/background.R:85-123). Normal-exponential background subtraction:
 * background (out-of-band + negative controls, probes in bgmask excluded) is
 * modelled Normal, true signal Exponential, and each channel's signal is
 * deconvolved via the inverse Mills ratio (D5), then shifted by offset (R's
 * default 15). combine_neg adds negative controls to the background, as in P.
 * If either channel has <100 positive background values, R returns the SigDF
 * unchanged; we do the same and set SESAME_STAT_NOOB_SKIPPED. D6: a huber MAD of
 * 0 sets SESAME_STAT_NOOB_MAD0 instead of aborting. */
int sesame_prep_noob(sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                     int combine_neg, double offset, sesame_err_t *err);

/* Beta values, one per probe, in index order. NaN for NA and for masked probes
 * when apply_mask is non-zero. out must hold sesame_index_nprobes() doubles. */
int sesame_get_betas(const sesame_sigdf_t *sdf, int apply_mask,
                     double *out, sesame_err_t *err);

/* ----------------------------------------------------------------- QC ---
 *
 * The sesameQC panel (R sesameQC_calcStats), computed from a *raw* SigDF. Each
 * metric is one field, tagged I (integer count) or D (real) for formatting. The
 * groups mirror R/QC.R: detection, numProbes, intensity, channel, dyeBias, betas.
 *
 * Two faithful-but-noted differences (see NUMERICS.md):
 *  - the detection group runs pOOBAH internally, so it carries the same mask
 *    lineage as `P`, and `num_dtna` counts probes with no signal in *either*
 *    channel (the D2 fix makes R's pOOBAH no longer return NA there);
 *  - the betas group is computed on a D->B->P-processed copy, exactly as
 *    sesameQC_calcStats_betas does (getBetas(pOOBAH(noob(dyeBiasNL(sdf))))).
 */
#define SESAME_QC_FIELDS(_) \
    _(I, num_dtna)      _(D, frac_dtna) \
    _(I, num_dt)        _(D, frac_dt) \
    _(I, num_dt_mk)     _(D, frac_dt_mk) \
    _(I, num_dt_cg)     _(D, frac_dt_cg) \
    _(I, num_dt_ch)     _(D, frac_dt_ch) \
    _(I, num_dt_rs)     _(D, frac_dt_rs) \
    _(I, num_probes)    _(I, num_probes_II) _(I, num_probes_IR) _(I, num_probes_IG) \
    _(I, num_probes_cg) _(I, num_probes_ch) _(I, num_probes_rs) \
    _(D, mean_intensity)   _(D, mean_intensity_MU) _(D, mean_ii) \
    _(D, mean_inb_grn)     _(D, mean_inb_red) \
    _(D, mean_oob_grn)     _(D, mean_oob_red) \
    _(I, na_intensity_M)   _(I, na_intensity_U) \
    _(I, na_intensity_ig)  _(I, na_intensity_ir)  _(I, na_intensity_ii) \
    _(I, InfI_switch_R2R)  _(I, InfI_switch_G2G) \
    _(I, InfI_switch_R2G)  _(I, InfI_switch_G2R) \
    _(D, medR) _(D, medG) _(D, topR) _(D, topG) _(D, RGratio) _(D, RGdistort) \
    _(D, mean_beta)    _(D, median_beta)    _(D, frac_unmeth)    _(D, frac_meth)    _(I, num_na)    _(D, frac_na) \
    _(D, mean_beta_cg) _(D, median_beta_cg) _(D, frac_unmeth_cg) _(D, frac_meth_cg) _(I, num_na_cg) _(D, frac_na_cg) \
    _(D, mean_beta_ch) _(D, median_beta_ch) _(D, frac_unmeth_ch) _(D, frac_meth_ch) _(I, num_na_ch) _(D, frac_na_ch) \
    _(D, mean_beta_rs) _(D, median_beta_rs) _(D, frac_unmeth_rs) _(D, frac_meth_rs) _(I, num_na_rs) _(D, frac_na_rs)

typedef struct {
#define SESAME_QC_DECL(t, nm) double nm;
    SESAME_QC_FIELDS(SESAME_QC_DECL)
#undef SESAME_QC_DECL
} sesame_qc_t;

/* Compute the panel from a raw SigDF. bgmask is the background mask (as for P/B);
 * required, since the detection and beta groups run pOOBAH internally. Does not
 * modify sdf. */
int sesame_qc_calc(const sesame_sigdf_t *sdf, const uint8_t *bgmask, int32_t bn,
                   sesame_qc_t *out, sesame_err_t *err);

/* Tab-separated metric-name header, no leading sample column and no trailing
 * tab. */
const char *sesame_qc_header(void);

/* Format the metric values tab-separated into buf (integer metrics as integers,
 * reals as %.10g, NaN as "NA"). Returns the length written, or -1 if truncated. */
int sesame_qc_format_row(const sesame_qc_t *q, char *buf, size_t n);

#ifdef __cplusplus
}
#endif

#endif /* SESAME_H */
