/* sha256.c -- SHA-256, now a shim over YAME's.
 *
 * This file used to carry its own FIPS 180-4 implementation. So did kycg's
 * src/digest.c and methscope-cli's src/digest.c: three copies of one hash
 * function, in three tools that all already link libyame.a. The implementation
 * moved to YAME/src/digest.c; what stays here is the symbol the rest of sesame
 * calls, so no caller had to change.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "internal.h"
#include "assets.h"

int sesame__sha256_file(const char *path, char out[65])
{
    return yame_assets_sha256_file(path, out);
}
