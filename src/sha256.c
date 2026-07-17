/* sha256.c -- FIPS 180-4 SHA-256.
 *
 * Self-contained on purpose: CommonCrypto is macOS-only and OpenSSL would be a
 * whole dependency for one hash. ~150 lines is cheaper than either.
 *
 * Used to verify downloaded index files against the pinned digest in the
 * registry -- a downloaded index silently changing would silently change
 * everyone's betas.
 *
 * SPDX-License-Identifier: AGPL-3.0-or-later
 *
 * Copyright (C) 2026-present Wanding Zhou
 * Part of sesame-cli, licensed under AGPL-3.0-or-later; see LICENSE.
 */
#include "internal.h"

#include <stdio.h>
#include <string.h>

typedef struct {
    uint32_t h[8];
    uint64_t len;      /* total message bytes */
    uint8_t  buf[64];
    size_t   n;        /* bytes currently in buf */
} sha256_ctx;

static const uint32_t K[64] = {
    0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,
    0x923f82a4u,0xab1c5ed5u,0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,
    0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,0xe49b69c1u,0xefbe4786u,
    0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
    0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,
    0x06ca6351u,0x14292967u,0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,
    0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,0xa2bfe8a1u,0xa81a664bu,
    0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
    0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,
    0x5b9cca4fu,0x682e6ff3u,0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,
    0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
};

#define ROR(x,n) (((x) >> (n)) | ((x) << (32 - (n))))
#define S0(x) (ROR(x, 2) ^ ROR(x,13) ^ ROR(x,22))
#define S1(x) (ROR(x, 6) ^ ROR(x,11) ^ ROR(x,25))
#define s0(x) (ROR(x, 7) ^ ROR(x,18) ^ ((x) >>  3))
#define s1(x) (ROR(x,17) ^ ROR(x,19) ^ ((x) >> 10))

static void sha256_block(sha256_ctx *c, const uint8_t *p)
{
    uint32_t w[64], a, b, cc, d, e, f, g, h;
    int i;

    for (i = 0; i < 16; i++)
        w[i] = ((uint32_t)p[i*4] << 24) | ((uint32_t)p[i*4+1] << 16) |
               ((uint32_t)p[i*4+2] << 8) | (uint32_t)p[i*4+3];
    for (i = 16; i < 64; i++)
        w[i] = s1(w[i-2]) + w[i-7] + s0(w[i-15]) + w[i-16];

    a=c->h[0]; b=c->h[1]; cc=c->h[2]; d=c->h[3];
    e=c->h[4]; f=c->h[5]; g=c->h[6];  h=c->h[7];

    for (i = 0; i < 64; i++) {
        uint32_t t1 = h + S1(e) + ((e & f) ^ (~e & g)) + K[i] + w[i];
        uint32_t t2 = S0(a) + ((a & b) ^ (a & cc) ^ (b & cc));
        h=g; g=f; f=e; e=d+t1; d=cc; cc=b; b=a; a=t1+t2;
    }
    c->h[0]+=a; c->h[1]+=b; c->h[2]+=cc; c->h[3]+=d;
    c->h[4]+=e; c->h[5]+=f; c->h[6]+=g;  c->h[7]+=h;
}

static void sha256_init(sha256_ctx *c)
{
    c->h[0]=0x6a09e667u; c->h[1]=0xbb67ae85u; c->h[2]=0x3c6ef372u; c->h[3]=0xa54ff53au;
    c->h[4]=0x510e527fu; c->h[5]=0x9b05688cu; c->h[6]=0x1f83d9abu; c->h[7]=0x5be0cd19u;
    c->len = 0; c->n = 0;
}

static void sha256_update(sha256_ctx *c, const void *data, size_t len)
{
    const uint8_t *p = (const uint8_t *)data;
    c->len += len;
    if (c->n) {
        size_t take = 64 - c->n;
        if (take > len) take = len;
        memcpy(c->buf + c->n, p, take);
        c->n += take; p += take; len -= take;
        if (c->n == 64) { sha256_block(c, c->buf); c->n = 0; }
    }
    while (len >= 64) { sha256_block(c, p); p += 64; len -= 64; }
    if (len) { memcpy(c->buf, p, len); c->n = len; }
}

static void sha256_final(sha256_ctx *c, uint8_t out[32])
{
    uint64_t bits = c->len * 8;
    uint8_t pad = 0x80;
    int i;
    sha256_update(c, &pad, 1);
    pad = 0;
    while (c->n != 56) sha256_update(c, &pad, 1);
    for (i = 7; i >= 0; i--) {
        uint8_t b = (uint8_t)((bits >> (i * 8)) & 0xff);
        sha256_update(c, &b, 1);
    }
    for (i = 0; i < 8; i++) {
        out[i*4]   = (uint8_t)(c->h[i] >> 24);
        out[i*4+1] = (uint8_t)(c->h[i] >> 16);
        out[i*4+2] = (uint8_t)(c->h[i] >> 8);
        out[i*4+3] = (uint8_t)(c->h[i]);
    }
}

/* Hex digest of a file. out must hold 65 bytes. Returns 0 on success. */
int sesame__sha256_file(const char *path, char out[65])
{
    static const char hex[] = "0123456789abcdef";
    uint8_t digest[32], buf[65536];
    sha256_ctx c;
    size_t got;
    int i;
    FILE *f = fopen(path, "rb");

    if (!f) return -1;
    sha256_init(&c);
    while ((got = fread(buf, 1, sizeof buf, f)) > 0)
        sha256_update(&c, buf, got);
    if (ferror(f)) { fclose(f); return -1; }
    fclose(f);
    sha256_final(&c, digest);

    for (i = 0; i < 32; i++) {
        out[i*2]   = hex[digest[i] >> 4];
        out[i*2+1] = hex[digest[i] & 0xf];
    }
    out[64] = '\0';
    return 0;
}
