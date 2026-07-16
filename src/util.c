/* util.c -- shared helpers.
 * SPDX-License-Identifier: MIT
 */
#include "internal.h"

#include <stdarg.h>
#include <stdio.h>

int sesame__fail(sesame_err_t *err, int code, const char *fmt, ...)
{
    if (err) {
        va_list ap;
        err->code = code;
        va_start(ap, fmt);
        vsnprintf(err->msg, sizeof(err->msg), fmt, ap);
        va_end(ap);
    }
    return code;
}

const char *sesame_strerror(int code)
{
    switch (code) {
    case SESAME_OK:              return "ok";
    case SESAME_ERR_IO:          return "I/O error";
    case SESAME_ERR_FORMAT:      return "malformed input";
    case SESAME_ERR_UNSUPPORTED: return "unsupported input";
    case SESAME_ERR_NOMEM:       return "out of memory";
    default:                     return "unknown error";
    }
}
