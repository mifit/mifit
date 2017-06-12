#ifndef mifit_util_utillib_h
#define mifit_util_utillib_h

#include <stdarg.h>
#include <stdio.h>
#include <string>

#include "io.h"
#include "system.h"

// Printf-style formatting of a std::string
std::string format(const char *fmt, ...);
std::string format_arg_list(const char *fmt, va_list args);

bool startsWith(const std::string &str, const std::string &substr);
bool startsWith(const std::string &str, const char *substr);

bool endsWith(const std::string &str, const std::string &substr);
bool endsWith(const std::string &str, const char *substr);

/**
 * Find the extension of a file (the characters after the last . in the filename).
 */
const char *file_extension(const char *file);

int mi_strncasecmp(const char *s1, const char *s2, unsigned int n);

#ifdef DEBUG
#define MI_ASSERT(expr) \
    ((expr) \
     ? static_cast<void>(0) \
     : throw ::format("assert failed: %s (%s:%d)", # expr, __FILE__, __LINE__))
#else
#define MI_ASSERT(expr)  static_cast<void>(0)
#endif

#endif // ifndef mifit_util_utillib_h
