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

/**
 * Returns the current directory.
 */
std::string getCurrentDirectory();

/**
 * Returns whether the given path is absolute.
 */
bool isAbsolutePath (const std::string &path);

/**
 * Converts the given path to an absolute path using the current directory.
 */
std::string toAbsolutePath(const char *path);

/**
 * Joins two paths.
 */
std::string joinPaths(const char *p1, const char *p2);

#ifdef DEBUG
#define MI_ASSERT(expr) \
    ((expr) \
     ? static_cast<void>(0) \
     : throw ::format("assert failed: %s (%s:%d)", # expr, __FILE__, __LINE__))
#else
#define MI_ASSERT(expr)  static_cast<void>(0)
#endif

#endif // ifndef mifit_util_utillib_h
