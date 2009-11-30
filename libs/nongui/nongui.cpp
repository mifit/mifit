#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string>
#include "nongui.h"
#include "util/utillib.h"

// things in this file provide non-gui implementations of Logger, emess,
// WaitCursor, etc

static std::string last_message("");

static void emess(const char *m)
{
    last_message = std::string(m);
    printf("%s\n", m);
}

void Logger::footer(const char *format, ...)
{
    char buf[1024];
    va_list argp;
    va_start(argp, format);
    vsprintf(buf, format, argp);
    if (std::string(buf) != last_message)
    {
        emess(buf);
    }
}

void Logger::log(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    emess(s.c_str());
    va_end(argp);
}

void Logger::message(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    emess(s.c_str());
    va_end(argp);
}

void Logger::debug(const char *format, ...)
{
#ifdef DEBUG
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    emess(s.c_str());
    va_end(argp);
#endif
}

WaitCursor::WaitCursor(const char *op)
{
    op_name = std::string(op);
    printf("Beginning Operation: %s\n", op);
}

WaitCursor::~WaitCursor()
{
    printf("\n");
}

bool WaitCursor::CheckForAbort()
{
    printf(".");
    return false;
}

