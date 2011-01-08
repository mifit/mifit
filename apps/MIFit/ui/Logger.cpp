#include "Logger.h"

#include <cstdarg>
#include <string>

#include <util/utillib.h>

#include "MIMainWindow.h"
#include <QMessageBox>

void Logger::log(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    va_end(argp);
    MIMainWindowLog(s);
}

void Logger::message(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    va_end(argp);
    QMessageBox::information(0, "Info", s.c_str());
}

void Logger::footer(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    va_end(argp);
    MIMainWindowLeftFooter(s);
}

void Logger::debug(const char *format, ...)
{
#ifdef DEBUG
    va_list argp;
    va_start(argp, format);
    std::string s = format_arg_list(format, argp);
    va_end(argp);
    MIMainWindowDebug(s);
#endif
}

