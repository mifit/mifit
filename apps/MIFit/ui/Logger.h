#ifndef mifit_ui_Logger_h
#define mifit_ui_Logger_h

#include <string>

class Logger
{
public:
    static void log(const char *format, ...);
    static void debug(const char *format, ...);

    static void message(const char *format, ...);
    static void footer(const char *format, ...);


    static void log(const std::string &str)
    {
        log(str.c_str());
    }

    static void debug(const std::string &str)
    {
        debug(str.c_str());
    }

    static void message(const std::string &str)
    {
        message(str.c_str());
    }

    static void footer(const std::string &str)
    {
        footer(str.c_str());
    }

};

#endif // mifit_ui_Logger_h
