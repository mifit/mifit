#ifndef util_io_h
#define util_io_h

#include <cstdio>
#include <cstdarg>
#include <string>
#include <memory>

class io
{

    static io *defaultIo_;

protected:

    static void setDefaultIo(io *defaultIo);

    io();

    virtual io *create() = 0;

public:
    virtual ~io();

    static io *defaultIo();
    virtual bool open(const char *file, const char *mode) = 0;
    virtual bool isOpen() = 0;
    virtual void close() = 0;
    virtual char *gets(char *s, size_t length) = 0;

    /**
     * Reads a line from given file. Note: the current implementation may not
     * be compatable with Mac formatted files.
     * @param line the string in which the line is stored; the line separator
     *          is not stored
     * @return the number of characters read including the line separator
     */
    virtual size_t readLine(std::string &line) = 0;

    virtual int printf(const char *format, ...) = 0;
    virtual int vprintf(const char *format, va_list argp) = 0;
    virtual int seek(long offset, int whence) = 0;
    virtual long tell() = 0;
    virtual void rewind() = 0;

    // For migration compatibility
    virtual void attach(FILE *fp) = 0;
    virtual FILE *fp() = 0;

};

#endif // ifndef util_io_h
