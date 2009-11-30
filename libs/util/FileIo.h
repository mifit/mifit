#ifndef util_FileIo_h
#define util_FileIo_h

#include <cstdio>
#include <cstdarg>
#include <string>
#include <map>

#include "io.h"

class FileIo : public io
{

    typedef std::map<FILE*, size_t> FilePointerRefCountMap;
    static FilePointerRefCountMap refCounts;

    FILE *filePointer_;
    size_t bufferSize_;

    virtual io *create();

public:

    static void setAsDefaultIo();

    FileIo();
    FileIo(const FileIo &fileIo);
    virtual ~FileIo();

    FileIo&operator=(const FileIo &fileIo);

    virtual bool open(const char *file, const char *mode);
    virtual bool isOpen();
    virtual void close();
    virtual char *gets(char *s, size_t length);
    virtual size_t readLine(std::string &line);
    virtual int printf(const char *format, ...);
    virtual int vprintf(const char *format, va_list argp);
    virtual int seek(long offset, int whence);
    virtual long tell();
    virtual void rewind();

    virtual void attach(FILE *fp);
    virtual FILE *fp();

};

#endif // ifndef util_FileIo_h
