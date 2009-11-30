#include "FileIo.h"
#include <string.h>

FileIo::FilePointerRefCountMap FileIo::refCounts;

void FileIo::setAsDefaultIo()
{
    io::setDefaultIo(new FileIo());
}

FileIo::FileIo()
    : filePointer_(NULL),
      bufferSize_(2048)
{
}

FileIo::~FileIo()
{
    close();
}

FileIo::FileIo(const FileIo &fileIo)
    : filePointer_(fileIo.filePointer_),
      bufferSize_(fileIo.bufferSize_)
{
    ++refCounts[filePointer_];
}

FileIo&FileIo::operator=(const FileIo &fileIo)
{
    filePointer_ = fileIo.filePointer_;
    bufferSize_ = fileIo.bufferSize_;
    ++refCounts[filePointer_];
    return *this;
}

io*FileIo::create()
{
    return new FileIo;
}

bool FileIo::open(const char *file, const char *mode)
{
    close();
    filePointer_ = fopen(file, mode);
    ++refCounts[filePointer_];
    return isOpen();
}

bool FileIo::isOpen()
{
    return filePointer_ != NULL;
}

void FileIo::close()
{
    if (filePointer_ != NULL)
    {
        if (refCounts.find(filePointer_) != refCounts.end())
        {
            --refCounts[filePointer_];
            if (refCounts[filePointer_] == 0)
            {
                refCounts.erase(filePointer_);
                fclose(filePointer_);
            }
        }
        filePointer_ = NULL;
    }
}

int FileIo::printf(const char *format, ...)
{
    va_list argp;
    va_start(argp, format);
    int result = vprintf(format, argp);
    va_end(argp);
    return result;
}

int FileIo::vprintf(const char *format, va_list argp)
{
    return vfprintf(filePointer_, format, argp);
}

char*FileIo::gets(char *s, size_t length)
{
    return fgets(s, (int)length, filePointer_);
}

size_t FileIo::readLine(std::string &line)
{
    char *buffer = new char[bufferSize_];
    size_t charCount = 0;
    line = "";
    while (1)
    {
        if (fgets(buffer, (int)bufferSize_, filePointer_) == NULL)
        {
            break;
        }
        size_t length = strlen(buffer);
        charCount += length;
        if (length == 0)
        {
            continue;
        }
        else if (buffer[length-1] != '\n')
        {
            line += buffer;
        }
        else
        {
            buffer[length-1] = '\0';
            if (length > 1 && buffer[length-2] == '\r')
            {
                buffer[length-2] = '\0';
            }
            line += buffer;
            break;
        }
    }
    delete[] buffer;
    return charCount;
}

int FileIo::seek(long offset, int whence)
{
    return fseek(filePointer_, offset, whence);
}

long FileIo::tell()
{
    return ftell(filePointer_);
}

void FileIo::rewind()
{
    ::rewind(filePointer_);
}

void FileIo::attach(FILE *fp)
{
    close();
    filePointer_ = fp;
    ++refCounts[filePointer_];
    // Increment reference count so it is never automatically
    // closed and deleted by FileIo
    ++refCounts[fp];

}

FILE*FileIo::fp()
{
    return filePointer_;
}

