#ifndef util_StringIo_h
#define util_StringIo_h

#include <cstdio>
#include <cstdarg>
#include <string>
#include <map>

#include "io.h"

/**
 * Intended for testing input and output without touching files.
 */
class StringIo : public io {

  struct StringIoData {
    std::string data_;
    size_t index_;
  };

  typedef std::map<StringIoData*, size_t> DataRefCountMap;
  static DataRefCountMap refCounts;

  StringIoData* data_;

  virtual io* create();

public:

  StringIo();
  StringIo(const char* data);
  StringIo(const StringIo& stringIo);
  virtual ~StringIo();
  
  StringIo& operator=(const StringIo& stringIo);

  void setAsDefaultIo();

  virtual bool open(const char* file, const char* mode);
  virtual bool isOpen();
  virtual void close();
  virtual char* gets(char* s, size_t length);
  virtual size_t readLine(std::string& line);
  virtual int printf(const char* format, ...);
  virtual int vprintf(const char* format, va_list argp);
  virtual int seek(long offset, int whence);
  virtual long tell();
  virtual void rewind();
  
  virtual void attach(FILE* fp);
  virtual FILE* fp();

};

#endif
