#include "StringIo.h"
#include <stdlib.h>
#include <string.h>

StringIo::DataRefCountMap StringIo::refCounts;

StringIo::StringIo() : data_(NULL) {
}

StringIo::StringIo(const char* data) : data_(NULL) {
  data_ = new StringIoData;
  data_->data_ = data;
  data_->index_ = 0;
  ++refCounts[data_];
}

StringIo::~StringIo() {
  close();
}

StringIo::StringIo(const StringIo& stringIo)
  : data_(stringIo.data_) {
  ++refCounts[data_];
}

StringIo& StringIo::operator=(const StringIo& stringIo) {
  data_ = stringIo.data_;
  ++refCounts[data_];
  return *this;
}

void StringIo::setAsDefaultIo() {
  io::setDefaultIo(new StringIo(*this));
}

io*  StringIo::create() {
  return new StringIo(*this);
}

bool StringIo::open(const char* /* file */, const char* mode) {
  bool opened = false;
  if (mode != NULL) {
    if (mode[0] == 'r' && data_ != NULL) {
      data_->index_ = 0;
      opened = true;
    }
  }
  return opened;
}

bool StringIo::isOpen() {
  return data_ != NULL;
}

void StringIo::close() {
  if (data_ != NULL) {
    --refCounts[data_];
    if (refCounts[data_] == 0) {
      delete data_;
    }
    data_ = NULL;
  }
}

int StringIo::printf(const char* format, ...) {
  va_list argp;
  va_start(argp, format);
  int result = vprintf(format, argp);
  va_end(argp);
  return result;
}

int StringIo::vprintf(const char* /* format */, va_list /* argp */) {
  return 0;
}

char* StringIo::gets(char* s, size_t length) {
  if (data_ == NULL || data_->index_ >= data_->data_.length()) {
    return NULL;
  }
  size_t start = data_->index_;
  size_t end = data_->index_;
  while (end < data_->data_.length() && (end-start) < length-1) {
    if (data_->data_[end] == '\n') {
      ++end;
      break;
    }
    ++end;
  }
  data_->index_ = end;
  size_t n = end - start;
  strncpy(s, data_->data_.c_str() + start, n);
  s[n] = '\0';
  return s;
}

size_t StringIo::readLine(std::string& line) {
  if (data_ == NULL || data_->index_ >= data_->data_.length()) {
    line = "";
    return 0;
  }
  size_t start = data_->index_;
  size_t end = data_->index_;
  while (end < data_->data_.length()) {
    if (data_->data_[end] == '\n') {
      ++end;
      break;
    }
    ++end;
  }
  size_t n = end - start;
  data_->index_ = end;
  if (end > 1 && data_->data_[end-1] == '\n') {
    --end;
    if (end > 1 && data_->data_[end-1] == '\r') {
      --end;
    }
  }
  line = data_->data_.substr(start, end - start);
  return n;
}

int StringIo::seek(long offset, int whence) {
  if (data_ == NULL) {
    return -1;
  }
  switch (whence) {
  default:
  case SEEK_SET:
    data_->index_ = offset;
    break;
  case SEEK_CUR:
    if (offset < 0 && ((size_t) -offset) > data_->index_) {
      data_->index_ = 0;
    } else {
      data_->index_ += offset;
    }
    break;
  case SEEK_END:
    if (offset < 0 && ((size_t) -offset) > data_->data_.length()) {
      data_->index_ = 0;
    } else {
      data_->index_ = data_->data_.length() + offset;
    }
    break;
  }
  return 0;
}

long StringIo::tell() {
  if (data_ == NULL) {
    return -1;
  }
  return (long) data_->index_;
}

void StringIo::rewind() {
  if (data_ != NULL) {
    data_->index_ = 0;
  }
}

void StringIo::attach(FILE* /* fp */) {
}

FILE* StringIo::fp() {
  return NULL;
}

