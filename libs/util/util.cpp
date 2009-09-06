#include "utillib.h"

#include <stdio.h>
#include <stdarg.h>
#include <string>
#include <string.h>

// Define used in format_arg_list
#ifdef _WIN32
#include <windows.h>
#undef PATH_MAX
const int PATH_MAX = MAX_PATH;
#define vsnprintf _vsnprintf
// vsnprintf returns negative value indicating truncation
#define VSNPRINTF_NEG_TRUNC_CODE
#else
#include <unistd.h>
#endif

// Variable argument list version of std::string format()
std::string format_arg_list(const char* fmt, va_list args) {
  if (!fmt) {
    return "";
  }
  int result = -1;
  int length = 256;
  char* buffer = 0;
  while (result < 0) {
    if (buffer) {
      delete[] buffer;
    }
    buffer = new char[length + 1];
    memset(buffer, 0, length + 1);
    result = vsnprintf(buffer, length, fmt, args);
#ifdef VSNPRINTF_NEG_TRUNC_CODE
    // vsnprintf returns negative value indicating truncation
    length *= 2;
#else
    // vsnprintf returns negative value indicating error
    // otherwise returns length written or length that would
    // have been written if truncation had not been required.
    if (result < 0) {
      buffer[0] = '\0';
      break;
    } else if (result >= length) {
      length = result + 1;
      result = -1;
    }
#endif
  }
  std::string s(buffer);
  delete[] buffer;
  return s;
}

std::string format(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  std::string s = format_arg_list(fmt, args);
  va_end(args);
  return s;
}

bool startsWith(const std::string& str, const std::string& substr) {
  return str.find(substr) == 0;
}

bool startsWith(const std::string& str, const char* substr) {
  return startsWith(str, std::string(substr));
}

bool endsWith(const std::string& str, const std::string& substr) {
  return str.rfind(substr) == (str.length() - substr.length());
}

bool endsWith(const std::string& str, const char* substr) {
  return endsWith(str, std::string(substr));
}

const char* file_extension(const char* file) {
  for (size_t i = strlen(file); i > 0; --i) {
    if (file[i-1] == '.') {
      return file+(i-1);
    }
  }

  // no extension, return null
  static char empty[2] = {0,0};
  return empty;
}

std::string getCurrentDirectory() {
  char dir[FILENAME_MAX + 1];
#ifdef _WIN32
  size_t len = GetCurrentDirectory(sizeof(dir), dir);
  if (len == 0) {
    throw std::string("Unable to get current directory");
  }
#else
  if (!getcwd(dir, sizeof(dir))) {
    throw std::string("Unable to get current directory");
  }
#endif
  return std::string(dir);
}

bool isAbsolutePath(const std::string& path) {
  if (path.size() > 0) {
    if (path[0] == '/') {
      return true;
    }
#ifdef _WIN32
    if ((path[0] == '\\') || (path[1] == ':')) {
      return true;
    }
#endif
  }
  return false;
}

std::string toAbsolutePath(const char* path) {
  std::string absPath = (path != NULL) ? path : "";
  if (!isAbsolutePath(path)) {
    absPath = joinPaths(getCurrentDirectory().c_str(), path);
  }
  return absPath;
}

std::string joinPaths(const char* p1, const char* p2) {
#ifdef _WIN32
  static const char* sep = "\\";
#else
  static const char* sep = "/";
#endif
  std::string path = (p1 != NULL) ? p1 : "";
  if (p2 != NULL) {
    path += sep;
    path += p2;
  }
  return path;
}
