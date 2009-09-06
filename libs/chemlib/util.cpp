#include "util.h"
#include <ctype.h>
#include <algorithm>
#include <functional>

namespace chemlib {


int strncasecmp(const char* s1, const char* s2, unsigned int n) {
  for (; 0 < n; ++s1, ++s2, --n) {
    if (toupper(*s1) != toupper(*s2) ) {
      return ((*(unsigned char*)s1 <
               *(unsigned char*)s2) ? -1 : +1);
    } else if (*s1 == '\0') {
      return (0);
    }
  }
  return (0);
}

} //namespace chemlib
