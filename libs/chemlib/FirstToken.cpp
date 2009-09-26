#include <cstdio>
#include <string>
#include "FirstToken.h"

static inline bool IsSmilesWhitespace(char c) {
  return (c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f');
}

namespace chemlib {

/////////////////////////////////////////////////////////////////////////////
// Function:    FirstToken
// Purpose:   Get the first (whitespace-delimited) token from a file
// Input:       Pointer to the (opened) file
// Output:    std::string with the first token
// Requires:
/////////////////////////////////////////////////////////////////////////////
std::string MIFirstToken(FILE* fp) {
  char c;
  while ((c = fgetc(fp)) != EOF) {
    if (!IsSmilesWhitespace(c)) {
      break;
    }
  }

  if (c == EOF) {
    return "";
  }

  std::string token;
  do {
    token += c;
  } while ((c = fgetc(fp)) != EOF && !IsSmilesWhitespace(c));

  return token;
}

} // namespace chemlib
