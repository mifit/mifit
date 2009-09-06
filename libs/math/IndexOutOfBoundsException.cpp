#include "IndexOutOfBoundsException.h"
#include <cstring>

namespace mi {
namespace math {

IndexOutOfBoundsException::IndexOutOfBoundsException()
  : detail(0) {
}

IndexOutOfBoundsException::IndexOutOfBoundsException(const char* detail) {

  this->detail = new char[strlen(detail)+1];
  strcpy(this->detail, detail);
}

IndexOutOfBoundsException::IndexOutOfBoundsException(
  const IndexOutOfBoundsException& ex) {

  detail = 0;
  if (ex.detail != 0) {
    detail = new char[strlen(ex.detail)+1];
    strcpy(detail, ex.detail);
  }
}

IndexOutOfBoundsException& IndexOutOfBoundsException::operator=(
  const IndexOutOfBoundsException& ex) {

  detail = 0;
  if (ex.detail != 0) {
    detail = new char[strlen(ex.detail)+1];
    strcpy(detail, ex.detail);
  }
  return *this;
}

IndexOutOfBoundsException::~IndexOutOfBoundsException() {
  if (detail != 0) {
    delete detail;
  }
}

const char* IndexOutOfBoundsException::what() const {
  return detail;
}

}
}
