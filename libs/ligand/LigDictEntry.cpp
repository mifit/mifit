#include "LigDictEntry.h"

#include <chemlib/RESIDUE_.h>

LigDictEntry::LigDictEntry(chemlib::RESIDUE* r)
: res(r) {
  
}

LigDictEntry::~LigDictEntry() {
  delete res;
}
