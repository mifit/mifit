#include "LigDictEntry.h"

LigDictEntry::LigDictEntry(chemlib::RESIDUE* r)
: res(r) {
  
}

LigDictEntry::~LigDictEntry() {
  delete res;
}