#include "LigDictEntry.h"

#include <chemlib/Monomer.h>

LigDictEntry::LigDictEntry(chemlib::Residue *r)
    : res(r)
{

}

LigDictEntry::~LigDictEntry()
{
    delete res;
}
