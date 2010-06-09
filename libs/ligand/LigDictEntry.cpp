#include "LigDictEntry.h"

#include <chemlib/RESIDUE_.h>

LigDictEntry::LigDictEntry(chemlib::Residue *r)
    : res(r)
{

}

LigDictEntry::~LigDictEntry()
{
    delete res;
}
