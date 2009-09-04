#include "ConfSaver.h"
#include "MIAtom.h"
#include <chemlib/RESIDUE_.h>

using namespace std;

namespace chemlib {

ConfSaver::ConfSaver(RESIDUE* res)
  : _res(res) {
  // the first position is empty. A value of 0 for a save token
  // indicates an error condition
  SaveSets.push_back(vector< list<SaveAtom>::iterator >() );
  _natoms = _res->atomCount();
}

ConfSaver::~ConfSaver() {
}

void ConfSaver::Save() {
  if (!Residue::isValid(_res) || (_res->atomCount() != _natoms) ) {
    throw "Residue has been modified in the middle of a conformation save";
  }

  SaveSets.push_back(vector< confAtomIter >());   //Add a vector of atoms for this conformer,
  SaveSets.back().reserve(_natoms);       //and make sure we have space for the atoms

  vector< confAtomIter >& current = SaveSets.back();   //Store a ref to the current conformer

  for (int i = 0; i < _natoms; i++) {
    if (_res->atom(i)->search_flag() || SaveSets.size() == 2) {
      current.push_back(
        atom_store.insert(atom_store.end(), SaveAtom(_res->atom(i))));
    } else {
      current.push_back(SaveSets[SaveSets.size()-2][i]);
    }
  }
}

void ConfSaver::Restore(unsigned int token) const {
  vector< list<SaveAtom>::iterator >::const_iterator i, e = SaveSets[token].end();
  for (i = SaveSets[token].begin(); i < e; ++i) {
    (*i)->Restore();
  }
}

void ConfSaver::RestoreLast() const {
  vector< list<SaveAtom>::iterator >::const_iterator i, e = SaveSets.back().end();
  for (i = SaveSets.back().begin(); i < e; ++i) {
    (*i)->Restore();
  }
}

void ConfSaver::ConvertToGeomSaver(GeomSaver& gs, MIMoleculeBase* model) {
  for (int i = 1; i <= NumberSets(); ++i) {
    Restore(i);
    gs.Save(_res, 1, model);
  }
}

int ConfSaver::NumberSets() const {
  return SaveSets.size() - 1;
}

const RESIDUE* ConfSaver::GetResidue() const {
  return _res;
}

} //namespace chemlib

