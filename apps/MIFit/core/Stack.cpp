#include <stdio.h>

#include "nonguilib.h"
#include "chemlib.h"
#include "RESIDUE_.h"

#include "Stack.h"
#include "RESIDUE.h"
#include "Molecule.h"

using namespace chemlib;

Stack::Stack() {
  changed = false;
  minimized = false;
}

Stack::~Stack() {
}

void Stack::moleculeToBeDeleted(MIMoleculeBase* mol) {
  Purge(mol);

  // Purge may not delete everything, if there was an atom pushed on the
  // stack where mol and/or res was null, the stack might still have a ref
  // to it, so we have to try harder
  std::vector<RESIDUE*> residues;
  for (MIIter<RESIDUE> res=mol->GetResidues(); res; ++res) {
    residues.push_back(res);
  }
  for (MIIter<RESIDUE> res=mol->GetSymmResidues(); res; ++res) {
    residues.push_back(res);
  }
  residuesToBeDeleted(mol,residues);
}

void Stack::residuesToBeDeleted(MIMoleculeBase* mol, std::vector<RESIDUE*>& residues) {
  std::vector<RESIDUE*>::iterator iter;
  for (iter = residues.begin(); iter != residues.end(); ++iter) {
    RESIDUE* residue = *iter;
    Purge(residue);
  }
  
  // Purge may not delete everything, if there was an atom pushed on the
  // stack where mol and/or res was null, the stack might still have a ref
  // to it, so we have to try harder
  MIAtomList atoms;
  for (size_t i=0; i < residues.size(); ++i) {
    atoms.insert(atoms.end(), residues[i]->atoms().begin(), residues[i]->atoms().end());
  }
  atomsToBeDeleted(mol,atoms);

}

void Stack::atomsToBeDeleted(MIMoleculeBase*, const MIAtomList& atoms) {
  for (unsigned int i=0;i<atoms.size(); ++i) {
    Purge(atoms[i]);
  }
}


void Stack::Push(MIAtom* atom, RESIDUE* res, Molecule* m) {
  if (!MIAtom::isValid(atom) || !Residue::isValid(res) || !MIMoleculeBase::isValid(m)) {
    return;
  }
  StackItem item;
  item.atom = atom;
  item.residue = res;
  item.molecule = m;
  bool wasEmpty = empty();
  data.push_back(item);
  changed = true;
  if (wasEmpty != empty()) {
    emptyChanged(empty());
  }
}

void Stack::Pop(MIAtom*& atom, RESIDUE*& res) {
  if (data.size() == 0) {
    atom = NULL;
    res = NULL;
    return;
  }
  StackItem item = data.back();
  dataPop();
  atom = item.atom;
  res = item.residue;
}

void Stack::Pop(MIAtom*& atom, RESIDUE*& res, Molecule*& m) {
  if (data.size() == 0) {
    atom = NULL;
    res = NULL;
    m = NULL;
    return;
  }
  StackItem item = data.back();
  dataPop();
  atom = item.atom;
  res = item.residue;
  m = item.molecule;
}

MIAtom* Stack::Pop() {
  if (data.size() == 0) {
    return NULL;
  }
  StackItem item = data.back();
  dataPop();
  return item.atom;
}

void Stack::dataPop() {
  bool wasEmpty = empty();
  data.pop_back();
  changed = true;
  if (wasEmpty != empty()) {
    emptyChanged(empty());
  }
}

void Stack::Peek(MIAtom*& atom, RESIDUE*& res, Molecule*& m) {
  if (data.size() == 0) {
    atom = NULL;
    res = NULL;
    m = NULL;
    return;
  }
  StackItem item = data.back();
  atom = item.atom;
  res = item.residue;
  m = item.molecule;
}

bool Stack::StackChanged() {
  return changed;
}

void Stack::ClearChanged() {
  changed = false;
}

void Stack::ExpandTopAllAtoms() {
  MIAtom* a;
  Molecule* m;
  RESIDUE* r;
  Pop(a, r, m);

  if (a && a->type() & AtomType::SYMMATOM) {
    Push(a, r, m);
    Logger::message("Can't expand symmetry atom(s) on stack");
    return;
  }

  for (int i = 0; i < r->atomCount(); i++) {
    Push(r->atom(i), r, m);
  }
}

void Stack::ExpandTop2AllAtoms() {
  MIAtom* a1;
  MIAtom* a2;
  Molecule* m1, * m2;
  RESIDUE* r1, * r2;
  Pop(a1, r1, m1);
  Pop(a2, r2, m2);
  if (m1 != m2) {
    Logger::message("Both atoms must be in same molecule");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
    return;
  }

  if ((a1 && a1->type() & AtomType::SYMMATOM) ||
      (a2 && a2->type() & AtomType::SYMMATOM)) {
    Logger::message("Can't expand symmetry atom(s) on stack");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
    return;
  }

  int i = 0;
  int i1 = -1;
  int i2 = -1;
  if (r1 == r2) {
    for (i = 0; i < r2->atomCount(); i++) {
      Push(r2->atom(i), r2, m2);
    }
    return;
  }
  MIIterBase<RESIDUE> * rip = 0;
  MIIter<RESIDUE> res = rip, ri_beg = rip, ri_end = rip;
  for (res = m1->GetResidues(); res; ++res) {
    if (res == r1) {
      ri_beg = res;
      i1 = i;
    }
    if (res == r2) {
      ri_end = res;
      i2 = i;
    }
    i++;
  }
  if (i1==-1 || i2==-1) {
    Logger::message("Error expanding stack: residue not found in model");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
  }

  if (i2 < i1) {
    res = ri_beg;
    ri_beg = ri_end;
    ri_end = res;
  }
  for (res = ri_beg; res; ++res) {
    for (i = 0; i < res->atomCount(); i++) {
      Push(res->atom(i), res, m2);
    }
    if ((RESIDUE*)res == (RESIDUE*)ri_end) {
      break;
    }
  }
}

void Stack::ExpandTop2Range() {
  MIAtom* a1, * a2, * a;
  Molecule* m1, * m2;
  RESIDUE* r1, * r2;
  Pop(a1, r1, m1);
  Pop(a2, r2, m2);
  if (m1 != m2) {
    Logger::message("Both atoms must be in same molecule");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
    return;
  }

  if ((a1 && a1->type() & AtomType::SYMMATOM) ||
      (a2 && a2->type() & AtomType::SYMMATOM)) {
    Logger::message("Can't expand symmetry atom(s) on stack");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
    return;
  }

  int i = 0, i1 = (-1), i2 = (-1);
  if (r1 == r2) {
    if ((a = atom_from_name("CA", *r2)) == NULL) {
      Push(a, r2, m2);
    } else {
      Push(r2->atom(0), r2, m2);
    }
    return;
  }

  MIIterBase<RESIDUE> * rip = 0;
  MIIter<RESIDUE> res = rip, ri_beg = rip, ri_end = rip;
  for (res = m1->GetResidues(); res; ++res) {
    if (res == r1) {
      ri_beg = res;
      i1 = i;
    }
    if (res == r2) {
      ri_end = res;
      i2 = i;
    }
    i++;
  }
  if (i1==-1 || i2==-1) {
    Logger::message("Error expanding stack: residue not found in model");
    Push(a2, r2, m2);
    Push(a1, r1, m1);
    return;
  }

  if (i2 < i1) {
    res = ri_beg;
    ri_beg = ri_end;
    ri_end = res;
  }
  for (res = ri_beg; res; ++res) {
    if ((a = atom_from_name("CA", *res)) != NULL) {
      Push(a, res, m2);
    } else {
      Push(res->atom(0), res, m2);
    }
    if ((RESIDUE*)res == (RESIDUE*)ri_end) {
      break;
    }
  }
}

bool Stack::InStack(RESIDUE* res) {
  DataContainer::iterator iter = data.begin();
  while (iter != data.end()) {
    StackItem item = *iter;
    if (item.residue == res) {
      return true;
    }
    ++iter;
  }
  return false;
}

void Stack::Purge(MIMoleculeBase* model) {
  DataContainer::iterator iter = data.begin();
  while (iter != data.end()) {
    StackItem item = *iter;
    if (item.molecule == model) {
      changed = true;
      data.erase(iter);
      iter = data.begin();
      continue;
    }
    ++iter;
  }
}

void Stack::Purge(RESIDUE* res) {
  DataContainer::iterator iter = data.begin();
  while (iter != data.end()) {
    StackItem item = *iter;
    if (item.residue == res) {
      changed = true;
      data.erase(iter);
      iter = data.begin();
      continue;
    }
    ++iter;
  }
}

void Stack::Purge(MIAtom* atom) {
  DataContainer::iterator iter = data.begin();
  while (iter != data.end()) {
    StackItem item = *iter;
    if (item.atom == atom) {
      changed = true;
      data.erase(iter);
      iter = data.begin();
      continue;
    }
    ++iter;
  }
}

