#include <chemlib/chemlib.h>
#include "HPhobe.h"

using namespace moldraw;
using namespace chemlib;
using namespace std;

namespace moldraw {
/////////////////////////////////////////////////////////////////////////////
// Function:    GetNonPolars
// Purpose:		Searches a sequence of residues for atoms capable of a hydrophobic
//				interaction. (Used in figure generation)
// Input:       A vector of residues to search
// Output:      A vector (passed by reference) of atoms that meet the criteria of an acceptor
// Requires:	The function IsNonPolar()
/////////////////////////////////////////////////////////////////////////////
void GetNonPolars(const std::vector <Residue*>& residues, std::vector <MIAtom*>& nonpolar_atoms) {

  std::vector<Residue*>::const_iterator res, resEnd = residues.end();
  MIAtom_const_iter atom, endAtom;
  for (res = residues.begin(); res != resEnd; ++res) {
    Residue& r = **res;
    const MIAtomList& atoms = r.atoms();
    endAtom = atoms.end();
    for (atom = atoms.begin(); atom != endAtom; ++atom) {
      if (IsNonPolar(**atom)) {
        nonpolar_atoms.push_back(*atom);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:	CheckHPhobe
// Purpose:		Checks whether a pair of atoms are close enough to form a
//				hydrophobic interaction
// Input:       A pair of atom pointers (order immaterial), and an HPHOBE struct
// Output:      Records the 2 atoms involved,  and the distance of the interaction
//				Returns true if the atoms meet the geometric criteria
// Requires:	MAX_HPHOBE_LENGTH to be defined
/////////////////////////////////////////////////////////////////////////////
bool CheckHPhobe(MIAtom* atom1, MIAtom* atom2, HPhobe& hphobe) {
  if (!IsNonPolar(*atom1) || !IsNonPolar(*atom2)) {
    return false;
  }

  hphobe.distance = AtomDist(*atom1, *atom2);
  if (hphobe.distance < MAX_HPHOBE_LENGTH) {
    hphobe.setAtom1(atom1);
    hphobe.setAtom2(atom2);
    return true;
  } else {
    return false;
  }
}

void RePointHphobe(HPhobe& hphobe, Ligand& site) {

  vector<Residue*>::iterator res;
  vector<MIAtom*>::const_iterator atm;

  for (res = site.residues.begin(); res != site.residues.end(); ++res) {
    for (atm = (*res)->atoms().begin(); atm != (*res)->atoms().end(); ++atm) {
      if (hphobe.getAtom1()->x() == (*atm)->x()
          && hphobe.getAtom1()->y() == (*atm)->y()
          && hphobe.getAtom1()->z() == (*atm)->z()) {
        hphobe.setAtom1(*atm);
        hphobe.res1 = *res;
      }
      if (hphobe.getAtom2()->x() == (*atm)->x()
          && hphobe.getAtom2()->y() == (*atm)->y()
          && hphobe.getAtom2()->z() == (*atm)->z()) {
        hphobe.setAtom2(*atm);
        hphobe.res2 = *res;
      }
    }
  }
  return;
}

} //namespace moldraw
