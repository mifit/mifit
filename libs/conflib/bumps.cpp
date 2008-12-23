#include "bumps.h"

#include "chemlib.h"
#include "RESIDUE_.h"

using namespace chemlib;

namespace conflib {

/////////////////////////////////////////////////////////////////////////////
// Function:    AssignBump
// Purpose:		Given two atoms, stores the SQUARE of the minimum distance of approach between
//				the two atoms
// Input:       Ptrs to the two atoms, an Bond to store the distance in, and a list of the bonds
//				in the residue or molecule
// Output:      Stores the SQUARE of the distance in the ideal_length field of the given Bond
//				Returns false, (and stores nothing) if the atoms won't bump in any conformation
// Requires:	The method FindRingSystems() has been run for the residue or molecule
/////////////////////////////////////////////////////////////////////////////
bool AssignBump(MIAtom* atom1, MIAtom* atom2, Bond& bump, const std::vector<Bond>& bonds) {

  if (AlreadyBonded(atom1, atom2, bonds)) {                         //Check if atoms are bonded
    return false;
  }

  MIAtomList nabors1, nabors2;
  GetNabors(atom1, bonds, nabors1);
  GetNabors(atom2, bonds, nabors2);

  MIAtom_iter nbr1;                          //Check if atoms share a neighbor
  MIAtom_iter nbr2;                          //(1-3 interactions)
  for (nbr1 = nabors1.begin(); nbr1 != nabors1.end(); ++nbr1) {
    for (nbr2 = nabors2.begin(); nbr2 != nabors2.end(); ++nbr2) {
      if (*nbr1 == *nbr2) {
        return false;
      }
    }
  }

  if (atom1->iscyclic()                                           //Don't track bumps between
      && atom2->iscyclic()                                        //atoms in the same ring system
      && (atom1->ring_system() == atom2->ring_system())) {
    return false;
  }

  bump.setAtom1(atom1);
  bump.setAtom2(atom2);
  bump.ideal_length = 1.6f * (CovalentRadius(bump.getAtom1()->atomicnumber())
                              + CovalentRadius(bump.getAtom2()->atomicnumber()));

  bump.ideal_length *= bump.ideal_length;                       //Stores the SQUARE of the minimum
  return true;                                                  //distance to speed comparisons
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GetBumps
// Purpose:		Loop to generate bump info for all pairs of atoms in residue
// Input:       Ptrs to the residue of interest, vector of Bonds to store the bumps,
//				vector of Bonds for all the bonds in the residue
// Output:      Stores the bump info in Bond structs
// Requires:	The method FindRingSystems() has been run for the residue
/////////////////////////////////////////////////////////////////////////////
void GetBumps(const RESIDUE* res, std::vector<Bond>& bumps, const std::vector<Bond>& bonds) {
  int i, j;
  Bond bump;
  for (i = 0; i < res->atomCount(); ++i) {
    for (j = i+1; j < res->atomCount(); ++j) {
      if (AssignBump(res->atom(i), res->atom(j), bump, bonds)) {
        bumps.push_back(bump);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    CheckBumps
// Purpose:		Given a residue, checks whether the atoms that have been moved
//				have resulted in any interatomic clashes
// Input:       Ptr to the RESIDUE, and a vector containing the bumps to check
// Output:      True if there are no bumps, false otherwise
// Requires:	The search_flag fields of the atoms have been set to a non-zero
//				value for atoms that have been moved since the last check.  (For
//				the first check, set all the search_flags.)
/////////////////////////////////////////////////////////////////////////////
bool CheckBumps(RESIDUE*, std::vector<Bond>& bumps) {
  std::vector<Bond>::iterator i = bumps.begin();
  std::vector<Bond>::iterator e = bumps.end();
  while (i < e) {

    //Only check atom pair for bumps if at least one has been moved
    if ((i->getAtom1()->search_flag() || i->getAtom2()->search_flag())
        && SquaredAtomDist(*i->getAtom1(), *i->getAtom2()) < i->ideal_length) {
      return false;
    }
    i++;
  }
  return true;
}

} //namespace conflib
