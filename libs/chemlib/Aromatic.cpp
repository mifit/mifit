#include <vector>
#include <sstream>
#include <fstream>
#include "Aromatic.h"
#include "Dictionary.h"
#include "mol_util.h"

namespace chemlib {


// Purpose:		Remove all data from the object
void Aromatic::Clear() {
  _atoms.clear();
  _bonds.clear();
}

// Purpose:		Include an existing atom in this aromatic system
void Aromatic::AddAtom(MIAtom* atom) {
  if (std::find(_atoms.begin(), _atoms.end(), atom) == _atoms.end()) {    //Don't duplicate
    _atoms.push_back(atom);
  }
}

// Purpose:		Include an existing bond in this aromatic system
void Aromatic::AddBond(Bond* bond) {
  if (std::find(_bonds.begin(), _bonds.end(), bond) == _bonds.end()) {    //Don't duplicate
    _bonds.push_back(bond);
  }
}

// Purpose:		Test whether an atom is in this aromatic system
// Output:      True if the atom is in this aromatic, false otherwise
// Requires:	The method DetectAromatics() has been run for the parent ring system
bool Aromatic::Contains(const MIAtom* query) const {
  return std::find(_atoms.begin(), _atoms.end(), query) != _atoms.end();
}

// Purpose:		Test whether a given bond is in this aromatic system
// Output:      True if the atom is in this aromatic, false otherwise
// Requires:	The method DetectAromatics() has been run for the parent ring system
bool Aromatic::Contains(const Bond* query) const {
  return std::find(_bonds.begin(), _bonds.end(), query) != _bonds.end();
}

// Purpose:		Get the # of atoms in this aromatic system
// Output:      # of atoms
// Requires:	The method DetectAromatics() has been run for the parent ring system
int Aromatic::NumAtoms() const {
  return _atoms.size();
}

/////////////////////////////////////////////////////////////////////////////
// Function:    NumBonds
// Purpose:		Get the # of bonds in this aromatic system
// Input:       Looked up from _bonds member
// Output:      # of bonds
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

int Aromatic::NumBonds() const {
  return _bonds.size();
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Print
// Purpose:		Reports summary information for the aromatic system
// Input:       Looked up from members of the Aromatic class
// Output:      Summary appended to a std::string
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

void Aromatic::Print(std::string& s) const {
  MIAtom_const_iter atm;
  std::vector<Bond*>::const_iterator bond;

  s += "Printing Aromatic...";
  s += "\n";

  for (atm = _atoms.begin(); atm != _atoms.end(); ++atm) {
    s += (*atm)->name();
    s += " ";
  }
  s += "\n";

  for (bond = _bonds.begin(); bond != _bonds.end(); ++bond) {
    s +=  (*bond)->getAtom1()->name();
    s += "-";
    s += (*bond)->getAtom2()->name();
    s += " ";
  }
  s += "\n";
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateConnTable
// Purpose:		Stores a matrix of booleans, where mat[i][j] is true if and only
//				atoms i and j are bonded
// Input:       Looked up from data stored in the object, and in the parent molecule
// Output:      Matrix written to _conn_table member
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

void Aromatic::GenerateConnTable() {
  _conn_table.newsize(NumAtoms(), NumAtoms());
  int i, j;

  for (i = 0; i < _conn_table.num_rows(); ++i) {             //Initialize the array with
    for (j = 0; j < _conn_table.num_cols(); ++j) {      //no connections
      _conn_table[i][j] = false;
    }
  }

  int xatom1, xatom2;
  std::vector<Bond*>::iterator bnd;
  for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd) {
    xatom1 = GetAtomIndex((*bnd)->getAtom1());              //Loop thru bonds, creating the
    xatom2 = GetAtomIndex((*bnd)->getAtom2());              //connections in a symmetric
    //matrix
    if (xatom1 != -1 && xatom2 != -1) {
      _conn_table[xatom1][xatom2] = true;
      _conn_table[xatom2][xatom1] = true;
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateSRData
// Purpose:		Stores the size of the smallest aromatic ring in which each atom
//				and bond is contained
// Input:       Looked up from data stored in the object, and in the parent molecule
// Output:      Ring sizes written to the MIAtom and Bond struct for each atom and bond
// Requires:	The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for the aromatic system
/////////////////////////////////////////////////////////////////////////////

void Aromatic::GenerateSRData() {
  MIAtom_const_iter atm;
  for (atm = _atoms.begin(); atm != _atoms.end(); ++atm) {
    (*atm)->set_smallest_aromatic_ring(SmallestRing(*atm));
  }

  std::vector<Bond*>::const_iterator bnd;
  for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd) {
    (*bnd)->smallest_aromatic_ring = SmallestRing(*bnd);
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GeneratePlane
// Purpose:		Gathers atoms in the plane of the aromatic system
// Input:       A reference to the dictionary for the parent molecule
// Output:      Stores the Plane object in the dictionary
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

void Aromatic::GeneratePlane(LigDictionary& dict) const {
  MIAtom_const_iter atm;
  MIAtom_const_iter nabor;

  Plane pln(dict);
  //	pln.SetResidue(dict.GetResidue());

  for (atm = _atoms.begin(); atm != _atoms.end(); ++atm) {
    pln.AddAtom(*atm);

    if ((*atm)->hybrid() == 2) {
      for (nabor = (*atm)->nabors().begin(); nabor != (*atm)->nabors().end(); ++nabor) {
        pln.AddAtom(*nabor);                    //Neighboring atoms are also
        //part of the ring plane
      }
    }
  }
  dict.AddPlane(pln);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    GenerateImpropers
// Purpose:		Enumerates all the torsions that are constrained by the aromatic
//				and assign the proper value
// Input:       A reference to the dictionary for the parent molecule
// Output:      Stores the Improper objects in the dictionary
// Requires:	The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

void Aromatic::GenerateImpropers(LigDictionary& dict) const {
  std::vector<Bond*>::const_iterator bnd;
  std::vector< MIAtomList > torsions;
  std::vector< MIAtomList >::iterator atms;
  Improper imp;
  int sr, x_tightest = 0, min_size;


  for (bnd = _bonds.begin(); bnd != _bonds.end(); ++bnd) {
    if (SmallestRing(*bnd) >= 10) {                         //For (very rare) large aromatic
      dict.UnspecDoubleBond(*bnd);                          //rings, allow torsion be either
    }                                                       //0 or 180 degrees (i.e. cis or trans)

    min_size = _atoms.size() + 1;

    torsions.clear();
    EnumerateTorsions(*bnd, torsions);

    for (atms = torsions.begin(); atms != torsions.end(); ++atms) {

      if (Contains((*atms)[0])                              //Check if the end atoms are in
          && Contains((*atms)[3])) {                        //the aromatic system,

        sr = SmallestRing((*atms)[0],
               (*atms)[1],
               (*atms)[2],
               (*atms)[3]);

        x_tightest = (sr < min_size) ?                          //Tracks the index of the
                     atms - torsions.begin() : x_tightest;      //torsion in smallest ring
        min_size = (sr < min_size) ? sr : min_size;
      }
    }

    for (atms = torsions.begin(); atms != torsions.end(); ++atms) {

      switch (HammingDistance(*atms, torsions[x_tightest])) {
        case 0: imp.ReInit(*atms, 0); break;                //the tightest one is set to 0
        case 1: imp.ReInit(*atms, 180); break;              //flips from tightest set to 180
        case 2: imp.ReInit(*atms, 0); break;                //angle opposite tightest is also 0
      }


      // This conditional addresses the rare case of sp3 atoms in aromatic rings.  (The only case
      // I know of is for cyclic sulfonamides, as in saccharin.)  In these cases, no
      // impropers are set involving the substituent atoms of the sp3 aromatic atom.
      if ((Contains((*atms)[0]) || (*atms)[1]->hybrid() == 2)
          && (Contains((*atms)[3]) || (*atms)[2]->hybrid() == 2)) {
        dict.AddImproper(imp, true);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SmallestRing
// Input:       A pointer to an atom
// Output:      Returns the size (in atoms) of the smallest aromatic ring
//				containing the atom
// Requires:	The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

int Aromatic::SmallestRing(const MIAtom* atom) const {
  int xatom;
  if ((xatom = GetAtomIndex(atom)) == -1) {         //Check that this atom is in the ring
    return -1;
  }

  int* path = new int[_atoms.size()];       //Get space for the arrays used in the
  bool* used = new bool[_atoms.size()];     //search
  InitializeArray(used, _atoms.size(), false);
  int min = _atoms.size() + 1;

  path[0] = xatom;                          //Initialize the search arrays, starting
  used[xatom] = true;                       //with the input atom

  ExtendPath(path, used, 1, min);           //Perform the recursive, depth-first search

  delete[] path;
  delete[] used;

  return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SmallestRing
// Purpose:		Perform a recursive, depth-first search to detect the smallest cycle
//				in a (molecular) graph
// Input:       A pointer to a bond
// Output:      Returns the size (in atoms) of the smallest aromatic ring
//				containing the bond
// Requires:	The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

int Aromatic::SmallestRing(const Bond* bond) const {
  int xatom1, xatom2;

  if ((xatom1 = GetAtomIndex(bond->getAtom1())) == -1) {    //Check that the 1st atom is in the ringsys
    return -1;
  }
  if ((xatom2 = GetAtomIndex(bond->getAtom2())) == -1) {    //Check that the 2nd atom is in the ringsys
    return -1;
  }

  int* path = new int[_atoms.size()];       //Get space for the arrays used in the
  bool* used = new bool[_atoms.size()];     //search
  InitializeArray(used, _atoms.size(), false);
  int min = _atoms.size();

  path[0] = xatom1;
  path[1] = xatom2;                         //Initialize the search arrays, starting
  used[xatom1] = true;                      //with the input atom
  used[xatom2] = true;

  ExtendPath(path, used, 2, min);           //Perform the recursive, depth-first search

  delete[] path;
  delete[] used;

  return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SmallestRing
// Input:       Pointers to three consecutive atoms.
// Output:      Returns the size (in atoms) of the smallest aromatic ring
//				containing the atoms
// Requires:	The input args be ordered correctly.  That is, atom1 must be bonded
//				to atom2 and atom2 must be bonded to atom3.
//				The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

int Aromatic::SmallestRing(const MIAtom* atom1,
                           const MIAtom* atom2,
                           const MIAtom* atom3) const {

  int xatom1, xatom2, xatom3;

  if ((xatom1 = GetAtomIndex(atom1)) == -1) {       //Check that the 1st atom is in the ringsys
    return -1;
  }
  if ((xatom2 = GetAtomIndex(atom2)) == -1) {       //Check that the 2nd atom is in the ringsys
    return -1;
  }
  if ((xatom3 = GetAtomIndex(atom3)) == -1) {       //Check that the 3rd atom is in the ringsys
    return -1;
  }

  int* path = new int[_atoms.size()];       //Get space for the arrays used in the
  bool* used = new bool[_atoms.size()];     //search
  InitializeArray(used, _atoms.size(), false);
  int min = _atoms.size();

  path[0] = xatom1;                         //Assumes the three atoms are in sequence
  path[1] = xatom2;
  path[2] = xatom3;                         //Initialize the search arrays, starting
  used[xatom1] = true;                      //with the input atoms
  used[xatom2] = true;
  used[xatom3] = true;

  ExtendPath(path, used, 3, min);           //Perform the recursive, depth-first search

  delete[] path;
  delete[] used;

  return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    SmallestRing
// Input:       Pointers to four consecutive atoms (in order)
// Output:      Returns the size (in atoms) of the smallest aromatic ring
//				containing the atoms
// Requires:	The input args be ordered correctly.  That is, atom1 must be bonded
//				to atom2 and atom2 must be bonded to atom3, etc.
//				The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

int Aromatic::SmallestRing(const MIAtom* atom1,
                           const MIAtom* atom2,
                           const MIAtom* atom3,
                           const MIAtom* atom4) const {

  int xatom1, xatom2, xatom3, xatom4;

  if ((xatom1 = GetAtomIndex(atom1)) == -1) {       //Check that the 1st atom is in the ringsys
    return -1;
  }
  if ((xatom2 = GetAtomIndex(atom2)) == -1) {       //Check that the 2nd atom is in the ringsys
    return -1;
  }
  if ((xatom3 = GetAtomIndex(atom3)) == -1) {       //Check that the 3rd atom is in the ringsys
    return -1;
  }
  if ((xatom4 = GetAtomIndex(atom4)) == -1) {       //Check that the 4th atom is in the ringsys
    return -1;
  }

  int* path = new int[_atoms.size()];       //Get space for the arrays used in the
  bool* used = new bool[_atoms.size()];     //search
  InitializeArray(used, _atoms.size(), false);
  int min = _atoms.size();

  path[0] = xatom1;                         //Assumes the four atoms are in sequence
  path[1] = xatom2;
  path[2] = xatom3;                         //Initialize the search arrays, starting
  path[3] = xatom4;                         //with the input atoms
  used[xatom1] = true;
  used[xatom2] = true;
  used[xatom3] = true;
  used[xatom4] = true;

  ExtendPath(path, used, 4, min);           //Perform the recursive, depth-first search

  delete[] path;
  delete[] used;

  return min;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ExtendPath
// Purpose:		Recursive, depth-first search for the smallest cycle in a molecule
// Input:       Sequence of indices of atoms representing a path thru the mol
//				Array of booleans tracking which atoms have been used in the path
//				Count of the number of atoms used in the path
// Output:      Stores the size of the smallest ring found so far in the integer
//				referred to in the final argument.
// Requires:	The method DetectAromatics() has been run for the parent ring system
//				The method GenerateConnTable() has been run for this object
/////////////////////////////////////////////////////////////////////////////

void Aromatic::ExtendPath(int* path,            //Ordered array of atom indices in the path
                          bool* used,           //Array that tracks which atoms have been used
                          int depth,            //Number of atoms in current path
                          int& min) const {      //Size of smallest ring found so far
  unsigned int i;
  for (i = 0; i < _atoms.size(); ++i) {              //Loop thru atoms in ringsys
    if (_conn_table[path[depth-1]][i] == false) {
      continue;                                 //Skip if this atom not connected to last
    }

    if (i == (unsigned int)path[0] && depth > 2) {          //Check if this completes a cycle
      min = depth;
    } else if (used[i] == false && depth + 1 < min) {
      path[depth] = i;
      used[i] = true;
      ExtendPath(path, used, depth + 1, min);
      used[i] = false;                          //Reset for next bond
    }
  }                                             //End loop over atoms
}

//////////////////////////////////////////////////////////////////////////////
// Function:    GetAtomIndex
// Purpose:		Get the index of an atom used in the internal storage of the
//				Aromatic object (i.e. _atoms, _conn_table, and the local variables
//				in the ring-finding & search functions.)
// Input:       Pointer to an atom contained in the aromatic
// Output:      Position of the atom in the _atoms vector
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

int Aromatic::GetAtomIndex(const MIAtom* query) const {
  MIAtom_const_iter atom_location;

  atom_location = std::find(_atoms.begin(), _atoms.end(), query);

  if (atom_location == _atoms.end()) {
    return -1;
  } else {
    return atom_location - _atoms.begin();
  }
}

}
