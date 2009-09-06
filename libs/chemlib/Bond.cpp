#include <vector>
#include "Bond.h"
#include "sequence_util.h"
#include "MIAtom.h"

namespace chemlib {

/////////////////////////////////////////////////////////////////////////////
// Function:    Default constructor
// Purpose:		Initializes data fields
// Input:       None
// Output:      None
// Note:
/////////////////////////////////////////////////////////////////////////////
Bond::Bond() : bondOrder(NORMALBOND), atom1(0), atom2(0), type(B_NORMAL),
  stereo(STEREO_NONE), ideal_length(-1.0f), tolerance(0.2f),
  dict_include(0), isaromatic(0), iscyclic(0), ring_system(0),
  smallest_ring_size(0), smallest_aromatic_ring(0) {
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Constructor
// Purpose:		Creates a bond object from an existing Bond struct
// Input:       A Bond object (passed by ref)
//				A vector of MIAtom structs containing the atoms in the bond
//				A vector of corresponding MIAtom objects to use for the new ptrs
// Output:      None
// Note:		Correspondence between MIAtom is inferred from their bondOrder
//				in the vector containers
/////////////////////////////////////////////////////////////////////////////
Bond::Bond(const Bond& bond,
           const MIAtomList& old_atoms,
           const MIAtomList& new_atoms) {

  std::string error_message;
  int i1 = GetIndex(static_cast<MIAtom*> (bond.getAtom1()), old_atoms);
  int i2 = GetIndex(static_cast<MIAtom*> (bond.getAtom2()), old_atoms);

  if (i1 < 0 || i2 < 0) {
    error_message = "Corrupted bond.";
    throw error_message;
  }

  bondOrder = bond.bondOrder;
  atom1 = new_atoms[i1];
  atom2 = new_atoms[i2];
  type = bond.type;
  stereo = bond.stereo;
  ideal_length = bond.ideal_length;
  tolerance = bond.tolerance;
  dict_include = bond.dict_include;
  isaromatic = bond.isaromatic;
  iscyclic = bond.iscyclic;
  ring_system = bond.ring_system;
  smallest_ring_size = bond.smallest_ring_size;
  smallest_aromatic_ring = bond.smallest_aromatic_ring;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Constructor
// Purpose:		Creates a bond object with the given parameters
// Input:       Two atom pointers for the two ends of the bond
//				A code for the bond bondOrder (see Xguicryst for codes)
//				A code for the stereochemistry (see Xguicryst for codes)
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
Bond::Bond(MIAtom* in_atom1,
           MIAtom* in_atom2,
           unsigned char in_bondOrder,
           char in_stereo) {
  Clear();

  atom1 = in_atom1;             //Input values
  atom2 = in_atom2;
  bondOrder = in_bondOrder;
  stereo = in_stereo;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ClearOrder
// Purpose:		Sets the bondOrder to a single bond, and sets aromaticity and stereo
//				flags to false.
// Input:       None
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Bond::ClearOrder() {
  bondOrder = SINGLEBOND;
  isaromatic = 0;
  stereo = 0;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    UpdateAtoms
// Purpose:		Translates the atom ptrs in a bond, usually when the atoms have
//				been copied to a new location
// Input:       Two vectors of atom ptrs that represent the old and new
//				locations of the atoms in the molecule
// Output:      True if the constituent atoms are in the "old_atoms" vector,
//				false otherwise
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool Bond::UpdateAtoms(const MIAtomList& old_atoms,
                       const MIAtomList& new_atoms) {

  int i1 = GetIndex(static_cast<MIAtom*> (atom1), old_atoms);
  int i2 = GetIndex(static_cast<MIAtom*> (atom2), old_atoms);

  if (i1 < 0 || i2 < 0) {
    return false;
  } else {
    atom1 = new_atoms[i1];
    atom2 = new_atoms[i2];
    return true;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    IsRotatable
// Purpose:		Eliminate cyclic and terminal bonds from consideration in
//				conformational searches.
// Input:       None
// Output:      True or false
// Requires:	The ring systems of the molecule must have been determined
/////////////////////////////////////////////////////////////////////////////

bool Bond::IsRotatable() const {
  if (iscyclic) {                                    //If this bond is in a ring, it is
    return false;                                   //not rotatable
  } else if (atom1->nabors().size() == 1) {              //If this bond is terminal, it is not
    return false;                                   //rotatable
  } else if (atom2->nabors().size() == 1) {
    return false;
  }

  return true;                                      //If we reach this line, the bond is
}                                                   //rotatable


/////////////////////////////////////////////////////////////////////////////
// Function:    Print
// Purpose:		Reports summary information for the bond
// Input:       None
// Output:      Summary returned as a std::string
// Requires:	The method DetectAromatics() has been run for the parent ring system
/////////////////////////////////////////////////////////////////////////////

std::string Bond::Print() const {

  std::string s;
  std::string bo;

  switch (bondOrder) {
    case 0: bo = "?"; break;
    case SINGLEBOND: bo = "single"; break;
    case DOUBLEBOND: bo = "double"; break;
    case TRIPLEBOND: bo = "triple"; break;
    case PARTIALDOUBLEBOND: bo = "aromatic"; break;
    case HYDROGENBOND: bo = "h-bond"; break;
    default: bo = "?";
  }

  s += atom1->name();
  s += "-";
  s += atom2->name();
  s += " ";
  s += bondOrder + '0';
  s += " ";
  s += bo;
  s += " ";
  s += isaromatic + '0';
  s += "\n";

  return s;
}

void Bond::copyBondOrders(const std::vector<Bond>& fromBonds, std::vector<Bond>& toBonds) {
  std::map<std::pair<MIAtom*, MIAtom*>, unsigned char> bondOrders;

  std::vector<Bond>::const_iterator constBondIter;
  for (constBondIter = fromBonds.begin(); constBondIter != fromBonds.end(); ++constBondIter) {
    Bond b = *constBondIter;
    if (b.getAtom1() < b.getAtom2()) {
      bondOrders[std::make_pair(b.getAtom1(), b.getAtom2())] = b.getOrder();
    } else {
      bondOrders[std::make_pair(b.getAtom2(), b.getAtom1())] = b.getOrder();
    }
  }

  std::vector<Bond>::iterator bondIter;
  for (bondIter = toBonds.begin(); bondIter != toBonds.end(); ++bondIter) {
    Bond b = *bondIter;
    if (b.getAtom1() < b.getAtom2()) {
      bondIter->setOrder(bondOrders[std::make_pair(b.getAtom1(), b.getAtom2())]);
    } else {
      bondIter->setOrder(bondOrders[std::make_pair(b.getAtom2(), b.getAtom1())]);
    }
  }

}

}
