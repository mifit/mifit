#include "mathlib.h"

#include "Constraint.h"
#include "Dictionary.h"
#include "MIAtom.h"

#include <algorithm>

namespace chemlib {

Bond BondLength::ToBond(const MIAtomList& new_atoms,
                        const std::vector<const MIAtom*>& old_atoms) {
  Bond bnd;
  std::string error_message;
  int i1 = GetIndex<const MIAtom*>(atom1, old_atoms);
  int i2 = GetIndex<const MIAtom*>(atom2, old_atoms);

  if (i1 < 0 || i2 < 0) {
    error_message = "Corrupted bond.";
    throw error_message;
  }

  bnd.setAtom1(new_atoms[i1]);
  bnd.setAtom2(new_atoms[i2]);
  bnd.ideal_length = ideal_dist;
  bnd.tolerance = tolerance;
  bnd.dict_include = dict_include;

  return bnd;
}

ANGLE Angle::ToANGLE(const MIAtomList& new_atoms,
                     const std::vector<const MIAtom*>& old_atoms) {
  ANGLE ang;

  std::string error_message;
  int i1 = GetIndex<const MIAtom*>(atom1, old_atoms);
  int i2 = GetIndex<const MIAtom*>(atom2, old_atoms);
  int i3 = GetIndex<const MIAtom*>(atom3, old_atoms);

  if (i1 < 0 || i2 < 0 || i3 < 0) {
    error_message = "Corrupted bond.";
    throw error_message;
  }

  ang.setAtom1(new_atoms[i1]);
  ang.setAtom2(new_atoms[i2]);
  ang.atom3 = new_atoms[i3];
  ang.ideal_angle = ideal_angle;
  ang.tolerance = tolerance;
  ang.iscyclic = iscyclic;
  ang.isaromatic = isaromatic;
  ang.ring_system = ring_system;
  ang.smallest_ring_size = smallest_ring_size;
  ang.smallest_aromatic_ring = smallest_aromatic_ring;

  return ang;
}

//Default Constructor
Plane::Plane() {
  _tol_set = false;
}

//Construct a plane with a given tolerance
Plane::Plane(float tol) {
  _tolerance = tol;
  _tol_set = true;
}

Plane::Plane(LigDictionary& dict) {
  _tolerance = dict.GetSigmaPlane();
  _tol_set = true;
}

void Plane::AddAtom(MIAtom* new_atom) {
  if (std::find(atoms.begin(), atoms.end(), new_atom) == atoms.end()) {   //Add new atom with
    atoms.push_back(new_atom);                                          //no duplication
  }
}

MIAtom* Plane::GetAtom(int x_atom) {
  if (x_atom >= 0 && x_atom < (int)atoms.size()) {                              //Make sure we get
    return atoms[x_atom];                                               //a valid index
  } else {
    return 0;                                                           //If invalid index,
  }                                                                     //return null ptr
}

int Plane::NumAtoms() {
  return atoms.size();
}

/*
   void Plane::Print() {
    std::vector <MIAtom *>::iterator atm;

    for(atm=atoms.begin();atm!=atoms.end();++atm) {						//Prints a list of
        std::cout << (*atm)->name << " ";									//atom names
    }
    std::cout << endl;
   }
 */
void Plane::Clear() {
  atoms.clear();
}

void Plane::Card(std::string& s) {
  char buf[10];                                                         //Appends the string with the
  sprintf(buf, " %d", atoms.size());                                    //number of atoms in the plane
  s.append(buf);

  MIAtom_iter atm;                                          //Appends the string
  for (atm = atoms.begin(); atm != atoms.end(); ++atm) {                       //with a space-separated
    s.append(" ");                                                     //list of atom names
    s.append((*atm)->name());
  }
}

bool Plane::operator==(const Plane& pln2) const {

  if (this->atoms.size() != pln2.atoms.size()) {
    return false;
  }

  MIAtom_const_iter atm;

  for (atm = this->atoms.begin(); atm != this->atoms.end(); ++atm) {
    if (std::find(pln2.atoms.begin(), pln2.atoms.end(), *atm) ==
        pln2.atoms.end()) {
      return false;
    }
  }

  return true;
}

PLANE Plane::ToPLANE(const MIAtomList& new_atoms,
                     const std::vector<const MIAtom*>& old_atoms,
                     const std::vector<RESIDUE*>& new_residues,
                     const std::vector<const Residue*>& old_residues) {
  PLANE pln;
  pln.natoms = atoms.size();

  std::string error_message;

  if ((pln.atoms = (MIAtom**)malloc(pln.natoms * sizeof(MIAtom*))) == 0) {
    error_message = "Out of memory in ToPLANE.";
    throw error_message;
  }

  int xAtom;
  for (int i = 0; i < pln.natoms; ++i) {
    xAtom = GetIndex<const MIAtom*>(atoms[i], old_atoms);

    if (xAtom < 0) {
      error_message = "Corrupted atom information in plane.";
      throw error_message;
    }

    pln.atoms[i] = new_atoms[xAtom];
  }

  pln.tolerance = _tolerance;

  int xRes = GetIndex<const Residue*>(_res, old_residues);

  if (xRes < 0) {
    error_message = "Corrupted residue information in plane.";
    throw error_message;
  }

  pln.res = new_residues[xRes];

  return pln;
}

void Torsion::AddAtom(MIAtom* new_atom) {
  if (std::find(atoms.begin(), atoms.end(), new_atom) == atoms.end()) {       //Add new atom with
    atoms.push_back(new_atom);                                          //no duplication
  }
}

/*
   void Torsion::Print() {
    std::vector <MIAtom *>::iterator atm;

    for(atm=atoms.begin();atm!=atoms.end();++atm) {
        std::cout << (*atm)->name  << " ";
    }
    std::cout << endl;
   }
 */
void Torsion::Clear() {
  atoms.clear();
  _torsnumber = 0;
}

void Torsion::SetIndex(int torsnumber) {
  _torsnumber = torsnumber;
}

void Torsion::Card(std::string& s) {
  char buf[6];
  sprintf(buf, "%d", _torsnumber);

  s.append(buf);
  MIAtom_iter atm;
  for (atm = atoms.begin(); atm != atoms.end(); ++atm) {
    s.append(" ");
    s.append((*atm)->name());
  }
}

TORSION Torsion::ToTORSION(const MIAtomList& new_atoms,
                           const std::vector<const MIAtom*>& old_atoms,
                           const std::vector<RESIDUE*>& new_residues,
                           const std::vector<const Residue*>& old_residues) {
  TORSION tors;

  std::string error_message;
  int i1 = GetIndex<const MIAtom*>(atoms[0], old_atoms);
  int i2 = GetIndex<const MIAtom*>(atoms[1], old_atoms);
  int i3 = GetIndex<const MIAtom*>(atoms[2], old_atoms);
  int i4 = GetIndex<const MIAtom*>(atoms[3], old_atoms);

  if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
    error_message = "Corrupted atom ptr in torsion.";
    throw error_message;
  }

  strcpy(tors.type, "CHI");
  sprintf(tors.type + 3, "%d", _torsnumber);

  tors.setAtom1(new_atoms[i1]);
  tors.setAtom2(new_atoms[i2]);
  tors.atom3 = new_atoms[i3];
  tors.atom4 = new_atoms[i4];
  tors.nideal = 0;

  int xRes = GetIndex<const Residue*>(_res, old_residues);

  if (xRes < 0) {
    error_message = "Corrupted residue information in torsion.";
    throw error_message;
  }

  tors.res = new_residues[xRes];

  return tors;
}

bool Torsion::operator==(const Torsion& tor2) const {

  //Return Equal(forward_compare) || Equal(reverse_compare)
  if (std::equal(this->atoms.begin(), this->atoms.end(), tor2.atoms.begin())) {
    return true;
  }
  if (std::equal(this->atoms.rbegin(), this->atoms.rend(), tor2.atoms.begin())) {
    return true;
  } else {
    return false;
  }
}

void Improper::AddAtom(MIAtom* new_atom) {
  if (std::find(atoms.begin(), atoms.end(), new_atom) == atoms.end()) {       //Add new atom with
    atoms.push_back(new_atom);                                          //no duplication
  }
}

MIAtom* Improper::GetAtom(int x_atom) {
  if (x_atom < (int)atoms.size() && x_atom >= 0) {                              //Make sure we get
    return atoms[x_atom];                                               //a valid index
  } else {
    return 0;                                                           //If invalid index,
  }                                                                     //return null ptr
}

const MIAtom* Improper::GetAtom(int x_atom) const {
  if (x_atom < (int)atoms.size() && x_atom >= 0) {                              //Make sure we get
    return atoms[x_atom];                                               //a valid index
  } else {
    return 0;                                                           //If invalid index,
  }                                                                     //return null ptr
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ReInit
// Purpose:		Reinitialize the atoms and specify an ideal angle of the improper
// Input:       A four-atom sequence to define the torsion and the value
//				of the ideal dihedral angle.
// Output:      Redefines "this" Improper object using the input atoms & angle
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Improper::ReInit(const MIAtomList& in_atoms,
                      double in_angle) {
  Clear();

  for (int i = 0; i < 4; ++i) {
    AddAtom(in_atoms[i]);
  }

  AddAngle(in_angle);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ReInit
// Purpose:		Reinitialize the atoms and ideal angles of the improper
// Input:       A four-atom sequence to define the torsion and a sequence of
//				acceptable values for the dihedral angle.
// Output:      Redefines "this" Improper object using the input atoms & angle
// Requires:
/////////////////////////////////////////////////////////////////////////////
void Improper::ReInit(const MIAtomList& in_atoms,
                      const std::vector<double>& in_angles) {
  unsigned int i;
  Clear();

  for (i = 0; i < 4; ++i) {                              //Add four atoms to the improper
    AddAtom(in_atoms[i]);
  }

  for (i = 0; i < in_angles.size(); ++i) {                   //Add all the angles to the improper
    AddAngle(in_angles[i]);
  }
}

//Don't need this--impropers should always have 4 atoms

//int Improper::NumAtoms() {
//	return atoms.size();
//}

void Improper::AddAngle(double new_angle) {
  _angles.push_back(new_angle);
}

/*
   void Improper::Print() {
    std::vector <MIAtom *>::iterator atm;

    for(atm=atoms.begin();atm!=atoms.end();++atm) {
        std::cout << (*atm)->name << " ";
    }

    std::vector <double>::iterator ang;
    for(ang=_angles.begin(); ang!=_angles.end(); ++ang) {
        std::cout << "\t" << *ang;
    }

    std::cout << endl;
   }
 */
//const std::vector<double> * Improper::GetAngles() {
//	return &(_angles);
//}

void Improper::Clear() {
  atoms.clear();
  _angles.clear();
  _impnumber = 0;
}

void Improper::SetIndex(int impnumber) {
  _impnumber = impnumber;
}

void Improper::Card(std::string& s) {
  char buf[6];
  sprintf(buf, "%d", _impnumber);
  s.append(buf);

  MIAtom_iter atm;
  for (atm = atoms.begin(); atm != atoms.end(); ++atm) {
    s.append(" ");
    s.append((*atm)->name());
  }

  std::vector<double>::iterator ang;
  for (ang = _angles.begin(); ang != _angles.end(); ++ang) {
    sprintf(buf, " %5.1f", *ang);
    s.append(buf);
  }
}

bool Improper::operator==(const Improper& imp2) const {

  if (this->atoms.size() != imp2.atoms.size()) {
    return false;
  }

  //Return Equal(forward_compare) || Equal(reverse_compare)
  if (std::equal(this->atoms.begin(), this->atoms.end(), imp2.atoms.begin())) {
    return true;
  }
  if (std::equal(this->atoms.rbegin(), this->atoms.rend(), imp2.atoms.begin())) {
    return true;
  } else {
    return false;
  }
}

/////////////////////////////////////////////////////////////////////////////
// Function:    ConvertToDistance
// Purpose:		Convert the torsion angle to a 1-4 distance
// Input:
//
// Output:
// Requires:
/////////////////////////////////////////////////////////////////////////////
double Improper::ConvertToDistance(const ConstraintList& geom) const {
  if (atoms.size() < 4) {
    return -1;
  }

  //Retrieve the geometric data for the constituent bonds and angles
  double a = geom.GetBondLength(atoms[0], atoms[1]);
  double b = geom.GetBondLength(atoms[1], atoms[2]);
  double c = geom.GetBondLength(atoms[2], atoms[3]);
  double alpha = geom.GetAngleLength(atoms[0], atoms[1], atoms[2]);
  double beta = geom.GetAngleLength(atoms[1], atoms[2], atoms[3]);

  if (a < 0 || b < 0 || c < 0 || alpha < 0 || beta < 0) {
    return -1;
  }
  //Get the cosine and sine of the bond angles
  //(applies the law of cosines: a^2 = b^2 + c^2 - 2bcCos(alpha))
  double cos_alpha = (a * a +
                      b * b -
                      alpha * alpha) /
                     (2 * a * b);
  double cos_beta =  (b * b +
                      c * c -
                      beta * beta) /
                    (2 * b * c);
  double sin_alpha = sqrt(1 - cos_alpha * cos_alpha);
  double sin_beta = sqrt(1 - cos_beta * cos_beta);

  //Calculate the 1-4 distance, applying the formula:
  //d = (D13)^2 + (D24)^2 - b^2 + 2acCos(alpha)Cos(beta) - 2acSin(alpha)Sin(Beta)Cos(Omega)
  //Where omega is the dihedral angle and D13 and D24 are the 1-3 distances associated
  //with the bond angles alpha and beta
  return sqrt(alpha * alpha +
           beta * beta -
           b * b +
           2 * a * c * cos_alpha * cos_beta -
           2 * a * c * sin_alpha * sin_beta * cos(DEG2RAD * _angles.front()));
}

TORSION Improper::ToTORSION(const MIAtomList& new_atoms,
                            const std::vector<const MIAtom*>& old_atoms,
                            const std::vector<RESIDUE*>& new_residues,
                            const std::vector<const Residue*>& old_residues) {
  TORSION imp;

  std::string error_message;
  int i1 = GetIndex<const MIAtom*>(atoms[0], old_atoms);
  int i2 = GetIndex<const MIAtom*>(atoms[1], old_atoms);
  int i3 = GetIndex<const MIAtom*>(atoms[2], old_atoms);
  int i4 = GetIndex<const MIAtom*>(atoms[3], old_atoms);

  if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0) {
    error_message = "Corrupted atom ptr in impion.";
    throw error_message;
  }

  strcpy(imp.type, "IMP");
  sprintf(imp.type + 3, "%d", _impnumber);

  imp.setAtom1(new_atoms[i1]);
  imp.setAtom2(new_atoms[i2]);
  imp.atom3 = new_atoms[i3];
  imp.atom4 = new_atoms[i4];

  imp.nideal = _angles.size();

  if (imp.nideal > 3) {
    imp.nideal = 3;
  }

  if (imp.nideal < 0) {
    imp.nideal = 0;
  }

  for (int i = 0; i < imp.nideal; ++i) {
    imp.ideal[i] = (float)_angles[i];
  }

  int xRes = GetIndex<const Residue*>(_res, old_residues);

  if (xRes < 0) {
    error_message = "Corrupted residue information in torsion.";
    throw error_message;
  }

  imp.res = new_residues[xRes];

  return imp;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddBond
// Purpose:		Adds a bond constraint to the list of constraints
// Input:       Bond object, specifying the distance between two atoms
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddBond(BondLength& bond) {
  Bonds.push_back(bond);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddAngle
// Purpose:		Adds an angle constraint to the list of constraints
// Input:       Angle object, specifying a 1-3 distance
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddAngle(Angle& angle) {
  Angles.push_back(angle);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddImproper
// Purpose:		Adds an improper constraint to the list of constraints
// Input:       Improper object, specifiying one or more torsion angles for a
//				sequence of four atoms
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddImproper(Improper& improper) {
  Impropers.push_back(improper);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddPlane
// Purpose:		Adds a planarity constraint to the list of constraints
// Input:       Plane object specifying the list of coplanar atoms
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddPlane(Plane& plane) {
  Planes.push_back(plane);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddBump
// Purpose:		Adds a bump constraint to the list of constraints
// Input:       Bump object which specifies two atoms and their minimum distance
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddBump(Bump& bump) {
  Bumps.push_back(bump);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddChiral
// Purpose:		Adds a chiraliy constraint to the list of constraints
// Input:       Chiral object specifying the atoms and their chirality
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddChiral(Chiral& chiral) {
  Chirals.push_back(chiral);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddTorsion
// Purpose:		Adds a torsion to the list of constraints
// Input:       Torsion object
// Output:      None
// Requires:
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::AddTorsion(Torsion& torsion) {
  Torsions.push_back(torsion);
}

/////////////////////////////////////////////////////////////////////////////
// Function:    Clear
// Purpose:		Reinitialize the object
// Input:       None
// Output:      None
// Requires:	Assumes the memory for atoms in the PLANE structs has been
//				allocated using operator new
/////////////////////////////////////////////////////////////////////////////
void ConstraintList::Clear() {
  Angles.clear();
  Bonds.clear();
  Impropers.clear();
  Planes.clear();
  Bumps.clear();
}

float ConstraintList::GetBondLength(const MIAtom* a1, const MIAtom* a2) const {
  std::vector<BondLength>::const_iterator bnd, end = Bonds.end();

  bnd = Bonds.begin();
  while (bnd != end) {
    if ((bnd->getAtom1() == a1 && bnd->getAtom2() == a2)
        || (bnd->getAtom2() == a1 && bnd->getAtom1() == a2)) {
      return bnd->ideal_dist;
    }
    ++bnd;
  }
  return -1.0F;
}

float ConstraintList::GetAngleLength(const MIAtom* a1, const MIAtom* a2, const MIAtom* a3) const {
  std::vector<Angle>::const_iterator ang, end = Angles.end();

  ang = Angles.begin();
  while (ang != end) {
    if ((ang->getAtom1() == a1 && ang->getAtom2() == a2 && ang->atom3 == a3)
        || (ang->atom3 == a1 && ang->getAtom2() == a2 && ang->getAtom1() == a3)) {
      return ang->ideal_angle;
    }
    ++ang;
  }
  return -1.0F;
}

}
