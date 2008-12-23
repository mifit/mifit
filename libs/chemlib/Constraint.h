#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <vector>
#include <fstream>
#include <ostream>

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue.h"
#include "Chiral.h"
#include "ANGLE.h"
#include "PLANE.h"
#include "TORSION.h"

namespace chemlib {

class ConstraintList;
class LigDictionary;

class BondLength {
public:
  MIAtom* atom1;
  MIAtom* atom2;
  float ideal_dist;
  float tolerance;
  int dict_include;


  MIAtom* getAtom1() const {
    return atom1;
  }

  void setAtom1(MIAtom* atom) {
    atom1 = atom;
  }

  MIAtom* getAtom2() const {
    return atom2;
  }

  void setAtom2(MIAtom* atom) {
    atom2 = atom;
  }

private:
  friend class Ligand;
  Bond ToBond(const MIAtomList&, const std::vector<const MIAtom*>&);
};

class Angle {
public:
  MIAtom* atom1;
  MIAtom* atom2;
  MIAtom* atom3;
  float ideal_angle;
  float tolerance;
  const Residue* res;                     //The residue the bond is from
  int iscyclic;                     //True if all three atoms are in the same ring system
  int isaromatic;                   //True if all three atoms are in the same aromatic system
  int ring_system;                      //Index of the ring system, -1 for acyclic angles
  int smallest_ring_size;               //Size of smallest ring in which the angle is contained
  int smallest_aromatic_ring;           //Size of smallest aromatic ring which has the angle

  MIAtom* getAtom1() const {
    return atom1;
  }

  void setAtom1(MIAtom* atom) {
    atom1 = atom;
  }

  MIAtom* getAtom2() const {
    return atom2;
  }

  void setAtom2(MIAtom* atom) {
    atom2 = atom;
  }

  MIAtom* getAtom3() const {
    return atom3;
  }

  void setAtom3(MIAtom* atom) {
    atom3 = atom;
  }

private:
  friend class Ligand;
  ANGLE ToANGLE(const MIAtomList&, const std::vector<const MIAtom*>&);
};

class Bump {
public:
  MIAtom* atom1;
  MIAtom* atom2;
  float min_d;
  float tolerance;

  MIAtom* getAtom1() const {
    return atom1;
  }

  void setAtom1(MIAtom* atom) {
    atom1 = atom;
  }

  MIAtom* getAtom2() const {
    return atom2;
  }

  void setAtom2(MIAtom* atom) {
    atom2 = atom;
  }

};

class Plane {
public:
  Plane();
  Plane(float);
  Plane(LigDictionary&);
  void AddAtom(MIAtom*);
  MIAtom* GetAtom(int);
  int NumAtoms();
  //		void Print();
  void Clear();
  void Card(std::string&);
  float GetTolerance() {
    return _tol_set ? _tolerance : 0;
  }

  void SetTolerance(float t) {
    _tolerance = t; _tol_set = true;
  }

  const Residue* GetResidue() {
    return _res;
  }

  void SetResidue(const Residue* r) {
    _res = r;
  }

  bool operator==(const Plane&) const;

protected:
  MIAtomList atoms;
  float _tolerance;
  bool _tol_set;
  const Residue* _res;

private:
  friend class Ligand;
  PLANE ToPLANE(const MIAtomList&, const std::vector<const MIAtom*>&,
                const std::vector<RESIDUE*>&, const std::vector<const Residue*>&);
};


class Torsion {
public:
  void AddAtom(MIAtom*);
  //		void Print();
  void Clear();
  void SetIndex(int);
  void Card(std::string&);

  const Residue* GetResidue() {
    return _res;
  }

  void SetResidue(const Residue* r) {
    _res = r;
  }

  bool operator==(const Torsion&) const;
protected:
  MIAtomList atoms;
  int _torsnumber;
  const Residue* _res;

private:
  friend class Ligand;
  TORSION ToTORSION(const MIAtomList&, const std::vector<const MIAtom*>&,
                    const std::vector<RESIDUE*>&, const std::vector<const Residue*>&);
};

class Improper {
public:
  void AddAtom(MIAtom*);
  MIAtom* GetAtom(int);
  const MIAtom* GetAtom(int) const;
  double GetAngle(int i) const {
    return _angles[i];
  }

  int NumAngles() {
    return _angles.size();
  }

  void AddAngle(double);
  //const std::vector<double> * GetAngles();
  void ReInit(const MIAtomList &, double);
  void ReInit(const MIAtomList&,
              const std::vector<double>&);
  //		void Print();
  void Clear();
  void SetIndex(int);
  void Card(std::string&);
  void IncludeInDict() {
    _dict_include = true;
  }

  void ExcludeFromDict() {
    _dict_include = false;
  }

  bool DictTest() {
    return _dict_include;
  }

  bool operator==(const Improper& imp2) const;
  double ConvertToDistance(const ConstraintList& geom) const;

  const Residue* GetResidue() {
    return _res;
  }

  void SetResidue(const Residue* r) {
    _res = r;
  }

protected:
  MIAtomList atoms;
  std::vector<double> _angles;                      //The allowable dihedral angle(s) (in degrees)
  int _impnumber;
  bool _dict_include;                           //False for imps used only in smiles2pdb, not in dict
  const Residue* _res;

private:
  friend class Ligand;
  TORSION ToTORSION(const MIAtomList&, const std::vector<const MIAtom*>&,
                    const std::vector<RESIDUE*>&, const std::vector<const Residue*>&);

};

class ConstraintList {
public:
  std::vector<BondLength> Bonds;
  std::vector<Angle> Angles;
  std::vector<Plane> Planes;
  std::vector<Improper> Impropers;
  std::vector<Bump> Bumps;
  std::vector<Chiral> Chirals;
  std::vector<Torsion> Torsions;
  void AddBond(BondLength&);
  void AddAngle(Angle&);
  void AddPlane(Plane&);
  void AddImproper(Improper&);
  void AddChiral(Chiral&);
  void AddBump(Bump&);
  void AddTorsion(Torsion&);
  void Clear();

  float GetBondLength(const MIAtom* a1, const MIAtom* a2) const;
  float GetAngleLength(const MIAtom* a1, const MIAtom* a2, const MIAtom* a3) const;
};

} //namespace chemlib

#endif //CONSTRAINT_H
