#ifndef LIGDICTIONARY_H
#define LIGDICTIONARY_H

#include "sequence_util.h"
#include "Chiral.h"
#include "Constraint.h"

#include <vector>
#include <algorithm>

namespace chemlib {

class MIAtom;
class Bond;
class Ligand;
class Residue;

class LigDictionary {
public:
  LigDictionary(Ligand*, const Residue*);
  void AnalyzeRings();
  void AddPlane(Plane&);
  void AddTorsion(Torsion& tor);
  void AddImproper(Improper& imp, bool report = true);
  //		void GenRingPlane(Cycle &);
  //		void GenRingImpropers(Cycle &);
  void GenDoubleBondImpropers();
  void FindAcycPlanes();
  void GenTorsions();
  void GenChirals();
  void AddChiral(const MIAtom *, int, int, int, int);
  //		double ImproperAngle(MIAtom *atom1, MIAtom *atom4, Cycle &rng);
  //		void Print();							//Just a debugging tool--prints info to cout
  void SetMolGeometry();
  //		void SetOneFours(Cycle &);
  bool CheckConsistency(const std::vector<int>&);
  float GetSigmaPlane() {
    return _sigmaplane;
  };
  const Residue* GetResidue() {
    return _res;
  };
  void UnspecDoubleBond(Bond*);
protected:
  Ligand* lig;
  const Residue* _res;
  std::vector<Plane> planes;
  std::vector<Torsion> torsions;
  std::vector<Improper> impropers;
  std::vector<Chiral> chirals;
  float _sigmaplane;
};

} //namespace chemlib


#endif //LIGDICTIONARY_H
