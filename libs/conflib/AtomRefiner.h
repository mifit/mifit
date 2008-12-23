#ifndef ATOM_REFINER_H
#define ATOM_REFINER_H

#include "chemlib.h"

namespace conflib {
class CoordGenerator;

class AtomRefiner {
public:
  AtomRefiner(CoordGenerator*);         //Constructor
  virtual ~AtomRefiner();                           //Destructor
  void SetNumberCycles(int n) {
    _nCycles = n;
  }

  int GetNumberCycles() {
    return _nCycles;
  }

  void RefineAtom(chemlib::MIAtom*, chemlib::MIAtom*, chemlib::MIAtom*);
  void SetConstraintList(chemlib::ConstraintList* cl) {
    cons = cl;
  }

private:
  void Min1D(chemlib::MIAtom*, float*, float&);
  void Bracket3DMin(chemlib::MIAtom*, float*, float&, float*, float*);
  void MoveAtom(chemlib::MIAtom*, float*);
  float ScoreAtom(chemlib::MIAtom*);
  float minimize_bonds(chemlib::MIAtom*);
  float minimize_angles(chemlib::MIAtom*);
  float minimize_planes(chemlib::MIAtom*);
  float minimize_torsions(chemlib::MIAtom*);
  float ScoreChirals(chemlib::MIAtom*);
  float ScoreBumps(chemlib::MIAtom*);
  float BondWeight, AngleWeight, PlaneWeight, TorsionWeight;
  int _nCycles;
  int _nTwists;
  chemlib::ConstraintList* cons;
  //		Ligand *_mol;
  CoordGenerator* _xyzGen;
};
} //namespace conflib

#endif //ATOM_REFINER_H
