#ifndef LIG_REFINER_H
#define LIG_REFINER_H

#include <chemlib/chemlib.h>

#include <vector>

namespace conflib {

class LigRefiner {
public:
  std::vector<chemlib::BondLength> RefiBonds;
  std::vector<chemlib::Angle> RefiAngles;
  std::vector<chemlib::Plane> RefiPlanes;
  std::vector<chemlib::Improper> RefiTorsions;
  std::vector<chemlib::Bump> RefiBumps;
  std::vector<chemlib::Residue*> RefiRes;

  float BondWeight, AngleWeight, PlaneWeight, TorsionWeight, BumpWeight;

  chemlib::Ligand* CurrentModel;

  int nCycles;
  bool RefiVerbose;

  LigRefiner();

  double Refine();
  double minimize_bonds(chemlib::BondLength*, unsigned int);
  double minimize_angles(chemlib::Angle*, unsigned int);
  double minimize_planes(chemlib::Plane*, unsigned int);
  double minimize_bumps(chemlib::Bump* bumps, unsigned int nbumps);
  double minimize_torsions();
  int resetderivatives();
  int applyderivatives();

  void AddRefiRes(chemlib::Residue* res);
  long SetRefiRes(std::vector<chemlib::Residue*>::iterator res1,
                  std::vector<chemlib::Residue*>::iterator res2);
  void SetModel(chemlib::Ligand* model);
  int BuildBumps();
  void Clear();

  int GetNumberCycles() {
    return nCycles;
  }

  void SetNumberCycles(int n) {
    nCycles = n;
  }

};

} //namespace conflib

#endif //LIG_REFINER_H
