#ifndef SUBSTITUENT_CLASS_H
#define SUBSTITUENT_CLASS_H

#include <vector>

#include "MIAtom_fwd.h"

namespace chemlib {

class Ligand;

class Substituent {
public:
  void Clear();
  void SetOrigin(MIAtom* origin) {
    _origin = origin;
  }

  void SetMolecule(const Ligand* mol) {
    _lig = mol;
  };
  void AddBranch(MIAtom* root);

  void Squash();
  void ZAxisRotate(double theta);
  std::vector<double> _direction;
private:
  MIAtom* _origin;
  MIAtomList _first_shell_atoms;
  const Ligand* _lig;
  MIAtomList _atoms;
  bool _issimple;
  bool _iscyclic;
  double _theta;

  void UpdateFirstShell();
  void UpdateFlags();
  void CalcDirection();
  void CalcTheta();

  void SingleSquash();
  void DoubleSquash();
  void CycleSquash();

  friend struct LeastTheta;             //Allow sorting on angle theta
};
} //namespace chemlib

#endif //SUBSTITUENT_H
