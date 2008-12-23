#ifndef COORD_GENERATOR_H
#define COORD_GENERATOR_H

#include "chemlib.h"

namespace conflib {
class sdgDistance;

void PositionAtom(chemlib::MIAtom* atom,
                  const chemlib::MIAtom* root1,
                  const chemlib::MIAtom* root2,
                  const chemlib::MIAtom* root3,
                  float length,
                  float angle,
                  float torsion);

void PositionAtom(chemlib::MIAtom* atom,
                  const chemlib::MIAtom* root1,
                  const chemlib::MIAtom* root2,
                  float length,
                  float angle);

class CoordGenerator {
public:
  CoordGenerator(chemlib::Ligand* mol, chemlib::Residue* res) {
    _lig = mol; _res = res;
  }

  double sdgRefine();
  void sdgGenMolStruct();
  double GenMolStruct(std::string& log);
  void PlaceAtom(chemlib::MIAtom*);
  void SetExtended(chemlib::MIAtom*);
  //	void OptimizeAtom(chemlib::MIAtom *);
  //	void OptimizeMolecule();
  //	void WritePDB(ofstream &);
private:
  void GatherAtomConstraints(chemlib::MIAtom*);
  void GetRootAtom(chemlib::MIAtom*);
  void GetRootVector();
  int GetRootNeighbor();
  void SumRootVectors();
  void CalcInitialCoordinates(chemlib::MIAtom*);
  int CountPlaced();
  int CountPlaced(chemlib::Plane*);
  void PartialPlane(chemlib::Plane*, chemlib::Plane*);
  chemlib::Ligand* _lig;
  chemlib::Residue* _res;
  chemlib::ConstraintList _master;
  chemlib::ConstraintList _current;
  int _nplaced;                     //Number of atoms placed thus far
  chemlib::MIAtom* _root;                       //Current atom will be attached to this atom
  chemlib::MIAtom* _rnabor;                 //A placed nabor of the root atom
  double _rvect[3];                 //Bond vector from _rnabor to _root
  int _nCycles;                     //Number of Cycles to refine atoms
  float _newbond[3];                //Unit vector for the next bond to be added
};

double GetTargetVolume(chemlib::Chiral& chiral, std::vector<sdgDistance>& dists);
sdgDistance FindDistance(chemlib::MIAtom* a1, chemlib::MIAtom* a2, std::vector<sdgDistance>& dists);

} //namespace conflib

#endif //COORD_GENERATOR_H
