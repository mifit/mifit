#ifndef mifit_legacy_MIMolOPt_h
#define mifit_legacy_MIMolOPt_h

#include <QObject>
#include <set>
#include <vector>
#include <list>
#include <map>

#include <math/mathlib.h>
#include <chemlib/chemlib.h>


class EMapBase;
class Molecule;
class InterpBox;

// called during ligand optimization and full optimization
class MIMolOptCheckPoint {
public:
  virtual ~MIMolOptCheckPoint() {
  }

  virtual bool operator()(chemlib::MIMoleculeBase*) = 0;
};


class MIMolOpt : public QObject {
    Q_OBJECT
public:
  MIMolOpt();
  virtual ~MIMolOpt();

  bool BuildMainchain(chemlib::RESIDUE* res, chemlib::MIMoleculeBase* model, EMapBase* emap, bool addAtomsToNextResidue = true);
  void RigidOptimize(std::vector<chemlib::MIAtom*>&, chemlib::MIMoleculeBase* fitmol, EMapBase* emap);
  void TorsionOptimize(std::vector<chemlib::MIAtom*>& CurrentAtoms, chemlib::MIMoleculeBase* fitmol, EMapBase* emap,
                       std::vector<chemlib::TORSION>& torsions, bool do_setup = true);
  void FullOptimize(std::vector<chemlib::MIAtom*>& CurrentAtoms, chemlib::MIMoleculeBase* fitmol, EMapBase* emap,
                    const float* center, InterpBox& box, unsigned int refine_level,
                    MIMolOptCheckPoint* checkpoint = 0);
  void LigandOptimize(std::vector<chemlib::MIAtom*>& CurrentAtoms, chemlib::MIMoleculeBase* fitmol, EMapBase* emap,
                      const float* center, InterpBox& box, unsigned int refine_level, chemlib::GeomSaver& conformations,
                      MIMolOptCheckPoint* checkpoint = 0);
  void MolecularReplace(chemlib::MIMoleculeBase* fitmol, EMapBase* emap);
  void RefiAllTorsions(chemlib::RESIDUE* reslist);

  void Purge(chemlib::MIMoleculeBase*);
  void Purge(chemlib::RESIDUE*);
  void Purge(chemlib::MIAtom*);
  void Purge(EMapBase*);

  void Do();
  bool CanUndo() {
    return SaveToken!=0 && geomsaver.NumberSets() > 1;
  }

  bool Undo();

  bool Redo();
  bool CanRedo() {
    return (int)SaveToken < geomsaver.NumberSets()-1;
  }

  void Reset();
  void Cancel();
  void Accept();
  void Refine();

  long SetRefiRes(chemlib::RESIDUE* res1, chemlib::RESIDUE* res2, chemlib::MIMoleculeBase* model, EMapBase* emap = NULL);

  void lockRefineTarget();
  void unlockRefineTarget();

  // checks to see is an atom is in those being currently refined
  bool IsBeingRefined(chemlib::MIAtom* atom);
  bool IsRefining() {
    return (RefiRes != NULL);
  }

  chemlib::MIMoleculeBase* GetCurrentModel() {
    return CurrentModel;
  }

  // not these
  void SetMap(EMapBase* map) {
    CurrentMap = map;
  }

  EMapBase* GetMap() {
    return CurrentMap;
  }

  int GetBondWeightI() {
    return ROUND(10.0F*BondWeight);
  }

  void SetBondWeightI(int v) {
    BondWeight = (float)v/10.0F;
  }

  int GetBumpWeightI() {
    return ROUND(10.0F*BumpWeight);
  }

  void SetBumpWeightI(int v) {
    BumpWeight = (float)v/10.0F;
  }

  int GetPlaneWeightI() {
    return ROUND(10.0F*PlaneWeight);
  }

  void SetPlaneWeightI(int v) {
    PlaneWeight = (float)v/10.0F;
  }

  int GetAngleWeightI() {
    return ROUND(10.0F*AngleWeight);
  }

  void SetAngleWeightI(int v) {
    AngleWeight = (float)v/10.0F;
  }

  int GetTorsionWeightI() {
    return ROUND(10.0F*TorsionWeight);
  }

  void SetTorsionWeightI(int v) {
    TorsionWeight = (float)v/10.0F;
  }

  void SetMapWeightI(int v) {
    MapWeight = (float)v/10.0F;
  }

  int GetMapWeightI() {
    return ROUND(10.0F*MapWeight);
  }

  int GetNumberCycles() {
    return nCycles;
  }

  void SetNumberCycles(int n) {
    nCycles = n;
  }

  bool GetRefineWhileFit() {
    return fit_while_refine;
  }

  void SetRefineWhileFit(bool v) {
    fit_while_refine = v;
  }

  bool GetVerbose() {
    return RefiVerbose;
  }

  void SetVerbose(bool v) {
    RefiVerbose = v;
  }

  unsigned int SaveToken;

  chemlib::MIMolDictionary dict;

Q_SIGNALS:
  void isRefiningChanged(bool);

protected:
  // these belong here
  int FindNeighbours(std::vector<chemlib::MIAtom*>& CurrentAtoms, std::vector<chemlib::MIAtom*>& Neighbours, chemlib::MIMoleculeBase* fitmol, float distance);
  EMapBase* CurrentMap;
  chemlib::MIMoleculeBase* CurrentModel;
  bool RefiVerbose;
  bool fit_while_refine;
  chemlib::GeomSaver geomsaver;
  chemlib::RESIDUE* RefiRes;
  chemlib::RESIDUE* ResActiveModel;
  bool refineTargetLocked;
  
  float BondWeight, AngleWeight, PlaneWeight, MapWeight, TorsionWeight, BumpWeight;
  int AutoFit;
  int nCycles;
  int nRefiRes;
  int nucleic;

  void internalSetRefiRes(chemlib::RESIDUE* residue, int nResidues);

  int getbonddist(chemlib::RESIDUE* res, chemlib::Bond* bond);

  int minimize_bonds(std::vector<chemlib::Bond>&, unsigned int);
  int minimize_angles(std::vector<chemlib::ANGLE>&, unsigned int);
  int minimize_planes(std::vector<chemlib::PLANE>&, unsigned int);
  int minimize_constraints(std::vector<chemlib::Bond>&, unsigned int nbonds);
  int minimize_bumps(std::vector<chemlib::Bond>&, unsigned int nbumps);
  int minimize_map();
  int minimize_phipsi();
  int minimize_torsions();

  int takestep(int seed);
  float scorestep(chemlib::Bond*, chemlib::ANGLE*, chemlib::PLANE*, chemlib::Bond*, chemlib::Bond*, std::vector<chemlib::MIAtom*>, int count);

  int resetatomderivatives(chemlib::MIAtom * atoms[], int natoms);
  int resetderivatives();
  int applyderivatives();

  void clearRefineTarget();

  float phipsi_energy(float phi, float psi);
  float StdevTorsions();
  float StdevPlanes();
  float StdevAngles();
  float StdevBonds();

public Q_SLOTS:
  void ConnectTo(chemlib::MIMoleculeBase *mol);
  void atomsToBeDeleted(chemlib::MIMoleculeBase* model, const std::vector<chemlib::MIAtom*> &atoms);
  void residuesToBeDeleted(chemlib::MIMoleculeBase* model, std::vector<chemlib::RESIDUE*> &res);
  void moleculeToBeDeleted(chemlib::MIMoleculeBase* molecule);
  void moleculeDeleted(chemlib::MIMoleculeBase* molecule);
};


//@{
// Returns the ideal length of a bond given two atoms.
//@}
float bond_dist(chemlib::MIAtom* atom1, chemlib::MIAtom* atom2);
//@{
// Returns the ideal bond radius given an atom.
//@}
float bond_radius(chemlib::MIAtom* atom);

struct TorsionMatch : public std::binary_function<chemlib::TORSION, chemlib::TORSDICT, bool> {
  bool operator()(const chemlib::TORSION& t, const chemlib::TORSDICT& td) const;
};

void lsqplane(chemlib::PLANE& plane);

#endif

