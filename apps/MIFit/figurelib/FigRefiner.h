#include <chemlib/chemlib.h>

#include "HBond.h"

#define SQ_ATOM_ATOM_THRESHOLD 6.25
#define SQ_BOND_ATOM_THRESHOLD 0.49
#define SQ_BOND_BOND_THRESHOLD 0.25
#define BOND_CLASH_WEIGHT 0.05
#define BOND_INTERSECT_PENALTY 10
#define SMALL_DISTANCE 0.01
#define MAX_FLIP_CYCLES 100

namespace moldraw {
void FigOpt(chemlib::Ligand& site, const std::vector<HBond>& hbonds);
double CalcFigFitness(const chemlib::Ligand& site, const std::vector<HBond>& hbonds);

double AtomClashTotal(const chemlib::Ligand& site);
double BondClashTotal(const chemlib::Ligand& site, const std::vector<HBond>& hbonds);

double AtomAtomClashScore(const chemlib::MIAtom& atom1, const chemlib::MIAtom& atom2);
double BondAtomClashScore(const HBond& bond, const chemlib::MIAtom& atom);
double BondBondClashScore(const HBond& bond1, const HBond& bond2);

inline double BondClashScore(double distance) {
  if (distance > SQ_BOND_ATOM_THRESHOLD) {
    return 0;
  } else if (distance > SMALL_DISTANCE) {
    return 1 / distance;
  } else {
    return 1 / SMALL_DISTANCE;              //Clamp to avoid infinite scores
  }
}

} //namespaece moldraw
