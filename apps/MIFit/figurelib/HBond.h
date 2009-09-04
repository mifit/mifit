#ifndef HBOND_H
#define HBOND_H

#include <chemlib/chemlib.h>

#include <vector>
#include <string>

#define MAX_HBOND_LENGTH 3.3F
#define MIN_HBOND_LENGTH 2.1F

namespace HBondLookup {
bool IsKnown(const chemlib::Residue& res);
bool IsKnownAcceptor(const chemlib::MIAtom& atom, const chemlib::Residue& res);
bool IsKnownDonor(const chemlib::MIAtom& atom, const chemlib::Residue& res);
}



namespace moldraw {

class HBond {
public:
  chemlib::MIAtom* donor;
  chemlib::MIAtom* acceptor;
  chemlib::Residue* donor_res;
  chemlib::Residue* acceptor_res;
  double distance;

  std::string Print() const;
};

inline std::string AppendHBond(std::string stringSoFar, const HBond& bond) {
  return stringSoFar + bond.Print();
}

void GetDonors(const std::vector<chemlib::Residue*>& residues, std::vector<chemlib::MIAtom*>& donor_atoms);
void GetAcceptors(const std::vector<chemlib::Residue*>& residues, std::vector<chemlib::MIAtom*>& acceptor_atoms);

bool CheckHBond(chemlib::MIAtom* donor, chemlib::MIAtom* acceptor, HBond& hbond);

inline bool IsDonor(const chemlib::MIAtom& atom, const chemlib::Residue& res) {

  if (HBondLookup::IsKnown(res)) {
    return HBondLookup::IsKnownDonor(atom, res);
  } else if ((atom.atomicnumber() == 7
              || atom.atomicnumber() == 8)
             && atom.hcount() > 0) {
    return true;
  } else {
    return false;
  }
}

inline bool IsAcceptor(const chemlib::MIAtom& atom, const chemlib::Residue& res) {
  if (HBondLookup::IsKnown(res)) {
    return HBondLookup::IsKnownAcceptor(atom, res);
  } else if (atom.atomicnumber() == 8) {
    return true;
  } else if (atom.atomicnumber() == 7
             && atom.nabors().size() == 2
             && atom.hybrid() == 2
             && atom.hcount() == 0) {
    return true;
  } else {
    return false;
  }
}

} //namespace moldraw
#endif //HBOND_H
