#ifndef REFMAC_DICT_H
#define REFMAC_DICT_H

#include "AtomTyper.h"

#include <vector>

namespace chemlib {

class RefmacAtomTyper : public AtomTyper {
public:
  RefmacAtomTyper(const RESIDUE& res, const std::vector<Bond>& bonds);

  virtual char* AtomType(const MIAtom* atom) const;

private:
  unsigned int Index(const char* name) const;

  char* Name(unsigned int index) const;

  char* TypeHydrogen(const MIAtom* atom) const;
  char* TypeCarbon(const MIAtom* atom) const;
  char* TypeNitrogen(const MIAtom* atom) const;
  char* TypeOxygen(const MIAtom* atom) const;
  char* TypeSilicon(const MIAtom* atom) const;
  char* TypePhosphorus(const MIAtom* atom) const;
  char* TypeSulfur(const MIAtom* atom) const;
  char* TypeGermanium(const MIAtom* atom) const;
  char* TypeArsenic(const MIAtom* atom) const;
  char* TypeOther(const MIAtom* atom) const;

};
}



#endif //REFMAC_DICT_H
