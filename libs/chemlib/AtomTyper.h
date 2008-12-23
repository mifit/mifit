#ifndef ATOM_TYPER_BASE_CLASS_H
#define ATOM_TYPER_BASE_CLASS_H

#include "Ligand.h"
#include "model.h"
#include "Bond.h"

//This may seem a bit clumsy: clients must first construct an AtomTyper using a
//RESIDUE and a vector of BondS, then they get the atom types out one-by-one.
//More streamlined would be to have one ftn that writes the type info into the
//MIAtom structs themselves or returns a list of atom types to the client.  But
//adding fields to the MIAtom struct has obvious drawbacks, particularly since
//we may eventually have many atom typing systems.  And returning a list seems
//kind of inflexible...

//The most obvious drawback is that the MIAtom names are again presumed to be
//unique identifiers within a residue.  This assumption is present in many parts
//of MIFit, though. We might consider having MIAtom contain a pointer
//to an MIAtom struct, so that we have a built-in cross-reference when we've
//created an Ligand that corresponds to an MIFit molecule. -KWB


namespace chemlib {

class AtomTyper {
public:
  virtual ~AtomTyper() {
  }

  virtual char* AtomType(const MIAtom* atom) const = 0;
  AtomTyper(const RESIDUE& res, const std::vector<Bond>& bonds);
protected:
  Ligand m_mol;
};
}

#endif //ATOM_TYPER_BASE_CLASS_H
