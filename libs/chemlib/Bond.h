#ifndef BOND_H
#define BOND_H

#include "MIAtom_fwd.h"
#include "model.h"

#include <vector>
#include <functional>

namespace chemlib {

class Bond {
private:
  unsigned char bondOrder;
  MIAtom* atom1;
  MIAtom* atom2;

public:
  
  static void copyBondOrders(const std::vector<Bond>& fromBonds, std::vector<Bond>& toBonds);
  
  Bond();                                               //Default constructor

  Bond(const Bond& bond,                                 //Construction from an existing bond
       const MIAtomList& old_atoms,          //(Need vectors to convert atm ptrs)
       const MIAtomList& new_atoms);

  Bond(MIAtom* in_atom1,                                //Construction with explicit parameters
       MIAtom* in_atom2,
       unsigned char in_bondOrder,
       char in_stereo);
  void Clear() {
    Bond b; *this = b;
  }                                   // doesn't free/delete atom pointers

  MIAtom* getAtom1() const {
    return atom1;
  }

  void setAtom1(MIAtom* atom) {
    atom1 = atom;
  }

  MIAtom* getAtom2() const {
    return atom2;
  }

  void setAtom2(MIAtom* atom) {
    atom2 = atom;
  }

  unsigned char getOrder() const {
    return bondOrder;
  }

  void setOrder(unsigned char order) {
    this->bondOrder = order;
  }

  void ClearOrder();
  bool UpdateAtoms(const MIAtomList& old_atoms,
                   const MIAtomList& new_atoms);


  inline void SetSingle() {
    bondOrder = SINGLEBOND;
  }

  inline void SetDouble() {
    bondOrder = DOUBLEBOND;
  }

  inline void SetTriple() {
    bondOrder = TRIPLEBOND;
  }

  inline void SetPartialDouble() {
    bondOrder = PARTIALDOUBLEBOND;
  }

  bool IsRotatable() const;

  inline void AssignIdealLength(float length) {
    ideal_length = length;
  }

  std::string Print() const;

  char type;
  char stereo;
  float ideal_length;
  float tolerance;
  int dict_include;
  int isaromatic;
  int iscyclic;
  int ring_system;                      //Index of the ring system, -1 for acyclic bonds
  int smallest_ring_size;               //Size of smallest ring in which the bond is contained
  int smallest_aromatic_ring;           //Size of smallest aromatic ring in which the bond is contained

private:
};

struct IsAromatic : public std::unary_function<Bond, bool> {
  bool operator()(const Bond& bond) const;
};

inline bool IsAromatic::operator ()(const Bond& bond) const {
  return bond.isaromatic != 0;
}

}   //namespace chemlib
#endif //BOND_H
