#ifndef MI_ANGLE_H
#define MI_ANGLE_H

#include "model.h"

namespace chemlib {
class MIAtom;
class RESIDUE;

class ANGLE {
public:
  ANGLE() : atom1(0), atom2(0), atom3(0), ideal_angle(0.0f),
    tolerance(0.0f), res(0), iscyclic(0), isaromatic(0),
    ring_system(0), smallest_ring_size(0),
    smallest_aromatic_ring(0) {
  }

  void Clear() {
    ANGLE a; *this = a;
  }                                    // note: doesn't free old mem, if any

  MIAtom* atom1;
  MIAtom* atom2;
  MIAtom* atom3;
  float ideal_angle;
  float tolerance;
  /**
   * the residue the angle is from
   */
  const RESIDUE* res;
  /**
   * True if all three atoms are in the same ring system
   */
  int iscyclic;
  /**
   * True if all three atoms are in the same aromatic system
   */
  int isaromatic;
  /**
   * Index of the ring system, -1 for acyclic angles
   */
  int ring_system;
  /**
   * Size of smallest ring in which the angle is contained
   */
  int smallest_ring_size;
  /**
   * Size of smallest aromatic ring which has the angle
   */
  int smallest_aromatic_ring;

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

  MIAtom* getAtom3() const {
    return atom3;
  }

  void setAtom3(MIAtom* atom) {
    atom3 = atom;
  }

};


}

#endif
