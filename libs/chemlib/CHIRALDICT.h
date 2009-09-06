#ifndef MI_CHIRAL_H
#define MI_CHIRAL_H

#include <cstring>
#include "model.h"

namespace chemlib {

class MIAtom;

#define CHIRAL_AUTO    0x0
#define CHIRAL_USERADD 0x1
#define CHIRAL_DELETED 0x2

class CHIRAL {
public:
  MIAtom* center;
  MIAtom* atom1;
  MIAtom* atom2;
  MIAtom* atom3;
  int order;
  /**
   * 0x0 = automatically generated
   * 0x1 = manually created
   * 0x2 = manually deleted and hence may NOT be regenerated
   */
  short flags;

  CHIRAL() : center(NULL), atom1(NULL), atom2(NULL), atom3(NULL), order(0), flags(0) {
  }

  void Clear() {
    CHIRAL c; *this = c;
  }                                     // doesn't clear memory

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

class CHIRALDICT {
public:
  CHIRALDICT() {
    memset(this, 0, sizeof(CHIRALDICT)); order = -1;
  }

  void Clear() {
    memset(this, 0, sizeof(CHIRALDICT)); order = -1;
  }

  bool AddAtom(const char* atom_name, int index) {
    if (index < 0 || index >= 3) {
      return false;
    } else {
      strncpy(name[index], atom_name, MAXATOMNAME);
      return true;
    }
  }

  char restype[MAXNAME];
  char center[MAXATOMNAME];
  char name[4][MAXATOMNAME];
  int order;
};

}

#endif
