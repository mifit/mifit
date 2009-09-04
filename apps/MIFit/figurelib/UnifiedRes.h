#ifndef UNIRES_H
#define UNIRES_H

#include <chemlib/chemlib.h>
#include "Drawing.h"

#include <vector>

#define FUZZ_SPAN 60

namespace moldraw {

class UnifiedRes {
public:
  std::string type;
  std::string name;
  unsigned short chain_id;
  float x;
  float y;
  std::vector<chemlib::MIAtom*> partners;               //Hydrophobic interaction partners

  void CalcDirection(double* direction, double* extent);        //direction in radians
  void Draw(Drawing* dp);
};
} //namespace moldraw

#endif //UNIRES_H
