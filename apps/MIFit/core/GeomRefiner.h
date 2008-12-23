#ifndef MI_GEOM_REFINER_H
#define MI_GEOM_REFINER_H

#include "moloptlib.h"

class Molecule;

class GeomRefiner : public MIMolOpt {
public:
  bool EditEntry(const char* type);
  void EditEntryCleanup(bool);
  unsigned long FindGeomErrors(Molecule* model, float error_threshold);
};

// get DictResList converted from vector of std::string to vector of std::string
std::vector<std::string> GetDictResList(GeomRefiner* gr);



#endif
