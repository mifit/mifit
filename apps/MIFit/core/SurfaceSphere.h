#ifndef mifit_model_SurfaceSphere_h
#define mifit_model_SurfaceSphere_h

#include <vector>
#include <math/mathlib.h>

class SurfaceSphere {

  float spacing;

  float radius;

  std::vector<APOINT> points;

public:

  SurfaceSphere();

  void build(float radius, float spacing);

  float getRadius();

  float getSpacing();

  std::vector<APOINT>& getPoints();

  void clearPoints();

};

#endif
