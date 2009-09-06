#ifndef mifit_figurelib_SingleBond_h
#define mifit_figurelib_SingleBond_h

#include "Shape.h"

namespace moldraw {

class SingleBond : Shape {
public:
  SingleBond(Drawing* dp,
             float x1,
             float y1,
             float x2,
             float y2,
             float width);
  virtual ~SingleBond();
  virtual void Draw();
private:
  float _x1;
  float _y1;
  float _x2;
  float _y2;
  float _width;
};

}

#endif
