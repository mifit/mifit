#ifndef mifit_figurelib_HydrogenBond_h
#define mifit_figurelib_HydrogenBond_h

#include "Shape.h"

namespace moldraw {

class HydrogenBond : Shape {
public:
  HydrogenBond(Drawing* dp,
               float x1,
               float y1,
               float x2,
               float y2,
               float width,
               float font_size,
               float dist);
  virtual ~HydrogenBond();
  virtual void Draw();
private:
  float _x1;
  float _y1;
  float _x2;
  float _y2;
  float _width;
  float _font_size;
  float _dist;
};

}

#endif
