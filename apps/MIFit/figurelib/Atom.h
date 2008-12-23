#ifndef mifit_figurelib_Atom_h
#define mifit_figurelib_Atom_h

#include "corelib.h"
#include "Shape.h"

namespace moldraw {

class Atom : Shape {
public:
  Atom(Drawing* dp,
       float x,
       float y,
       float r,
       const PaletteColor &color);
  virtual ~Atom();
  virtual void Draw();
private:
  float _x;
  float _y;
  float _r;
  PaletteColor _color;
};
}

#endif
