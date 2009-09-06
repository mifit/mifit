#ifndef mifit_figurelib_Label_h
#define mifit_figurelib_Label_h

#include <string>
#include "Shape.h"

namespace moldraw {

class Label : Shape {
public:
  Label(Drawing* dp,
        float x,
        float y,
        float off,
        float off_x,
        float off_y,
        float font_size,
        const std::string& name);
  virtual ~Label();
  virtual void Draw();
private:
  float _x;
  float _y;
  float _off;
  float _off_x;
  float _off_y;
  float _font_size;
  std::string _name;
};

}

#endif
