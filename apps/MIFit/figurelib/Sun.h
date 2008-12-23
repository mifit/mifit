#ifndef mifit_figurelib_Sun_h
#define mifit_figurelib_Sun_h

#include <string>
#include "Shape.h"

namespace moldraw {

class Sun : Shape {
public:
  Sun(Drawing* dp,
      float x,
      float y,
      float r,
      float x_dir,
      float y_dir,
      float extent,
      float font_size,
      float width,
      const std::string& label);
  virtual ~Sun();
  virtual void Draw();
private:
  float _x;
  float _y;
  float _r;
  float _x_dir;
  float _y_dir;
  float _extent;
  float _font_size;
  float _width;
  std::string _label;
};

}
#endif
