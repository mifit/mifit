#include "Sun.h"

#include "Drawing.h"

namespace moldraw {

Sun::Sun(Drawing* dp,
         float x,
         float y,
         float r,
         float x_dir,
         float y_dir,
         float extent,
         float font_size,
         float width,
         const std::string& label) : Shape(dp) {
  _x = x;
  _y = y;
  _r = r;
  _x_dir = x_dir;
  _y_dir = y_dir;
  _extent = extent;
  _font_size = font_size;
  _width = width;
  _label = label;
}

Sun::~Sun() {
}

void Sun::Draw() {
  _dp->DrawSun(_x, _y, _r, _width, _font_size, _label);
}

}
