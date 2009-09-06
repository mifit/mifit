#include "Label.h"

#include "Drawing.h"

namespace moldraw {

Label::Label(Drawing* dp,
             float x,
             float y,
             float off,
             float off_x,
             float off_y,
             float font_size,
             const std::string& name) : Shape(dp) {
  _x = x;
  _y = y;
  _off = off;
  _off_x = off_x;
  _off_y = off_y;
  _font_size = font_size;
  _name = name;
}

Label::~Label() {
}

void Label::Draw() {
  _dp->DrawLabel(_x, _y, _off, _off_x, _off_y, _font_size, _name);
}

}
