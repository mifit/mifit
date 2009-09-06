#include "HydrogenBond.h"

#include "Drawing.h"

namespace moldraw {

HydrogenBond::HydrogenBond(Drawing* dp,
                           float x1,
                           float y1,
                           float x2,
                           float y2,
                           float width,
                           float font_size,
                           float dist) : Shape(dp) {
  _x1 = x1;
  _y1 = y1;
  _x2 = x2;
  _y2 = y2;
  _width = width;
  _font_size = font_size;
  _dist = dist;
}

HydrogenBond::~HydrogenBond() {
}

void HydrogenBond::Draw() {
  _dp->DrawHydrogenBond(_x1, _y1, _x2, _y2, _width, _font_size, _dist);
}

}
