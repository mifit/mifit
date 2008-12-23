#include "DoubleBond.h"

#include "Drawing.h"

namespace moldraw {

DoubleBond::DoubleBond(Drawing* dp,
                       float x1,
                       float y1,
                       float x2,
                       float y2,
                       float width) : Shape(dp) {
  _x1 = x1;
  _y1 = y1;
  _x2 = x2;
  _y2 = y2;
  _width = width;
}

DoubleBond::~DoubleBond() {
}

void DoubleBond::Draw() {
  _dp->DrawDoubleBond(_x1, _y1, _x2, _y2, _width);
}

}
