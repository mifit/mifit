#include "Spokes.h"

#include "Drawing.h"

namespace moldraw
{

Spokes::Spokes(Drawing *dp,
               float x1,
               float y1,
               float r1,
               float x2,
               float y2,
               float r2,
               float extent,
               float width1,
               float width2)
    : Shape(dp)
{
    _x1 = x1;
    _y1 = y1;
    _r1 = r1;
    _x2 = x2;
    _y2 = y2;
    _r2 = r2;
    _extent = extent;
    _width1 = width1;
    _width2 = width2;
}

void Spokes::Draw()
{
    _dp->DrawSpokes(_x1, _y1, _r1, _x2, _y2, _r2, _extent, _width1, _width2);
}

}
