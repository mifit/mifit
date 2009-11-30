#include "Atom.h"

#include "Drawing.h"

namespace moldraw
{

Atom::Atom(Drawing *dp,
           float x,
           float y,
           float r,
           const PaletteColor &color)
    : Shape(dp)
{
    _x = x;
    _y = y;
    _r = r;
    _color = color;
}

Atom::~Atom()
{
}

void Atom::Draw()
{
    _dp->DrawAtom(_x, _y, _r, _color);
}

}
