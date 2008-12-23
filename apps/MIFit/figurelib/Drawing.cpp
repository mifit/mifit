#include "Drawing.h"
#include <cmath>

namespace moldraw {

// note: these must all be different
const float Drawing::ANNOTATION_LABEL_SIZE = 20.0f;
const float Drawing::HBOND_LABEL_SIZE = 11.0f;
const float Drawing::ATOM_LABEL_SIZE = 12.0f;
const float Drawing::SUN_LABEL_SIZE = 13.0;
const float Drawing::RESIDUE_LABEL_SIZE = 14.0f;


const float Drawing::ATOM_LABEL_OFFSET = 1.2f; // actual offset is ATOM_LABEL_OFFSET * ATOM_RADIUS
const float Drawing::ATOM_RADIUS = 0.256f;   // 1/6 avg C-C bond distance
const float Drawing::UNIRES_RADIUS = 0.77f; // 1/2 avc C-C bond distance
const float Drawing::SPOKE_EXTENT = 0.385f; // 1/4 avg C-C bond distance

// these are in pixels
const float Drawing::ATOM_SPOKE_WIDTH = 1.0f;
const float Drawing::SUN_BORDER_WIDTH = 2.0f;
const float Drawing::STD_BOND_WIDTH = 1.0f;
const float Drawing::LIG_BOND_WIDTH = 2.0f;


const int Drawing::MAX_SUN_CYCLES = 600;
const int Drawing::MAX_RESLABEL_CYCLES = 600;


void Drawing::DrawAtom(float x, float y, float r, const PaletteColor& color)
{
  draw_item i;
  i.type=draw_item::Atom;
  i.x1=x;
  i.y1=y;
  i.r=r;
  i.color=color;
  _items.push_back(i);
}

void Drawing::DrawSingleBond(float x1, float y1, float x2, float y2, float width) {
  draw_item i;
  i.type=draw_item::SingleBond;
  i.x1=x1;
  i.y1=y1;
  i.x2=x2;
  i.y2=y2;
  i.width=width;
  _items.push_back(i);
}

void Drawing::DrawDoubleBond(float x1, float y1, float x2, float y2, float width)
{
  draw_item i;
  i.type=draw_item::DoubleBond;
  i.x1=x1;
  i.y1=y1;
  i.x2=x2;
  i.y2=y2;
  i.width=width;
  _items.push_back(i);
}

void Drawing::DrawHydrogenBond(float x1, float y1, float x2, float y2, float width, float font_size, float dist)
{
  draw_item i;
  i.type=draw_item::HydrogenBond;
  i.x1=x1;
  i.y1=y1;
  i.x2=x2;
  i.y2=y2;
  i.width=width;
  i.font_size=font_size;
  i.dist=dist;
  _items.push_back(i);
}

void Drawing::DrawSun(float x, float y, float r, float width, float font_size, const std::string& name)
{
  draw_item i;
  i.type=draw_item::Sun;
  i.x1=x;
  i.y1=y;
  i.r=r;
  i.width=width;
  i.font_size=font_size;
  i.text=name;
  _items.push_back(i);

}

void Drawing::DrawSpokes(float x1, float y1, float r1, float x2, float y2, float r2, float, float width1, float width2)
{
  draw_item i;
  i.type=draw_item::Spokes;
  i.x1=x1;
  i.y1=y1;
  i.r=r1;
  i.x2=x2;
  i.y2=y2;
  i.r2=r2;
  i.width=width1;
  i.width2=width2;
  _items.push_back(i);
}

void Drawing::DrawLabel(float x, float y, float off, float off_x, float off_y, float font_size, const std::string& content)
{
  draw_item i;
  i.type=draw_item::Label;
  i.x1=x;
  i.y1=y;
  i.off=off;
  i.off_x=off_x;
  i.off_y=off_y;
  i.font_size=font_size;
  i.text=content;
  _items.push_back(i);
}




}
