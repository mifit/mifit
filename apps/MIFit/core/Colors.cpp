#include "Colors.h"

const int Colors::PDEPTH = 10;

const int Colors::BLACK = 0;
const int Colors::YELLOW = 1;
const int Colors::RED = 2;
const int Colors::BLUE = 3;
const int Colors::CYAN = 4;
const int Colors::GREEN = 5;
const int Colors::MAGENTA = 6;
const int Colors::PINK = 7;
const int Colors::ORANGE = 8;
const int Colors::BROWN = 9;
const int Colors::WHITE = 10;
const int Colors::CUSTOM1 = 11;
const int Colors::CUSTOM2 = 12;
const int Colors::CUSTOM3 = 13;
const int Colors::CUSTOM4 = 14;
const int Colors::CUSTOM5 = 15;
const int Colors::CUSTOM6 = 16;
const int Colors::CUSTOM7 = 17;
const int Colors::CUSTOM8 = 18;
const int Colors::CUSTOM9 = 19;
const int Colors::CUSTOM10 = 20;
const int Colors::MAP1 = 21;
const int Colors::MAP2 = 22;
const int Colors::MAP3 = 23;
const int Colors::MAP4 = 24;
const int Colors::MAP5 = 25;
const int Colors::MAP6 = 26;
const int Colors::MAP7 = 27;
const int Colors::MAP8 = 28;
const int Colors::MAP9 = 29;
const int Colors::MAP10 = 30;
const int Colors::CONTOUR1 = 31;
const int Colors::CONTOUR2 = 32;
const int Colors::CONTOUR3 = 33;
const int Colors::CONTOUR4 = 34;
const int Colors::CONTOUR5 = 35;
const int Colors::CONTOUR6 = 36;
const int Colors::CONTOUR7 = 37;
const int Colors::CONTOUR8 = 38;
const int Colors::CONTOUR9 = 39;
const int Colors::CONTOUR10 = 40;

// positions of color start in map
const int Colors::PBLACK = 0;
const int Colors::PYELLOW = 10;
const int Colors::PRED = 20;
const int Colors::PBLUE = 30;
const int Colors::PCYAN = 40;
const int Colors::PGREEN = 50;
const int Colors::PMAGENTA = 60;
const int Colors::PPINK = 70;
const int Colors::PORANGE = 80;
const int Colors::PBROWN = 90;
const int Colors::PWHITE = 100;
const int Colors::PCUSTOM1 = 110;
const int Colors::PCUSTOM2 = 120;
const int Colors::PCUSTOM3 = 130;
const int Colors::PCUSTOM4 = 140;
const int Colors::PCUSTOM5 = 150;
const int Colors::PCUSTOM6 = 160;
const int Colors::PCUSTOM7 = 170;
const int Colors::PCUSTOM8 = 180;
const int Colors::PCUSTOM9 = 190;
const int Colors::PCUSTOM10 = 200;
const int Colors::PMAP1 = 210;
const int Colors::PMAP2 = 220;
const int Colors::PMAP3 = 230;
const int Colors::PMAP4 = 240;
const int Colors::PMAP5 = 250;
const int Colors::PMAP6 = 260;
const int Colors::PMAP7 = 270;
const int Colors::PMAP8 = 280;
const int Colors::PMAP9 = 290;
const int Colors::PMAP10 = 300;
const int Colors::PCONTOUR1 = 310;
const int Colors::PCONTOUR2 = 320;
const int Colors::PCONTOUR3 = 330;
const int Colors::PCONTOUR4 = 340;
const int Colors::PCONTOUR5 = 350;
const int Colors::PCONTOUR6 = 360;
const int Colors::PCONTOUR7 = 370;
const int Colors::PCONTOUR8 = 380;
const int Colors::PCONTOUR9 = 390;
const int Colors::PCONTOUR10 = 400;

// coloring methods- number = value returned from color dialog
const int Colors::COLORC = 0;
const int Colors::COLORALL = 1;
const int Colors::COLORSECOND = 2;
const int Colors::COLORBVALUE = 3;
const int Colors::COLORTYPE = 4;
const int Colors::COLORCHARGE = 5;
const int Colors::COLORHYDRO = 6;
const int Colors::COLORSHAPELY = 7;
const int Colors::COLOROFF = (-1);
const int Colors::COLORASIS = 8;
const int Colors::COLORON = 9;

// secondary structure colors
const int Colors::HELIX = 0;
const int Colors::SHEET = 1;
const int Colors::COIL = 2;
const int Colors::TURN = 3;

char Colors::colornames[Colors_NUMBERCOLORS][8] = {
  "Black", "Yellow", "Red", "Blue", "Cyan", "Green", "Magenta",
  "Pink", "Orange", "Brown", "White",
  "User 1", "User 2", "User 3", "User 4", "User 5", "User 6", "User 7", "User 8", "User 9", "User 10",
  "Map 1", "Map 2", "Map 3", "Map 4", "Map 5", "Map 6", "Map 7", "Map 8", "Map 9", "Map 10",
  "Level 1", "Level 2", "Level 3", "Level 4", "Level 5", "Level 6", "Level 7", "Level 8", "Level 9", "Level10"
};

std::vector<std::string> Colors::atomnames;
std::vector<std::string> Colors::atomcolors;

// these are declared as global in drawline.cpp
// the position of each color in the palette after
// realizing the palette to the screen
int Colors::palindex[Colors_NUMBERPALETTE];
unsigned char Colors::RPallette[Colors_NUMBERPALETTE];
unsigned char Colors::GPallette[Colors_NUMBERPALETTE];
unsigned char Colors::BPallette[Colors_NUMBERPALETTE];

// secondary structure colors
int Colors::sec_colors[] = {MAGENTA, RED, BLUE, GREEN};

//BValue colors
int Colors::BValueColors[Colors_NumBValueColors] = {
  PCONTOUR1, PCONTOUR2, PCONTOUR3, PCONTOUR4, PCONTOUR5, PCONTOUR6, PCONTOUR7, PCONTOUR8, PCONTOUR9, PCONTOUR10
};

int Colors::BValueRanges[Colors_NumBValueColors] = {
  300, 600, 1000, 1500, 2000, 2500, 3500, 5000, 7500, 10000
};

