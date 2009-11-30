#ifndef COLORS_H
#define COLORS_H

#include <string>
#include <vector>

static const int Colors_NCOLORS = 41;
static const int Colors_NUMBERCOLORS = 41;
static const int Colors_NUMBERPALETTE = 410;
const int Colors_NumBValueColors = 10;

class Colors
{
public:

    // 2-6-02 added more colors - dem
    // each color has 10 depths (decreasing brightness level)
    static const int PDEPTH;

    static const int BLACK;
    static const int YELLOW;
    static const int RED;
    static const int BLUE;
    static const int CYAN;
    static const int GREEN;
    static const int MAGENTA;
    static const int PINK;
    static const int ORANGE;
    static const int BROWN;
    static const int WHITE;
    static const int CUSTOM1;
    static const int CUSTOM2;
    static const int CUSTOM3;
    static const int CUSTOM4;
    static const int CUSTOM5;
    static const int CUSTOM6;
    static const int CUSTOM7;
    static const int CUSTOM8;
    static const int CUSTOM9;
    static const int CUSTOM10;
    static const int MAP1;
    static const int MAP2;
    static const int MAP3;
    static const int MAP4;
    static const int MAP5;
    static const int MAP6;
    static const int MAP7;
    static const int MAP8;
    static const int MAP9;
    static const int MAP10;
    static const int CONTOUR1;
    static const int CONTOUR2;
    static const int CONTOUR3;
    static const int CONTOUR4;
    static const int CONTOUR5;
    static const int CONTOUR6;
    static const int CONTOUR7;
    static const int CONTOUR8;
    static const int CONTOUR9;
    static const int CONTOUR10;

    // positions of color start in map
    static const int PBLACK;
    static const int PYELLOW;
    static const int PRED;
    static const int PBLUE;
    static const int PCYAN;
    static const int PGREEN;
    static const int PMAGENTA;
    static const int PPINK;
    static const int PORANGE;
    static const int PBROWN;
    static const int PWHITE;
    static const int PCUSTOM1;
    static const int PCUSTOM2;
    static const int PCUSTOM3;
    static const int PCUSTOM4;
    static const int PCUSTOM5;
    static const int PCUSTOM6;
    static const int PCUSTOM7;
    static const int PCUSTOM8;
    static const int PCUSTOM9;
    static const int PCUSTOM10;
    static const int PMAP1;
    static const int PMAP2;
    static const int PMAP3;
    static const int PMAP4;
    static const int PMAP5;
    static const int PMAP6;
    static const int PMAP7;
    static const int PMAP8;
    static const int PMAP9;
    static const int PMAP10;
    static const int PCONTOUR1;
    static const int PCONTOUR2;
    static const int PCONTOUR3;
    static const int PCONTOUR4;
    static const int PCONTOUR5;
    static const int PCONTOUR6;
    static const int PCONTOUR7;
    static const int PCONTOUR8;
    static const int PCONTOUR9;
    static const int PCONTOUR10;

    // coloring methods- number = value returned from color dialog
    static const int COLORC;
    static const int COLORALL;
    static const int COLORSECOND;
    static const int COLORBVALUE;
    static const int COLORTYPE;
    static const int COLORCHARGE;
    static const int COLORHYDRO;
    static const int COLORSHAPELY;
    static const int COLOROFF;
    static const int COLORASIS;
    static const int COLORON;

    // secondary structure colors
    static const int HELIX;
    static const int SHEET;
    static const int COIL;
    static const int TURN;

    static char colornames[Colors_NUMBERCOLORS][8];

    static std::vector<std::string> atomnames;
    static std::vector<std::string> atomcolors;
    static int findColorNumber(const std::string &name);

    // these are declared as global in drawline.cpp
    // the position of each color in the palette after
    // realizing the palette to the screen
    static int palindex[Colors_NUMBERPALETTE];
    static unsigned char RPallette[Colors_NUMBERPALETTE];
    static unsigned char GPallette[Colors_NUMBERPALETTE];
    static unsigned char BPallette[Colors_NUMBERPALETTE];

    // secondary structure colors
    static int sec_colors[];

    //BValue colors
    static int BValueColors[Colors_NumBValueColors];
    static int BValueRanges[Colors_NumBValueColors];
};

//@{
// Returns a palette color index given a color index
// @param color a color such as WHITE or BLACK. See Colors.h
//@}
inline int PaletteIndex(int color)
{
    return color*10;
}

class PaletteColor
{
public:
    PaletteColor(unsigned char r = 0, unsigned char g = 0, unsigned char b = 0) : red(r), green(g), blue(b)
    {
    }
    unsigned char red;
    unsigned char green;
    unsigned char blue;
};

class MIPalette
{
public:
    std::vector<PaletteColor> colors;
};

#endif /* COLORS_H*/
