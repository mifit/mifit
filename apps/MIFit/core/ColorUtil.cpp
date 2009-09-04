#include <string>

#include <chemlib/chemlib.h>

#include <math/mathlib.h>

#include "Colors.h"
#include "ColorUtil.h"


#if !defined (_WIN32)
static
int
strnicmp(const char* s1, const char* s2, int n) {
  std::string c1, c2;
  c1 = (char*)s1;
  c2 = (char*)s2;
  c1=MIToUpper(c1);
  c2=MIToUpper(c2);
  return strncmp(c1.c_str(), c2.c_str(), n);
}

#endif

void init_colornames() {
  Colors::atomcolors.clear();
  Colors::atomnames.clear();
  Colors::atomnames.push_back("");
  Colors::atomnames.push_back("C");
  Colors::atomnames.push_back("O");
  Colors::atomnames.push_back("N");
  Colors::atomnames.push_back("H");
  Colors::atomnames.push_back("S");
  Colors::atomnames.push_back("F");
  Colors::atomnames.push_back("1");
  Colors::atomnames.push_back("2");
  Colors::atomnames.push_back("3");
  Colors::atomnames.push_back("4");
  Colors::atomnames.push_back("P");
  Colors::atomnames.push_back("CU");
  Colors::atomnames.push_back("ZN");

  Colors::atomcolors.push_back(Colors::colornames[Colors::GREEN]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::YELLOW]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::RED]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::BLUE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::CYAN]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::MAGENTA]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::MAGENTA]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::CUSTOM5]);
  Colors::atomcolors.push_back(Colors::colornames[Colors::WHITE]);
}



short color_by_name(const char* name) {
  short color = Colors::WHITE; // default atom color
  for (unsigned int i = 0; i < Colors::atomnames.size(); i++) {
    if (!strnicmp(name, Colors::atomnames[i].c_str(), Colors::atomnames[i].size())) {
      color = Colors::findColorNumber(Colors::atomcolors[i]);
    }
  }
  if (color < 0) {
    color = 0;
  }
  if (color >= Colors_NUMBERCOLORS) {
    color = Colors_NUMBERCOLORS-1;
  }
  return (color);
}

short color_by_name(const char* name, char altloc) {
  short color = color_by_name(name);
  if (!(altloc == ' ' || altloc == 'A' || name[0] != 'C')) {
    color = Colors::CONTOUR6 + altloc - 'A';
    color %= Colors_NUMBERCOLORS;
  }
  if (color == Colors::BLACK) {
    color = Colors::WHITE;
  }
  return color;
}

short secstrcolor(char ss) {
  if (ss == 'H') {
    return Colors::sec_colors[Colors::HELIX];
  }
  if (ss == 'S') {
    return Colors::sec_colors[Colors::SHEET];
  }
  // more types here
  // default return coil
  return (Colors::sec_colors[Colors::COIL]);
}

short bvaluecolor(int b) {
  for (int i = 0; i < Colors_NumBValueColors; i++) {
    if (b < Colors::BValueRanges[i]) {
      return Colors::BValueColors[i];
    }
  }
  return Colors::BValueColors[Colors_NumBValueColors-1];
}

short bvaluecolor(float bvalue) {
  int b = ROUND(bvalue*100.0F);
  for (int i = 0; i < Colors_NumBValueColors; i++) {
    if (b < Colors::BValueRanges[i]) {
      return Colors::BValueColors[i];
    }
  }
  return Colors::BValueColors[Colors_NumBValueColors-1];
}

int Colors::findColorNumber(const std::string &name)
{
  for (int i=0; i < Colors_NCOLORS; ++i) {
    std::string cn(colornames[i]);
    if (cn == name) {
      return i;
    }
  }
  return -1;
}


