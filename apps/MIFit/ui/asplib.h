#ifndef ASPLIB_H_
#define ASPLIB_H_

#include "chemlib.h"

class asplib {
public:

  /* flag values - 0 = connect to next and no mark - default */
  /* magic number meaning autoscale */
  static const float GR_AUTO;
  static const int GR_BREAK;
  static const int GR_CROSS;
  static const int GR_BOX;
  static const int GR_CIRCLE;
  static const int GR_DASH;
  static const int GR_LABEL;
  static const int GR_SMALLCROSS;
  static const int GR_SMALLBOX;
  static const int GR_SMALLCIRCLE;
  static const int GR_POINT;

  /* stuff for graphstyle */
  /* convert x axis from (sin(theta)/lambda)**3 to resolution in Angstroms */
  static const int GR_CUBEDTORES;

};

#endif //ASPLIB_H_
