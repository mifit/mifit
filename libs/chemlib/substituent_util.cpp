#include <vector>
#include <algorithm>
#include <cmath>

#include <math/mathlib.h>

#include "substituent_util.h"


namespace chemlib {
void SortSubsCounterClockwise(std::vector<Substituent>& subs) {

  std::sort(subs.begin(), subs.end(), LeastTheta() );
}

void EqualizeSpacing(std::vector < Substituent >& subs) {
  //Calculate the angle for uniform spacing
  double spacing_angle = DEG2RAD * 360.0 / (subs.size() + 1);

  double target_arc, current_arc;

  //Loop over substituents
  unsigned int i;
  for (i = 0; i < subs.size(); ++i) {
    target_arc = (i + 1) * spacing_angle;

    current_arc = atan2(subs[i]._direction[1], subs[i]._direction[0]);

    //Rotate around the Z-axis to correct the direction of this substituent
    subs[i].ZAxisRotate(target_arc - current_arc);
  }
}

} //namespace chemlib
