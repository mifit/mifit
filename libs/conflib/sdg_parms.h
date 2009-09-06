#ifndef SDG_PARAMETER_DEFS_H
#define SDG_PARAMETER_DEFS_H

namespace conflib {

namespace sdg_params {
//"Temperature" factors determines step size for correcting
//distances and volumes
const float SDG_DIST_TEMP = 1.0F;
const float SDG_VOL_TEMP = 1.0F;
const int SDG_NCYCLES = 50;


//Multiply by nAtoms to get # of steps per cycle
const int SDG_NSTEP_FACTOR = 100;
//Rate of temperature decrease
const float SDG_ANNEAL_DIST = (0.9F / SDG_NCYCLES);
//Rate of temperature decrease
const float SDG_ANNEAL_VOL = (0.9F / SDG_NCYCLES);
//Maximum chance of adjusting a volume
const float SDG_MIN_VOL_ODDS = 0.5;
//Prevents division by zero
const double SDG_NONZERO = 1E-10;

const int SDG_OPENEND_NONE = 0;
const int SDG_OPENEND_LOW = 1;
const int SDG_OPENEND_HIGH = 2;
}
}
#endif

