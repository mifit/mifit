#ifndef VALENCE_LOOKUP_H
#define VALENCE_LOOKUP_H

#include <vector>

namespace chemlib {

namespace ValenceType {
const unsigned int vtH = 1;
const unsigned int vtO = 2;
const unsigned int vtB = 3;
const unsigned int vtC = 4;
const unsigned int vtS = 5;
const unsigned int vtN = 6;
const unsigned int vtI = 7;
}

std::vector<int> GetValenceStates(int atomicnumber);

unsigned int GetValenceType(int atomicnumber);
}

#endif //VALENCE_LOOKUP_H
