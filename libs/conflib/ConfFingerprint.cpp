#include <chemlib/chemlib.h>

#include <climits>
#include <algorithm>
#include <cstdio>

#include "ConfFingerprint.h"

//#include <mifit/legacy/mifit_algorithm.h>

using namespace std;
using namespace conflib;
using namespace chemlib;

ConfFingerprint::ConfFingerprint(const vector<Bond>& distances) {
  vector<Bond>::const_iterator i, e = distances.end();
  for (i = distances.begin(); i != e; ++i) {
    AddCode(*i);
  }
  sort(_distanceCodes.begin(), _distanceCodes.end(), greater<int>());
#if 0
  FILE* f = fopen("./test.log", "a");
  for (unsigned int n = 0; n < _distanceCodes.size(); ++n) {
    fprintf(f, "%d,", _distanceCodes[n]);
  }
  fprintf(f, "\n");
  fclose(f);
#endif
}

bool ConfFingerprint::operator ==(const ConfFingerprint& rhs) {
  return _distanceCodes == rhs._distanceCodes;
}

#ifdef ROUND
#undef ROUND
#endif
#define ROUND(a) ((a) > 0 ? (int)((a)+0.5) : (int)((a)-0.5))

int ConfFingerprint::CodeDistance(const Bond& dist) {
  float d = (float)(2.0 * SquaredAtomDist(*dist.getAtom1(), *dist.getAtom2()));

  return std::min(dist.getAtom1()->atomicnumber() - 1, 53) +
         54 * std::min(dist.getAtom1()->hybrid() - 1, 2) +
         162 * std::min(dist.getAtom2()->atomicnumber() - 1, 53) +
         8748 * std::min(dist.getAtom2()->hybrid() - 1, 2) +
         26244 * std::min(ROUND(d), 81000);
}

