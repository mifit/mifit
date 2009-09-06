#ifndef CONFLIB_CONFORMATION_FINGERPRINT_H
#define CONFLIB_CONFORMATION_FINGERPRINT_H

#include <vector>

namespace chemlib {
class Bond;
}

namespace conflib {

class ConfFingerprint {
public:
  ConfFingerprint(const std::vector<chemlib::Bond>& distances);
  bool operator==(const ConfFingerprint& rhs);

  void AddCode(const chemlib::Bond& dist) {
    _distanceCodes.push_back(CodeDistance(dist));
  }

protected:
  static int CodeDistance(const chemlib::Bond& dist);

  std::vector<int> _distanceCodes;
};

}

#endif //CONFLIB_CONFORMATION_FINGERPRINT_H
