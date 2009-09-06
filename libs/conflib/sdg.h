#include <chemlib/chemlib.h>

namespace conflib {

class sdgDistance {
private:
  chemlib::MIAtom& _a1;
  chemlib::MIAtom& _a2;
  bool _isRange;
  bool _isOpenEnd;                  //For distances with no upper limit
  double _lower;
  double _upper;
public:
  //	sdgDistance() : _a1(MIAtom a1()), _a2(MIAtom a2()) {_lower = 0.0F; _upper = 0.0F;}
  sdgDistance(chemlib::MIAtom& a1, chemlib::MIAtom& a2, double lower, double upper, bool isRange);
  sdgDistance(chemlib::MIAtom& a1, chemlib::MIAtom& a2, double ideal_dist, bool isRange);

  sdgDistance(const sdgDistance&);
  sdgDistance& operator=(const sdgDistance&);
  bool operator==(const sdgDistance& d2) const;
  double Measure();
  double Score();
  double GetIdeal();
  void Tweak(float pace);
};

class sdgVolume {
  chemlib::MIAtom& _center;
  chemlib::MIAtom& _a1;
  chemlib::MIAtom& _a2;
  chemlib::MIAtom& _a3;
  bool _isRange;
  int _openEnd;
  double _lower;
  double _upper;
public:
  sdgVolume(chemlib::MIAtom& center, chemlib::MIAtom& a1, chemlib::MIAtom& a2, chemlib::MIAtom& a3, double lower, double upper, int openEnd);
  sdgVolume(chemlib::MIAtom& center, chemlib::MIAtom& a1, chemlib::MIAtom& a2, chemlib::MIAtom& a3, double ideal_volume, int openEnd);
  sdgVolume(chemlib::MIAtom& center, chemlib::MIAtom& a1, chemlib::MIAtom& a2, chemlib::MIAtom& a3, double ideal_volume);

  sdgVolume(const sdgVolume&);
  sdgVolume& operator=(const sdgVolume&);
  bool operator==(const sdgVolume& vol2) const;
  double Measure();
  double Score();
  void Tweak(float pace);
};

class sdgEngine {
private:
  std::vector<sdgDistance>& _distances;
  std::vector<sdgVolume>& _volumes;
  int _nVol;                    //Vector sizes included for efficiency
  int _nDist;
  int _nSteps;                  //Steps per cycle (proportional to # of atoms)
  double _vol_odds;             //The chance, for each step, of tweaking a volume constraint (set in constructor)
  const std::vector<chemlib::MIAtom*>& _atoms;

public:
  sdgEngine(std::vector<sdgDistance>& dists, std::vector<sdgVolume>& vols, const std::vector<chemlib::MIAtom*>& atoms);
  double DoOptimize();
  void DoCycle(float, float);
  void DoStep(float, float);
  void Explode(float factor);
};

/*
   bool MatchDist(sdgDistance &d1, sdgDistance &d2);
    struct MatchDist : public std::binary_function<const sdgDistance, const sdgDistance, bool> {
        inline bool operator()(const sdgDistance, const sdgDistance) const {
            return sdgDistance._a1
    };
 */

}
