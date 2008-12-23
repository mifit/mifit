#ifndef CONFLIB_CONFORMATION_ITERATOR_H
#define CONFLIB_CONFORMATION_ITERATOR_H

#include "chemlib.h"

#include "FlexTorsion.h"
#include "ConfFingerprint.h"

namespace conflib {

class ConfIterator {
public:
  //Constructor
  //	ConfIterator(chemlib::RESIDUE *res, vector<chemlib::Bond> &bonds);

  //Destructor
  virtual ~ConfIterator() {
  }

  //Iteration
  virtual bool Next() = 0;

  //Batch Generation
  virtual int GenerateConfs(chemlib::ConfSaver& confs, int max) = 0;

  //Ensemble size info
  virtual bool FiniteConfs() const = 0;
  virtual int NumberConfs() const = 0;
};

class ConfEnumerator : public ConfIterator {
public:
  //Constructor
  ConfEnumerator(chemlib::RESIDUE* res, std::vector<chemlib::Bond>& bonds, std::vector<chemlib::TORSION>& torsions);

  //Destructor
  virtual ~ConfEnumerator();

  //Iteration
  virtual bool Next();

  //Batch Generation
  virtual int GenerateConfs(chemlib::ConfSaver& confs, int max = 1000);

  //Ensemble size info
  virtual bool FiniteConfs() const {
    return true;
  }

  virtual int NumberConfs() const {
    return _nTheory;
  }

  chemlib::RESIDUE* _res;

private:
  int _nTry;
  int _nTheory;
  bool _init;
  std::vector<FlexTorsion> _flexors;
  std::vector<chemlib::Bond> _bumps;
  std::vector<ConfFingerprint> _prints;
  //	vector< vector< int > > _distanceHashKeys;

  bool CheckConf();

  //Declare these private to prevent use
  //Copy Constructor
  ConfEnumerator(const ConfEnumerator& rhs);
  //Assignment Operator
  ConfEnumerator& operator=(const ConfEnumerator& rhs);
};

//Class to sample conformations by randomly positioning torsions, largely
//used for big molecules in lieu of exhaustive enumeration

class ConfSampler : public ConfIterator {
public:
  //Constructor
  ConfSampler(chemlib::RESIDUE* res, std::vector<chemlib::Bond>& bonds, std::vector<chemlib::TORSION>& torsions);

  //Destructor
  virtual ~ConfSampler();

  //Iteration
  virtual bool Next();

  //Batch Generation
  virtual int GenerateConfs(chemlib::ConfSaver& confs, int max);

  //Ensemble size info
  virtual bool FiniteConfs() const {
    return false;
  }

  virtual int NumberConfs() const {
    return -1;
  }

  chemlib::RESIDUE* _res;

private:
  std::vector<FlexTorsion> _flexors;
  std::vector<chemlib::Bond> _bumps;
  std::vector<ConfFingerprint> _prints;         //Tracks which conformations have been generated already

  bool CheckConf();                         //Checks for bumps and redundancy

  //Declare these private to prevent use
  //Copy Constructor
  ConfSampler(const ConfSampler& rhs);
  //Assignment Operator
  ConfSampler& operator=(const ConfSampler& rhs);
};

} //namespace conflib
#endif //CONFLIB_CONFORMATOIN_ITERATOR_H
