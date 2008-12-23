#ifndef MMCIF_DATA_CLASSES_H
#define MMCIF_DATA_CLASSES_H

#include <vector>
#include <map>
#include <string>

/////////////////////////////////////////////////////////////////////////////
// Class:        CifLoop
// Purpose:     Container to store the column-names and columns from a loop
//				in an mmCIF file
/////////////////////////////////////////////////////////////////////////////
class CifLoop {
public:
  std::string _category;                //e.g. "chem_comp", "chem_comp_atom", "chem_comp_bond"
  std::vector<std::string> _names;              //array of all the keywords, e.g. "chem_comp_atom.atom_id"
  std::string _values;                  //The concatenation of all the values

  CifLoop();
  void Clear();
  void PrintNames(std::string& out);
};

/////////////////////////////////////////////////////////////////////////////
// Class:        CifDataBlock
// Purpose:     Container to take a string with a block from an mmCIF file and
//				parse it into the name-value pairs and loops.
/////////////////////////////////////////////////////////////////////////////
class CifDataBlock {
private:
  bool _isParsed;

  void Parse();
public:
  std::map<std::string, std::string> _items;        //Name-Value pairs defined outside of loops
  std::vector<CifLoop> _loops;                  //Loop data
  std::string _block;
  CifDataBlock();

  void SetInput(const char* input);
  void Clear();

  bool FindLoop(std::string category, CifLoop& loop);
};

/////////////////////////////////////////////////////////////////////////////
// Class:       HeaderKeyIndices
// Purpose:     Container to associate the column numbers in an mmCIF header
//				loop with the data fields
/////////////////////////////////////////////////////////////////////////////
class HeaderKeyIndices {
public:
  int xDescLevel;

  HeaderKeyIndices(std::vector<std::string>& names);
};

/////////////////////////////////////////////////////////////////////////////
// Class:       AtomKeyIndices
// Purpose:     Container to associate the column numbers in an mmCIF atom
//				loop with the data fields
/////////////////////////////////////////////////////////////////////////////
class AtomKeyIndices {
public:
  int xRes;
  int xName;
  int xSymbol;
  int xCharge;
  int xX;
  int xY;
  int xZ;

  AtomKeyIndices(std::vector<std::string>& names);
};

class BondKeyIndices {
public:
  int xRes;
  int xAtom1;
  int xAtom2;
  int xOrder;
  int xLength;
  int xTolerance;

  BondKeyIndices(std::vector<std::string>& names);
};

class AngleKeyIndices {
public:
  int xRes;
  int xAtom1;
  int xAtom2;
  int xAtom3;
  int xAngle;
  int xDist;
  int xTolerance;
  int xDistTolerance;

  AngleKeyIndices(std::vector<std::string>& names);
};

class TorsionKeyIndices {
public:
  int xRes;
  int xName;
  int xAtom1;
  int xAtom2;
  int xAtom3;
  int xAtom4;
  int xAngle;
  int xTolerance;
  int xPeriod;

  TorsionKeyIndices(const std::vector<std::string>& names);
};

class TorsionValueKeyIndices {
public:
  int xRes;
  int xName;
  int xAngle;
  int xTolerance;
  int xDist;
  int xDistTolerance;

  TorsionValueKeyIndices(std::vector<std::string>& names);
};

class PlaneKeyIndices {
public:
  int xRes;
  int xName;
  int xAtom;
  int xTolerance;
  int xAtomCount;
  int xHvyAtomCount;

  PlaneKeyIndices(std::vector<std::string>& names);
};

class PlaneAtomKeyIndices {
public:
  int xRes;
  int xName;
  int xAtom;
  int xTolerance;

  PlaneAtomKeyIndices(std::vector<std::string>& names);
};

class ChiralKeyIndices {
public:
  int xRes;
  int xName;
  int xCenter;
  int xAtom1;
  int xAtom2;
  int xAtom3;
  int xVolume;
  int xConfig;
  int xDegree;                              //# of atoms attached to chrl center
  int xHvyDegree;                           //# of attached non-hydrogen atoms

  ChiralKeyIndices(std::vector<std::string>& names);
};

class ChiralAtomKeyIndices {
public:
  int xName;
  int xRes;
  int xAtom;

  ChiralAtomKeyIndices(std::vector<std::string>& names);
};

std::string GetCategory(const std::vector<std::string>& names);
void ClipCategory(std::vector<std::string>& names);

unsigned char DecodeCifBondOrder(const std::string& value);


#endif //CIF_DATA_CLASSES_H
