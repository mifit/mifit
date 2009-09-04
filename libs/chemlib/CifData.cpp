#include <functional>
#include <algorithm>

#include "CifData.h"
#include "CifParser.h"
#include "model.h"
#include <util/system.h>

using namespace std;

CifLoop::CifLoop() {
}

void CifLoop::Clear() {
  _values.clear();
  _category.clear();
  _names.clear();
}

void CifLoop::PrintNames(std::string& out) {
  for (unsigned int i = 0; i < _names.size(); ++i) {
    out += _names[i];
    out += "\n";
  }
}

CifDataBlock::CifDataBlock() {
  _isParsed = false;
}

void CifDataBlock::SetInput(const char* input) {
  _block = input;
  //	CifLoop loop;
  //	loop._values = new char[strlen(input)] ;
  //	strcpy(loop._values, input);
  //	_loops.push_back(loop);
}

void CifDataBlock::Clear() {
  _items.clear();
  _loops.clear();
  _block.clear();
  _isParsed = false;
}

bool CifDataBlock::FindLoop(std::string category, CifLoop& loop) {
  if (!_isParsed) {
    Parse();
  }
  std::vector<CifLoop>::iterator i, e = _loops.end();

  for (i = _loops.begin(); i != e; ++i) {
    if (i->_category == category) {
      loop = *i;
      return true;
    }
  }

  return false;
}

#define CTP_WANTNAME 0
#define CTP_WANTVALUE 1
#define CTP_LOOP 2

void CifDataBlock::Parse() {
  CifTokenizer toker(_block.c_str());
  int state = CTP_WANTNAME;
  std::string name;
  CifLoop loop;

  std::string t;
  while (toker.GetToken(t)) {
    if (t[0] == '#') {
      continue;
    }

    switch (state) {
      case CTP_WANTNAME:
        if (t[0] == '_') {
          name = t;
          state = CTP_WANTVALUE;
        } else if (strncmp(t.c_str(), "loop_", 5) == 0) {
          loop.Clear();
          toker.SlurpNames(loop._names);
          toker.SlurpValues(loop._values);
          loop._category = GetCategory(loop._names);
          if (!loop._category.empty()) {
            ClipCategory(loop._names);
          }
          _loops.push_back(loop);
          state = CTP_WANTNAME;
        }
        break;

      case CTP_WANTVALUE:
        _items[name] = t;
        state = CTP_WANTNAME;
        break;
    }
  }

  _block.clear();
  _isParsed = true;
}

std::string GetCategory(const std::vector<std::string>& names) {
  std::string category, prefix;

  for (unsigned int i = 0; i < names.size(); ++i) {
    prefix = MIBeforeFirst(names[i], '.');
    if (prefix.empty() || prefix == category) {
      continue;
    } else if (category.empty()) {
      category = prefix;
    } else {
      return std::string("");
    }
  }

  return category;
}

void ClipCategory(std::vector<std::string>& names) {
  std::string extension;                        //the portion of the name that follows the "." char

  for (unsigned int i = 0; i < names.size(); ++i) {
    if ((extension = MIAfterFirst(names[i], '.')) != std::string("")) {        //compare to empty
      names[i] = extension;                                         //string to skip
    }                                                                   //names w/o periods
  }
}

HeaderKeyIndices::HeaderKeyIndices(std::vector<std::string>& names) {
  xDescLevel = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "desc_level") {
      xDescLevel = i;
    }
  }
}

AtomKeyIndices::AtomKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xSymbol = -1;
  xCharge = -1;
  xX = -1;
  xY = -1;
  xZ = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "id"
               || names[i] == "atom_id"
               || names[i] == "label_atom_id"
               || names[i] == "label_alt_id") {
      xName = i;
    } else if (names[i] == "type_symbol") {
      xSymbol = i;
    } else if (names[i] == "charge"
               || names[i] == "partial_charge") {
      xCharge = i;
    } else if (names[i] == "x"
               || names[i] == "cartn_x"
               || names[i] == "model_Cartn_x") {
      xX = i;
    } else if (names[i] == "y"
               || names[i] == "cartn_y"
               || names[i] == "model_Cartn_y") {
      xY = i;
    } else if (names[i] == "z"
               || names[i] == "cartn_z"
               || names[i] == "model_Cartn_z") {
      xZ = i;
    }
  }
}

BondKeyIndices::BondKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xAtom1 = xAtom2 = -1;
  xOrder = -1;
  xLength = -1;
  xTolerance = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    }
    if (names[i] == "atom_id_1") {
      xAtom1 = i;
    }
    if (names[i] == "atom_id_2") {
      xAtom2 = i;
    }
    if (names[i] == "type"
        || names[i] == "value_order") {
      xOrder = i;
    }
    if (names[i] == "value_dist") {
      xLength = i;
    }
    if (names[i] == "value_dist_esd") {
      xTolerance = i;
    }
  }
}

AngleKeyIndices::AngleKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xAtom1 = xAtom2 = xAtom3 = -1;
  xAngle = -1;
  xDist = -1;
  xTolerance = -1;
  xDistTolerance = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "atom_id_1"
               || names[i] == "atom_site_id_1"
               || names[i] == "atom_site_label_id_1") {
      xAtom1 = i;
    } else if (names[i] == "atom_id_2"
               || names[i] == "atom_site_id_2"
               || names[i] == "atom_site_label_id_2") {
      xAtom2 = i;
    } else if (names[i] == "atom_id_3"
               || names[i] == "atom_site_id_3"
               || names[i] == "atom_site_label_id_3") {
      xAtom3 = i;
    } else if (names[i] == "value"
               || names[i] == "value_angle") {
      xAngle = i;
    } else if (names[i] == "value_dist") {
      xDist = i;
    } else if (names[i] == "value_esd"
               || names[i] == "value_angle_esd") {
      xTolerance = i;
    } else if (names[i] == "value_dist_esd") {
      xDistTolerance = i;
    }
  }
}

TorsionKeyIndices::TorsionKeyIndices(const std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xAtom1 = -1;
  xAtom2 = -1;
  xAtom3 = -1;
  xAtom4 = -1;
  xAngle = -1;
  xTolerance = -1;
  xPeriod = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "id") {
      xName = i;
    } else if (names[i] == "atom_id_1"
               || names[i] == "atom_site_id_1"
               || names[i] == "atom_site_label_id_1") {
      xAtom1 = i;
    } else if (names[i] == "atom_id_2"
               || names[i] == "atom_site_id_2"
               || names[i] == "atom_site_label_id_2") {
      xAtom2 = i;
    } else if (names[i] == "atom_id_3"
               || names[i] == "atom_site_id_3"
               || names[i] == "atom_site_label_id_3") {
      xAtom3 = i;
    } else if (names[i] == "atom_id_4"
               || names[i] == "atom_site_id_4"
               || names[i] == "atom_site_label_id_4") {
      xAtom4 = i;
    } else if (names[i] == "value"
               || names[i] == "value_angle") {
      xAngle = i;
    } else if (names[i] == "value_esd"
               || names[i] == "value_angle_esd") {
      xTolerance = i;
    } else if (names[i] == "period") {
      xPeriod = i;
    }
  }
}

TorsionValueKeyIndices::TorsionValueKeyIndices(std::vector<std::string>& names) {
  xName = -1;
  xRes = -1;
  xAngle = xTolerance = -1;
  xDist = xDistTolerance = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "tor_id") {
      xName = i;
    } else if (names[i] == "angle") {
      xAngle = i;
    } else if (names[i] == "dist") {
      xDist = i;
    } else if (names[i] == "angle_esd") {
      xTolerance = i;
    } else if (names[i] == "dist_esd") {
      xDistTolerance = i;
    }
  }
}

PlaneKeyIndices::PlaneKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xAtom = -1;
  xTolerance = -1;
  xAtomCount = -1;
  xHvyAtomCount = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "plane_id"
               || names[i] == "id") {
      xName = i;
    } else if (names[i] == "atom_id") {
      xAtom = i;
    } else if (names[i] == "dist_esd") {
      xTolerance = i;
    } else if (names[i] == "number_atoms_all") {
      xAtomCount = i;
    } else if (names[i] == "number_atoms_nh") {
      xHvyAtomCount = i;
    }
  }
}

PlaneAtomKeyIndices::PlaneAtomKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xAtom = -1;
  xTolerance = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "plane_id") {
      xName = i;
    } else if (names[i] == "atom_id") {
      xAtom = i;
    } else if (names[i] == "dist_esd") {
      xTolerance = i;
    }
  }
}

ChiralKeyIndices::ChiralKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xCenter = -1;
  xAtom1 = xAtom2 = xAtom3 = -1;
  xVolume = -1;
  xConfig = -1;
  xDegree = xHvyDegree = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "id") {
      xName = i;
    } else if (names[i] == "atom_id"
               || names[i] == "atom_id_centre") {
      xCenter = i;
    } else if (names[i] == "atom_id_1"
               || names[i] == "atom_site_id_1"
               || names[i] == "atom_site_label_id_1") {
      xAtom1 = i;
    } else if (names[i] == "atom_id_2"
               || names[i] == "atom_site_id_2"
               || names[i] == "atom_site_label_id_2") {
      xAtom2 = i;
    } else if (names[i] == "atom_id_3"
               || names[i] == "atom_site_id_3"
               || names[i] == "atom_site_label_id_3") {
      xAtom3 = i;
    } else if (names[i] == "volume"
               || names[i] == "volume_sign") {
      xVolume = i;
    } else if (names[i] == "atom_config") {
      xConfig = i;
    } else if (names[i] == "number_atoms_all") {
      xDegree = i;
    } else if (names[i] == "number_atoms_nh") {
      xHvyDegree = i;
    }
  }
}

ChiralAtomKeyIndices::ChiralAtomKeyIndices(std::vector<std::string>& names) {
  xRes = -1;
  xName = -1;
  xAtom = -1;

  for (unsigned int i = 0; i < names.size(); ++i) {
    if (names[i] == "comp_id"
        || names[i] == "label_comp_id") {
      xRes = i;
    } else if (names[i] == "chir_id") {
      xName = i;
    } else if (names[i] == "atom_id") {
      xAtom = i;
    }
  }
}

unsigned char DecodeCifBondOrder(const std::string& value) {
  if (strcasecmp(value.c_str(), "single")
      || strcasecmp(value.c_str(), "sing")) {
    return SINGLEBOND;
  } else if (strcasecmp(value.c_str(), "double")
             || strcasecmp(value.c_str(), "doub")) {
    return DOUBLEBOND;
  } else if (strcasecmp(value.c_str(), "triple")
             || strcasecmp(value.c_str(), "trip")) {
    return TRIPLEBOND;
  } else if (strcasecmp(value.c_str(), "aromatic")
             || strcasecmp(value.c_str(), "arom")
             || strcasecmp(value.c_str(), "delo")) {
    return PARTIALDOUBLEBOND;
  }

  return NORMALBOND;
}

