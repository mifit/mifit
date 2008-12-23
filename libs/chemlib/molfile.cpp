#include "molfile.h"
#include "MIMolIOBase.h"
#include "RESIDUE_.h"
#include "system.h"
#include "atom_util.h"
#include "mol_util.h"


#include <functional>

#define MAX_MOLFILE_LINE 180

#include <time.h>

using namespace std;

namespace chemlib {

molfile::molfile() {
}

molfile::~molfile() {
}

bool molfile::SquirtAtoms(FILE* fp, const RESIDUE& res, const vector<Bond>& bonds,
                          const vector<CHIRAL>& chirals) {
  int chirality;
  for (int i = 0; i < res.atomCount(); ++i) {

    if (SearchChirals(res.atom(i), chirals) != 0) {
      chirality = ChiralCode(res.atom(i), res, bonds);
    } else {
      chirality = 0;
    }

    if (!fprintf(fp, "%10.4f%10.4f%10.4f %-2s %2d%3d%3d%3d%3d%3d\n",
          res.atom(i)->x(),
          res.atom(i)->y(),
          res.atom(i)->z(),
          Left_Atomic_Name(res.atom(i)->atomicnumber()),
          0,                                                //default isotope atomic mass
          ChargeCode(res.atom(i)->formal_charge()),
          chirality,
          0,
          0,
          0)) {
      return false;
    }
  }
  return true;
}

bool molfile::SquirtBonds(FILE* fp, const RESIDUE& res, const vector<Bond>& bonds) {
  int x1, x2;

  for (unsigned int i = 0; i < bonds.size(); ++i) {
    if (bonds[i].getOrder() == HYDROGENBOND || bonds[i].getOrder() == IONICBOND) {
      continue;
    }
    if ((x1 = GetAtomIndex(bonds[i].getAtom1(), res)) < 0) {
      continue;
    }
    if ((x2 = GetAtomIndex(bonds[i].getAtom2(), res)) < 0) {
      continue;
    }
    if (!fprintf(fp, "%3d%3d%3d%3d  0  0\n",
          x1 + 1,
          x2 + 1,
          BondCode(bonds[i].getOrder()),
          0)) {
      return false;
    }
  }
  return true;
}

//bool molfile::GetAtomIndex(int &index, const MIAtom *patom, const RESIDUE &res) {
//	index = patom - res.atom(0);				//Subtract ptrs to get index
//
//	return (index >= 0) && (index < res.natoms);
//}

int molfile::ChargeCode(int charge) {
  switch (charge) {
    case 1:
      return 3;
    case 2:
      return 2;
    case 3:
      return 1;
    case -1:
      return 5;
    case -2:
      return 6;
    case -3:
      return 7;
    default:
      return 0;
  }
}

int molfile::ChiralCode(const MIAtom* patom, const RESIDUE& res, const vector<Bond>& bonds) {
  if (patom->chiral_class() != CH_TETRAHEDRAL) {
    return 0;
  }

  vector<MIAtom*> nabors;
  for (int i = 0; i < res.atomCount(); ++i) {
    if (patom != res.atom(i) && AlreadyBonded(patom, res.atom(i), bonds)) {
      nabors.push_back(res.atom(i));
    }
  }

  if (nabors.size() < 3) {
    return 0;
  }

  if (SignedAtomVolume(*patom, *nabors[0], *nabors[1], *nabors[2]) < 0.0) {
    return 1;
  } else {
    return 2;
  }
}

int molfile::BondCode(unsigned char order) {
  if (order == SINGLEBOND) {
    return 1;
  } else if (order == DOUBLEBOND) {
    return 2;
  } else if (order == TRIPLEBOND) {
    return 3;
  } else if (order == PARTIALDOUBLEBOND) {          //Not strictly valid under the MOLfile spec
    return 4;                                   //should kekulize to either single or double bnd
  } else {                                          //Use single bond as default
    return 1;
  }
}

int molfile::ProcessCharge(int charge_code) {
  switch (charge_code) {
    case 1:
      return 3;
    case 2:
      return 2;
    case 3:
      return 1;
    case 4:
      return 0;                 //Specifies "doublet radical"
    case 5:
      return -1;
    case 6:
      return -2;
    case 7:
      return -3;
    default:
      return 0;
  }
}

void molfile::ProcessChiral(MIAtom* patom, int chiral_code) {

  if (chiral_code == 0) {
    patom->chiral_class(CH_NONE);
  } else {
    patom->chiral_class(CH_TETRAHEDRAL);
    patom->chiral_order(chiral_code);
  }
}

unsigned char molfile::ProcessOrder(int order_code) {
  switch (order_code) {
    case 1:
      return SINGLEBOND;
    case 2:
      return DOUBLEBOND;
    case 3:
      return TRIPLEBOND;
    case 4:
      return PARTIALDOUBLEBOND;
    default:
      return SINGLEBOND;
  }
}

char molfile::ProcessStereo(int stereo_code, int bond_order) {
  if (bond_order == SINGLEBOND) {
    switch (stereo_code) {
      case 1:
        return STEREO_WEDGE;
      case 6:
        return STEREO_DASHED;
      default:
        return STEREO_NONE;
    }
  } else if (bond_order == DOUBLEBOND) {
    switch (stereo_code) {
      case 0:
        return STEREO_INFER;
      case 3:
        return STEREO_EITHER;
      default:
        return STEREO_NONE;
    }
  }
  return STEREO_NONE;
}

bool molfile::Write(FILE* fp, MIMolInfo& mol) {
  time_t t;
  time(&t);
  struct tm* foo = localtime(&t);

  while (Residue::isValid(mol.res)) {
    int ischiral = mol.chirals.empty() ?  0 : 1;
    //int nbonds = CountResBonds(*mol.res, mol.bonds);

    if (!fprintf(fp, "%3s\n", mol.res->type().c_str())) {                                   //Molecule/residue name
      return false;
    }
    if (!fprintf(fp, "  MIFit   %02d%02d%02d%02d%02d3D                              \n",
          foo->tm_mon, foo->tm_mday, foo->tm_year%100, foo->tm_hour, foo->tm_min)) {
      return false;
    }
    if (!fprintf(fp, "\n")) {                                 //Comments line
      return false;
    }
    if (!fprintf(fp, "%3d%3d  0  0  %1d  0              1 V2000\n",
          mol.res->atomCount(), mol.bonds.size(), ischiral)) {
      return false;
    }

    if (!SquirtAtoms(fp, *mol.res, mol.bonds, mol.chirals)) {
      return false;
    }

    if (!SquirtBonds(fp, *mol.res, mol.bonds)) {
      return false;
    }

    if (!fprintf(fp, "M  END\n")) {
      return false;
    }
    mol.res = mol.res->next();
  }
  return true;
}

bool molfile::SlurpAtoms(FILE* fp, RESIDUE* res, int natoms) {
  char line[MAX_MOLFILE_LINE], symbol[3];
  int charge_code, chiral_code;

  for (int i = 0; i < natoms; ++i) {
    res->addAtom(new MIAtom);

    if (!fgets(line, MAX_MOLFILE_LINE, fp)) {
      return false;
    }
    res->atom(i)->setAtomnumber(i+1);

    //		sscanf(line, "%10f%10f%10f%4s%2d%3d%3d%3d%3d%3d",
    //			&res->atom(i)->x, &res->atom(i)->y, &res->atom(i)->z,
    //			&symbol, &mass_diff, &charge_code, &chiral_code,
    //			&hcount, &stereo_care, &valence);

    res->atom(i)->setPosition((float)atof(std::string(line, 10).c_str()),
        (float)atof(std::string(line + 10, 10).c_str()),
        (float)atof(std::string(line + 20, 10).c_str()));
    strncpy(symbol, std::string(line + 31, 2).c_str(), 3);
    charge_code = atoi(std::string(line + 36, 3).c_str());
    chiral_code = atoi(std::string(line + 39, 3).c_str());
    //Be sure the element symbol is upper case and right-justified
    if (symbol[1] == ' ') {
      symbol[1] = toupper(symbol[0]);
      symbol[0] = ' ';
    } else {
      symbol[0] = toupper(symbol[0]);
      symbol[1] = toupper(symbol[1]);
    }

    //Lookup the atomic number
    res->atom(i)->setAtomicnumber(Atomic_Number(symbol));
    //Give the atom a unique name
    res->atom(i)->setName(format("%s%d", Left_Atomic_Name(res->atom(i)->atomicnumber()),
      res->atom(i)->atomnumber()).c_str());

    //Give the atom a color
    if (MIGetColorSetter()) {
      (*MIGetColorSetter())(res->atom(i));
    }

    //Lookup the charge
    if (charge_code != 0) {
      res->atom(i)->set_formal_charge(ProcessCharge(charge_code));
    }

    if (chiral_code != 0) {
      ProcessChiral(res->atom(i), chiral_code);
    }
  }
  return true;
}

bool molfile::SlurpBonds(FILE* fp, vector<Bond>& bonds, const RESIDUE* res, int nbonds) {
  char line[MAX_MOLFILE_LINE];
  int xAtom1, xAtom2, order, stereo;
  Bond bond;

  for (int i = 0; i < nbonds; ++i) {
    bond.Clear();
    if (!fgets(line, MAX_MOLFILE_LINE, fp)) {
      return false;
    }

    sscanf(line, "%3d%3d%3d%3d", &xAtom1, &xAtom2, &order, &stereo);

    bond.setAtom1(res->atom(xAtom1 - 1));
    bond.setAtom2(res->atom(xAtom2 - 1));

    //Lookup the bond order
    bond.setOrder(ProcessOrder(order));
    bond.stereo = ProcessStereo(stereo, bond.getOrder());

    bonds.push_back(bond);
  }
  return true;
}

bool molfile::Read(FILE* fp, MIMolInfo& mol) {
  char line[MAX_MOLFILE_LINE];

  // clear current mol
  MIMolInfo foo;
  mol = foo;
  mol.res = new RESIDUE();

  fgets(line, MAX_MOLFILE_LINE, fp);        //mol header lines
  fgets(line, MAX_MOLFILE_LINE, fp);
  if (strlen(line) >= 22) {
    char type[3];
    type[0] = line[20];
    type[1] = line[21];
    type[2] = '\0';
    mol.res->setType(type);
  }

  fgets(line, MAX_MOLFILE_LINE, fp);        //ctab header lines
  fgets(line, MAX_MOLFILE_LINE, fp);
  std::string natoms(line, 3);
  std::string nbonds(line + 3, 3);

  SlurpAtoms(fp, mol.res, atoi(natoms.c_str()));
  SlurpBonds(fp, mol.bonds, mol.res, atoi(nbonds.c_str()));
  SlurpProperties(fp, mol.res);
  return true;
}

bool molfile::SlurpProperties(FILE* fp, RESIDUE* res) {
  bool first_charge = true;
  int xAtom, charge;
  char line[MAX_MOLFILE_LINE];

  while (fgets(line, MAX_MOLFILE_LINE, fp)) {
    if (!strncmp(line, "M  END", 6)) {
      return true;
    }
    if (!strncmp(line, "M  CHG", 6)) {
      if (first_charge) {
        ClearCharges(res);
        first_charge = false;
      }
      MolPropertyLine mpl(line);
      int n = mpl.NumEntries();
      while (n > 0) {
        mpl.GetEntry(xAtom, charge);
        res->atom(xAtom)->set_formal_charge(charge);
        --n;
      }
    }
  }
  return false;
}

const CHIRAL* molfile::SearchChirals(const MIAtom* atom, const vector<CHIRAL>& chirals) {
  vector<CHIRAL>::const_iterator i, e;
  //Are there any chirals?
  if (chirals.empty() ) {
    return NULL;
  }

  i = chirals.begin();
  e = chirals.end();
  for (; i != e; i++) {
    if (i->center == atom
        && i->flags != CHIRAL_DELETED) {
      return &*i;
    }
  }
  return NULL;
}

MolPropertyLine::MolPropertyLine(const char* input) {
  _line = input;
  _length = _line.length();                     //Cache this since it won't change
  _cur_pos = 10;
}

int MolPropertyLine::NumEntries() {
  int stated_n = atoi(_line.substr(6, 3).c_str());
  int n = 0;
  size_t pos = 10;
  long tmp = 0l;
  std::string xAtom, value;

  int i = 0;
  while (i < stated_n && pos + 7 < _length) {
    xAtom = _line.substr(pos, 3);
    value = _line.substr(pos+4, 3);
    if (MIStringToNumber(xAtom, tmp)
        && MIStringToNumber(value, tmp)) {
      n++;
    }
    pos += 8;
  }
  return (n < stated_n) ? n : stated_n;
}

bool MolPropertyLine::GetEntry(int& xAtom, int& value) {
  long xAtom_lng = 0l, value_lng = 0l;
  if (MIStringToNumber(_line.substr(_cur_pos, 3), xAtom_lng) == false) {
    return false;
  }
  if (MIStringToNumber(_line.substr(_cur_pos+4, 3), value_lng) == false) {
    return false;
  }

  xAtom = (int) xAtom_lng;
  value = (int) value_lng;
  return true;
}

}
