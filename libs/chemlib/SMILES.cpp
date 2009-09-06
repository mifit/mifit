#include "SMILES.h"
#include "MIMolIOBase.h"

#include <chemlib/RESIDUE_.h>
#include "FirstToken.h"
#include "SmilesReader.h"

using namespace std;

namespace chemlib {

SMILES::SMILES() {
}

SMILES::~SMILES() {
}

bool SMILES::Read(FILE* fp, MIMolInfo& mol) {
  std::string error, smiles = MIFirstToken(fp);
  return Read(smiles, mol);
}

bool SMILES::Read(const std::string& smiles, MIMolInfo& mol) {
  std::string error;

  // clear current mol
  MIMolInfo foo;
  mol = foo;
  mol.res = SmilesToMol(smiles, mol.bonds, error);
  if (mol.res) {
    mol.res->setSecstr('X');
    mol.res->set_chain_id(' ');
  }
  return mol.res != 0;
}

}
