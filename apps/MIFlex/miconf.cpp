#include "nongui.h"
#include "chemlib.h"
#include "conflib.h"
#include "moloptlib.h"
#include "ligandlib.h"
#include <utillib.h>

#include "common.h"

using namespace chemlib;
using namespace conflib;

bool MILigandConfs(const std::string& ligand_filename) {
  // create geomrefiner and dictionary
  MIMolOpt geomrefiner;
  char* MOLIMAGEHOME = getenv("MOLIMAGEHOME");
  if (!MOLIMAGEHOME) {
    Logger::log("MOLIMAGEHOME not defined, aborting\n");
    return false;
  }

  // load dictionary
  std::string newDict(MOLIMAGEHOME);
#ifndef _WIN32
  newDict += "/data/dict.noh.pdb";
#else
  newDict += "\\data\\dict.noh.pdb";
#endif
  geomrefiner.dict.LoadDefaultDictionary(newDict.c_str(), MOLIMAGEHOME);
  MISetDictionary(&geomrefiner.dict);

  // load molecule, model
  MIMoleculeBase* fitmol = LoadLigand(ligand_filename.c_str(), &geomrefiner);
  if (!fitmol) {
    Logger::log("Couldn't load ligand, aborting.");
    return false;
  }
  // get and load emap
  GeomSaver confs;
  int retval = MIGenConfs(confs, fitmol, &geomrefiner, true);

  // clean up
  delete fitmol;
  return retval;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    Logger::log("Usage: %s model_file", argv[0]);
    return -1;
  }

  return (int)MILigandConfs(argv[1]);
}

