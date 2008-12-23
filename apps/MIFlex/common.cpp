#include "nongui.h"
#include "chemlib.h"
#include "conflib.h"
#include "maplib.h"
#include "moloptlib.h"
#include "ligandlib.h"
#include <utillib.h>

#include "common.h"

using namespace chemlib;
using namespace conflib;

bool MIGenConfs(GeomSaver confs, MIMoleculeBase* fitmol, MIMolOpt* opt, bool print) {
  RESIDUE *fitres=fitmol->getResidues();
  if (fitres != 0 && fitres->natoms() > 60) {
    Logger::log("Too many atoms for the refiner (more than 60)"
      "- action canceled");
    return false;
  }

  GeomSaver entryConfs;
  const char* type = fitres->type.c_str();
  GetConfs(entryConfs, opt->dict.GetDictResidue(type, 0), &opt->dict, fitmol);
  Logger::debug("%d conformations before generation", entryConfs.NumberSets());
  if (entryConfs.NumberSets() <= 2) {
    int numberOfConformations =
      conflib::GenerateEnsemble(opt->dict.GetDictResidue(type, 0),
        *(opt->dict.GetDictBonds(type, 0)),
        &opt->dict, true);
    Logger::log("Generated %d confirmations", numberOfConformations);
  }

  if (!GetConfs(confs, fitres, &opt->dict, fitmol)) {
    Logger::log("Error getting conformers");
    return false;
  } 

  Logger::log("Generated %d conformations for residue %s.",
              confs.NumberSets(), fitres->type.c_str());

  // Aaaaaa...nd, print the results
  for (int i=0; print && i < confs.NumberSets(); ++i) {
    char buf[5];
    sprintf(buf,"%4d",i);
    fitres->name=buf;
    confs.Restore(i);
    SavePDB(stdout, fitres, NULL, 0, false);
  }

  return true;
}

MIMoleculeBase* LoadMol(const char* fname) {
  MIMoleculeBase* mol;
  RESIDUE* res;
  std::vector<Bond> connects;
  Bond* edges = NULL;

  FILE* fp = 0;
  if ((fp = fopen(fname, "r")) == NULL) {
    Logger::log("Can't open file %s\n", fname);
    return 0;
  }
  if ((res = LoadPDB(fp, &connects)) == NULL) {
    fclose(fp);
    Logger::log("Can't read PDB file from %s\n", fname);
    return 0;
  }
  fclose(fp);
  if (connects.size() > 0) {
    edges = &connects[0];
  }
  mol = new MIMoleculeBase(res, "", edges, connects.size());
  Logger::log("Loaded molecule with %d residues\n", mol->getnresidues());
  return mol;
}

MIMoleculeBase* EditEntry(MIMolOpt* opt, const char* type) {
  if (opt->IsRefining()) {
    Logger::message("Can not edit dictionary entry while refining\nCancel Refine and start again");
    return 0;
  }
  RESIDUE* dictres = opt->dict.GetDictResidue(type, 0);
  //RESIDUE * reslist = new RESIDUE;
  RESIDUE* reslist = new RESIDUE(*dictres);
  MIMoleculeBase* model = new MIMoleculeBase(reslist, "Dictionary", NULL, 0);

  //Store state of constraint prefs to restore and
  bool tmp_ca = opt->dict.GetConstrainCA();
  bool tmp_ends = opt->dict.GetConstrainEnds();

  //Don't add constraints to refinement in the editor
  opt->dict.SetConstrainCA(false);
  opt->dict.SetConstrainEnds(false);
  opt->SetRefiRes(reslist, reslist, model);
  opt->dict.SetConstrainCA(tmp_ca);
  opt->dict.SetConstrainEnds(tmp_ends);
  opt->RefiAllTorsions(reslist);
  return model;
}

MIMoleculeBase* LoadLigand(const std::string& filename,
                           MIMolOpt* opt,
                           const std::string& code) {
  std::string smi_err;                          //Error message returned from smiles library
  std::string err_report;                   //Report of smiles error to user

  // get MI data
  MIMolInfo mi;

  MIMolIOBase fio;
  int sel = fio.getReaderIndex(filename);
  if (!fio.Read(mi, filename.c_str(), sel)) {
    return 0;
  }

  LigDictEntry entry(mi.res);
  entry.bonds = mi.bonds;
  entry.angles = mi.angles;
  entry.torsions = mi.tordict;
  entry.planes = mi.planedict;
  entry.chirals = mi.chiralsdict;

  if (entry.res->natoms() == 0) {
    Logger::log("No atoms were found in the file %s. Please check the file format.", filename.c_str());
    delete mi.res;
    mi.res = 0;
    return 0;
  }

  if (code.size() > 0) {
    entry.res->type = code;
    entry.res->name = "1";
  }


  if (DupeAtomNames(entry.res) != 0) {
    entry.res->renameAtomsToUnique();
  }

  std::string wildcard("*");
  wildcard+=file_extension(filename.c_str());
  LigPostProcessor lpp(entry, wildcard.c_str());
  lpp.Process();
  entry.res->prefbonds = entry.bonds;
  entry.res->prefangles = entry.angles;

  unsigned int level;
  if (!opt->dict.DictHCheck(entry.res, level)) {
    return 0;
    // FIXME: enable code to load new dictionary!
    //     std::string path, name, ext;
    //     wxSplitPath(Application::instance()->getDictionary().c_str(), &path, &name, &ext);
    //     path += "/";
    //     if(level == DictionaryHLevel::NoHydrogens)
    //       path += "dict.noh.pdb";
    //     if(level == DictionaryHLevel::Polar)
    //       path += "dict.polarh.pdb";
    //     if(level == DictionaryHLevel::All)
    //       path += "dict.allh.pdb";
    //     opt->dict.LoadDictionary(path, false, true, level);
  }


  //This code substitutes for the "LoadDictionary" function
  opt->dict.LoadRes(entry.res, true, true);
  for (std::vector<TORSDICT>::iterator tor = entry.torsions.begin(); tor != entry.torsions.end(); ++tor) {
    opt->dict.AddTorsion(*tor);
  }
  for (std::vector<PLANEDICT>::iterator pln = entry.planes.begin(); pln != entry.planes.end(); ++pln) {
    opt->dict.AddPlane(*pln);
  }
  for (std::vector<CHIRALDICT>::iterator chrl = entry.chirals.begin(); chrl != entry.chirals.end(); ++chrl) {
    opt->dict.AddChiral(*chrl);
  }
  return EditEntry(opt, entry.res->type.c_str());
}

