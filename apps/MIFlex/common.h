#ifndef MIFLEX_COMMON_H
#define MIFLEX_COMMON_H

#include <string>

#include "chemlib.h"
#include "moloptlib.h"

bool MIGenConfs(chemlib::GeomSaver confs, chemlib::MIMoleculeBase *fitmol, MIMolOpt *opt, bool print = false);
chemlib::MIMoleculeBase *LoadMol(const char *fname);
chemlib::MIMoleculeBase *EditEntry(MIMolOpt *opt, const char *type);
chemlib::MIMoleculeBase *LoadLigand(const std::string &filename,
                                    MIMolOpt *opt,
                                    const std::string &code = "LIG");

#endif // MIFLEX_COMMON_H
