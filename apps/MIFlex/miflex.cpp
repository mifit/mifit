#include "nongui.h"
#include "chemlib.h"
#include "conflib.h"
#include "maplib.h"
#include "moloptlib.h"
#include "ligandlib.h"
#include <utillib.h>

#include "common.h"

#ifdef _WIN32
#define _MVS
#define i386
#include <mmtzlib.h>
#undef _MVS
#define strncasecmp strnicmp
#else
#include <mmtzlib.h>
#endif

using namespace chemlib;
using namespace conflib;

static void TranslateAtomsToCenter(MIMoleculeBase *mol, std::vector<MIAtom*> &atoms,
                                   float *center)
{
    float cx = center[0];
    float cy = center[1];
    float cz = center[2];
    float dx = 0, dy = 0, dz = 0;
    unsigned int i;
    for (i = 0; i < atoms.size(); i++)
    {
        dx += atoms[i]->x-cx;
        dy += atoms[i]->y-cy;
        dz += atoms[i]->z-cz;
    }
    dx /= (float)atoms.size();
    dy /= (float)atoms.size();
    dz /= (float)atoms.size();
    for (i = 0; i < atoms.size(); i++)
    {
        atoms[i]->x -= dx;
        atoms[i]->y -= dy;
        atoms[i]->z -= dz;
    }
    mol->SetCoordsChanged(true);
}

// NOTE: CurrentAtoms.size() must be non-zero before calling this
bool MIFlexLigandFit(EMapBase *currentmap,
                     MIMoleculeBase *model,
                     MIMoleculeBase *fitmol,
                     MIMolOpt *opt,
                     std::vector<MIAtom*> &CurrentAtoms,
                     float *center,
                     MIMolOptCheckPoint *ckpt = 0,
                     InterpBox *BoundingBox = 0)
{
    RESIDUE *fitres = fitmol->getResidues();
    if (!currentmap)
    {
        Logger::log("Sorry there is no map to refine against - action canceled");
        return false;
    }
    if (!currentmap->HasDensity() )
    {
        Logger::log("Sorry there is no map density to refine against "
                    "- action canceled");
        return false;
    }
    if (currentmap->mapheader->resmin > 2.5)
    {
        Logger::log("Map resolution is not high enough for the refiner - action canceled");
        return false;
    }

    if (fitres != 0 && fitres->natoms() > 60)
    {
        Logger::log("Too many atoms for the refiner (more than 60) - action canceled");
        return false;
    }

    GeomSaver confs;
    if (!MIGenConfs(confs, fitmol, opt))
    {
        Logger::log("Couldn't generate conformers");
        return false;
    }

    // refine into position
    WaitCursor wait("Refine Ligand");
    TranslateAtomsToCenter(fitmol, CurrentAtoms, center);
    if (!BoundingBox)
    {
        BoundingBox = new InterpBox(CurrentAtoms, currentmap);
    }
    BoundingBox->ZeroModel(model->getResidues());
    opt->LigandOptimize(
        CurrentAtoms, fitmol, currentmap,
        center, *BoundingBox, Refine_Level_Thorough, confs, ckpt);
    //    opt->FullOptimize(CurrentAtoms, fitmol, currentmap,
    //      viewpoint, this, *BoundingBox, Refine_Level_Thorough);
    //    opt->FullOptimize(CurrentAtoms, fitmol, currentmap,
    //      viewpoint, this, *BoundingBox, Refine_Level_Optimize);
    delete BoundingBox;
    BoundingBox = NULL;


    FILE *fil = fopen("miflex_out.pdb", "w");
    SavePDB(fil, fitres, NULL, 0, false);
    fclose(fil);

    return true;
}

bool MIFlexLigandFit(EMapBase *currentmap,
                     MIMoleculeBase *model,
                     MIMoleculeBase *fitmol,
                     MIMolOpt *opt,
                     float *center,
                     MIMolOptCheckPoint *ckpt = 0,
                     InterpBox *BoundingBox = 0)
{
    RESIDUE *fitres = fitmol->getResidues();
    std::vector<MIAtom*> CurrentAtoms;
    for (int i = 0; i < fitres->natoms(); ++i)
    {
        CurrentAtoms.push_back(fitres->atoms[i]);
    }
    return MIFlexLigandFit(currentmap, model, fitmol, opt, CurrentAtoms,
                           center, ckpt, BoundingBox);
}


class EMap
    : public EMapBase
{
public:
    bool PromptForColumnLabels(unsigned int /* num_cols */,
                               mmtz_column_ *col,
                               int &foindex, int &fcindex, int &fomindex,
                               int &phsindex, int &sigfindex, int &freeRindex)
    {
        Logger::log("Using column labels: FO=%s, FC=%s, FOM=%s, PHIC=%s, SIGF=%s, R_FREE=%s\n",
                    (foindex != -1    ? col[foindex].label    : ""),
                    (fcindex != -1    ? col[fcindex].label    : ""),
                    (fomindex != -1   ? col[fomindex].label   : ""),
                    (phsindex != -1   ? col[phsindex].label   : ""),
                    (sigfindex != -1  ? col[sigfindex].label  : ""),
                    (freeRindex != -1 ? col[freeRindex].label : ""));
        return true;
    }
};



bool MIFlexLigandFit(const std::string &model_filename,
                     const std::string &ligand_filename,
                     const std::string &map_filename,
                     float center[3])
{

    // create geomrefiner and dictionary
    MIMolOpt geomrefiner;
    char *MOLIMAGEHOME = getenv("MOLIMAGEHOME");
    if (!MOLIMAGEHOME)
    {
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

    /* set default scattering factors */
    Logger::log("Read in %d scattering factors",
                MIMapInitializeScatteringFactorTables("", MOLIMAGEHOME));


    // load molecule, model
    MIMoleculeBase *model = LoadMol(model_filename.c_str());
    MIMoleculeBase *fitmol = LoadLigand(ligand_filename.c_str(), &geomrefiner);
    if (!model || !fitmol)
    {
        Logger::log("Missing model or fitmol, aborting.");
        return false;
    }

    // get and load emap
    std::string ext = file_extension(map_filename.c_str());
    EMap *map = new EMap;

    //FIXME: support pre-loading of map header from model if non-mtz file
    //map->mapheader->set(model->GetMapHeader());

    if (ext == ".map")
    {
        map->LoadMapFile(map_filename.c_str());
    }
    else   // assume it's a phase file
    {
        //FIXME: set column names, maptype from arguments instead of hard-coded values here
        std::string maptypestr = "Direct FFT";
        std::string fo = "FWT";
        std::string fc = "";
        std::string fom = "";
        std::string phi = "PHWT";
        map->UseColumnLabels(fo, fc, fom, phi, "", "FreeR_flag");
        unsigned int maptype = MapTypeForString(maptypestr);

        map->LoadMapPhaseFile(map_filename.c_str());
        map->SFCalc(model->getResidues());
        map->FFTMap(maptype);
    }

    int retval = MIFlexLigandFit(map, model, fitmol, &geomrefiner, center);

    // clean up
    delete map;
    delete model;
    delete fitmol;

    return retval;
}

int main(int argc, char **argv)
{
    if (argc < 7)
    {
        Logger::log("Usage: %s model_file ligand_file map_file x y z", argv[0]);
        return -1;
    }

    // get center
    float center[3];
    char *err[3];
    center[0] = strtod(argv[4], &err[0]);
    center[1] = strtod(argv[5], &err[1]);
    center[2] = strtod(argv[6], &err[2]);

    if (err[0] == argv[4]
        || err[1] == argv[5]
        || err[2] == argv[6])
    {
        Logger::log("Error in center parameters\n");
        return 0;
    }

    return (int)MIFlexLigandFit(argv[1], argv[2], argv[3], center);
}

