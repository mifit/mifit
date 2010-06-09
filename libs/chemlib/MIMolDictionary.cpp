#include <vector>
#include <algorithm>

#include <nongui/nonguilib.h>
#include <math/mathlib.h>

#include <util/system.h>
#include "Residue.h"
#include <chemlib/Residue.h>
#include "MIMolDictionary.h"
#include "mol_util.h"
#include "mol_util_private.h"
#include "DictResidue.h"
#include "mmCIF.h"

namespace chemlib
{

//@{
// A temp structure to hold explicit bond information stored
// in a dictionary file
//@}
class BondDICT
{
public:
    char resname[MAXNAME];
    char aname1[MAXATOMNAME];
    char aname2[MAXATOMNAME];
    float ideal_length;
};


//@{
// A temp structure to hold explicit bond information stored
// in a dictionary file
//@}
class ANGLEDICT
{
public:
    char resname[MAXNAME];
    char aname1[MAXATOMNAME];
    char aname2[MAXATOMNAME];
    char aname3[MAXATOMNAME];
    float ideal_length;
};

bool ScanBond(std::multimap<std::string, BondDICT> &bondMap, const char *line);
bool CreateDictBonds(Residue *res, const std::multimap<std::string, BondDICT> &bondMap);
bool ScanAngle(std::multimap<std::string, ANGLEDICT> &angleMap, const char *line);
bool CreateDictAngles(Residue *res, const std::multimap<std::string, ANGLEDICT> &angleMap);

using namespace chemlib;
using namespace std;


static const float IDEAL_DIST = 3.80F;

// NOTE: for historical purposes, this file came from MIMolOpt.cpp.  Diff
// this file against a comtemporary version of MIMolOpt to see changes
//

//@{
// A synonym, another name for the same atom, dictionary entry.
//@}
class SYNONYM
{
public:
    char resname[MAXNAME];
    char aname1[MAXATOMNAME];
    char aname2[MAXATOMNAME];
};
static std::vector<SYNONYM> Synonym;

// atom_from_name, including synonyms
MIAtom *MIAtomFromNameIncludingSynonyms(const char *name, const Residue *residue)
{
    MIAtom *a = atom_from_name(name, *residue);
    if (a)
    {
        return a;
    }

    int nSynonyms = Synonym.size();
    if (nSynonyms == 0)
    {
        return 0;
    }
    for (int j = 0; j < nSynonyms; j++)
    {
        if (strcmp(residue->type().c_str(), Synonym[j].resname) == 0)
        {
            if ((strcmp(Synonym[j].aname1, name) == 0)
                || (strcmp(Synonym[j].aname2, name) == 0))
            {
                for (int i = 0; i < residue->atomCount(); i++)
                {
                    if ((strcmp(residue->atom(i)->name(), Synonym[j].aname1) == 0)
                        || (strcmp(residue->atom(i)->name(), Synonym[j].aname2) == 0))
                    {
                        return residue->atom(i);
                    }
                }
            }
        }
    }

    return 0;
}

using namespace chemlib;
using namespace std;

MIMolDictionary::MIMolDictionary()
    : cres(new Residue)
{
    RefiHBonds = false;
    RefiSecStruct = false;
    constrain_CA = true;
    constrain_Ends = true;
    modified = false;
    nResDict = 0;
    ResDict = 0;
    Alpha2 = 0;
    Beta2 = 0;
    sigmaangle = 0.03F;
    sigmabond = 0.015F;
    sigmabump = 0.15F;
    sigmaplane = 0.005F;
    sigmatorsion = 6.0F;
    CYCLEIZE = 0;
    HLevel = DictionaryHLevel::Unknown;

    cres->clearAtoms();
    cres->setType(std::string("CONST"));
    cres->setName(std::string("RAINT"));
}

bool MIMolDictionary::LoadDefaultDictionary(const std::string &dictionary,
                                            const std::string &homedir)
{
    FILE *fp = NULL;
    if (dictionary.size() != 0)
    {
        fp = fopen(dictionary.c_str(), "r");
        if (fp)
        {
            std::string path, name, ext;
            MISplitPath(dictionary, &path, &name, &ext);
            name = MIToLower(name);
            if (name.find("allh") != std::string::npos)
            {
                HLevel = DictionaryHLevel::All;
            }
            else if (name.find("polarh") != std::string::npos)
            {
                HLevel = DictionaryHLevel::All;
            }
            else if (name.find("noh") != std::string::npos)
            {
                HLevel = DictionaryHLevel::NoHydrogens;
            }
        }
        else
        {
            Logger::message("Unable to find dictionary - Insert, replace and refine functions will not work\n"
                            "You can set the dictionary using File/Preferences...");
            return false;
        }
    }
    if (fp)
    {
        if (RefiPlanes.size() > 0)
        {
            for (unsigned int i = 0; i < RefiPlanes.size(); i++)
            {
                free(RefiPlanes[i].atoms);
            }
        }
        RefiPlanes.clear();
        PlaneDict.clear();
        Synonym.clear();
        TorsDict.clear();
        ChiralDict.clear();
        TorsNames.clear();
        if (ResDict)
        {
            FreeResidueList(ResDict);
        }
        Logger::log("Loading dictionary:");
        Logger::log(dictionary.c_str());
        LoadDict(fp, false);
    }

    // load fragments
    vector<Bond> connects;
    std::string frag_dir(homedir.c_str());
    frag_dir += "/data/fragment_library";
    std::string file = frag_dir;
    file += "/beta2.pdb";
    if (fp)
    {
        fclose(fp);
        fp = NULL;
    }
    fp = fopen(file.c_str(), "r");
    if (fp)
    {
        if (Beta2)
        {
            FreeResidueList(Beta2);
        }
        Beta2 = LoadPDB(fp, &connects);
        Logger::log("Loaded beta fragment");
        fclose(fp);
        fp = NULL;
    }
    else
    {
        Logger::log("Unable to load beta fragment from %s", file.c_str());
    }
    connects.clear();
    file = frag_dir;
    file += "/alpha2.pdb";
    if (fp)
    {
        fclose(fp);
        fp = NULL;
    }
    fp = fopen(file.c_str(), "r");
    if (fp)
    {
        if (Alpha2)
        {
            FreeResidueList(Alpha2);
        }
        Alpha2 = LoadPDB(fp, &connects);
        Logger::log("Loaded alpha fragment");
        fclose(fp);
        fp = NULL;
    }
    else
    {
        Logger::log("Unable to load alpha fragment from %s", file.c_str());
    }
    return true;
}

MIMolDictionary::~MIMolDictionary()
{
    unsigned int i;
    if (ResDict)
    {
        FreeResidueList(ResDict);
    }
    if (Alpha2)
    {
        FreeResidueList(Alpha2);
    }
    if (Beta2)
    {
        FreeResidueList(Beta2);
    }
    if (RefiPlanes.size() > 0)
    {
        for (i = 0; i < RefiPlanes.size(); i++)
        {
            free(RefiPlanes[i].atoms);
        }
    }
    delete cres;
}

void MIMolDictionary::build_map()
{
    DictMap.clear();
    dict_map::iterator p;
    std::string s;
    for (Residue *res = ResDict; res != NULL; res = res->next())
    {
        if (!Monomer::isValid(res))
        {
            continue;
        }
        if (res->atomCount() > 500)
        {
            Logger::message("%d atoms found in dictionary residue %s!  Skipping residue...",
                            res->atomCount(), res->type().c_str());
            continue;
        }

        p = DictMap.insert(dict_map::value_type(res->type().c_str(), DictResidue(res)));
        //printf("Inserted residue %s into DictMap with key %s\n",
        //res->type().c_str(),p->first.c_str());
    }

    std::set<std::string> residue_set;
    for (dict_map::iterator p = DictMap.begin(); p != DictMap.end(); ++p)
        residue_set.insert(p->first);

    Logger::log("Inserted %d residues into the dictionary, with %d total conformers", residue_set.size(), DictMap.size());
}

int MIMolDictionary::CountConformers(const std::string &type) const
{
    int n = 0;
    dict_map::const_iterator p = DictMap.begin();
    while (p != DictMap.end())
    {
        if (type == p->second.residue()->type())
        {
            n++;
        }
        p++;
    }
    return n;
}

void MIMolDictionary::RestrainEnds(Residue *RefiRes, int nRefiRes)
{
    Residue *res = RefiRes;
    int n = 0, i;
    if (RefiRes && nRefiRes > 0)
    {
        /* find N terminus */
        for (i = 0; i < res->atomCount(); i++)
        {
            if (!strcmp("N", res->atom(i)->name()))
            {
                AddConstraint(res->atom(i), "0.20");
            }
        }
        /* find C terminus */
        while (Monomer::isValid(res) && n < nRefiRes-1)
        {
            n++;
            res = res->next();
        }
        if (Monomer::isValid(res))
        {
            for (i = 0; i < res->atomCount(); i++)
            {
                if (!strcmp("C", res->atom(i)->name()))
                {
                    AddConstraint(res->atom(i), "0.20");
                }
            }
        }
    }
}

int MIMolDictionary::AreBonded(const char *restype, MIAtom *atom1, MIAtom *atom2, Bond &bond)
{
    // returns 0 if not bonded, -1 if not in the dictionary but bonded by distance
    // and 1 if bonded in dictionary
    dict_map::iterator p;
    vector<Bond> *dbonds;
    float dlimit, d;
    const char *b1, *b2;
    size_t i;

    /* check to see if OK to bond */
    if (atom1->altloc() != ' ' && atom2->altloc() != ' ' && atom1->altloc() != atom2->altloc())
    {
        return 0;
    }
    p = DictMap.find(restype);
    if (p != DictMap.end())
    {
        dbonds = p->second.Bonds();
        for (i = 0; i < dbonds->size(); i++)
        {
            b1 = (*dbonds)[i].getAtom1()->name();
            b2 = (*dbonds)[i].getAtom2()->name();
            if ((strcmp(b1, atom1->name()) == 0 && strcmp(b2, atom2->name()) == 0)
                || (strcmp(b1, atom2->name()) == 0 && strcmp(b2, atom1->name()) == 0))
            {
                // check to see if deleted by checking for a negative tolerance
                d = (*dbonds)[i].tolerance;
                if (d < 0.0)
                {
                    return 0;
                }
                else
                {
                    bond = (*dbonds)[i];                //Copy the bond information
                    return 1;
                }
            }
        }
    }
    else
    {
        // check to see if distance bondable
        dlimit = (float)(BondLimit(atom1->name())+ BondLimit(atom2->name()));
        if ((d = (float)AtomDist(*atom1, *atom2)) < dlimit)
        {
            return -1;
        }
    }
    return 0;
}

/*
 * build an array of bonds and angles for a list of
 * residues given a dictionary of guide residues with perfect geometry
 * as an example to follow
 */
int MIMolDictionary::FindGeom(Residue *reslist, int nres, Residue *ResActiveModel)
{
    if (EmptyDictCheck() == false)
    {
        return 0;
    }
    ANGLE angle;
    Bond bond;
    PLANE plane;
    int i, j;
    MIAtom *a1, *a2, *a3;
    MIAtom *b1, *b2, *b3;
    Residue *res2;
    Residue *prev = NULL, *next, *start;
    TORSION *chi = NULL;
    std::string mess;
    int n = 0, na, it;
    static char TorsionNames[][MAXNAME] = { "CHI1", "CHI2", "CHI3", "CHI4",
                                            "CHI5", "CHI6", "CHI7", "CHI8", "CHI9", "CHI10", "CHI11", "CHI12",
                                            "IMP1", "IMP2", "IMP3", "IMP4", "IMP5", "IMP6", "IMP7",
                                            "IMP8", "IMP9", "IMP10", "IMP11", "IMP12", "IMP13", "IMP14", "IMP15",
                                            "IMP16", "IMP17", "IMP18", "IMP19", "IMP20", "IMP21", "IMP22", "IMP23",
                                            "IMP24", "IMP25", "IMP26", "IMP27", "IMP28", "IMP29", "IMP30", "IMP31",
                                            "IMP32", "IMP33", "IMP34", "IMP35", "IMP36", "IMP37", "IMP38", "IMP39",
                                            "IMP40", "IMP41", "IMP42", "IMP43", "IMP44", "IMP45", "IMP46", "IMP47"};

    Clear();
    start = reslist;
    dict_map::iterator p;
    vector<Bond> *dbonds;
    vector<ANGLE> *dangles;
    bond.setOrder(NORMALBOND); //Just to make sure we have it initilized.
    while (reslist != NULL && n < nres)
    {
        p = DictMap.find(reslist->type().c_str());
        if (p == DictMap.end())
        {
            Logger::log("residue type %s not found in dictionary",
                        reslist->type().c_str());
        }
        else
        {
            dbonds = p->second.Bonds();
            dangles = p->second.Angles();

            next = reslist->next();
            /* no next if at last  residue */
            if (n == nres-1)
            {
                next = NULL;
            }

            /* get bonds */
            for (i = 0; (unsigned int)i < (*dbonds).size(); i++)
            {
                a1 = (*dbonds)[i].getAtom1();
                a2 = (*dbonds)[i].getAtom2();
                bond.setAtom1(MIAtomFromNameIncludingSynonyms(a1->name(), reslist));
                bond.setAtom2(MIAtomFromNameIncludingSynonyms(a2->name(), reslist));
                if (bond.getAtom1() == NULL || bond.getAtom2() == NULL)
                {
                    /*
                       printf("Can't find bond %s %s in %s %s\n",a1->name,a2->name,reslist->type,reslist->name,res->type());
                     */
                }
                else
                {
                    bond.ideal_length = (*dbonds)[i].ideal_length;
                    bond.tolerance = sigmabond;
                    RefiBonds.push_back(bond);
                }
            }
            /* get angles */
            for (i = 0; (unsigned int)i < (*dangles).size(); i++)
            {
                a1 = (*dangles)[i].getAtom1();
                a2 = (*dangles)[i].getAtom2();
                a3 = (*dangles)[i].atom3;
                angle.setAtom1(MIAtomFromNameIncludingSynonyms(a1->name(), reslist));
                angle.setAtom2(MIAtomFromNameIncludingSynonyms(a2->name(), reslist));
                angle.atom3 = MIAtomFromNameIncludingSynonyms(a3->name(), reslist);
                if (angle.getAtom1() == NULL || angle.getAtom2() == NULL || angle.atom3 == NULL)
                {
                    /*
                       printf("Can't find angle %s %s in %s %s\n",a1->name,a3->name,reslist->type,reslist->name,res->type());
                     */
                }
                else
                {
                    angle.ideal_angle = (*dangles)[i].ideal_angle;
                    angle.tolerance = sigmaangle;
                    angle.res = reslist;
                    RefiAngles.push_back(angle);
                }
            }
            /* riding hydrogens due to non-standard shelx names, we need to check these hydrogens explicitly*/
            if ((strcmp(reslist->type().c_str(), "VAL") == 0 || strcmp(reslist->type().c_str(), "ILE") == 0) && MIAtomFromNameIncludingSynonyms("HG1", reslist))
            {
                MIAtom *CG1 = MIAtomFromNameIncludingSynonyms("CG1", reslist);
                MIAtom *CG2 = MIAtomFromNameIncludingSynonyms("CG2", reslist);
                MIAtom *CD1 = MIAtomFromNameIncludingSynonyms("CD1", reslist);
                MIAtom *CB = MIAtomFromNameIncludingSynonyms("CB", reslist);
                MIAtom *HG11 = 0, *HG12 = 0, *HG13 = 0;
                MIAtom *HG21 = 0, *HG22 = 0, *HG23 = 0;
                MIAtom *HD11 = 0, *HD12 = 0, *HD13 = 0;
                if (CG1 && CG2)
                {
                    for (i = 0; i < reslist->atomCount(); i++)
                    {
                        MIAtom *a = reslist->atom(i);
                        // HG1 is ambiguous - there are three on CG1 and one on CG2, use distance to determine
                        if (strcmp(a->name(), "HG1") == 0)
                        {
                            if (HG11 == 0)
                            {
                                HG11 = a;
                            }
                            else if (HG12 == 0)
                            {
                                HG12 = a;
                            }
                            else
                            {
                                HG13 = a;
                            }
                        }
                        if (strcmp(a->name(), "HG2") == 0)
                        {
                            if (HG21 == 0)
                            {
                                HG21 = a;
                            }
                            else if (HG22 == 0)
                            {
                                HG22 = a;
                            }
                            else
                            {
                                HG23 = a;
                            }
                        }
                        if (strcmp(a->name(), "HD1") == 0)
                        {
                            if (HD11 == 0)
                            {
                                HD11 = a;
                            }
                            else if (HD12 == 0)
                            {
                                HD12 = a;
                            }
                            else
                            {
                                HD13 = a;
                            }
                        }
                    }
                    if (HD11 && CD1)
                    {
                        bond.setAtom1(HD11);
                        bond.setAtom2(CD1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD12 && CD1)
                    {
                        bond.setAtom1(HD12);
                        bond.setAtom2(CD1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD13 && CD1)
                    {
                        bond.setAtom1(HD13);
                        bond.setAtom2(CD1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }

                    if (HD11 && CG1)
                    {
                        bond.setAtom1(HD11);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD12 && CG1)
                    {
                        bond.setAtom1(HD12);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD13 && CG1)
                    {
                        bond.setAtom1(HD13);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }

                    if (HD11 && HD12)
                    {
                        bond.setAtom1(HD11);
                        bond.setAtom2(HD12);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD11 && HD13)
                    {
                        bond.setAtom1(HD11);
                        bond.setAtom2(HD13);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD13 && HD12)
                    {
                        bond.setAtom1(HD13);
                        bond.setAtom2(HD12);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG11 && CG1)
                    {
                        bond.setAtom1(HG11);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG12 && CG1)
                    {
                        bond.setAtom1(HG12);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG13 && CG1)
                    {
                        bond.setAtom1(HG13);
                        bond.setAtom2(CG1);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG11 && HG13)
                    {
                        bond.setAtom1(HG11);
                        bond.setAtom2(HG13);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG11 && HG12)
                    {
                        bond.setAtom1(HG11);
                        bond.setAtom2(HG12);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG12 && HG13)
                    {
                        bond.setAtom1(HG12);
                        bond.setAtom2(HG13);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG21 && CG2)
                    {
                        bond.setAtom1(HG21);
                        bond.setAtom2(CG2);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG22 && CG2)
                    {
                        bond.setAtom1(HG22);
                        bond.setAtom2(CG2);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG23 && CG2)
                    {
                        bond.setAtom1(HG23);
                        bond.setAtom2(CG2);
                        bond.ideal_length = 0.960F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG21 && HG23)
                    {
                        bond.setAtom1(HG21);
                        bond.setAtom2(HG23);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG21 && HG22)
                    {
                        bond.setAtom1(HG21);
                        bond.setAtom2(HG22);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HG22 && HG23)
                    {
                        bond.setAtom1(HG22);
                        bond.setAtom2(HG23);
                        bond.ideal_length = 1.576F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG11)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG11);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG12)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG12);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG13)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG13);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG21)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG21);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG22)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG22);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (CB && HG23)
                    {
                        bond.setAtom1(CB);
                        bond.setAtom2(HG23);
                        bond.ideal_length = 2.063F;
                        bond.tolerance = sigmaangle /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                }
            }

            if (strcmp(reslist->type().c_str(), "GLN") == 0 && MIAtomFromNameIncludingSynonyms("HE2", reslist))
            {
                MIAtom *NE2 = MIAtomFromNameIncludingSynonyms("NE2", reslist);
                MIAtom *CD = MIAtomFromNameIncludingSynonyms("CD", reslist);
                MIAtom *OE1 = MIAtomFromNameIncludingSynonyms("OE1", reslist);
                MIAtom *CG = MIAtomFromNameIncludingSynonyms("CG", reslist);
                MIAtom *HE21 = 0, *HE22 = 0;
                if (NE2)
                {
                    for (i = 0; i < reslist->atomCount(); i++)
                    {
                        MIAtom *a = reslist->atom(i);
                        if (strcmp(a->name(), "HE2") == 0)
                        {
                            if (HE21 == 0)
                            {
                                HE21 = a;
                            }
                            else
                            {
                                HE22 = a;
                            }
                        }
                    }
                    if (HE21 && NE2)
                    {
                        bond.setAtom1(HE21);
                        bond.setAtom2(NE2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HE22 && NE2)
                    {
                        bond.setAtom1(HE22);
                        bond.setAtom2(NE2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HE22 && CD)
                    {
                        bond.setAtom1(HE22);
                        bond.setAtom2(CD);
                        bond.ideal_length = 1.912F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HE21 && CD)
                    {
                        bond.setAtom1(HE21);
                        bond.setAtom2(CD);
                        bond.ideal_length = 1.912F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HE21 && HE22)
                    {
                        bond.setAtom1(HE21);
                        bond.setAtom2(HE22);
                        bond.ideal_length = 1.490F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HE21 && HE22 && NE2 && CD && CG && OE1)
                    {
                        plane.natoms = 6;
                        if ((plane.atoms = (MIAtom**)malloc(plane.natoms*sizeof(MIAtom*))) != NULL)
                        {
                            plane.atoms[0] = HE21;
                            plane.atoms[1] = HE22;
                            plane.atoms[2] = NE2;
                            plane.atoms[3] = OE1;
                            plane.atoms[4] = CD;
                            plane.atoms[5] = CG;
                            plane.tolerance = sigmaplane;
                            plane.res = reslist;
                            RefiPlanes.push_back(plane);
                        }
                    }
                }
            }

            if (strcmp(reslist->type().c_str(), "ASN") == 0 && MIAtomFromNameIncludingSynonyms("HD2", reslist))
            {
                MIAtom *ND2 = MIAtomFromNameIncludingSynonyms("ND2", reslist);
                MIAtom *CB = MIAtomFromNameIncludingSynonyms("CB", reslist);
                MIAtom *OD1 = MIAtomFromNameIncludingSynonyms("OD1", reslist);
                MIAtom *CG = MIAtomFromNameIncludingSynonyms("CG", reslist);
                MIAtom *HD21 = 0, *HD22 = 0;
                if (ND2)
                {
                    for (i = 0; i < reslist->atomCount(); i++)
                    {
                        MIAtom *a = reslist->atom(i);
                        if (strcmp(a->name(), "HD2") == 0)
                        {
                            if (HD21 == 0)
                            {
                                HD21 = a;
                            }
                            else
                            {
                                HD22 = a;
                            }
                        }
                    }
                    if (HD21 && ND2)
                    {
                        bond.setAtom1(HD21);
                        bond.setAtom2(ND2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD22 && ND2)
                    {
                        bond.setAtom1(HD22);
                        bond.setAtom2(ND2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD22 && CG)
                    {
                        bond.setAtom1(HD22);
                        bond.setAtom2(CG);
                        bond.ideal_length = 1.912F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD21 && CG)
                    {
                        bond.setAtom1(HD21);
                        bond.setAtom2(CG);
                        bond.ideal_length = 1.912F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD21 && HD22)
                    {
                        bond.setAtom1(HD21);
                        bond.setAtom2(HD22);
                        bond.ideal_length = 1.490F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HD21 && HD22 && ND2 && CB && CG && OD1)
                    {
                        plane.natoms = 6;
                        if ((plane.atoms = (MIAtom**)malloc(plane.natoms*sizeof(MIAtom*))) != NULL)
                        {
                            plane.atoms[0] = HD21;
                            plane.atoms[1] = HD22;
                            plane.atoms[2] = ND2;
                            plane.atoms[3] = OD1;
                            plane.atoms[4] = CB;
                            plane.atoms[5] = CG;
                            plane.tolerance = sigmaplane;
                            plane.res = reslist;
                            RefiPlanes.push_back(plane);
                        }
                    }
                }
            }

            if (strcmp(reslist->type().c_str(), "ARG") == 0 && MIAtomFromNameIncludingSynonyms("HH1", reslist))
            {
                MIAtom *NH1 = MIAtomFromNameIncludingSynonyms("NH1", reslist);
                MIAtom *NH2 = MIAtomFromNameIncludingSynonyms("NH2", reslist);
                MIAtom *NE = MIAtomFromNameIncludingSynonyms("NE", reslist);
                MIAtom *CZ = MIAtomFromNameIncludingSynonyms("CZ", reslist);
                MIAtom *CD = MIAtomFromNameIncludingSynonyms("CD", reslist);
                MIAtom *HH11 = 0, *HH12 = 0;
                MIAtom *HH21 = 0, *HH22 = 0;
                if (NH1 && NH2)
                {
                    for (i = 0; i < reslist->atomCount(); i++)
                    {
                        MIAtom *a = reslist->atom(i);
                        if (strcmp(a->name(), "HH1") == 0)
                        {
                            if (HH11 == 0)
                            {
                                HH11 = a;
                            }
                            else
                            {
                                HH12 = a;
                            }
                        }
                        if (strcmp(a->name(), "HH2") == 0)
                        {
                            if (HH21 == 0)
                            {
                                HH21 = a;
                            }
                            else
                            {
                                HH22 = a;
                            }
                        }
                    }
                    if (HH11 && NH1)
                    {
                        bond.setAtom1(HH11);
                        bond.setAtom2(NH1);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH12 && NH1)
                    {
                        bond.setAtom1(HH12);
                        bond.setAtom2(NH1);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH21 && NH2)
                    {
                        bond.setAtom1(HH21);
                        bond.setAtom2(NH2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH22 && NH2)
                    {
                        bond.setAtom1(HH22);
                        bond.setAtom2(NH2);
                        bond.ideal_length = 0.860F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH11 && CZ)
                    {
                        bond.setAtom1(HH11);
                        bond.setAtom2(CZ);
                        bond.ideal_length = 1.908F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH12 && CZ)
                    {
                        bond.setAtom1(HH12);
                        bond.setAtom2(CZ);
                        bond.ideal_length = 1.908F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH21 && CZ)
                    {
                        bond.setAtom1(HH21);
                        bond.setAtom2(CZ);
                        bond.ideal_length = 1.908F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH22 && CZ)
                    {
                        bond.setAtom1(HH22);
                        bond.setAtom2(CZ);
                        bond.ideal_length = 1.908F;
                        bond.tolerance = sigmabond /*/20*/;
                        RefiBonds.push_back(bond);
                    }
                    if (HH11 && HH12)
                    {
                        bond.setAtom1(HH11);
                        bond.setAtom2(HH12);
                        bond.ideal_length = 1.490F;
                        bond.tolerance = sigmabond/20;
                        RefiBonds.push_back(bond);
                    }
                    if (HH21 && HH22)
                    {
                        bond.setAtom1(HH21);
                        bond.setAtom2(HH22);
                        bond.ideal_length = 1.490F;
                        bond.tolerance = sigmabond/20;
                        RefiBonds.push_back(bond);
                    }
                    if (HH11 && HH12 && HH21 && HH22 && NH2 && NH1 && CZ && NE && CD)
                    {
                        plane.natoms = 9;
                        if ((plane.atoms = (MIAtom**)malloc(plane.natoms*sizeof(MIAtom*))) != NULL)
                        {
                            plane.atoms[0] = HH11;
                            plane.atoms[1] = HH12;
                            plane.atoms[2] = HH21;
                            plane.atoms[3] = HH22;
                            plane.atoms[4] = NH2;
                            plane.atoms[5] = NH1;
                            plane.atoms[6] = CZ;
                            plane.atoms[7] = NE;
                            plane.atoms[8] = CD;
                            plane.tolerance = sigmaplane;
                            plane.res = reslist;
                            RefiPlanes.push_back(plane);
                        }
                    }
                }
            }

            /* find planes */
            for (i = 0; (unsigned int)i < PlaneDict.size(); i++)
            {
                if (strcmp(PlaneDict[i].restype, reslist->type().c_str()))
                {
                    continue;
                }
                /* found a match */
                if ((plane.atoms = (MIAtom**)malloc(PlaneDict[i].natoms*sizeof(MIAtom*))) == NULL)
                {
                    continue;
                }
                na = 0;
                for (j = 0; j < PlaneDict[i].natoms; j++)
                {
                    if ((plane.atoms[na] = MIAtomFromNameIncludingSynonyms(PlaneDict[i].name[j], reslist)) != NULL)
                    {
                        na++;
                    }
                }
                if (na < 3)
                {
                    free(plane.atoms);
                    continue;
                }
                plane.natoms = na;
                plane.tolerance = sigmaplane;
                plane.res = reslist;
                RefiPlanes.push_back(plane);
            }
            CHIRAL chiral;
            for (unsigned int i = 0; i < ChiralDict.size(); ++i)
            {
                if (strcmp(ChiralDict[i].restype, reslist->type().c_str()))
                {
                    continue;
                }
                chiral.center = MIAtomFromNameIncludingSynonyms(ChiralDict[i].center, reslist);
                chiral.setAtom1(MIAtomFromNameIncludingSynonyms(ChiralDict[i].name[0], reslist));
                chiral.setAtom2(MIAtomFromNameIncludingSynonyms(ChiralDict[i].name[1], reslist));
                chiral.atom3 = MIAtomFromNameIncludingSynonyms(ChiralDict[i].name[2], reslist);
                chiral.order = ChiralDict[i].order;
                chiral.flags = CHIRAL_AUTO;

                if (chiral.center != 0
                    && chiral.getAtom1() != 0
                    && chiral.getAtom2() != 0
                    && chiral.atom3 != 0)
                {
                    RefiChirals.push_back(chiral);
                }
            }

            /* peptide angles  */
            if ((next != NULL && !(reslist->linkage_type()&CTERMINUS)
                 && (reslist->linkage_type()&PEPTIDE)
                 && (next->linkage_type()&PEPTIDE))
                || (next == NULL && (reslist->linkage_type()&CTERMINUS)
                    && CYCLEIZE == 1))
            {
                a1 = a2 = a3 = b1 = b2 = b3 = NULL;
                if (next == NULL)
                {
                    next = start;
                }
                /* satisfied peptide criteria */
                a1 = MIAtomFromNameIncludingSynonyms("C", reslist);
                a2 = MIAtomFromNameIncludingSynonyms("O", reslist);
                a3 = MIAtomFromNameIncludingSynonyms("CA", reslist);
                b1 = MIAtomFromNameIncludingSynonyms("N", next);
                b2 = MIAtomFromNameIncludingSynonyms("CA", next);
                b3 = MIAtomFromNameIncludingSynonyms("H", next);
                if (!b3)
                {
                    b3 = MIAtomFromNameIncludingSynonyms("HN", next);
                }
                if (!b3)
                {
                    b3 = MIAtomFromNameIncludingSynonyms("H0", next);
                }
                /* break the chain if greater than 4.5 A between
                 * O and N atoms
                 */
                if (a1 && a2 && a3 && b1 && b2 && (AtomDist(*a1, *b1) < 4.5))
                {
                    na = 5;
                    if (b3)
                    {
                        na = 6;
                    }
                    if ((plane.atoms = (MIAtom**)malloc(na*sizeof(MIAtom*))) != NULL)
                    {
                        plane.natoms = na;
                        plane.tolerance = sigmaplane *5.0f;
                        plane.atoms[0] = a1;
                        plane.atoms[1] = a2;
                        plane.atoms[2] = a3;
                        plane.atoms[3] = b1;
                        plane.atoms[4] = b2;
                        if (b3)
                        {
                            plane.atoms[5] = b3;
                        }
                        plane.res = reslist;
                        RefiPlanes.push_back(plane);
                    }
                    /* add O-C-N angle */
                    angle.setAtom1(a2);
                    angle.setAtom2(a1);
                    angle.atom3 = b1;
                    angle.ideal_angle = 2.27F;
                    angle.tolerance = sigmaangle;
                    angle.res = reslist;
                    RefiAngles.push_back(angle);
                    /* add CA-C-N */
                    angle.setAtom1(a3);
                    angle.setAtom2(a1);
                    angle.atom3 = b1;
                    angle.ideal_angle = 2.39F;

                    angle.tolerance = sigmaangle;
                    angle.res = reslist;
                    RefiAngles.push_back(angle);
                    /* add CA'-N-C */
                    angle.setAtom1(b2);
                    angle.setAtom2(b1);
                    angle.atom3 = a1;
                    angle.ideal_angle = 2.45F;
                    angle.tolerance = sigmaangle;
                    angle.res = reslist;
                    RefiAngles.push_back(angle);
                    /* C - N - H */
                    if (b3)
                    {
                        angle.setAtom1(b3);
                        angle.setAtom2(b1);
                        angle.atom3 = a1;
                        angle.ideal_angle = 1.93F;
                        angle.tolerance = sigmaangle;
                        angle.res = reslist;
                        RefiAngles.push_back(angle);
                        TORSION t;
                        // H-N-C-CA
                        t.setAtom1(b3);
                        t.setAtom2(b1);
                        t.atom3 = a1;
                        t.atom4 = a3;
                        t.ideal[0] = 0.0;
                        t.nideal = 1;
                        t.tolerance = 5.0F;
                        t.res = reslist;
                        strcpy(t.type, "HNCCA");
                        RefiTorsions.push_back(t);
                        // H-N-C-O
                        t.setAtom1(b3);
                        t.setAtom2(b1);
                        t.atom3 = a1;
                        t.atom4 = a2;
                        t.ideal[0] = 180.0;
                        t.nideal = 1;
                        t.tolerance = 5.0F;
                        t.res = reslist;
                        strcpy(t.type, "HNCO");
                        RefiTorsions.push_back(t);
                    }
                    /* peptide bond C-N */
                    bond.setAtom1(a1);
                    bond.setAtom2(b1);
                    bond.ideal_length = 1.32F;
                    bond.tolerance = sigmabond;
                    RefiBonds.push_back(bond);
                    /* link CA's  */
                    bond.setAtom1(a3);
                    bond.setAtom2(b2);
                    bond.ideal_length = IDEAL_DIST;
                    bond.tolerance = 3.0F*sigmabond;
                    RefiBonds.push_back(bond);
                }
            }
            /* if just a CA atom link it to before and aft */
            if (next != NULL
                && (reslist->atomCount() == 1 || next->atomCount() == 1)
                && (a1 = MIAtomFromNameIncludingSynonyms("CA", reslist)) != NULL
                && (b1 = MIAtomFromNameIncludingSynonyms("CA", next)) != NULL
                && AtomDist(*a1, *b1) < 4.5)
            {
                /* link CA's  */
                bond.setAtom1(a1);
                bond.setAtom2(b1);
                bond.ideal_length = IDEAL_DIST;
                bond.tolerance = 5.0F*sigmabond;
                RefiBonds.push_back(bond);
            }
            /* add phi-psi and omega torsion */
            if (next != NULL && prev != NULL
                && !(reslist->linkage_type()&NTERMINUS)
                && !(reslist->linkage_type()&CTERMINUS)
                && (reslist->linkage_type()&PEPTIDE)
                && (next->linkage_type()&PEPTIDE)
                && (prev->linkage_type()&PEPTIDE) )
            {
                TORSION *phi = getTORSION(reslist, "PHI", prev);
                TORSION *psi = getTORSION(reslist, "PSI", prev);
                /* use omega of next residue */
                TORSION *omega = getTORSION(next, "OMEGA", reslist);
                a1 = MIAtomFromNameIncludingSynonyms("O", reslist);
                if (phi && psi && omega && a1)
                {
                    strcpy(phi->type, "PHI");
                    RefiPhiPsis.push_back(*phi);
                    strcpy(psi->type, "PSI");
                    RefiPhiPsis.push_back(*psi);
                    // added both values to allow cis and trans planes
                    omega->ideal[0] = 0.0;
                    omega->ideal[1] = 180.0;
                    omega->nideal = 2;
                    strcpy(omega->type, "OMEGA");
                    RefiPhiPsis.push_back(*omega);
                    //omega->ideal[0] = 0.0;
                    omega->setAtom1(a1);
                    strcpy(omega->type, "OMEGA\'");
                    RefiPhiPsis.push_back(*omega);
                }
                if (psi)
                {
                    free(psi);
                }
                if (phi)
                {
                    free(phi);
                }
                if (omega)
                {
                    free(omega);
                }
            }
            /* chi1 torsions and others */
            for (it = 0; (unsigned int)it < (sizeof(TorsionNames)/MAXNAME); it++)
            {
                if ((chi = getTORSION(reslist, TorsionNames[it])) != NULL)
                {
                    if (chi->nideal > 0)
                    {
                        strcpy(chi->type, TorsionNames[it]);
                        RefiTorsions.push_back(*chi);
                    }
                }
                if (chi)
                {
                    free(chi);
                    chi = NULL;
                }
            }
            if (RefiSecStruct)
            {
                if (reslist->secstr()  == 'H' || reslist->secstr()  == 'S')
                {
                    TORSION *phi = getTORSION(reslist, "PHI", prev);
                    if (phi)
                    {
                        phi->nideal = 1;
                        if (reslist->secstr()  == 'H')
                        {
                            phi->ideal[0] = 300.0;
                            strcpy(phi->type, "PHI-H");
                        }
                        else
                        {
                            phi->ideal[0] = 260.0;
                            strcpy(phi->type, "PHI-S");
                        }
                        RefiTorsions.push_back(*phi);
                        free(phi);
                    }
                    TORSION *psi = getTORSION(reslist, "PSI", prev);
                    if (psi)
                    {
                        psi->nideal = 1;
                        if (reslist->secstr()  == 'H')
                        {
                            psi->ideal[0] = 315.0;
                            strcpy(psi->type, "PHI-H");
                        }
                        else
                        {
                            psi->ideal[0] = 135.0;
                            strcpy(psi->type, "PHI-S");
                        }
                        RefiTorsions.push_back(*psi);
                        free(psi);
                    }
                }
            }
            /* main-chain hydrogen bonds */
            if (RefiHBonds)
            {
                if (reslist->linkage_type()&PEPTIDE)
                {
                    a1 = MIAtomFromNameIncludingSynonyms("O", reslist);
                    a2 = MIAtomFromNameIncludingSynonyms("N", reslist);
                    if (a1)
                    {
                        /// search thru model for N atoms
                        res2 = ResActiveModel;
                        while (res2 != NULL)
                        {
                            if ((res2 != reslist) && (res2->linkage_type()&PEPTIDE))
                            {
                                if ((b1 = MIAtomFromNameIncludingSynonyms("N", res2)) != NULL)
                                {
                                    if (hbondable(*a1, *b1, *reslist, *res2))
                                    {
                                        AddConstraint(a1, b1, "2.85", "0.3");
                                    }
                                }
                            }
                            res2 = res2->next();
                        }
                    }
                    if (a2)
                    {
                        // search thru model for O atoms
                        res2 = ResActiveModel;
                        while (res2 != NULL)
                        {
                            if (res2 != reslist && res2->linkage_type()&PEPTIDE)
                            {
                                if ((b2 = MIAtomFromNameIncludingSynonyms("O", res2)) != NULL)
                                {
                                    if (hbondable(*a2, *b2, *reslist, *res2))
                                    {
                                        AddConstraint(a2, b2, "2.85", "0.3");
                                    }
                                }
                            }
                            res2 = res2->next();
                        }
                    }
                }
            }
        }

        // termini
        {
            MIAtom *OT1 = MIAtomFromNameIncludingSynonyms("OT1", reslist);
            MIAtom *OT2 = MIAtomFromNameIncludingSynonyms("OT2", reslist);
            MIAtom *CA = MIAtomFromNameIncludingSynonyms("CA", reslist);
            MIAtom *C = MIAtomFromNameIncludingSynonyms("C", reslist);
            if (OT1 && OT2 && CA && C)
            {
                bond.setAtom1(OT1);
                bond.setAtom2(C);
                bond.ideal_length = 1.251F;
                bond.tolerance = sigmabond;
                RefiBonds.push_back(bond);
                bond.setAtom1(OT2);
                bond.setAtom2(C);
                bond.ideal_length = 1.251F;
                bond.tolerance = sigmabond;
                RefiBonds.push_back(bond);
                bond.setAtom1(OT2);
                bond.setAtom2(OT1);
                bond.ideal_length = 2.201F;
                bond.tolerance = sigmaangle;
                RefiBonds.push_back(bond);
                bond.setAtom1(OT1);
                bond.setAtom2(CA);
                bond.ideal_length = 2.367F;
                bond.tolerance = sigmaangle;
                RefiBonds.push_back(bond);
                bond.setAtom1(OT2);
                bond.setAtom2(CA);
                bond.ideal_length = 2.367F;
                bond.tolerance = sigmaangle;
                RefiBonds.push_back(bond);
                plane.natoms = 4;
                if ((plane.atoms = (MIAtom**)malloc(plane.natoms*sizeof(MIAtom*))) != NULL)
                {
                    plane.atoms[0] = OT1;
                    plane.atoms[1] = OT2;
                    plane.atoms[2] = C;
                    plane.atoms[3] = CA;
                    plane.tolerance = sigmaplane;
                    plane.res = reslist;
                    RefiPlanes.push_back(plane);
                }
            }
            MIAtom *H0A = MIAtomFromNameIncludingSynonyms("H0A", reslist);
            MIAtom *H0B = MIAtomFromNameIncludingSynonyms("H0B", reslist);
            MIAtom *H0C = MIAtomFromNameIncludingSynonyms("H0C", reslist);
            MIAtom *N = MIAtomFromNameIncludingSynonyms("N", reslist);
            CA = MIAtomFromNameIncludingSynonyms("CA", reslist);
            if (H0A && H0B && H0C && CA && N)
            {
                bond.setAtom1(H0A);
                bond.setAtom2(N);
                bond.ideal_length = 0.89F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0B);
                bond.setAtom2(N);
                bond.ideal_length = 0.89F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0C);
                bond.setAtom2(N);
                bond.ideal_length = 0.89F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0A);
                bond.setAtom2(H0B);
                bond.ideal_length = 1.454F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0A);
                bond.setAtom2(H0C);
                bond.ideal_length = 1.454F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0C);
                bond.setAtom2(H0B);
                bond.ideal_length = 1.454F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0A);
                bond.setAtom2(CA);
                bond.ideal_length = 1.949F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0B);
                bond.setAtom2(CA);
                bond.ideal_length = 1.949F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
                bond.setAtom1(H0C);
                bond.setAtom2(CA);
                bond.ideal_length = 1.949F;
                bond.tolerance = sigmabond/10;
                RefiBonds.push_back(bond);
            }
        }
        n++;
        prev = reslist;

        reslist = reslist->next();
    }
    return 1;
}

void MIMolDictionary::ConstrainCalpha(Residue *RefiRes, int nRefiRes)
{
    Residue *res = RefiRes;
    int n = 0, i;
    if (RefiRes && nRefiRes > 0)
    {
        while (Monomer::isValid(res) && n < nRefiRes)
        {
            for (i = 0; i < res->atomCount(); i++)
            {
                if (!strcmp("CA", res->atom(i)->name()))
                {
                    AddConstraint(res->atom(i), "0.7");
                }
            }
            n++;
            res = res->next();
        }
    }
}

int MIMolDictionary::AddConstraint(MIAtom *a1, const char *sigma)
{
    float tol, d;
    Bond bond;

    if (a1 == NULL)
    {
        return (0);
    }

    MIAtom *a2 = new MIAtom;
    a2->copyShallow(*a1);
    a2->setName("SAME");
    a2->addType(AtomType::FREEZE | AtomType::DUMMYATOM);
    d = 0.0001F;

    bond.setAtom1(a1);
    bond.setAtom2(a2);
    tol = (float)atof(sigma);
    if (tol <= .0001f)
    {
        tol = .0001f;
    }
    if (d <= .0001f)
    {
        d = .0001f;
    }
    bond.tolerance = tol;
    bond.ideal_length = d;
    RefiConstraints.push_back(bond);
    return 1;
}

int MIMolDictionary::AddConstraint(MIAtom *a1, MIAtom *a2, const char *target, const char *sigma)
{
    float tol, d;
    Bond bond;

    if (a1 == NULL)
    {
        return (0);
    }
    if (!strncmp(target, "As Is", 5))
    {
        d = (float)AtomDist(*a1, *a2);
    }
    else
    {
        d = (float)atof(target);
    }
    bond.setAtom1(a1);
    bond.setAtom2(a2);
    tol = (float)atof(sigma);
    if (tol <= .0001F)
    {
        tol = .0001F;
    }
    if (d <= .0001F)
    {
        d = .0001F;
    }
    bond.tolerance = tol;
    bond.ideal_length = d;
    RefiConstraints.push_back(bond);
    return 1;
}

int MIMolDictionary::BuildBumps(Residue *RefiRes, int nRefiRes)
{
    Residue *res, *res2;
    MIAtom *a1, *a2;
    Bond bond;
    int n2, n;
    int i, j, k;
    int found;
    float d;
    int nres = nRefiRes;
    Residue *reslist = RefiRes;
    res = reslist;

    // build some maps first to keep from having to search the entire RefiBonds/ RefiAtoms each time!
    std::map<MIAtom*, std::set<MIAtom*> > bond_map;
    for (k = 0; (unsigned int)k < RefiBonds.size(); k++)
    {
        a1 = RefiBonds[k].getAtom1();
        a2 = RefiBonds[k].getAtom2();
        bond_map[a1].insert(a2);
        bond_map[a2].insert(a1);
    }

    std::map<MIAtom*, std::set<MIAtom*> > angle_map;
    for (k = 0; (unsigned int)k < RefiAngles.size(); k++)
    {
        a1 = RefiAngles[k].getAtom1();
        a2 = RefiAngles[k].atom3;
        angle_map[a1].insert(a2);
        angle_map[a2].insert(a1);
    }

    /* loop thru every atom in the reslist and add a bump if
     * reasonably close */
    n = 0;
    RefiBumps.clear();
    while (Monomer::isValid(res) && n < nres)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            res2 = reslist;
            n2 = 0;
            a1 = res->atom(i);
            while (res2 != NULL && n2 < nres)
            {
                for (j = 0; j < res2->atomCount(); j++)
                {
                    a2 = res2->atom(j);
                    if (a1 < a2
                        && !MIAtom::MIIsHydrogen(a1) && !MIAtom::MIIsHydrogen(a2)
                        && (d = (float)AtomDist(*a1, *res2->atom(j))) < 4.3f)
                    {
                        /* are they bonded ? */
                        found = 0;
                        if (res == res2 && d < 2.95f && !MIAtom::MIIsMainChainAtom(a1) && !MIAtom::MIIsMainChainAtom(a2))
                        {
                            found = 1;
                        }
                        if (res->next())
                        {
                            if (res->next() == res2 && d < 3.5f)
                            {
                                if (MIAtom::MIIsMainChainAtom(a1) && MIAtom::MIIsMainChainAtom(a2))
                                {
                                    found = 1;
                                }
                            }
                        }

                        if (res == res2 && a1->altloc() != a2->altloc())
                        {
                            found = 1;
                        }

                        if (found == 0)
                        {
                            std::set<MIAtom*> *nbors = &bond_map[a1];
                            if (nbors->find(a2) != nbors->end())
                            {
                                found = 1;
                            }
                        }
                        if (found == 0)
                        {
                            std::set<MIAtom*> *nbors = &angle_map[a1];
                            if (nbors->find(a2) != nbors->end())
                            {
                                found = 1;
                            }
                        }

                        // angle due to PRO being cyclical
                        if (strcmp(res->type().c_str(), "PRO") == 0 || strcmp(res2->type().c_str(), "PRO") == 0)
                        {
                            if ((strcmp(a1->name(), "CD") == 0 && strcmp(a2->name(), "C") == 0)
                                || (strcmp(a2->name(), "CD") == 0 && strcmp(a1->name(), "C") == 0))
                            {
                                found = 1;
                            }
                        }
                        if (found == 0)
                        {
                            bond.setAtom1(a1);
                            bond.setAtom2(a2);
                            if ((a1->name()[0] == 'O' && a2->name()[0] == 'N')
                                || (a1->name()[0] == 'O' && a2->name()[0] == 'O')
                                || (a1->name()[0] == 'N' && a2->name()[0] == 'N')
                                || (a2->name()[0] == 'O' && a1->name()[0] == 'N'))
                            {
                                bond.ideal_length = 2.75F;
                                bond.tolerance = sigmabump*3.0F;
                            }
                            else
                            {
                                bond.ideal_length = 3.1F;
                                bond.tolerance = sigmabump;
                            }
                            RefiBumps.push_back(bond);
                        }
                    }
                }
                n2++;
                res2 = res2->next();
            }
        }
        n++;
        res = res->next();
    }
    return 1;
}

int MIMolDictionary::GetResidueTorsions(Residue *res, vector<TORSION> &torsions)
{
    int nFound = 0;
    unsigned int i, in;
    char name[10];
    MIAtom *a[4];
    TORSION torsion;

    for (i = 0; i < (unsigned int)TorsDict.size(); i++)
    {
        if (strcmp(res->type().c_str(), TorsDict[i].restype) == 0)
        {
            a[0] = a[1] = a[2] = a[3] = NULL;
            //  we only want free torsions - ones with nideal > 0 are impropers
            if (TorsDict[i].nideal > 0)
            {
                continue;
            }
            for (in = 0; in < 4; in++)
            {
                strcpy(name, TorsDict[i].name[in]);
                // don't want ones that involve the next or previous residue
                if (strchr(name, '-') != NULL)
                {
                    continue;
                }
                if (strchr(name, '+') != NULL)
                {
                    continue;
                }
                a[in] = MIAtomFromNameIncludingSynonyms(name, res);
            }
            if (a[0] && a[1] && a[2] && a[3])
            {
                torsion.setAtom1(a[0]);
                torsion.setAtom2(a[1]);
                torsion.atom3 = a[2];
                torsion.atom4 = a[3];
                torsion.nideal = TorsDict[i].nideal;
                for (int in = 0; in < torsion.nideal; in++)
                {
                    torsion.ideal[in] = TorsDict[i].ideal[in];
                }
                torsion.tolerance = 5.0F;
                torsion.res = res;
                torsions.push_back(torsion);
                nFound++;
            }
        }
    }
    // if nFound is 0 check to see if because there is no dictionary for this type
    if (nFound == 0)
    {
        if (GetDictResidue(res->type().c_str()) == NULL)
        {
            Logger::log("Warning: No dictionary entry found for type %s", res->type().c_str());
        }
    }
    return nFound;
}

TORSION*MIMolDictionary::getTORSION(Residue *res, const char *type, Residue *prev) const
{
    /* search for a torsion in dict that matches
     * res type and torsion type and then
     * find the 4 atoms that match  and return a TORSION structure pointer
     */
    char name[10];
    TORSION *torsion;
    unsigned int i, in;
    Residue *next = NULL, *current;
    MIAtom *a[4];
    a[0] = a[1] = a[2] = a[3] = NULL;
    current = res;
    if (res == NULL || type == NULL)
    {
        return NULL;
    }
    if (Monomer::isValid(res->next()))
    {
        next = res->next();
    }
    //if(res->prev_res!=NULL)prev = res->prev_res;
    for (i = 0; i < TorsDict.size(); i++)
    {
        if (!strcmp(type, TorsDict[i].type)
            && !strcmp(res->type().c_str(), TorsDict[i].restype))
        {
            break;
        }
    }
    if (i >= TorsDict.size())
    {
        return NULL;
    }
    for (in = 0; in < 4; in++)
    {
        strcpy(name, TorsDict[i].name[in]);
        res = current;
        if (strchr(name, '-') != NULL)
        {
            res = prev;
            name[strlen(name)-1] = '\0';
        }
        if (strchr(name, '+') != NULL)
        {
            res = next;
            name[strlen(name)-1] = '\0';
        }
        if (Monomer::isValid(res))
        {
            a[in] = MIAtomFromNameIncludingSynonyms(name, res);
        }
    }
    if (!a[0] || !a[1] || !a[2] || !a[3])
    {
        return NULL;
    }
    torsion = (TORSION*)calloc(1, sizeof(TORSION));
    /* detect out of memory condition */
    if (!torsion)
    {
        return NULL;
    }
    torsion->setAtom1(a[0]);
    torsion->setAtom2(a[1]);
    torsion->atom3 = a[2];
    torsion->atom4 = a[3];
    torsion->nideal = TorsDict[i].nideal;
    torsion->tolerance = 5.0F;
    for (int in = 0; in < torsion->nideal; in++)
    {
        torsion->ideal[in] = TorsDict[i].ideal[in];
    }
    strcpy(torsion->type, TorsDict[i].type);
    torsion->res = current;
    return torsion;
}

bool MIMolDictionary::LoadDictionary(const char *path, bool append, bool replace,
                                     unsigned int new_h_level)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        return false;
    }
    Logger::log("Loading Dictionary:");
    Logger::log(path);
    if (replace && GetModified())
    {
        // give user a chance to save if modified?
    }
    if (append)
    {
        SetModified(true);
    }
    if (!append && replace)
    {
        HLevel = DictionaryHLevel::Unknown;
    }
    if (new_h_level != DictionaryHLevel::Unknown)
    {
        HLevel = new_h_level;
    }
    bool result = LoadDict(fp, append, replace) != 0;
    fclose(fp);
    return result;
}

static bool in_to_add(vector<std::string> &to_add, char *type)
{
    for (size_t i = 0; i < to_add.size(); i++)
    {
        if (strcmp(to_add[i].c_str(), type) == 0)
        {
            return true;
        }
    }
    return false;
}

void MIMolDictionary::DeleteConformers(const std::string &type)
{
    Residue *next, *oldres;

    oldres = ResDict;
    while (oldres != NULL)
    {
        next = oldres->next();
        if (strcmp(oldres->type().c_str(), type.c_str()) == 0)
        {

            oldres->removeFromList();
            if (oldres == ResDict)
            {
                ResDict = next;
            }
            FreeResidueList(oldres);
        }
        oldres = next;
    }
}

static bool ContainsAdjacentAtoms(const TORSDICT &t)
{
    if (strchr(t.name[0], '-') != NULL)
    {
        return true;
    }
    if (strchr(t.name[0], '+') != NULL)
    {
        return true;
    }
    if (strchr(t.name[1], '-') != NULL)
    {
        return true;
    }
    if (strchr(t.name[1], '+') != NULL)
    {
        return true;
    }
    if (strchr(t.name[2], '-') != NULL)
    {
        return true;
    }
    if (strchr(t.name[2], '+') != NULL)
    {
        return true;
    }
    if (strchr(t.name[3], '-') != NULL)
    {
        return true;
    }
    if (strchr(t.name[3], '+') != NULL)
    {
        return true;
    }
    return false;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    LoadRes
// Purpose:		Adds a list of residues to the dictionary
// Input:       A pointer to the first residue in the list
//				A flag for appending rather than discarding the existing dict.
//				A flag for replacing existing residues rather than discarding
//				input residues that duplicate existing codes.
// Output:      Rebuilds the dictionary & returns the number of residues in the
//				new dictionary
// Requires:
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::LoadRes(Residue *respdb, bool append, bool replace_old, bool rplc_backbone_tor)
{
    Residue *res;
    unsigned int i, j;
    vector<std::string> to_add;
    Residue *newres;
    Residue *prev, *next, *oldres;

    //Can't append to non-existent dictionary!
    if (ResDict == NULL)
    {
        append = false;
    }

    if (replace_old && append)
    {
        char oldtype[MAXNAME];
        newres = respdb;
        while (Monomer::isValid(newres))
        {
            oldres = ResDict;
            while (oldres != NULL)
            {
                next = oldres->next();
                if (strcmp(oldres->type().c_str(), newres->type().c_str()) == 0)
                {
                    next = oldres->next();
                    strcpy(oldtype, oldres->type().c_str());
                    oldres->removeFromList();
                    if (oldres == ResDict)
                    {
                        ResDict = next;
                    }
                    FreeResidueList(oldres);
                    j = 0;
                    while (j < PlaneDict.size())
                    {
                        if (strcmp(PlaneDict[j].restype, oldtype) == 0)
                        {
                            PlaneDict.erase(PlaneDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                    j = 0;
                    while (j < TorsDict.size())
                    {
                        if (strcmp(TorsDict[j].restype, oldtype) == 0
                            && (rplc_backbone_tor
                                || !ContainsAdjacentAtoms(TorsDict[j])))
                        {
                            TorsDict.erase(TorsDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                    j = 0;
                    while (j < ChiralDict.size())
                    {
                        if (strcmp(ChiralDict[j].restype, oldtype) == 0)
                        {
                            ChiralDict.erase(ChiralDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                }
                oldres = next;
            }
            newres = newres->next();
        }
    }

    if (!replace_old && append)
    {
        //		char oldtype[MAXNAME];
        newres = ResDict;
        while (newres != NULL)
        {
            oldres = respdb;
            while (Monomer::isValid(oldres))
            {
                next = oldres->next();
                if (strcmp(oldres->type().c_str(), newres->type().c_str()) == 0)
                {
                    if (oldres != respdb)
                    {
                        prev = oldres->prev();
                    }
                    else
                    {
                        prev = NULL;
                    }
                    next = oldres->next();
                    if (prev != NULL)
                    {
                        prev->insertResidue(next);
                    }
                    else
                    {
                        respdb = next;
                    }
                }
                oldres = next;
            }
            newres = newres->next();
        }
    }

    newres = respdb;
    while (Monomer::isValid(newres))
    {
        to_add.push_back(newres->type());
        newres = newres->next();
    }

    //Skip the rest if we don't have any residues to add, to protect against
    //bad input
    if (to_add.empty())
    {
        return nResDict;
    }

    // get rid of old ones if not replacing and a dictionary list already exists
    if (append && !replace_old && ResDict != NULL)
    {
        oldres = ResDict;
        while (oldres != NULL)
        {
            i = 0;
            while (i < to_add.size())
            {
                if (strcmp(to_add[i].c_str(), oldres->type().c_str()) == 0)
                {
                    to_add.erase(to_add.begin() + i);
                }
                else
                {
                    i++;
                }
            }
            oldres = oldres->next();
        }
    }

    //Make copies of the residues for the dictionary
    newres = CopyResList(respdb);

    // add to list or replace
    if (!append)
    {
        if (RefiPlanes.size() > 0)
        {
            for (i = 0; i < RefiPlanes.size(); i++)
            {
                free(RefiPlanes[i].atoms);
            }
        }
        RefiPlanes.clear();
        PlaneDict.clear();
        Synonym.clear();
        TorsDict.clear();
        ChiralDict.clear();
        if (ResDict)
        {
            FreeResidueList(ResDict);
        }
        ResDict = newres;
    }
    else
    {
        res = ResDict;
        if (res)
        {
            while (res->next() != NULL)
            {
                res = res->next();
            }
            res->insertResidue(newres);
        }
    }


    /* count the length of the list */
    if ((res = ResDict) == NULL)
    {
        nResDict = 0;
    }
    else
    {
        nResDict = 1;
        while ((res = res->next()) != NULL)
        {
            nResDict++;
        }
    }

    for (res = ResDict; res != NULL; res = res->next())
    {
        res->setName1(singleletter(res->type().c_str()));
    }

    build_map();
    SetModified(true);
    return (nResDict);
}

bool MIMolDictionary::AddConfs(Residue *res, const std::string res_type)
{
    if (DictContains(res_type))
    {
        Residue *last = ResDict;
        while (last->next())
        {
            last = last->next();
        }
        last->insertResidue(res);
        build_map();
        SetModified(true);
        return true;
    }
    else
    {
        return false;
    }
}

bool MIMolDictionary::AddConfs(const ConfSaver &confs, bool replace)
{
    int n = std::min(confs.NumberSets(), 999999);
    Residue *dictres = GetDictResidue(confs.GetResidue()->type().c_str(), 0);
    if (dictres != confs.GetResidue())
    {
        return false;
    }


    ConfSaver orig(dictres);
    orig.Save();

    Residue *first_new, *current_new, *prev_new;
    for (int i = 0; i < n; ++i)
    {
        char buf[MAXNAME];
        confs.Restore(i+1);
        current_new = new Residue(*dictres);
        sprintf(buf, "%d", i+1);
        current_new->setName(std::string(buf));

        if (i == 0)
        {
            first_new = current_new;
        }
        else
        {
            prev_new->insertResidue(current_new);
        }
        prev_new = current_new;
    }
    if (replace)
    {
        DeleteConformers(dictres->type());
    }
    else
    {
        orig.Restore(1);
    }

    Residue *last = ResDict;
    if (last)
    {
        while (last->next())
        {
            last = last->next();
        }
        last->insertResidue(first_new);
    }
    build_map();
    SetModified(true);
    return true;
}

/////////////////////////////////////////////////////////////////////////////
// Function:	GetConfs
// Purpose:		Gets all the conformers from the dictionary for a given residue,
//				and stores them in a GeomSaver, associating them with the Atom
//				objects in the input residue
// Input:       A GeomSaver object for storing the conformers
//				A residue ptr which is the source of the atoms and the residue type
//				A MIMolOpt which is the source of the dict w/conformers
//				A MIMoleculeBase ptr, required for the GeomSaver
// Output:		Writes new coordinate sets to the GeomSaver argument
// Requires:	That the atom names be the same in the model and the dictionary
/////////////////////////////////////////////////////////////////////////////

bool GetConfs(GeomSaver &confs, Residue *res, MIMolDictionary *dict, MIMoleculeBase *model, unsigned int max)
{
    //First save the current coords of res, to restore at the end of this ftn
    GeomSaver tmp;
    tmp.Save(res, 1, model);

    //Now loop through the conformers in the dictionary, temporarily setting the coords
    //of res to that conformer and then saving it to the geomsaver
    int i, nFound = 0;
    int nConf = std::min(dict->GetNumberInDict(res->type().c_str()), max);
    Residue *conf;

    for (i = 0; i < nConf; ++i)
    {
        conf = dict->GetDictResidue(res->type().c_str(), i);
        if (MICopyCoords(res, conf))
        {
            confs.Save(res, 1, model);
            nFound++;
        }
    }

    tmp.Restore(1);
    return nFound > 0;
}

Residue *ExpandConfs(const Residue *single, const GeomSaver &confs)
{
    int n = std::min(confs.NumberSets() - 1, 999999);           //Cap at 1 million confs
    if (n <= 0)
    {
        Residue *current = new Residue(*single);
        return current;
    }

    Residue *first, *current, *prev;
    for (int i = 0; i < n; ++i)
    {
        char buf[MAXNAME];
        confs.Restore(i+1);
        current = new Residue(*single);
        sprintf(buf, "%d", i+1);
        current->setName(std::string(buf));

        if (i == 0)
        {
            first = current;
        }
        else
        {
            prev->insertResidue(current);
        }
        prev = current;
    }

    return first;
}


/*  a routine to read in the dictionary */
int MIMolDictionary::LoadDict(FILE *fp, bool append, bool replace_old)
{
    float value;
    Residue *res, *respdb;
    char buf[2000];
    vector<Bond> connects;
    unsigned int i, j;
    vector<std::string> to_add;
    Residue *newres;
    Residue *next, *oldres;

    respdb = LoadPDB(fp, &connects);

    if (replace_old && append)
    {
        char oldtype[MAXNAME];
        newres = respdb;
        while (newres != NULL)
        {
            oldres = ResDict;
            while (oldres != NULL)
            {
                next = oldres->next();
                if (strcmp(oldres->type().c_str(), newres->type().c_str()) == 0)
                {

                    strcpy(oldtype, oldres->type().c_str());
                    oldres->removeFromList();
                    if (oldres == ResDict)
                    {
                        ResDict = next;
                    }
                    FreeResidueList(oldres);

                    j = 0;
                    while (j < PlaneDict.size())
                    {
                        if (strcmp(PlaneDict[j].restype, oldtype) == 0)
                        {
                            PlaneDict.erase(PlaneDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                    j = 0;
                    while (j < TorsDict.size())
                    {
                        if (strcmp(TorsDict[j].restype, oldtype) == 0)
                        {
                            TorsDict.erase(TorsDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                    j = 0;
                    while (j < ChiralDict.size())
                    {
                        if (strcmp(ChiralDict[j].restype, oldtype) == 0)
                        {
                            ChiralDict.erase(ChiralDict.begin()+j);
                        }
                        else
                        {
                            j++;
                        }
                    }
                }
                oldres = next;
            }
            newres = newres->next();
        }
    }

    if (!replace_old && append)
    {
        char oldtype[MAXNAME];
        newres = ResDict;
        while (newres != NULL)
        {
            oldres = respdb;
            while (oldres != NULL)
            {
                next = oldres->next();
                if (strcmp(oldres->type().c_str(), newres->type().c_str()) == 0)
                {

                    strcpy(oldtype, oldres->type().c_str());
                    oldres->removeFromList();
                    if (oldres == respdb)
                    {
                        respdb = next;
                    }
                    FreeResidueList(oldres);
                }
                oldres = next;
            }
            newres = newres->next();
        }
    }

    newres = respdb;
    while (newres != NULL)
    {
        to_add.push_back(newres->type());
        newres = newres->next();
    }

    // get rid of old ones if not replacing and a dictionary list already exists
    if (append && !replace_old && ResDict != NULL)
    {
        oldres = ResDict;
        while (oldres != NULL)
        {
            i = 0;
            while (i < to_add.size())
            {
                if (strcmp(to_add[i].c_str(), oldres->type().c_str()) == 0)
                {
                    to_add.erase(to_add.begin() + i);
                }
                else
                {
                    i++;
                }
            }
            oldres = oldres->next();
        }
    }

    // add to list or replace
    if (!append)
    {
        if (RefiPlanes.size() > 0)
        {
            for (i = 0; i < RefiPlanes.size(); i++)
            {
                free(RefiPlanes[i].atoms);
            }
        }
        RefiPlanes.clear();
        PlaneDict.clear();
        Synonym.clear();
        TorsDict.clear();
        ChiralDict.clear();
        if (ResDict)
        {
            FreeResidueList(ResDict);
        }
        ResDict = respdb;
    }
    else
    {
        res = ResDict;
        if (res)
        {
            while (res->next() != NULL)
            {
                res = res->next();
            }
            res->insertResidue(respdb);
        }
    }


    /* count the length of the list */
    if ((res = ResDict) == NULL)
    {
        nResDict = 0;
    }
    else
    {
        nResDict = 1;
        while ((res = res->next()) != NULL)
        {
            nResDict++;
        }
    }

    for (res = ResDict; res != NULL; res = res->next())
    {
        res->setName1(singleletter(res->type().c_str()));
    }

    rewind(fp);
    TORSDICT tdict;
    PLANEDICT pdict;
    SYNONYM sdict;
    /* read in torsions */
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (strncmp(buf, "TORSION", 7))
        {
            continue;
        }
        if (sscanf(buf, "%*s%s%s%s%s%s%s",
                   tdict.restype,
                   tdict.type,
                   tdict.name[0],
                   tdict.name[1],
                   tdict.name[2],
                   tdict.name[3]) == 6)
        {
            if (!in_to_add(to_add, tdict.restype))
            {
                continue;
            }
            tdict.nideal = 0;
            /* load preferred values if any - max of 3  */
            if (sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%f", &value) == 1)
            {
                while (value < 0.0)
                {
                    value += 360.0;
                }
                while (value >= 360.0)
                {
                    value -= 360.0;
                }
                tdict.ideal[0] = value;
                tdict.nideal++;
                if (sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*f%f", &value) == 1)
                {
                    while (value < 0.0)
                    {
                        value += 360.0;
                    }
                    while (value >= 360.0)
                    {
                        value -= 360.0;
                    }
                    tdict.ideal[1] = value;
                    tdict.nideal++;
                    if (sscanf(buf, "%*s%*s%*s%*s%*s%*s%*s%*f%*f%f", &value) == 1)
                    {
                        while (value < 0.0)
                        {
                            value += 360.0;
                        }
                        while (value >= 360.0)
                        {
                            value -= 360.0;
                        }
                        tdict.ideal[2] = value;
                        tdict.nideal++;
                    }
                }
            }
            TorsNames.insert(std::string(tdict.type));
            TorsDict.push_back(tdict);
        }
    }
    Logger::log("%d torsions in dictionary", (int)TorsDict.size());

    rewind(fp);
    /* read in PlaneDict */
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (strncmp(buf, "PLANE", 5))
        {
            continue;
        }
        memset(&pdict, 0, sizeof(PLANEDICT));
        if (sscanf(buf, "%*s%s%d", pdict.restype, &pdict.natoms) != 2)
        {
            continue;
        }
        if (!in_to_add(to_add, pdict.restype))
        {
            continue;
        }
        if (pdict.natoms > MAXPLANE)
        {
            pdict.natoms = MAXPLANE;
        }
        sscanf(buf, "%*s%*s%*d%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s",
               pdict.name[0],
               pdict.name[1],
               pdict.name[2],
               pdict.name[3],
               pdict.name[4],
               pdict.name[5],
               pdict.name[6],
               pdict.name[7],
               pdict.name[8],
               pdict.name[9],
               pdict.name[10],
               pdict.name[11],
               pdict.name[12],
               pdict.name[13],
               pdict.name[14],
               pdict.name[15],
               pdict.name[16],
               pdict.name[17],
               pdict.name[18],
               pdict.name[19]
               );
        PlaneDict.push_back(pdict);
    }
    Logger::log("%d planes in dictionary", (int)PlaneDict.size());

    /* read in synonyms */
    rewind(fp);
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (strncmp(buf, "SYNON", 5))
        {
            continue;
        }
        if (3 == sscanf(buf, "%*s%s%s%s",
                        sdict.resname,
                        sdict.aname1,
                        sdict.aname2))
        {
            Synonym.push_back(sdict);
        }
    }
    Logger::log("%d synonyms in dictionary", (int)Synonym.size());

    /* read in bonds & angles, storing the info in maps for quick lookup by residue type */
    rewind(fp);
    multimap<std::string, BondDICT> bondMap;
    multimap<std::string, ANGLEDICT> angleMap;
    while (fgets(buf, sizeof buf, fp) != NULL)
    {
        if (!strncmp(buf, "ANGLE", 5))
        {
            ScanAngle(angleMap, buf);
        }
        if (!strncmp(buf, "BOND", 4))
        {
            ScanBond(bondMap, buf);
        }
    }

    //Create the "prefbonds" and "prefangles" information for each residue,
    //looking up the angles and bonds from the maps and looking up atom ptrs
    //using the atom names.
    for (res = ResDict; res != NULL; res = res->next())
    {
        if (bondMap.count(res->type()) != 0)
        {
            CreateDictBonds(res, bondMap);
        }
        if (angleMap.count(res->type()) != 0)
        {
            CreateDictAngles(res, angleMap);
        }
    }


    build_map();

    return (nResDict);

}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddTorsion
// Purpose:		Adds a torsion (or improper) to the dictionary
// Input:       A reference to the torsion struct to be added
// Output:      Stores the torsion object in the dictionary and
//				returns 1; or returns 0 if the torsion could not be added
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddTorsion(const TORSION &torsion)
{
    TORSDICT td;

    if (torsion.res == 0                 //Check that all the pointers are valid
        || torsion.getAtom1() == 0
        || torsion.getAtom2() == 0
        || torsion.atom3 == 0
        || torsion.atom4 == 0)
    {
        return 0;
    }

    strcpy(td.restype, torsion.res->type().c_str());
    strcpy(td.type, torsion.type);

    strcpy(td.name[0], torsion.getAtom1()->name());
    strcpy(td.name[1], torsion.getAtom2()->name());
    strcpy(td.name[2], torsion.atom3->name());
    strcpy(td.name[3], torsion.atom4->name());

    if (torsion.nideal < 0              //Check that we have a valid number of ideal angles
        || torsion.nideal > 3)
    {
        return 0;
    }

    td.nideal = torsion.nideal;

    for (int i = 0; i < td.nideal; ++i)
    {
        td.ideal[i] = torsion.ideal[i];
    }

    TorsDict.push_back(td);

    return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddPlane
// Purpose:		Adds a torsion (or improper) to the dictionary
// Input:       A reference to the torsion struct to be added
// Output:      Stores the plane in the dictionary and
//				returns 1; or returns 0 if the plane could not be added
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddPlane(const PLANE &plane)
{
    PLANEDICT pd;

    if (plane.natoms <= 0                    //Check that we received a valid plane
        || plane.res == 0
        || plane.atoms == 0)
    {
        return 0;
    }

    if (plane.natoms > MAXPLANE)             //Check that we have enough space to hold the
    {
        return 0;                           //atom names
    }

    strcpy(pd.restype, plane.res->type().c_str());

    pd.natoms = plane.natoms;

    for (int i = 0; i < pd.natoms; ++i)
    {
        if (plane.atoms[i] == 0)            //Check each atom ptr
        {
            return 0;
        }
        strcpy(pd.name[i], plane.atoms[i]->name());
    }

    PlaneDict.push_back(pd);

    return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddChiral
// Purpose:		Adds a chiral center to the dictionary
// Input:       A reference to the chiral to be added, and the name of the residue
// Output:      Stores the plane in the dictionary and
//				returns 1; or returns 0 if the chiral wasn't valid
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddChiral(const CHIRAL &chiral, const char *res_type)
{
    CHIRALDICT cd;

    if (chiral.center == 0                   //Check that we received a valid chiral
        || chiral.getAtom1() == 0
        || chiral.getAtom2() == 0
        || chiral.atom3 == 0
        || res_type == 0)
    {
        return 0;
    }


    strncpy(cd.restype, res_type, MAXNAME);
    strncpy(cd.center, chiral.center->name(), MAXATOMNAME);
    strncpy(cd.name[0], chiral.getAtom1()->name(), MAXATOMNAME);
    strncpy(cd.name[1], chiral.getAtom2()->name(), MAXATOMNAME);
    strncpy(cd.name[2], chiral.atom3->name(), MAXATOMNAME);

    cd.order = chiral.order;

    ChiralDict.push_back(cd);

    return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddTorsion
// Purpose:		Adds a torsion (or improper) to the dictionary
// Input:       A reference to the torsion struct to be added
// Output:      Stores the torsion object in the dictionary and
//				returns 1; or returns 0 if the torsion could not be added
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddTorsion(const TORSDICT &torsion)
{
    dict_map::iterator p;
    if ((p = DictMap.find(torsion.restype)) == DictMap.end())
    {
        return 0;
    }

    for (int i = 0; i < 4; ++i)
    {
        if (CountAtomsByName(torsion.name[i], p->second.residue()) != 1)
        {
            return 0;
        }
    }

    TorsDict.push_back(torsion);
    return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddPlane
// Purpose:		Adds a torsion (or improper) to the dictionary
// Input:       A reference to the torsion struct to be added
// Output:      Stores the plane in the dictionary and
//				returns 1; or returns 0 if the plane could not be added
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddPlane(const PLANEDICT &plane)
{

    dict_map::iterator p;
    if ((p = DictMap.find(plane.restype)) == DictMap.end())
    {
        return 0;
    }

    if (plane.natoms <= 0)
    {
        return 0;
    }

    for (int i = 0; i < plane.natoms; ++i)
    {
        if (CountAtomsByName(plane.name[i], p->second.residue()) != 1)
        {
            return 0;
        }
    }

    PlaneDict.push_back(plane);
    return 1;
}

/////////////////////////////////////////////////////////////////////////////
// Function:    AddChiral
// Purpose:		Adds a chiral center to the dictionary
// Input:       A reference to the chiral to be added, and the name of the residue
// Output:      Stores the plane in the dictionary and
//				returns 1; or returns 0 if the chiral wasn't valid
// Requires:	Doesn't check for duplication
/////////////////////////////////////////////////////////////////////////////
int MIMolDictionary::AddChiral(const CHIRALDICT &chiral)
{
    if (strlen(chiral.restype) ==  0)
    {
        return 0;
    }

    dict_map::iterator p;
    if ((p = DictMap.find(chiral.restype)) == DictMap.end())
    {
        return 0;
    }

    for (int i = 0; i < 4; ++i)
    {
        if (CountAtomsByName(chiral.name[i], p->second.residue()) != 1)
        {
            return 0;
        }
    }

    ChiralDict.push_back(chiral);
    return 1;
}

void MIMolDictionary::Clear()
{
    RefiAngles.clear();
    RefiBonds.clear();
    for (size_t i = 0; i < RefiPlanes.size(); i++)
    {
        free(RefiPlanes[i].atoms);
    }
    RefiPlanes.clear();
    RefiPhiPsis.clear();
    RefiTorsions.clear();
    RefiChirals.clear();
    RefiConstraints.clear();
    RefiBumps.clear();
}

unsigned int MIMolDictionary::GetNumberInDict(const char *type)
{
    unsigned int n = 0;
    dict_map::iterator p1 = DictMap.lower_bound(type);
    if (p1 == DictMap.end())
    {
        return n;
    }
    dict_map::iterator p2 = DictMap.upper_bound(type);
    while (p1 != p2 && p1 != DictMap.end())
    {
        n++;
        p1++;
    }
    return n;
}

vector<std::string> MIMolDictionary::GetDictResList()
{
    dict_map::iterator p = DictMap.begin();
    std::vector<std::string> choices;
    choices.clear();
    if (p != DictMap.end())
    {
        choices.push_back(p->second.residue()->type());
    }
    while (p != DictMap.end())
    {
        if (p->second.residue()->type() != choices.back())
        {
            choices.push_back(p->second.residue()->type());
        }
        p++;
    }
    return choices;
}

Residue*MIMolDictionary::GetDictResidue(const char *type, int nconformer)
{
    if (EmptyDictCheck() == false)
    {
        return NULL;
    }
    dict_map::iterator p = GetDictEntry(type, nconformer);
    if (p == DictMap.end())
    {
        return NULL;
    }
    else
    {
        return p->second.residue();
    }
}

Residue*MIMolDictionary::GetDictResidue(const char single, int nconformer)
{
    const char *type = tripleletter(single);
    if (type)
    {
        return GetDictResidue(type, nconformer);
    }
    else
    {
        return NULL;
    }
}

vector<Bond>*MIMolDictionary::GetDictBonds(const char *type, int nconformer)
{
    dict_map::iterator p = GetDictEntry(type, nconformer);
    if (p == DictMap.end())
    {
        return NULL;
    }
    else
    {
        return p->second.Bonds();
    }
}

vector<ANGLE>*MIMolDictionary::GetDictAngles(const char *type, int nconformer)
{
    dict_map::iterator p = GetDictEntry(type, nconformer);
    if (p == DictMap.end())
    {
        return NULL;
    }
    else
    {
        return p->second.Angles();
    }
}

// Chirals not in DictResidue; doesn't work
//
// vector<CHIRAL> * MIMolDictionary::GetDictChirals(const char *type, int nconformer)
// {
//      dict_map::iterator p = GetDictEntry(type, nconformer);
//      if(p == DictMap.end()) return NULL;
//      else return p->second.Chirals();
// }

dict_map::iterator MIMolDictionary::GetDictEntry(const char *type, int nconformer)
{
    dict_map::iterator p = DictMap.find(type);
    if (p != DictMap.end() && nconformer > 0)
    {
        int ntype = GetNumberInDict(type);
        nconformer %= ntype;
        while (nconformer > 0)
        {
            p++;
            nconformer--;
        }
    }
    return p;
}

unsigned int MIMolDictionary::BuildInternalBumpBonds(MIAtomList &CurrentAtoms, vector<Bond> &bonds)
{
    size_t i, j;
    Bond bond;
    MIAtom *a1, *a2;
    for (i = 0; i < CurrentAtoms.size()-1; i++)
    {
        a1 = CurrentAtoms[i];
        for (j = i+1; j < CurrentAtoms.size(); j++)
        {
            a2 = CurrentAtoms[j];
            if (AtomDist(*a1, *a2) > 3.30)
            {
                bond.setAtom1(a1);
                bond.setAtom2(a2);
                if (MIAtom::MIIsHBondable(a1, a2))
                {
                    bond.ideal_length = 2.3F;
                }
                else
                {
                    bond.ideal_length = 2.7F;
                }
                bonds.push_back(bond);
            }

        }
    }
    return bonds.size();
}

//FIXME: move Stddev functions into this class from MIMolOpt?

bool MIMolDictionary::DictHCheck(Residue *res, unsigned int &level_return)
{
    level_return = HLevel;

    if (EmptyDictCheck() == false)
    {
        return false;
    }
    // this routine checks that the H-level (none, polar or all) is
    // the same as the residue to be replaced.  If not the user is given
    // an oppurtunity to load the proper dictionary

    // if the level of the dictionary is unknown (the user loaded her own dict for example)
    // just go with it
    if (HLevel == DictionaryHLevel::Unknown)
    {
        return true;
    }

    int nall = 0;
    int npolar = 0;
    unsigned int level = DictionaryHLevel::NoHydrogens;
    bool peptide = IsPeptide(*res);
    bool nucleic = IsNucleic(res);

    // only check petides or nucleic acids
    if (!(peptide || nucleic))
    {
        return true;
    }

    for (int i = 0; i < res->atomCount(); i++)
    {
        if (MIAtom::MIIsHydrogen(res->atom(i)))
        {
            nall++;
        }
        if (IsPolarH(res->atom(i), res) )
        {
            npolar++;
        }
    }
    if (nall > 0)
    {
        level = DictionaryHLevel::All;
    }
    if (npolar > 0 && npolar == nall)
    {
        level = DictionaryHLevel::Polar;
    }

    if (level != HLevel)
    {
        level_return = level;
        return false;
    }
    return true;
}

bool MIMolDictionary::EmptyDictCheck()
{
    if (nResDict > 0)
    {
        return true;
    }

    Logger::message("Serious Error! There is no dictionary loaded.  Set a dictionary in File/Preferences...");
    return false;
}

bool MIMolDictionary::DictContains(const std::string &res) const
{
    return DictMap.find(res.c_str()) != DictMap.end();
}

bool MIMolDictionary::SaveDictionary(const char *pathname, const char *res_type)
{
    FILE *fp;
    if ((fp = fopen(pathname, "w")) == NULL)
    {
        Logger::message("Unable to open file - write protected?");
        return false;
    }

    if (res_type != NULL)
        return fwriteDictEntry(fp, res_type);

    // save the entire dictionary
    std::set<std::string> residue_set;

    for (dict_map::iterator p = DictMap.begin(); p != DictMap.end(); ++p)
        residue_set.insert(p->first);
    for (std::set<std::string>::iterator i = residue_set.begin();
         i!=residue_set.end(); ++i)
    {
        fwriteDictEntry(fp, i->c_str());
    }

    fclose(fp);
    SetModified(false);
    return true;
}

bool MIMolDictionary::fwriteDictEntry(FILE *fp, const char *res_type)
{
    Residue *res = GetDictResidue(res_type);
    size_t i;
    if (!res)
    {
        return false;
    }

    // iterate through all the confomers and save each residue in PDB format
    unsigned int nconf = GetNumberInDict(res_type);
    for (i = 0; i < nconf; i++)
    {
        res = GetDictResidue(res_type, i);
        if (!res)
            return false;
        Residue *tmp = res->next();
        res->setNext(NULL);
        bool result = SavePDB(fp, res, NULL, 0, false);
        res->setNext(tmp);
        if (!result)
            return false;
    }

    vector<Bond> *bonds = GetDictBonds(res_type);
    vector<ANGLE> *angles = GetDictAngles(res_type);

    for (i = 0; i < bonds->size(); i++)
    {
        fprintf(fp, "BOND %s %s %s %f\n", res_type, (*bonds)[i].getAtom1()->name(), (*bonds)[i].getAtom2()->name(), (*bonds)[i].ideal_length);
    }
    for (i = 0; i < angles->size(); i++)
    {
        fprintf(fp, "ANGLE %s %s %s %f %s\n", res_type, (*angles)[i].getAtom1()->name(), (*angles)[i].getAtom2()->name(),
                (*angles)[i].ideal_angle, (*angles)[i].atom3->name());
    }

    // find planes
    int j;
    for (i = 0; i < PlaneDict.size(); i++)
    {
        if (strcmp(PlaneDict[i].restype, res_type))
        {
            continue;
        }
        /* found a match */

        //PLANE PHE 6 CG CD1 CD2 CE1 CE2 CZ
        fprintf(fp, "PLANE %s %d", PlaneDict[i].restype, PlaneDict[i].natoms);
        for (j = 0; j < PlaneDict[i].natoms; j++)
        {
            fprintf(fp, " %s", PlaneDict[i].name[j]);
        }
        fprintf(fp, "\n");
    }

    // find torsions
    for (i = 0; i < TorsDict.size(); i++)
    {
        if (strcmp(TorsDict[i].restype, res_type))
        {
            continue;
        }
        //TORSION TYR	IMP5   CG   CD1  CE1  CZ     0.0
        fprintf(fp, "TORSION %s %s", TorsDict[i].restype, TorsDict[i].type);
        for (j = 0; j < 4; j++)
        {
            fprintf(fp, " %s", TorsDict[i].name[j]);
        }
        for (j = 0; j < TorsDict[i].nideal; j++)
        {
            fprintf(fp, " %8.1f", TorsDict[i].ideal[j]);
        }
        fprintf(fp, "\n");
    }

    return true;
}

bool MIMolDictionary::fwriteDictEntry_mmCIF(FILE *fp, const char *res_type)
{
    bool retval;

    //Lets grab the residue and populate all the necessary vectors and whatnot
    Residue *dictres = GetDictResidue(res_type);
    if (!dictres)
    {
        return false;
    }
    Residue *reslist = new Residue(*dictres);

    //Note: this was changed from SetRefiRes (a fn in MIMolOpt) to FindGeom,
    //which simply creates the data structures we need w/o enabling
    //refinement, thus there's no need to Cancel() below, just Clear()
    FindGeom(reslist, 1, 0);

    //Give all this off to the fios::mmCIF object and tell it to go
    mmCIF fio;
    fio.WriteTorsions(false);
    MIMolInfo mi;
    mi.res = reslist;
    mi.bonds = RefiBonds;
    mi.angles = RefiAngles;
    mi.torsions = RefiTorsions;
    mi.planes = RefiPlanes;
    mi.chirals = RefiChirals;

    retval = fio.Write(fp, mi);

    // clean up
    delete reslist;
    Clear();
    return retval;
}

void MIMolDictionary::RemoveConstraints()
{
    vector<Bond>::iterator i, e;
    i = RefiConstraints.begin();
    e = RefiConstraints.end();
    for (; i != e; i++)
    {
        if (i->getAtom1()->type() & (AtomType::FREEZE | AtomType::DUMMYATOM))
        {
            delete i->getAtom1();
        }
        if (i->getAtom2()->type() & (AtomType::FREEZE | AtomType::DUMMYATOM))
        {
            delete i->getAtom2();
        }
    }
    RefiConstraints.clear();
}

int MIMolDictionary::GetFlexibleTorsions(vector <TORSION> &torsions, Residue *res) const
{
    int n = 0;

    //Lookup the flexible torsions from the dictionary
    vector<TORSDICT>::const_iterator td = TBegin();
    vector<TORSDICT>::const_iterator e = TEnd();
    TORSION *t;
    while (td != e)
    {
        if (strcmp(td->restype, res->type().c_str()) == 0             //3 disqualifiers:
            && !ContainsAdjacentAtoms(*td)              //	1) From different residue type
            && td->nideal != 1)                         //	2) Defined using atoms from neighboring residues
        {
            t = getTORSION(res, td->type);                  //  3) Improper  i.e. not flexible
            if (t != 0)
            {
                torsions.push_back(*t);
                free(t);
                n++;
            }
        }
        td++;
    }

    return n++;
}

bool ScanBond(multimap<std::string, BondDICT> &bondMap, const char *line)
{
    BondDICT bond;
    if (4 == sscanf(line, "%*s%s%s%s%f", bond.resname,
                    bond.aname1,
                    bond.aname2,
                    &bond.ideal_length))
    {
        bondMap.insert(pair<std::string, BondDICT>(bond.resname, bond));
        return true;
    }
    else
    {
        return false;
    }
}

bool CreateDictBonds(Residue *res, const multimap<std::string, BondDICT> &bondMap)
{
    multimap<std::string, BondDICT>::const_iterator b;
    pair<multimap<std::string, BondDICT>::const_iterator,
         multimap<std::string, BondDICT>::const_iterator> bounds;

    res->clearPrefBonds();
    bounds = bondMap.equal_range(res->type());
    for (b = bounds.first; b != bounds.second; b++)
    {
        Bond bond;
        bond.ideal_length = b->second.ideal_length;
        bond.setAtom1(MIAtomFromNameIncludingSynonyms(b->second.aname1, res));
        bond.setAtom2(MIAtomFromNameIncludingSynonyms(b->second.aname2, res));
        res->addPrefBond(bond);
    }
    return true;
}

bool ScanAngle(multimap<std::string, ANGLEDICT> &angleMap, const char *line)
{
    ANGLEDICT angle;
    if (5 == sscanf(line, "%*s%s%s%s%f%s", angle.resname,
                    angle.aname1,
                    angle.aname2,
                    &angle.ideal_length,
                    angle.aname3))
    {
        angleMap.insert(pair<std::string, ANGLEDICT>(angle.resname, angle));
        return true;
    }
    else
    {
        return false;
    }
}

bool CreateDictAngles(Residue *res, const multimap<std::string, ANGLEDICT> &angleMap)
{
    multimap<std::string, ANGLEDICT>::const_iterator a;
    pair<multimap<std::string, ANGLEDICT>::const_iterator,
         multimap<std::string, ANGLEDICT>::const_iterator> bounds;

    res->clearPrefAngles();
    bounds = angleMap.equal_range(res->type());
    for (a = bounds.first; a != bounds.second; a++)
    {
        ANGLE angle;
        angle.ideal_angle = a->second.ideal_length;
        angle.setAtom1(MIAtomFromNameIncludingSynonyms(a->second.aname1, res));
        angle.setAtom2(MIAtomFromNameIncludingSynonyms(a->second.aname2, res));
        angle.atom3 = MIAtomFromNameIncludingSynonyms(a->second.aname3, res);
        res->addPrefAngle(angle);
    }
    return true;
}

static MIMolDictionary *DICTIONARY = NULL;
MIMolDictionary *MIGetDictionary()
{
    return DICTIONARY;
}

void MISetDictionary(MIMolDictionary *d)
{
    DICTIONARY = d;
}

}
