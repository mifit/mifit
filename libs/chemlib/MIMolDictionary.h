#ifndef mifit_legacy_MIMolDictionary_h
#define mifit_legacy_MIMolDictionary_h

#include <set>
#include <vector>
#include <list>
#include <map>

#include "MIAtom_fwd.h"
#include "RESIDUE_fwd.h"
#include "PLANE.h"
#include "Bond.h"
#include "TORSION.h"
#include "ANGLE.h"
#include "CHIRALDICT.h"
#include "ConfSaver.h"
#include "GeomSaver.h"
#include "PDB.h"


class DictEditCanvas;
class MIMolOpt;

namespace chemlib
{

    class DictResidue;

    typedef std::multimap<std::string, DictResidue> dict_map;


// Hydrogen level in dictionary
    namespace DictionaryHLevel
    {
        const unsigned int Unknown = 0;
        const unsigned int NoHydrogens = 1;
        const unsigned int Polar = 2;
        const unsigned int All = 3;
    }

// refine options for OnFullOptimize
    namespace Refine_Level
    {
        const unsigned int None = 0;
        const unsigned int Quick = 1;
        const unsigned int Thorough = 2;
        const unsigned int Optimize = 3;
    }


    class MIMolDictionary
    {
    public:
        friend class ::MIMolOpt;

        MIMolDictionary();
        ~MIMolDictionary();

        bool LoadDefaultDictionary(const std::string &dictionary,
                                   const std::string &homedir);
        bool LoadDictionary(const char *path, bool append = false, bool replace = false,
                            unsigned int new_h_level = DictionaryHLevel::Unknown);
        int LoadRes(RESIDUE *respdb, bool append = false, bool replace_old = true,
                    bool rplc_backbone_tor = true);
        //	int LoadRes(RESIDUE *respdb, std::vector<ANGLE> &angles

        bool SaveDictionary(const char *pathname, const char *res_type = NULL);
        bool fwriteDictEntry(FILE *fp, const char *res_type);
        bool fwriteDictEntry_mmCIF(FILE *fp, const char *res_type);

        bool EmptyDictCheck();
        bool DictHCheck(RESIDUE *res, unsigned int &level);
        bool DictContains(const std::string &res) const;

        unsigned int GetNumberInDict(const char *type);

        RESIDUE *GetBeta2()
        {
            return Beta2;
        }

        RESIDUE *GetAlpha2()
        {
            return Alpha2;
        }

        std::vector<Bond> *GetDictBonds(const char *type, int nconformer = 0);
        std::vector<ANGLE> *GetDictAngles(const char *type, int nconformer = 0);
        //std::vector<CHIRAL> * GetDictChirals(const char *type, int nconformer=0);

        RESIDUE *GetDictResidue(const char *type, int nconformer = 0);
        RESIDUE *GetDictResidue(const char single, int nconformer = 0);
        std::vector<std::string> GetDictResList();
        dict_map::iterator GetDictEntry(const char *type, int nconformer = 0);
        RESIDUE *GetResdict()
        {
            return ResDict;
        }

        int GetFlexibleTorsions(std::vector <TORSION> &torsions, RESIDUE *res) const;
        TORSION *getTORSION(RESIDUE *res, const char *type, RESIDUE *prev = NULL) const;
        int GetResidueTorsions(RESIDUE *res, std::vector<TORSION> &torsions);
        int AreBonded(const char *restype, MIAtom *atom1, MIAtom *atom2, Bond &bond);

        int CountConformers(const std::string &type) const;
        void DeleteConformers(const std::string &type);

        bool AddConfs(RESIDUE *res, const std::string res_type);
        bool AddConfs(const ConfSaver &confs, bool replace = false);

        int AddTorsion(const TORSION &torsion);
        int AddPlane(const PLANE &plane);
        int AddChiral(const CHIRAL &chiral, const char *res_type);
        int AddTorsion(const TORSDICT &torsion);
        int AddPlane(const PLANEDICT &plane);
        int AddChiral(const CHIRALDICT &chiral);

        int AddConstraint(MIAtom *a1, const char *sigma);
        int AddConstraint(MIAtom *a1, MIAtom *a2, const char *d, const char *s);
        void ConstrainCalpha(RESIDUE*, int);
        void RestrainEnds(RESIDUE*, int);
        void RemoveConstraints();

        std::vector<TORSDICT>::const_iterator TBegin() const
        {
            return TorsDict.begin();
        }

        std::vector<TORSDICT>::const_iterator TEnd() const
        {
            return TorsDict.end();
        }

        bool GetModified()
        {
            return modified;
        }

        void SetModified(bool m = true)
        {
            modified = m;
        }

        bool GetConstrainCA()
        {
            return constrain_CA;
        }

        void SetConstrainCA(bool v)
        {
            constrain_CA = v;
        }

        bool GetConstrainEnds()
        {
            return constrain_Ends;
        }

        void SetConstrainEnds(bool v)
        {
            constrain_Ends = v;
        }

        float GetSigmaBond()
        {
            return sigmabond;
        }

        void SetSigmaBond(float v)
        {
            sigmabond = v;
        }

        float GetSigmaAngle()
        {
            return sigmaangle;
        }

        void SetSigmaAngle(float v)
        {
            sigmaangle = v;
        }

        float GetSigmaPlane()
        {
            return sigmaplane;
        }

        void SetSigmaPlane(float v)
        {
            sigmaplane = v;
        }

        float GetSigmaBump()
        {
            return sigmabump;
        }

        void SetSigmaBump(float v)
        {
            sigmabump = v;
        }

        float GetSigmaTorsion()
        {
            return sigmatorsion;
        }

        void SetSigmaTorsion(float v)
        {
            sigmatorsion = v;
        }

        void SetRefiSecStruct(bool v)
        {
            RefiSecStruct = v;
        }

        bool GetRefiSecStruct()
        {
            return RefiSecStruct;
        }

        std::vector<Bond> RefiBonds;
        std::vector<ANGLE> RefiAngles;
        std::vector<PLANE> RefiPlanes;
        std::vector<TORSION> RefiPhiPsis;
        std::vector<TORSION> RefiTorsions;
        std::vector<Bond> RefiConstraints;
        std::vector<Bond> RefiBumps;
        std::vector<CHIRAL> RefiChirals;

        std::set<std::string> TorsNames;

    private:
        friend class ::DictEditCanvas;

        int BuildBumps(RESIDUE *RefiRes, int nRefiRes);
        int FindGeom(RESIDUE *reslist, int nres, RESIDUE *ResActiveModel);
        int LoadDict(FILE *fp, bool append = false, bool replace = false);
        unsigned int BuildInternalBumpBonds(MIAtomList &CurrentAtoms, std::vector<Bond> &bonds);
        void Clear();
        void build_map();

        bool RefiHBonds;
        bool RefiSecStruct;
        bool constrain_CA;
        bool constrain_Ends;
        bool modified;

        int nResDict;
        RESIDUE *ResDict;

        RESIDUE *Alpha2;
        RESIDUE *Beta2;
        RESIDUE *cres;

        dict_map DictMap;

        float sigmaangle;
        float sigmabond;
        float sigmabump;
        float sigmaplane;
        float sigmatorsion;
        int CYCLEIZE; // currently, there's no way to set this
        unsigned int HLevel;

        std::vector<CHIRALDICT> ChiralDict;
        std::vector<PLANEDICT> PlaneDict;
        std::vector<TORSDICT> TorsDict;
    };

    MIAtom *MIAtomFromNameIncludingSynonyms(const char *name, const RESIDUE *residue);

//Retrieves conformations from the dictionary, and creates a geomsaver associated
//with the given residue and molecule
    bool GetConfs(GeomSaver &confs, RESIDUE *res, MIMolDictionary *dict, MIMoleculeBase *model, unsigned int max = 10000);

//Converts conformation data from a GeomSaver to a list of residues
    RESIDUE *ExpandConfs(const RESIDUE *single, const GeomSaver &confs);



    MIMolDictionary *MIGetDictionary();
    void MISetDictionary(MIMolDictionary*);

}
#endif // mifit_legacy_MIMolDictionary_h

