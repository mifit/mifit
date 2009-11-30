#ifndef mifit_fios_mmCIF_h
#define mifit_fios_mmCIF_h

#include <cstdio>
#include <string>
#include <map>
#include <vector>


#include "CifData.h"
#include "MIMolIOBase.h"

namespace chemlib
{
    class RefmacAtomTyper;

    class mmCIF : public Reader, public Writer
    {

        int _squirt_torsions; // 0 = false, 1 = true, 2 = ask user
        RefmacAtomTyper *atyper;

        std::string GetBondOrder(unsigned char);
        void GetBonds(chemlib::MIAtom *atom, std::vector<chemlib::Bond*>&, std::vector<chemlib::Bond>&);
        void SquirtTreeFindNext(chemlib::MIAtom *prev, chemlib::MIAtom *curr, chemlib::MIAtom **next,
                                std::map<chemlib::MIAtom*, bool> &alreadyprinted, std::vector<chemlib::Bond*> &sorted_neighbors,
                                std::map<chemlib::Bond*, bool> &bond_map, std::vector<chemlib::Bond*>::iterator &i);
        bool AtomName(chemlib::MIAtom*, std::map<chemlib::MIAtom*, int>&, std::string&);
        bool SquirtAtoms(FILE*, std::map<chemlib::MIAtom*, int>&, chemlib::RESIDUE*);
        bool SquirtTree(FILE *fp, std::string &resname,
                        std::map<chemlib::MIAtom*, int> &atom_map, chemlib::RESIDUE *res, std::vector<chemlib::Bond> &bonds);
        bool SquirtTreeAtom(FILE *fp, std::string &resname, chemlib::MIAtom *patom,
                            chemlib::MIAtom *catom, chemlib::MIAtom *natom, std::map<chemlib::MIAtom*, int> &atom_map, bool);
        bool SquirtBonds(FILE*, std::string&, std::map<chemlib::MIAtom*, int>&, std::vector<chemlib::Bond>&);
        bool SquirtAngles(FILE*, std::string&, std::map<chemlib::MIAtom*, int>&,
                          std::vector<chemlib::ANGLE>&);
        bool SquirtTorsions(FILE*, std::string&, std::map<chemlib::MIAtom*, int>&,
                            std::vector<chemlib::TORSION>&);
        bool SquirtChirality(FILE *fp, std::string &resname,
                             std::map<chemlib::MIAtom*, int> &atom_map, chemlib::RESIDUE *res, std::vector<chemlib::Bond> &bonds);
        bool SquirtPlanes(FILE*, std::string&, std::map<chemlib::MIAtom*, int>&,
                          std::vector<chemlib::PLANE>&);
        bool CreateMonomerTree(FILE *fp, std::string &resname, chemlib::MIAtom *prev,
                               chemlib::MIAtom *curr, std::map<chemlib::MIAtom*, int> &atom_map,
                               std::map<chemlib::MIAtom*, bool> &alreadyprinted, std::vector<chemlib::Bond> &bonds,
                               std::map<chemlib::Bond*, bool> &bond_map, bool squirtend);

        bool SlurpAtoms(CifLoop &loop, chemlib::RESIDUE *res);
        bool SlurpBonds(CifLoop &loop, std::map<std::string, chemlib::MIAtom* > &atom_map, std::vector<chemlib::Bond> &bonds);
        bool SlurpAngles(CifLoop &loop, std::map<std::string, chemlib::MIAtom* > &atom_map, std::vector<chemlib::ANGLE> &angles, chemlib::RESIDUE *res);
        bool SlurpTorsions(CifLoop &loop, std::map<std::string, chemlib::TORSDICT> &torsion_map);
        bool SlurpPlanes(CifLoop &loop, std::map<std::string, chemlib::PLANEDICT> &plane_map);
        bool SlurpChirals(CifLoop &loop, std::map<std::string, chemlib::CHIRALDICT> &chiral_map);
        bool SlurpHeader(CifLoop &loop, chemlib::RESIDUE *res);

        //Some helper functions to gather data from subordinate loops
        bool SlurpTorsionValues(CifLoop &loop, std::map<std::string, chemlib::TORSDICT > &torsion_map);
        bool SlurpPlaneAtoms(CifLoop &loop, std::map<std::string, chemlib::PLANEDICT > &plane_map);
        bool SlurpChiralAtoms(CifLoop &loop, std::map<std::string, chemlib::CHIRALDICT > &chiral_map);

    public:
        mmCIF();
        virtual ~mmCIF();
        void WriteTorsions(bool);
        virtual std::string getDescription()
        {
            return "Refmac mmCIF Dictionary File (*.cif)";
        }

        virtual std::string getExtension()
        {
            return "*.cif";
        }

        virtual bool Write(FILE *fp, MIMolInfo &mol);
        virtual bool Read(FILE *fp, MIMolInfo &mol);
    };

    class MITorsionWritePrompt
    {
    public:
        virtual ~MITorsionWritePrompt()
        {
        }

        virtual bool operator()() = 0;
    };

    MITorsionWritePrompt *MIGetTorsionWritePrompt();
    void MISetTorsionWritePrompt(MITorsionWritePrompt *p);

}

#endif // ifndef mifit_fios_mmCIF_h

