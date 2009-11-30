#ifndef HPHOBE_H
#define HPHOBE_H

#include <chemlib/chemlib.h>


#include <vector>

#define MAX_HPHOBE_LENGTH 3.8F

namespace moldraw
{

    class HPhobe
    {
    public:
        chemlib::MIAtom *atom1;
        chemlib::MIAtom *atom2;
        chemlib::Residue *res1;
        chemlib::Residue *res2;
        double distance;

        chemlib::MIAtom *getAtom1() const
        {
            return atom1;
        }

        void setAtom1(chemlib::MIAtom *atom)
        {
            atom1 = atom;
        }

        chemlib::MIAtom *getAtom2() const
        {
            return atom2;
        }

        void setAtom2(chemlib::MIAtom *atom)
        {
            atom2 = atom;
        }

    };

    void GetNonPolars(const std::vector<chemlib::Residue*> &residues, std::vector<chemlib::MIAtom*> &nonpolar_atoms);
    bool CheckHPhobe(chemlib::MIAtom *atom1, chemlib::MIAtom *atom2, HPhobe &hphobe);

    void RePointHphobe(HPhobe &hphobe, chemlib::Ligand &lig);

    inline bool IsNonPolar(const chemlib::MIAtom &atom)
    {
        if (atom.atomicnumber() == 6)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

} //namespace moldraw
#endif //HPHOBE_H
