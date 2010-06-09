#ifndef FLEX_TORSION_H
#define FLEX_TORSION_H

#include <vector>

namespace chemlib
{
    class MIAtom;
    class Bond;
    class Residue;
    class TORSION;
}
//Torsion class to: determine ideal angles, store the atoms
//under the torsion, and rotate the atoms to the next angle.  These are needed to
//generate ensembles of conformations for the new fitting algorithms.  -KWB 10/6/5

namespace conflib
{

    class FlexTorsion
    {
    public:
        void Set();                 //Set torsion to current ideal angle
        void Pick(int index);       //Set torsion to Nth ideal angle (random access)

        FlexTorsion(const chemlib::TORSION &t, const std::vector<chemlib::Bond> &bonds, bool set = false);

        FlexTorsion(const FlexTorsion &rhs);
        FlexTorsion&operator=(const FlexTorsion &rhs);
        //		void GatherFlexAtoms(const AtomGraphIndex *distances, const chemlib::RESIDUE *res);
        void SetIdeal(int n, const float *angles);
        void GuessIdeal();
        bool Advance();
        float Measure();
        int NumAngles() const
        {
            return _dihedrals.size();
        }

    protected:
        chemlib::MIAtom *_a1;
        chemlib::MIAtom *_a2;
        chemlib::MIAtom *_a3;
        chemlib::MIAtom *_a4;
        std::vector<float> _dihedrals;
        std::vector<float>::iterator _dihedral_index;
        std::vector<chemlib::MIAtom*> _flex_atoms;
        std::vector<chemlib::Bond>::const_iterator _bond;
    };

    chemlib::MIAtom *SetUpTorsion(const chemlib::Residue *res,
                                  const std::vector<chemlib::Bond> &bonds,
                                  chemlib::MIAtom *atom1,
                                  chemlib::MIAtom *atom2,
                                  std::vector<chemlib::MIAtom*> &flex_atoms);
} //namespace conflib

#endif //FLEX_TORSION_H
