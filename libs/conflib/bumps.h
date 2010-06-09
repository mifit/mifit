#ifndef CHECK_MANAGE_BUMP_H
#define CHECK_MANAGE_BUMP_H

#include <vector>

namespace chemlib
{
    class MIAtom;
    class Bond;
    class Residue;
}

namespace conflib
{
//
    bool AssignBump(chemlib::MIAtom *atom1, chemlib::MIAtom *atom2, chemlib::Bond &bump, const std::vector<chemlib::Bond> &bonds);
    void GetBumps(const chemlib::Residue *res, std::vector<chemlib::Bond> &bumps, const std::vector<chemlib::Bond> &bonds);
    bool CheckBumps(chemlib::Residue *res, std::vector<chemlib::Bond> &bumps);
}

#endif //CHECK_MANAGE_BUMP_H
