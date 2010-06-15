#ifndef ATOM_UTIL_H
#define ATOM_UTIL_H

#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <utility>

#include "MIAtom_fwd.h"
#include "Bond.h"
#include "Residue_fwd.h"
#include "mol_util.h"


namespace chemlib
{

//Queries an atom as to whether its name matches the given string
    struct MatchesAtomName : public std::binary_function<MIAtom*, std::string, bool>
    {
        bool operator()(const MIAtom *atom, const std::string &query) const;
    };

    inline bool MatchesAtomName::operator ()(const MIAtom *atom, const std::string &query) const
    {
        //		return query.CmpNoCase(atom.name) != 0;
        return atom->name() == query;
    }

//Queries a bond as to whether it contains the given atom
    struct ContainsAtom : public std::binary_function<Bond, const MIAtom*, bool>
    {
        bool operator()(const Bond &bond, const MIAtom *atom) const;
    };

    inline bool ContainsAtom::operator ()(const Bond &bond, const MIAtom *atom) const
    {
        return bond.getAtom1() == atom || bond.getAtom2() == atom;
    }

//Queries a bond as to whether it contains the given atom pair
    struct ContainsAtoms : public std::binary_function<Bond, std::pair<const MIAtom*, MIAtom*>, bool>
    {
        bool operator()(const Bond &bond, std::pair<const MIAtom*, const MIAtom*> atom_pair) const;
    };

    inline bool ContainsAtoms::operator ()(const Bond &bond,
                                           std::pair<const MIAtom*, const MIAtom*> atom_pair) const
    {
        return (bond.getAtom1() == atom_pair.first && bond.getAtom2() == atom_pair.second)
               || (bond.getAtom1() == atom_pair.second && bond.getAtom2() == atom_pair.first);
    }

    struct ContainsCoAtoms : public std::binary_function<Bond, std::pair<const MIAtom*, const MIAtom*>, bool>
    {
        bool operator()(const Bond &bond, std::pair<const MIAtom*, const MIAtom*> atom_pair) const;
    };

    inline bool ContainsCoAtoms::operator ()(const Bond &bond,
                                             std::pair<const MIAtom*, const MIAtom*> atom_pair) const
    {
        return (bond.getAtom1() == atom_pair.first && bond.getAtom2() == atom_pair.second)
               || (bond.getAtom1() == atom_pair.second && bond.getAtom2() == atom_pair.first);
    }

//Queries a bond as to whether it contains the given atom
    struct ContainsNewAtom : public std::binary_function<Bond, const MIAtom*, bool>
    {
        bool operator()(const Bond &bond, const MIAtom *atom) const;
    };

    inline bool ContainsNewAtom::operator ()(const Bond &bond, const MIAtom *atom) const
    {
        return bond.getAtom1() == atom || bond.getAtom2() == atom;
    }

//Queries a bond as to whether it contains the given atom & another heavy atom
    struct ContainsHvyNabor : public std::binary_function<Bond, const MIAtom*, bool>
    {
        inline bool operator()(const Bond &bond, const MIAtom *atom) const
        {
            return (bond.getAtom1() == atom && bond.getAtom2()->atomicnumber() != 1)
                   || (bond.getAtom2() == atom && bond.getAtom1()->atomicnumber() != 1);
        }

    };

//Returns a vector of atoms that are bonded to the give atom
    void GetNabors(const MIAtom *atom,
                   const std::vector<Bond> &bonds,
                   MIAtomList &nabors);

    struct IsHydrogen : std::unary_function<const MIAtom*, bool>
    {
        inline bool operator()(const MIAtom *atom) const
        {
            return atom->atomicnumber() == 1;
        }

    };

    inline int Degree(const MIAtom *atom)
    {
        return atom->nabors().size() + atom->hcount();
    }

    inline int Degree(const MIAtom *atom, const std::vector<Bond> bonds)
    {
        return std::count_if(bonds.begin(), bonds.end(), std::bind2nd(ContainsAtom(), atom));
    }

    inline int HeavyDegree(const MIAtom *atom)
    {
        return std::count_if(atom->nabors().begin(),
                             atom->nabors().end(),
                             std::not1(IsHydrogen()));
    }

    inline int HeavyDegree(const MIAtom *atom, const std::vector<Bond> &bonds)
    {
        return std::count_if(bonds.begin(),
                             bonds.end(),
                             std::bind2nd(ContainsHvyNabor(), atom));
    }

    bool AlreadyBonded(const MIAtom *atom1,
                       const MIAtom *atom2,
                       const std::vector<Bond> &bonds);

    inline std::vector<Bond>::iterator
    GetBond(const MIAtom *atom1,
            const MIAtom *atom2,
            std::vector<Bond> &bonds)
    {
        return std::find_if(bonds.begin(),
                            bonds.end(),
                            std::bind2nd(ContainsCoAtoms(), std::pair<const MIAtom*, const MIAtom*>(atom1, atom2)));
    }

    inline std::vector<Bond>::const_iterator
    GetBond(const MIAtom *atom1,
            const MIAtom *atom2,
            const std::vector<Bond> &bonds)
    {
        return std::find_if(bonds.begin(),
                            bonds.end(),
                            std::bind2nd(ContainsCoAtoms(), std::pair<const MIAtom*, const MIAtom*>(atom1, atom2)));
    }

    float AngleFromGeom(int geometry);


    float CalcBestHDistance(const MIAtom *donor,
                            const MIAtom *acceptor,
                            const MIAtom *nabor,
                            double bondlength);


    inline double*
    SumBondVectors(double *v, const MIAtom *atom)
    {
        v[0] += atom->x();
        v[1] += atom->y();
        v[2] += atom->z();

        return v;
    }

    bool DirectNextBond(const MIAtom *atom, const std::vector<Bond> &bonds, double *v);

    inline int GetNumHydrogens(const MIAtom &atom)
    {
        return std::count_if(atom.nabors().begin(), atom.nabors().end(), IsHydrogen())
               + atom.hcount();
    }

//	CanTakeDouble(const MIAtom &atom, const std::vector<Bond> &bonds);
//	CanTakeTriple(const MIAtom &atom, const std
    MIAtom *GetDoublePartner(const MIAtom *atom, std::vector<Bond> &bonds);
    MIAtom *GetTriplePartner(const MIAtom *atom, std::vector<Bond> &bonds);

    int PositionHydrogens(const MIAtom *atom, MIAtomList &atoms);
    void Add2Tetrahedrals(const MIAtom *atom, double length, MIAtomList &atoms);

    int GetAtomIndex(const MIAtom *patm, const Residue &res);

//Queries a residue as to whether it contains the given bond
    struct ResContainsBond : public std::binary_function<Bond, const Residue, bool>
    {
        bool operator()(const Bond &bond, const Residue &res) const;
    };


    void ClearCharges(const Residue *res);

    void CenterOfMass(const MIAtomList &atoms, double *com);

    double SignedAtomVolume(const MIAtom &a1, const MIAtom &a2, const MIAtom &a3, const MIAtom &a4);

    inline double VolAtomGradNorm(MIAtom &a1, MIAtom &a2, MIAtom &a3, MIAtom &a4)
    {
        double x1 = (a2.y()*a3.z() - a2.z()*a3.y() - (a2.y()*a4.z() - a2.z()*a4.y()) + a3.y()*a4.z() - a3.z()*a4.y());
        double y1 = (a2.x()*a3.z() - a2.z()*a3.x() - (a2.x()*a4.z() - a2.z()*a4.x()) + a3.x()*a4.z() - a3.z()*a4.x());
        double z1 = (a2.x()*a3.y() - a2.y()*a3.x() - (a2.x()*a4.y() - a2.y()*a4.x()) + a3.x()*a4.y() - a3.y()*a4.x());
        double x2 = (a1.y()*a3.z() - a1.z()*a3.y() - (a1.y()*a4.z() - a1.z()*a4.y()) + a3.y()*a4.z() - a3.z()*a4.y());
        double y2 = (a1.x()*a3.z() - a1.z()*a3.x() - (a1.x()*a4.z() - a1.z()*a4.x()) + a3.x()*a4.z() - a3.z()*a4.x());
        double z2 = (a1.x()*a3.y() - a1.y()*a3.x() - (a1.x()*a4.y() - a1.y()*a4.x()) + a3.x()*a4.y() - a3.y()*a4.x());
        double x3 = (a1.y()*a2.z() - a1.z()*a2.y() - (a1.y()*a4.z() - a1.z()*a4.y()) + a2.y()*a4.z() - a2.z()*a4.y());
        double y3 = (a1.x()*a2.z() - a1.z()*a2.x() - (a1.x()*a4.z() - a1.z()*a4.x()) + a2.x()*a4.z() - a2.z()*a4.x());
        double z3 = (a1.x()*a2.y() - a1.y()*a2.x() - (a1.x()*a4.y() - a1.y()*a4.x()) + a2.x()*a4.y() - a2.y()*a4.x());
        double x4 = (a1.y()*a2.z() - a1.z()*a2.y() - (a1.y()*a3.z() - a1.z()*a3.y()) + a2.y()*a3.z() - a2.z()*a3.y());
        double y4 = (a1.x()*a2.z() - a1.z()*a2.x() - (a1.x()*a3.z() - a1.z()*a3.x()) + a2.x()*a3.z() - a2.z()*a3.x());
        double z4 = (a1.x()*a2.y() - a1.y()*a2.x() - (a1.x()*a3.y() - a1.y()*a3.x()) + a2.x()*a3.y() - a2.y()*a3.x());

        return (x1 * x1 + y1 * y1 + z1 * z1 + x2 * x2 + y2 * y2 + z2 * z2
                +x3 * x3 + y3 * y3 + z3 * z3 + x4 * x4 + y4 * y4 + z4 * z4) / 36.0;
    }

    float ZByName(const char *name);

    void renameResidueAtomsToUnique(const Monomer *res);

} //namespace chemlib

#endif //ATOM_UTIL_H
