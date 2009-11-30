#include <math/mathlib.h>

#include "CovalentGeom.h"
#include "atom_util.h"
#include "math_util.h"
#include "mol_util.h"
#include <chemlib/RESIDUE_.h>
#include "Matrix.h"

using namespace std;

namespace chemlib
{

int GetAtomIndex(const MIAtom *patm, const RESIDUE &res)
{
    int i = 0;

    while (i < res.atomCount())
    {
        if (patm == res.atom(i))
        {
            return i;
        }
        ++i;
    }
    return -1;
}

void GetNabors(const MIAtom *atom,
               const std::vector<Bond> &bonds,
               MIAtomList &nabors)
{

    std::vector<Bond>::const_iterator bnd;

    bnd = bonds.begin();
    while (bnd != bonds.end())
    {
        if (bnd->getAtom1() == atom)
        {
            nabors.push_back(bnd->getAtom2());
        }
        if (bnd->getAtom2() == atom)
        {
            nabors.push_back(bnd->getAtom1());
        }
        ++bnd;
    }
}

bool AlreadyBonded(const MIAtom *atom1,
                   const MIAtom *atom2,
                   const std::vector<Bond> &bonds)
{
    vector<Bond>::const_iterator i = bonds.begin();
    vector<Bond>::const_iterator e = bonds.end();

    while (i != e)
    {
        if ((i->getAtom1() == atom1 && i->getAtom2() == atom2)
            || (i->getAtom1() == atom2 && i->getAtom2() == atom1))
        {
            return true;
        }
        ++i;
    }
    return false;
}

float AngleFromGeom(int geometry)
{
    if (geometry == LINEAR)
    {
        return 180.0F;
    }
    else if (geometry == TRIGONAL_PLANAR)
    {
        return 120.0F;
    }
    else if (geometry == TETRAHEDRAL)
    {
        return 109.47F;
    }
    else if (geometry == SQUARE_PLANAR)
    {
        return 90.0F;
    }
    else if (geometry == TRIGONAL_BIPYRAMIDAL)
    {
        return 109.47F;
    }
    else if (geometry == OCTAHEDRAL)
    {
        return 90.0F;
    }
    return 0;
}

float CalcBestHDistance(const MIAtom *donor,
                        const MIAtom *acceptor,
                        const MIAtom *nabor,
                        double c)                                   //Length of h-bond

{
    double b = AtomDist(*donor, *acceptor);

    //		double c = conflib::IdealBondLength(donor->atomicnumber, 1, 1);		//Length of single
    //bond to Hydrogen

    double angle = AngleFromGeom(TETRAHEDRAL)
                   -CalcAtomAngle(*nabor, *donor, *acceptor);                  //Can be neg, since
                                                                               //we take the cos
    return (float)sqrt(b*b + c*c - 2*b*c*cos(DEG2RAD*angle));
}

//Given an atom, returns a vector with the idealized position for an additional bond to be
//added.  If the vector returned is very short, this indicates that the valences have already
//been filled.  (e.g. a trigonal planar atom that already has 3 bonds, each 120 degrees apart.)
bool DirectNextBond(const MIAtom *atom, const std::vector<Bond> &bonds, double *v)
{
    if (v == 0)
    {
        return false;
    }
    v[0] = v[1] = v[2] = 0.0;
    MIAtomList nabors;
    GetNabors(atom, bonds, nabors);
    MIAtom_const_iter nbr;
    //		double bv[3];
    for (nbr = nabors.begin(); nbr != nabors.end(); ++nbr)
    {
        v[0] += atom->x() - (*nbr)->x();
        v[1] += atom->y() - (*nbr)->y();
        v[2] += atom->z() - (*nbr)->z();
    }


    //		std::accumulate(nabors.begin(), nabors.end(), v, SumBondVectors);
    //		ScaleVect(v, 1.0/nabors.size());

    return true;
}

//Given an atom, returns a vector with the idealized position for an additional bond to be
//added.  If the vector returned is very short, this indicates that the valences have already
//been filled.  (e.g. a trigonal planar atom that already has 3 bonds, each 120 degrees apart.)
bool DirectNextBond(const MIAtom *atom, double *v)
{
    if (v == 0)
    {
        return false;
    }
    v[0] = v[1] = v[2] = 0.0;
    MIAtom_const_iter nbr;

    for (nbr = atom->nabors().begin(); nbr != atom->nabors().end(); ++nbr)
    {
        v[0] += atom->x() - (*nbr)->x();
        v[1] += atom->y() - (*nbr)->y();
        v[2] += atom->z() - (*nbr)->z();
    }


    //		std::accumulate(nabors.begin(), nabors.end(), v, SumBondVectors);
    //		ScaleVect(v, 1.0/nabors.size());

    return true;
}

MIAtom *GetDoublePartner(const MIAtom *atom, std::vector<Bond> &bonds)
{

    float max_elneg = -1.0f;
    float min_bond_length = 2000.0F;
    MIAtom *partner = 0;

    for (unsigned int i = 0; i < atom->nabors().size(); ++i)
    {
        Bond *bond = &bonds[atom->bondnumbers()[i]];
        MIAtom *nabor = atom->nabors()[i];

        float elneg = (float)Electroneg(nabor->atomicnumber());
        float bond_length = (float)AtomDist(*atom, *nabor);

        if (bond->getOrder() != SINGLEBOND)             //skip if we've already set this bond order
        {
            continue;
        }
        if (nabor->hybrid() == 3
            && nabor->nabors().size() != 1)       //can't form double bond with sp3 atom,
        {
            continue;                           //but terminal atoms might not *really* be sp3
        }
        if (UnusedValences(*nabor, bonds) == 0)         //skip if there are no more available valences
        {
            continue;
        }
        if (CurrentValence(*nabor, bonds) != (int)nabor->nabors().size()
            && nabor->hybrid() == 2)              //skip if we've already set a double bnd for an sp2
        {
            continue;                               //neighbor, to avoid allenes
        }
        if (elneg < max_elneg)                      //skip if we've already found a better partner
        {
            continue;
        }
        if (elneg == max_elneg                      //break ties by picking the shorter bond
            && bond_length > min_bond_length)
        {
            continue;
        }

        partner = nabor;                            //If we reach here, choose this atom for
        max_elneg = elneg;                          //now and store "fitness" values for comparison
        min_bond_length = bond_length;
    }

    return partner;
}

MIAtom *GetTriplePartner(const MIAtom *atom, std::vector<Bond> &bonds)
{

    float max_elneg = -1.0f;
    float min_bond_length = 2000.0F;
    MIAtom *partner = 0;

    for (unsigned int i = 0; i < atom->nabors().size(); ++i)
    {
        Bond *bond = &bonds[atom->bondnumbers()[i]];
        MIAtom *nabor = atom->nabors()[i];

        float elneg = (float)Electroneg(nabor->atomicnumber());
        float bond_length = (float)AtomDist(*atom, *nabor);

        if (bond->getOrder() != SINGLEBOND)             //skip if we've already set this bond order
        {
            continue;
        }
        if (nabor->hybrid() != 1
            && nabor->nabors().size() != 1)       //can't form triple bond with anything but an sp,
        {
            continue;                           //but any terminal atom could be an sp
        }
        if (UnusedValences(*nabor, bonds) < 2)          //skip if there aren't enough available valences
        {
            continue;
        }
        if (elneg < max_elneg)                      //skip if we've already found a better partner
        {
            continue;
        }
        if (elneg == max_elneg                      //break ties by picking the shorter bond
            && bond_length > min_bond_length)
        {
            continue;
        }

        partner = nabor;                            //If we reach here, choose this atom for
        max_elneg = elneg;                          //now and store "fitness" values for comparison
        min_bond_length = bond_length;
    }

    return partner;
}

int PositionHydrogens(const MIAtom *atom, MIAtomList &atoms)
{
    double v[3];
    double length = IdealBondLength(atom->atomicnumber(), 1, 1);
    double angle;

    MIAtom *h = new MIAtom;
    h->copyPosition(*atom);

    switch (atom->nabors().size())
    {
    case 0:
        return 0;
    case 1:                                     //Shouldn't get here for hydrogen bonding calcs
        DirectNextBond(atom, v);
        h->translate((float)(length * v[0]),
                     (float)(length * v[1]),
                     (float)(length * v[2]));
        atoms.push_back(h);
        return 1;
    case 2:
        angle = CalcAtomAngle(*atom->nabors()[0],
                              *atom,
                              *atom->nabors()[1]);
        if (angle > 155)
        {
            return 0;
        }
        else if (angle > 116)
        {
            DirectNextBond(atom, v);
            h->translate((float)(length * v[0]),
                         (float)(length * v[1]),
                         (float)(length * v[2]));
            atoms.push_back(h);
            return 1;
        }
        else
        {
            Add2Tetrahedrals(atom, length, atoms);
            return 2;
        }
    case 3:
        DirectNextBond(atom, v);
        h->translate((float)(length * v[0]),
                     (float)(length * v[1]),
                     (float)(length * v[2]));
        atoms.push_back(h);
        return 1;
    default:
        return 0;
    }

}

void Add2Tetrahedrals(const MIAtom *atom, double length, MIAtomList &atoms)
{
    double b1[3], b2[3], h1[3], h2[3];
    double neg_sum[3];
    double cross[3];
    double neg_cross[3];

    b1[0] = atom->nabors()[0]->x() - atom->x();
    b1[1] = atom->nabors()[0]->y() - atom->y();
    b1[2] = atom->nabors()[0]->z() - atom->z();

    b2[0] = atom->nabors()[1]->x() - atom->x();
    b2[1] = atom->nabors()[1]->y() - atom->y();
    b2[2] = atom->nabors()[1]->z() - atom->z();

    NormVect(b1);
    NormVect(b2);

    neg_sum[0] = -(b1[0] + b2[0]);
    neg_sum[1] = -(b1[0] + b2[0]);
    neg_sum[2] = -(b1[0] + b2[0]);

    CrossVects(b1, b2, cross);
    CrossVects(b2, b1, neg_cross);

    NormVect(cross);
    NormVect(neg_cross);

    h1[0] = neg_sum[0] + cross[0];
    h1[1] = neg_sum[1] + cross[1];
    h1[2] = neg_sum[2] + cross[2];

    h2[0] = neg_sum[0] + neg_cross[0];
    h2[1] = neg_sum[1] + neg_cross[1];
    h2[2] = neg_sum[2] + neg_cross[2];

    ScaleVect(h1, length / VectLength(h1));
    ScaleVect(h2, length / VectLength(h2));

    MIAtom *h = new MIAtom;
    h->setPosition((float)(atom->x() + h1[0]),
                   (float)(atom->y() + h1[1]),
                   (float)(atom->z() + h1[2]));
    atoms.push_back(h);

    h = new MIAtom;
    h->setPosition((float)(atom->x() + h2[0]),
                   (float)(atom->y() + h2[1]),
                   (float)(atom->z() + h2[2]));
    atoms.push_back(h);
}

void ClearCharges(const RESIDUE *res)
{
    //		if (Residue::isValid(res)) {
    if (res != NULL)
    {
        const MIAtomList &atoms = res->atoms();
        for (MIAtom_const_iter atomIter = atoms.begin(); atomIter != atoms.end(); ++atomIter)
        {
            (*atomIter)->set_formal_charge(0);
        }
    }
}

void CenterOfMass(const MIAtomList &atoms, double *com)
{
    MIAtom_const_iter i, e = atoms.end();

    com[0] = com[1] = com[2] = 0.0;
    for (i = atoms.begin(); i != e; ++i)
    {
        MIAtom &a = **i;
        com[0] += a.x();
        com[1] += a.y();
        com[2] += a.z();
    }

    com[0] /= (double) atoms.size();
    com[1] /= (double) atoms.size();
    com[2] /= (double) atoms.size();
}

double SignedAtomVolume(const MIAtom &a1, const MIAtom &a2, const MIAtom &a3, const MIAtom &a4)
{
    double x1 = a1.x();
    double y1 = a1.y();
    double z1 = a1.z();
    double x2 = a2.x();
    double y2 = a2.y();
    double z2 = a2.z();
    double x3 = a3.x();
    double y3 = a3.y();
    double z3 = a3.z();
    double x4 = a4.x();
    double y4 = a4.y();
    double z4 = a4.z();

    double det;
    det =  x1*y2*z3 + y1*z2*x3 + z1*x2*y3
          -z1*y2*x3 - y1*x2*z3 - x1*z2*y3
          -x1*y2*z4 - y1*z2*x4 - z1*x2*y4
          +z1*y2*x4 + y1*x2*z4 + x1*z2*y4
          +x1*y3*z4 + y1*z3*x4 + z1*x3*y4
          -z1*y3*x4 - y1*x3*z4 - x1*z3*y4
          -x2*y3*z4 - y2*z3*x4 - z2*x3*y4
          +z2*y3*x4 + y2*x3*z4 + x2*z3*y4;

    return det / 6.0;
}


bool ResContainsBond::operator ()(const Bond &bond, const RESIDUE &res) const
{
    return (find(res.atoms().begin(), res.atoms().end(), bond.getAtom1()) != res.atoms().end())
           && (find(res.atoms().begin(), res.atoms().end(), bond.getAtom2()) != res.atoms().end());
}

float ZByName(const char *name)
{
    float z = 6.7f; /*average for proteins */
    // note a switch is faster than a tree of if statements, since it is implemented as a single jump statement

    switch (name[0])
    {
    case 'C': return (6.0);
        break;
    case 'N': return (7.0);
        break;
    case 'O': return (7.0);
        break;
    case 'S': if (name[1] == 'E')
        {
            return (32.0);
        }
        return (16.0);
        break;
    case 'H': return (1.0);
        break;
    case 'P': if (name[1] == 'T')
        {
            return (78.0);
        }
        return (15.0);
        break;
    case '1': return (1.0);
        break;
    case '2': return (1.0);
        break;
    case '3': return (1.0);
        break;
    case '4': return (1.0);
        break;
    case '5': return (1.0);
        break;
    case 'F': if (name[1] == 'E')
        {
            return (26.0);
        }
        break;
    case 'M': if (name[1] == 'N')
        {
            return (25.0);
        }
        break;
    }
    return (z);
}

void renameResidueAtomsToUnique(const Residue *res)
{
    for (int i = 0; i < res->atomCount(); ++i)
    {
        MIAtom *atom = res->atom(i);
        const char *atomicSymbol = chemlib::Left_Atomic_Name(atom->atomicnumber());
        atom->setName(format("%s%d", atomicSymbol, i+1).c_str());
    }
}

} //namespace chemlib
