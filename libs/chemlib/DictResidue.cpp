#include <chemlib/RESIDUE_.h>
#include "Bond.h"
#include "ANGLE.h"
#include "DictResidue.h"
#include "atom_util.h"

using namespace chemlib;
using namespace std;

DictResidue::DictResidue(RESIDUE *r)
{
    residue = r; // just keep a reference to the residue, not a copy

    if (residue->prefBonds().size() == 0 && residue->prefAngles().size() == 0)
    {
        Build();
    }
    else if (residue->prefAngles().size() == 0)
    {
        bonds.clear();
        bonds.assign(r->prefBonds().begin(), r->prefBonds().end());
        BuildAngles();
    }
    else if (residue->prefBonds().size() == 0)
    {
        angles.clear();
        angles.assign(r->prefAngles().begin(), r->prefAngles().end());
        BuildBonds();
    }
    else
    {
        bonds.clear();
        bonds.assign(r->prefBonds().begin(), r->prefBonds().end());
        angles.clear();
        angles.assign(r->prefAngles().begin(), r->prefAngles().end());
    }
}

DictResidue::~DictResidue()
{
}

void DictResidue::BuildAngles()
{

    if (!Residue::isValid(residue))
    {
        return;
    }

    angles.clear();

    MIAtom *a1, *a2, *a3;
    MIAtom *b1, *b2, *b3, *b4;
    ANGLE angle;
    /* search for angles - ends of two joined bonds */
    for (unsigned int i = 0; i < bonds.size(); i++)
    {
        b1 = bonds[i].getAtom1();
        b2 = bonds[i].getAtom2();
        for (unsigned int j = i+1; j < bonds.size(); j++)
        {
            b3 = bonds[j].getAtom1();
            b4 = bonds[j].getAtom2();
            a1 = a2 = a3 = NULL;
            /* see if these two bonds share a
             * common atom */
            if (b1 == b3)
            {
                a2 = b1;
                a1 = b2;
                a3 = b4;
            }
            else if (b1 == b4)
            {
                a2 = b1;
                a1 = b2;
                a3 = b3;
            }
            else if (b2 == b3)
            {
                a2 = b2;
                a1 = b1;
                a3 = b4;
            }
            else if (b2 == b4)
            {
                a2 = b2;
                a1 = b1;
                a3 = b3;
            }
            if (!a1 || !a2 || !a3)
            {
                continue;
            }
            /* find these 3 atoms in reslist */
            angle.setAtom1(a1);
            angle.setAtom2(a2);
            angle.atom3 = a3;
            angle.ideal_angle = (float)AtomDist(*a1, *a3);
            angle.tolerance = 0.03F;
            angle.res = residue;
            angles.push_back(angle);
        }
    }
}

void DictResidue::BuildBonds()
{
    Bond bond;
    unsigned int i, j;
    double dlimit, d;
    MIAtom *atom1, *atom2;

    bonds.clear();
    if (!Residue::isValid(residue))
    {
        return;
    }
    for (i = 0; i < (unsigned int)residue->atomCount(); i++)
    {
        atom1 = residue->atom(i);
        for (j = i+1; j < (unsigned int)residue->atomCount(); j++)
        {
            atom2 = residue->atom(j);
            // check to see if distance bondable
            dlimit = BondLimit(atom1->name())+ BondLimit(atom2->name());
            if ((d = AtomDist(*atom1, *atom2)) < dlimit)
            {
                bond.setAtom1(atom1);
                bond.setAtom2(atom2);
                bond.ideal_length = (float)d;
                bond.tolerance = 0.015F;
                bonds.push_back(bond);
            }
        }
    }
}

bool DictResidue::Build()
{
    Bond bond;
    unsigned int i, j;
    double dlimit, d;
    MIAtom *atom1, *atom2;

    bonds.clear();
    angles.clear();
    if (!Residue::isValid(residue))
    {
        return false;
    }
    for (i = 0; i < (unsigned int)residue->atomCount(); i++)
    {
        atom1 = residue->atom(i);
        for (j = i+1; j < (unsigned int) residue->atomCount(); j++)
        {
            atom2 = residue->atom(j);
            // check to see if distance bondable
            dlimit = BondLimit(atom1->name())+ BondLimit(atom2->name());
            if ((d = AtomDist(*atom1, *atom2)) < dlimit)
            {
                bond.setAtom1(atom1);
                bond.setAtom2(atom2);
                bond.ideal_length = (float)d;
                bond.tolerance = 0.015F;
                bonds.push_back(bond);
            }
        }
    }
    MIAtom *a1, *a2, *a3;
    MIAtom *b1, *b2, *b3, *b4;
    ANGLE angle;
    /* search for angles - ends of two joined bonds */
    for (i = 0; i < bonds.size(); i++)
    {
        b1 = bonds[i].getAtom1();
        b2 = bonds[i].getAtom2();
        for (j = i+1; j < bonds.size(); j++)
        {
            b3 = bonds[j].getAtom1();
            b4 = bonds[j].getAtom2();
            a1 = a2 = a3 = NULL;
            /* see if these two bonds share a
             * common atom */
            if (b1 == b3)
            {
                a2 = b1;
                a1 = b2;
                a3 = b4;
            }
            else if (b1 == b4)
            {
                a2 = b1;
                a1 = b2;
                a3 = b3;
            }
            else if (b2 == b3)
            {
                a2 = b2;
                a1 = b1;
                a3 = b4;
            }
            else if (b2 == b4)
            {
                a2 = b2;
                a1 = b1;
                a3 = b3;
            }
            if (!a1 || !a2 || !a3)
            {
                continue;
            }
            /* find these 3 atoms in reslist */
            angle.setAtom1(a1);
            angle.setAtom2(a2);
            angle.atom3 = a3;
            angle.ideal_angle = (float)AtomDist(*a1, *a3);
            angle.tolerance = 0.03F;
            angle.res = residue;
            angles.push_back(angle);
        }
    }
    return true;
}

RESIDUE*DictResidue::Residue()
{
    return residue;
}

const RESIDUE*DictResidue::Residue() const
{
    return residue;
}

std::vector<Bond>*DictResidue::Bonds()
{
    return &bonds;
}

std::vector<ANGLE>*DictResidue::Angles()
{
    return &angles;
}

