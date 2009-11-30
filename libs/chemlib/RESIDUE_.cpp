#include <nongui/nonguilib.h>
#include <math/mathlib.h>
#include "model.h"

#include "atom_util.h"
#include "mol_util.h"
#include <chemlib/RESIDUE_.h>

namespace chemlib
{

RESIDUE*RESIDUE::insertResidue(RESIDUE *residue)
{
    MI_ASSERT(residue != NULL);
    if (residue == NULL)
    {
        return NULL;
    }
    MI_ASSERT(residue->prev_res == NULL);
    RESIDUE *formerNext = next_res;
    next_res = residue;
    residue->prev_res = this;
    RESIDUE *residueTail = residue;
    while (residueTail->next_res != NULL)
    {
        residueTail = residueTail->next_res;
    }
    if (formerNext != NULL)
    {
        formerNext->prev_res = residueTail;
    }
    residueTail->next_res = formerNext;
    return residueTail;
}

RESIDUE*RESIDUE::removeFromList(RESIDUE *toResidue)
{
    if (toResidue == NULL)
    {
        toResidue = this;
    }
    else
    {
#if DEBUG
        // Ensure that toResidue follows this residue in list
        RESIDUE *res = this;
        while (res != NULL && res != toResidue)
        {
            res = res->next_res;
        }
        MI_ASSERT(res == toResidue);
#endif
    }
    if (prev_res != NULL)
    {
        prev_res->next_res = toResidue->next_res;
    }
    if (toResidue->next_res != NULL)
    {
        toResidue->next_res->prev_res = prev_res;
    }
    RESIDUE *result = toResidue->next_res;
    prev_res = NULL;
    toResidue->next_res = NULL;
    return result;
}

void FreeResidueList(RESIDUE *reslist)
{
    RESIDUE *res;
    while (Residue::isValid(reslist))
    {
        res = reslist->next();
        delete reslist;
        reslist = res;
    }
}

int CountAtomsByName(const char *atom_name, const RESIDUE *res)
{
    if (!Residue::isValid(res))
    {
        return 0;
    }
    int n = 0;
    for (int i = 0; i < res->atomCount(); ++i)
    {
        if (strncmp(atom_name, res->atom(i)->name(), chemlib::MAXATOMNAME) == 0)
        {
            n++;
        }
    }
    return n;
}

int DupeAtomNames(const RESIDUE *res)
{
    if (!Residue::isValid(res))
    {
        return 0;
    }
    int n = 0;
    for (int i = 0; i < res->atomCount(); ++i)
    {
        if (CountAtomsByName(res->atom(i)->name(), res) != 1)
        {
            n++;
        }
    }
    return n;
}

RESIDUE *CopyResList(const RESIDUE *oldres)
{
    RESIDUE *current = NULL;
    RESIDUE *prev = NULL;
    RESIDUE *start = NULL;

    while (Residue::isValid(oldres))
    {
        current = new RESIDUE(*oldres);
        if (prev)
        {
            prev->setNext(current);
        }
        else
        {
            start = current;
        }
        oldres = oldres->next();
        prev = current;
    }

    if (current)
    {
        current->setNext(NULL);
    }
    return start;
}

RESIDUE::RESIDUE()
    : Residue(),
      next_res(NULL),
      prev_res(NULL)
{
}

RESIDUE::RESIDUE(const RESIDUE &lhs)
    : Residue(lhs),
      next_res(NULL),
      prev_res(NULL)
{

    // Cache atoms from both residues for migrating prefbonds and prefangles
    MIAtomList old_atoms;
    MIAtomList new_atoms;
    for (int i = 0; i < lhs.atomCount(); i++)
    {
        old_atoms.push_back(lhs.atoms_[i]);
        new_atoms.push_back(atoms_[i]);
    }

    // copy and update prefbonds and prefangles if necessary
    for (unsigned int i = 0; i < lhs.prefbonds.size(); ++i)
    {
        Bond b(lhs.prefbonds[i], old_atoms, new_atoms);
        prefbonds.push_back(b);
    }

    for (unsigned int i = 0; i < lhs.prefangles.size(); ++i)
    {
        int i1 = GetIndex<MIAtom*>(lhs.prefangles[i].getAtom1(), old_atoms);
        int i2 = GetIndex<MIAtom*>(lhs.prefangles[i].getAtom2(), old_atoms);
        int i3 = GetIndex<MIAtom*>(lhs.prefangles[i].getAtom3(), old_atoms);
        ANGLE ang = lhs.prefangles[i];

        if (i1 < 0 || i2 < 0 || i3 < 0)
        {
            std::string error_message = "Corrupted angle.";
            throw error_message;
        }
        ang.atom1 = new_atoms[i1];
        ang.atom2 = new_atoms[i2];
        ang.atom3 = new_atoms[i3];
        prefangles.push_back(ang);
    }
}

RESIDUE&RESIDUE::operator=(RESIDUE other)
{
    using std::swap;

    swap(*this, other);
    return *this;
}

RESIDUE::~RESIDUE()
{
}

RESIDUE::RESIDUE(const Residue &r)
    : Residue(r),
      next_res(NULL),
      prev_res(NULL)
{
}


/////////////////////////////////////////////////////////////////////////////
// Function:  AtomVectMatchesRes
// Purpose:   Checks if the atom pointers in a vector correspond exactly
//        to those in a given RESIDUE object
// Input:       A vector atom pointers
//        A residue ptr
// Output:    True if the pointers in the vector and residue point to the
//        same atoms
// Requires:
/////////////////////////////////////////////////////////////////////////////
bool AtomVectMatchesRes(const MIAtomList &ptrs, const RESIDUE *res)
{
    if (!Residue::isValid(res))
    {
        return false;
    }
    if (ptrs.size() != (unsigned int)res->atomCount())
    {
        return false;
    }

    int i;
    for (i = 0; i < res->atomCount(); ++i)
    {
        if (std::count(ptrs.begin(), ptrs.end(), res->atom(i)) != 1)
        {
            return false;
        }
    }
    return true;
}

void PrintResidueDebugInfo(RESIDUE *res1)
{
    if (!Residue::isValid(res1))
    {
        return;
    }
    printf("Res1 (%s) contains %d atoms\n", res1->name().c_str(), res1->atomCount());
    for (int i = 0; i < res1->atomCount(); ++i)
    {
        printf("\tatom %d name %s pointer %p\n", i, res1->atom(i)->name(), res1->atom(i));
    }
    printf("prefbonds size %d prefangles size %d\n", res1->prefBonds().size(), res1->prefAngles().size());
    for (unsigned int i = 0; i < res1->prefBonds().size(); ++i)
    {
        printf("bond %d is between %p and %p\n",
               i, res1->prefBonds()[i].getAtom1(), res1->prefBonds()[i].getAtom2());
    }
}

chemlib::MIAtom *atom_default(const chemlib::RESIDUE *res)
{
    if (!Residue::isValid(res))
    {
        return NULL;
    }
    chemlib::MIAtom *a = atom_from_name("CA", *res);
    if (a)
    {
        return a;
    }
    else
    {
        return res->atom(0);
    }
}

static void fixname(char *name)
{
    // remove white space in name strings
    size_t j, l;
    while (isspace(name[0]) )
    {
        l = strlen(name);
        for (j = 0; j < l; j++)
        {
            name[j] = name[j+1];
        }
    }
    while (isspace(name[strlen(name)-1]))
    {
        name[strlen(name)-1] = '\0';
    }
}

static void fixname(std::string &name)
{
    char buf[MAXNAME];
    strncpy(buf, name.c_str(), MAXNAME);
    fixname(buf);
    name = std::string(buf);
}

void RESIDUE::fixnames(RESIDUE *res)
{
    fixname(res->name_);
    fixname(res->type_);
    for (int i = 0; i < res->atomCount(); i++)
    {
        res->atom(i)->fixname();
    }
}

const std::string RESIDUE::liststring(RESIDUE *res)
{
    char chainid;
    int chainno;
    std::string buff;
    chainid = (char)(res->chain_id() & 255);
    if (chainid == ' ')
    {
        chainid = '_';
    }
    chainno = res->chain_id()/256;
    buff = format("%s %s %c %d",
                  res->type().c_str(), res->name().c_str(), chainid, chainno);
    if ((res->linkage_type()&NTERMINUS) && (res->linkage_type()&PEPTIDE))
    {
        buff += " Nter";
    }
    if (res->linkage_type()&CTERMINUS)
    {
        if (res->linkage_type()&PEPTIDE)
        {
            buff += " Cter";
        }
        else
        {
            buff += " Ter";
        }
    }
    return (buff);
}

}
