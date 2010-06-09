#ifndef MI_RESIDUE_H
#define MI_RESIDUE_H

#include <string>

#include "RESIDUE_fwd.h"
#include "MIAtom_fwd.h"
#include "Bond.h"
#include "ANGLE.h"
#include "Residue.h"
#include "iterator.h"
#include <util/utillib.h>

namespace chemlib
{

#define MAXNAME 6


    class Residue : public Monomer
    {
    public:
        Residue();
        ~Residue();

        //copy ctor, only copies single residue (i.e. doesn't follow next())
        Residue(const Residue &lhs);
        Residue&operator=(Residue other);

        Residue(const Monomer &lhs);

        Residue *next() const
        {
            return next_res;
        }
        void setNext(Residue *r)
        {
            next_res = r;
        }
        Residue *prev() const
        {
            return prev_res;
        }
        void setPrev(Residue *r)
        {
            prev_res = r;
        }


        const std::vector<Bond>&prefBonds() const
        {
            return prefbonds;
        }
        void addPrefBond(const Bond &b)
        {
            prefbonds.push_back(b);
        }
        void setPrefBonds(const std::vector<Bond> &bonds)
        {
            prefbonds = bonds;
        }
        void clearPrefBonds()
        {
            prefbonds.clear();
        }

        const std::vector<ANGLE>&prefAngles() const
        {
            return prefangles;
        }
        void addPrefAngle(const ANGLE &a)
        {
            prefangles.push_back(a);
        }
        void setPrefAngles(const std::vector<ANGLE> &angles)
        {
            prefangles = angles;
        }
        void clearPrefAngles()
        {
            prefangles.clear();
        }

        static MIIterBase<Residue> *getIter(Residue *start)
        {
            return new MIDoublyLinkedListIter<Residue>(start);
        }

        /**
         * Inserts a residue after this residue.
         * Inserting a NULL is ignored in order to prevent dangling lists.
         * Inserting multiple residues is supported, if the residue passed
         * is the head of the list. If the residue contains a link to a
         * previous residue (namely is not the head of the list), an
         * MI_ASSERT failure occurs.
         *
         * @return residue or the tail of residue if a list
         */
        Residue *insertResidue(Residue *residue);

        /**
         * Removes this residue to the given residue from the list.
         * If the toResidue is NULL, removes just this residue.
         * In DEBUG, an MI_ASSERT failure occurs if the toResidue is not
         * found in the list after this residue.
         *
         * @return the next residue in the list
         */
        Residue *removeFromList(Residue *toResidue = NULL);

        static void fixnames(Residue *res);

        static const std::string liststring(Residue *res);

    private:
        Residue *next_res;
        Residue *prev_res;

        std::vector<Bond> prefbonds;
        std::vector<ANGLE> prefangles;

    };


//@{
// Copy a list of residues, returning a newly allocated list of residues
//@}
    Residue *CopyResList(const Residue *oldres);

//@{
// Finds atoms in a residue with the given name, returns the # of atoms
// with that name (should always be zero or one, if the names are unique!)
//@}
    int CountAtomsByName(const char *atom_name, const Residue *res);

//@{
// Finds atoms in a residue with the given name, returns the # of atoms
// with that name (should always be zero or one, if the names are unique!)
//@}
    int DupeAtomNames(const chemlib::Residue *res);


//@{
// Delete the memory associated with a residue list.
//@}
    void FreeResidueList(Residue *reslist);

//@{
// Checks that a vector of atom pointers matches up exactly with those
// in the given residue
//@}
    bool AtomVectMatchesRes(const MIAtomList &ptrs, const Residue *res);

/**
 * Returns the default a pointer to an atom given a pointer to a residue.
 * For amino acids this is the CA atom, otherwise it returns the first atom.
 */
    chemlib::MIAtom *atom_default(const chemlib::Residue *res);


}



#endif // ifndef MI_RESIDUE_H
