#include <vector>
#include <sstream>
#include <fstream>

#include "Residue.h"
#include "mol_util.h"
#include "mol_util_private.h"
#include "atom_util.h"

using namespace std;
namespace chemlib
{

//#define MEMORY_CORRUPTION_DEBUG
#if defined(MEMORY_CORRUPTION_DEBUG)
static unsigned int magic_delete = UINT_MAX;
static unsigned int DELETION_COUNT = 0;
static std::map<Residue*, unsigned int> DELETED_RESIDUES;
#endif

Monomer::MonomerRefCountMap Monomer::refCounts;

bool Monomer::isValid(const Monomer *res)
{
    bool result = false;
    if (res != NULL && refCounts.find(res) != refCounts.end())
    {
        result = refCounts[res] > 0;
    }
#ifdef MEMORY_CORRUPTION_DEBUG
    if (!result && res)
    {
        printf("Invalid residue %p queried, was deletion #%d\n", res, DELETED_RESIDUES[res]);
    }
#endif

    return result;
}

Monomer::Monomer()
    : linkage_type_(0),
      chain_id_(' '),
      secstr_('U'),
      name1_(0),
      seqpos_(0),
      flags_(0),
      confomer_(0),
      x_(0.0f),
      y_(0.0f)
{
    ++refCounts[this];
    type_ = "";
    name_ = "";
}

Monomer::Monomer(const Monomer &rhs)
{
    ++refCounts[this];
    *this = rhs;

    // must do a deep copy of the atoms
    atoms_.clear();
    for (unsigned int i = 0; i < rhs.atoms_.size(); ++i)
    {
        MIAtom *newAtom = new MIAtom;
        newAtom->copyShallow(*rhs.atoms_[i]);
        atoms_.push_back(newAtom);
    }
}

Monomer&Monomer::operator=(const Monomer &rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    // delete old atoms and copy new ones
    for (unsigned int i = 0; i < atoms_.size(); ++i)
    {
        delete atoms_[i];
    }
    atoms_.clear();
    for (unsigned int i = 0; i < rhs.atoms_.size(); ++i)
    {
        MIAtom *newAtom = new MIAtom;
        newAtom->copyShallow(*rhs.atoms_[i]);
        atoms_.push_back(newAtom);
    }

    type_ = rhs.type_;
    name_ = rhs.name_;
    linkage_type_ = rhs.linkage_type_;
    chain_id_ = rhs.chain_id_;
    secstr_ = rhs.secstr_;
    name1_ = rhs.name1_;
    seqpos_ = rhs.seqpos_;
    flags_ = rhs.flags_;
    confomer_ = rhs.confomer_;
    x_ = rhs.x_;
    y_ = rhs.y_;

    return *this;
}


Monomer::~Monomer()
{
    clear_residue_from_atom_cache(this);

    for (unsigned int i = 0; i < atoms_.size(); ++i)
    {
        delete atoms_[i];
    }
    --refCounts[this];
    if (refCounts[this] == 0)
    {
        refCounts.erase(this);
    }
#if defined(MEMORY_CORRUPTION_DEBUG)
    if (DELETED_RESIDUES.size()==magic_delete)
        printf("Set magic_delete in debugger and set breakpoint here.\n");
    DELETED_RESIDUES[this] = DELETION_COUNT++;
#endif
}

void Monomer::setAtoms(const MIAtomList &a)
{
    MIAtomList::iterator iter = atoms_.begin();
    for (; iter != atoms_.end(); ++iter)
        delete *iter;
    atoms_.clear();
    atoms_ = a;
}

const MIAtom*Monomer::atomByName(const std::string &name) const
{
    MIAtom_const_iter i;

    i = find_if(atoms_.begin(), atoms_.end(), std::bind2nd(MatchesAtomName(), name));

    if (i != atoms_.end())
    {
        return *i;                              //Converts iterator to ptr
    }
    else
    {
        return 0;
    }
}

int Monomer::indexOfAtom(const MIAtom *patm) const
{
    for (unsigned int i = 0; i < atoms_.size(); ++i)
    {
        if (atoms_[i] == patm)
        {
            return i;
        }
    }
    return -1;
}

}
