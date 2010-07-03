#ifndef MONOMER_H
#define MONOMER_H

#include <algorithm>
#include <map>
#include "MIAtom_fwd.h"
#include "Bond.h"

namespace chemlib
{

    class Monomer
    {

        typedef std::map<const Monomer*, size_t> MonomerRefCountMap;

        /**
         * Stores reference counts to objects of this class. Currently,
         * only increments and decrements in constructor and destructor.
         * Used to check if object has been deleted with isValid method.
         */
        static MonomerRefCountMap refCounts;

    public:

        /**
         * Returns whether the given residue is still valid (not been deleted).
         */
        static bool isValid(const Monomer *res);

        Monomer();
        Monomer(const Monomer&);
        Monomer&operator=(const Monomer&);

        virtual ~Monomer();

        const MIAtomList&atoms() const
        {
            return atoms_;
        }
        int atomCount() const
        {
            return (int)atoms_.size();
        }

        MIAtom *atom(size_t index) const
        {
            return atoms_[index];
        }
        const MIAtom *atomByName(const std::string &name) const;
        int indexOfAtom(const MIAtom*) const;

        void setAtoms(const MIAtomList &a);
        void reserveAtoms(int natoms)
        {
            atoms_.reserve(natoms);
        }
        void addAtom(MIAtom *atom)
        {
            atoms_.push_back(atom);
        }

        // FIXME: clearAtoms does not delete MIAtom objects. In some places it is used
        // it probably should.
        void clearAtoms()
        {
            atoms_.clear();
        }
        bool removeAtom(MIAtom *a)
        {
            bool result = false;
            MIAtom_iter ai
                = std::find(atoms_.begin(), atoms_.end(), a);
            if (ai != atoms_.end())
            {
                atoms_.erase(ai);
                result = true;
            }
            return result;
        }

        const std::string type() const
        {
            return type_;
        }
        void setType(const std::string &t)
        {
            type_ = t;
        }
        const std::string name() const
        {
            return name_;
        }
        void setName(const std::string &n)
        {
            name_ = n;
        }
        unsigned short linkage_type() const
        {
            return linkage_type_;
        }
        void set_linkage_type(unsigned short i)
        {
            linkage_type_ = i;
        }
        void add_linkage_type(unsigned short i)
        {
            linkage_type_ |= i;
        }
        unsigned short chain_id() const
        {
            return chain_id_;
        }
        char getChainId() const
        {
            return (char)(chain_id_ & 255);
        }
        unsigned short getSegmentNumber() const
        {
            return chain_id_ / 256;
        }
        void set_chain_id(unsigned short i)
        {
            chain_id_ = i;
        }
        char secstr() const
        {
            return secstr_;
        }
        void setSecstr(char c)
        {
            secstr_ = c;
        }
        char name1() const
        {
            return name1_;
        }
        void setName1(char c)
        {
            name1_ = c;
        }
        unsigned int seqpos() const
        {
            return seqpos_;
        }
        void setSeqpos(unsigned int i)
        {
            seqpos_ = i;
        }
        unsigned short flags() const
        {
            return flags_;
        }
        void setFlags(unsigned short i)
        {
            flags_ = i;
        }
        unsigned short confomer() const
        {
            return confomer_;
        }
        void setConfomer(unsigned short i)
        {
            confomer_ = i;
        }

        float x() const
        {
            return x_;
        }
        float y() const
        {
            return y_;
        }
        void setPosition(float x, float y)
        {
            x_ = x;
            y_ = y;
        }


    protected:
        MIAtomList atoms_;
        std::string type_;
        std::string name_;
        unsigned short linkage_type_;
        unsigned short chain_id_;
        char secstr_;
        char name1_;
        unsigned int seqpos_;
        unsigned short flags_;
        unsigned short confomer_;

        float x_;
        float y_;

    };

}   //namespace chemlib
#endif //MONOMER_H
