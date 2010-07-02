#ifndef ATOM_H
#define ATOM_H

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <util/utillib.h>
#include <math/mathlib.h>
#include <QObject>

#include "MIAtom_fwd.h"

namespace chemlib
{

    class MIAtom : public QObject
    {
        typedef std::map<const MIAtom*, size_t> AtomRefCountMap;

        /**
         * Stores reference counts to objects of this class. Currently,
         * only increments and decrements in constructor and destructor.
         * Used to check if object has been deleted with isValid method.
         */
        static AtomRefCountMap refCounts;

    public:

        /**
         * Returns whether the given atom is still valid (not been deleted).
         */
        static bool isValid(const MIAtom *atom);

        MIAtom();
        ~MIAtom();

    private:
        // Atoms are not copyable or assignable because not all data will properly
        // transfer to another atom.
        // For copying pertinent data, use copy...() member functions.
        MIAtom(const MIAtom&);
        MIAtom&operator=(const MIAtom&);

    public:

        void copyShallow(const MIAtom &atom);

        int CountCyclicBonds() const;
        void ClearHybrid();

        static std::string chiralClassToString(int chiralClass);
        int chiral_class() const
        {
            return chiral_class_;
        }

        void chiral_class(int value)
        {
            chiral_class_ = value;
        }

        int chiral_order() const
        {
            return chiral_order_;
        }

        void chiral_order(int value)
        {
            chiral_order_ = value;
        }

        // set the position of this atom
        void setPosition(float x, float y, float z)
        {
            this->x_ = x;
            this->y_ = y;
            this->z_ = z;
        }

        //@{
        // copy the x,y,z values from atom into the array a[3].
        //@}
        void getPosition(float a[3]) const;
        void getPosition(double a[3]) const;

        void copyPosition(const MIAtom &from);

        // Returns the distance between this atom and atom b.
        double distance(MIAtom *b) const
        {
            return sqrt((double)((x_ - b->x_) * (x_ - b->x_) + (y_ - b->y_) * (y_ - b->y_) + (z_ - b->z_) * (z_ - b->z_)));
        }

        bool isHidden() const
        {
            return color_ < 0;
        }

        const char *name() const
        {
            return name_;
        }
        void setName(const char *name)
        {
            strncpy(name_, name, MAXATOMNAME);
            name_[MAXATOMNAME-1] = 0;
        }
        void fixname();
        unsigned int type() const
        {
            return type_;
        }
        void setType(unsigned int type)
        {
            type_ = type;
        }
        void addType(unsigned int type)
        {
            type_ |= type;
        }
        void removeType(unsigned int type)
        {
            type_ &= ~type;
        }
        float BValue() const
        {
            return BValue_;
        }
        void setBValue(float bvalue)
        {
            BValue_ = bvalue;
        }
        float occ() const
        {
            return occ_;
        }
        void setOcc(float occ)
        {
            occ_ = occ;
        }
        short color() const
        {
            return color_;
        }
        void setColor(short c)
        {
            color_ = c;
        }
        float x() const
        {
            return x_;
        }
        void setX(float x)
        {
            x_ = x;
        }
        float y() const
        {
            return y_;
        }
        void setY(float y)
        {
            y_ = y;
        }
        float z() const
        {
            return z_;
        }
        void setZ(float z)
        {
            z_ = z;
        }

        float dx() const
        {
            return dx_;
        }
        float dy() const
        {
            return dy_;
        }
        float dz() const
        {
            return dz_;
        }
        void addDelta(float dx, float dy, float dz)
        {
            dx_ += dx;
            dy_ += dy;
            dz_ += dz;
        }
        void resetDelta()
        {
            dx_ = 0.0f;
            dy_ = 0.0f;
            dz_ = 0.0f;
        }
        float weight() const
        {
            return weight_;
        }
        void addWeight(float w)
        {
            weight_ += w;
        }
        void resetWeight()
        {
            weight_ = 0.0f;
        }
        unsigned char radius_type() const
        {
            return radius_type_;
        }
        void set_radius_type(unsigned char c)
        {
            radius_type_ = c;
        }

        int symmop() const
        {
            return symmop_;
        }
        void setSymmop(int i)
        {
            symmop_ = i;
        }

        void copyChemicalData(const MIAtom &source);

        char altloc() const
        {
            return altloc_;
        }
        void setAltloc(char c)
        {
            altloc_ = c;
        }
        int atomnumber() const
        {
            return atomnumber_;
        }
        void setAtomnumber(int n)
        {
            atomnumber_ = n;
        }
        short atomicnumber() const
        {
            return atomicnumber_;
        }
        void setAtomicnumber(short n)
        {
            atomicnumber_ = n;
        }
        int isaromatic() const
        {
            return isaromatic_;
        }
        void setIsaromatic(int i)
        {
            isaromatic_ = i;
        }
        int hybrid() const
        {
            return hybrid_;
        }
        void setHybrid(int i)
        {
            hybrid_ = i;
        }
        int geom() const
        {
            return geom_;
        }
        void setGeom(int i)
        {
            geom_ = i;
        }
        int mass() const
        {
            return mass_;
        }
        void setMass(int m)
        {
            mass_ = m;
        }
        int formal_charge() const
        {
            return formal_charge_;
        }
        void set_formal_charge(int i)
        {
            formal_charge_ = i;
        }
        float charge() const
        {
            return charge_;
        }
        void setCharge(float c)
        {
            charge_ = c;
        }
        int hcount() const
        {
            return hcount_;
        }
        void setHcount(int i)
        {
            hcount_ = i;
        }
        int iscyclic() const
        {
            return iscyclic_;
        }
        void setIscyclic(int i)
        {
            iscyclic_ = i;
        }
        int ring_system() const
        {
            return ring_system_;
        }
        void set_ring_system(int i)
        {
            ring_system_ = i;
        }
        int smallest_aromatic_ring() const
        {
            return smallest_aromatic_ring_;
        }
        void set_smallest_aromatic_ring(int i)
        {
            smallest_aromatic_ring_ = i;
        }
        int search_flag() const
        {
            return search_flag_;
        }
        void set_search_flag(int i)
        {
            search_flag_ = i;
        }

        const MIAtomList&nabors() const
        {
            return nabors_;
        }
        void addNabor(MIAtom *atom)
        {
            nabors_.push_back(atom);
        }
        const std::vector<int>&bondnumbers() const
        {
            return bondnumbers_;
        }
        void addBondnumber(int bondNumber)
        {
            bondnumbers_.push_back(bondNumber);
        }

        float U(size_t index) const;
        float U(size_t index, float value);
        bool hasAnisotropicity() const;
        void newAnisotropicity();
        void deleteAnisotropicity();

        void translate(float x, float y, float z)
        {
            x_ += x;
            y_ += y;
            z_ += z;
        }

        void scale(float s)
        {
            x_ *= s;
            y_ *= s;
            z_ *= s;
        }


        inline bool IsHydrogen() const
        {
            return atomicnumber_ == 1;
        }

        inline bool IsCarbon() const
        {
            return atomicnumber_ == 6;
        }

        inline bool IsNitrogen() const
        {
            return atomicnumber_ == 7;
        }

        inline bool IsOxygen() const
        {
            return atomicnumber_ == 8;
        }

        inline bool IsSulfur() const
        {
            return atomicnumber_ == 16;
        }

        // return the default atom radius for an atom of the given type
        static float MIAtomRadiusForType(int i);

        float getRadius() const;

        //@{
        // Returns an internal atom type given a name.
        //@}
        static int MIGetAtomTypeFromName(const char *name);

        static bool MIIsMainChainAtom(const MIAtom*);
        static bool MIIsSideChainAtom(const MIAtom*);
        static bool MIIsMainChainDNAAtom(const MIAtom*);
        static bool MIIsSideChainDNAAtom(const MIAtom*);
        static bool MIIsHydrogen(const MIAtom*);
        static bool MIIsHBondable(const MIAtom *a, const MIAtom *b);

        static std::string liststring(const MIAtom *atom);

        static const MIAtom *GetHeavyNeighbor(const MIAtom *atom);

    private:
        char name_[MAXATOMNAME];
        unsigned int type_;
        float BValue_;
        float occ_;
        short color_;
        float x_;
        float y_;
        float z_;

        float dx_;
        float dy_;
        float dz_;
        float weight_;
        unsigned char radius_type_;

        int symmop_;

        char altloc_; // identifies split residues (alternate location). column 17
        int atomnumber_;
        short atomicnumber_;
        int isaromatic_;
        int hybrid_;
        int geom_;
        int mass_;
        int formal_charge_;
        float charge_;
        int hcount_;
        int iscyclic_;
        int ring_system_; // Index of the ring system, -1 for acyclic atoms
        int smallest_aromatic_ring_; // Size of smallest aromatic ring in which the atom is contained
        int search_flag_; // Flag used in ring detection and aromaticity searches

        MIAtomList nabors_;
        std::vector<int> bondnumbers_;

        int chiral_class_;
        int chiral_order_;

        float *U_;

    };

    struct IsHeavy : std::unary_function<MIAtom*, bool>
    {
        inline bool operator()(const MIAtom *atom) const
        {
            return atom->atomicnumber() != 1;
        }
    };

    inline const MIAtom*MIAtom::GetHeavyNeighbor(const MIAtom *atom)
    {
        MIAtom_const_iter atm;

        atm = std::find_if(atom->nabors().begin(),
                           atom->nabors().end(),
                           IsHeavy());
        return *atm;
    }

    inline void MIAtom::fixname()
    {
        // remove white space in name strings
        size_t j, l;
        while (isspace(name_[0]) )
        {
            l = strlen(name_);
            for (j = 0; j < l; j++)
            {
                name_[j] = name_[j+1];
            }
        }
        while (isspace(name_[strlen(name_)-1]))
        {
            name_[strlen(name_)-1] = '\0';
        }
    }

//#define CATCH_FREE_RESIDUE_TEMPLATE_CODE
#ifdef CATCH_FREE_RESIDUE_TEMPLATE_CODE
    class TORSION;
    template <typename ARG> void free(ARG *a)
    {
        foo((void*)a, -1);
    }

    template <> void free<void>(void *a)
    {
        free(a);
    }

    template <> void free<char>(char *a)
    {
        free(a);
    }

    template <> void free<char*>(char **a)
    {
        free(a);
    }

    template <> void free<float>(float *a)
    {
        free(a);
    }

    template <> void free<double>(double *a)
    {
        free(a);
    }

    template <> void free<int>(int *a)
    {
        free(a);
    }

    template <> void free<MIAtom>(MIAtom *a)
    {
        free(a);
    }

    template <> void free<MIAtom*>(MIAtom **a)
    {
        free(a);
    }

    template <> void free<LINE>(LINE *a)
    {
        free(a);
    }

    template <> void free<TORSION>(TORSION *a)
    {
        free(a);
    }

#endif // ifdef CATCH_FREE_RESIDUE_TEMPLATE_CODE

}   //namespace chemlib
#endif //ATOM_H
