#ifndef MI_TORSION_H
#define MI_TORSION_H

#include "model.h"

namespace chemlib
{
    class MIAtom;
    class Residue;

    class TORSDICT
    {
    public:
        TORSDICT()
        {
            memset(this, 0, sizeof(TORSDICT));
        }

        void Clear()
        {
            memset(this, 0, sizeof(TORSDICT));
        }

        bool AddAngle(float angle)
        {
            if (nideal == 3)
            {
                return false;
            }
            else
            {
                ideal[nideal] = angle;
                nideal++;
                return true;
            }
        }

        char restype[MAXNAME];
        char type[11];
        char name[4][MAXATOMNAME]; /* atom names */
        float ideal[3];
        int nideal;
    };

    class TORSION
    {
    public:
        TORSION() : atom1(0), atom2(0), atom3(0), atom4(0), nideal(0), tolerance(0.0f),
            res(0), cyclic(0), aromatic(0), ring_system(0), smallest_ring_size(0),
            smallest_aromatic_ring(0), smallest_saturated_ring(0)
        {
            type[0] = '\0';
            ideal[0] = ideal[1] = ideal[2] = 0.0f;
        }

        void Clear()
        {
            TORSION t;
            *this = t;
        }

        MIAtom *atom1;
        MIAtom *atom2;
        MIAtom *atom3;
        MIAtom *atom4;
        float ideal[3];
        int nideal;
        float tolerance;
        /**
         * the residue the torsion is from
         */
        const Residue *res;
        char type[11];
        /**
         * True if the central bond (atom2->atom3) is cyclic
         */
        int cyclic;
        /**
         * True if the central bond is aromatic
         */
        int aromatic;
        /**
         * Index of the ring system of the central bond, -1 if acyclic
         */
        int ring_system;
        /**
         * Size of smallest ring in which the central bnd is contained
         */
        int smallest_ring_size;
        /**
         * Size of smallest aromatic ring which has the central bnd
         */
        int smallest_aromatic_ring;
        /**
         * Size of smallest saturated (all sp3) ring which has central bnd
         */
        int smallest_saturated_ring;

        MIAtom *getAtom1() const
        {
            return atom1;
        }

        void setAtom1(MIAtom *atom)
        {
            atom1 = atom;
        }

        MIAtom *getAtom2() const
        {
            return atom2;
        }

        void setAtom2(MIAtom *atom)
        {
            atom2 = atom;
        }

        MIAtom *getAtom3() const
        {
            return atom3;
        }

        void setAtom3(MIAtom *atom)
        {
            atom3 = atom;
        }

        MIAtom *getAtom4() const
        {
            return atom4;
        }

        void setAtom4(MIAtom *atom)
        {
            atom4 = atom;
        }

    };


}
#endif // ifndef MI_TORSION_H
