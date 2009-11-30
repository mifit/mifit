#ifndef MI_PLANE_H
#define MI_PLANE_H

#include <cstring>
#include "model.h"

namespace chemlib
{
    class MIAtom;
    class RESIDUE;

    class PLANE
    {
    public:
        PLANE() : atoms(NULL), natoms(0), tolerance(0.005f), res(NULL)
        {
        }

        void Clear()
        {
            PLANE p;
            *this = p;
        }                              // note: doesn't free old mem, if any

        MIAtom **atoms;
        int natoms;
        float tolerance;
        /**
         * normal vector to plane.
         */
        float vm[3];
        /**
         * distance from origin of normal vector.
         */
        float d;
        /**
         * the residue the bond is from.
         */
        const RESIDUE *res;

    };
#define MAXPLANE 20
//@{
// A plane dictionary entry.
//@}
    class PLANEDICT
    {
    public:
        PLANEDICT()
        {
            memset(this, 0, sizeof(PLANEDICT));
        }

        void Clear()
        {
            memset(this, 0, sizeof(PLANEDICT));
        }

        bool AddAtom(const char *atom_name)
        {
            if (natoms + 1 <= MAXPLANE)
            {
                strncpy(name[natoms], atom_name, MAXNAME);
                natoms++;
                return true;
            }
            else
            {
                return false;
            }
        }

        char restype[MAXNAME];
        int natoms;
        char name[MAXPLANE][MAXNAME]; /* if more than 20 use overlapping planes*/
    };

}


#endif // ifndef MI_PLANE_H
