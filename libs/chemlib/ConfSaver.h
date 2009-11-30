#ifndef mifit_ConfSaver_h
#define mifit_ConfSaver_h

#include <vector>
#include <list>

#include "GeomSaver.h"
#include "SaveAtom.h"

namespace chemlib
{

    class RESIDUE;
    class MIMoleculeBase;

    typedef  std::list<SaveAtom>::iterator confAtomIter;

    class ConfSaver
    {
        RESIDUE *_res;
        int _natoms;
        ConfSaver(const ConfSaver &rhs);
        ConfSaver&operator=(const ConfSaver &rhs);

        std::list<SaveAtom> atom_store;
        std::vector < std::vector < confAtomIter > > SaveSets;

    public:

        ConfSaver(RESIDUE *res);
        virtual ~ConfSaver();
        void Save();
        void Restore(unsigned int set) const;
        void RestoreLast() const;
        int NumberSets() const;
        const RESIDUE *GetResidue() const;
        void ConvertToGeomSaver(GeomSaver &gs, MIMoleculeBase *model);
    };
}
#endif // ifndef mifit_ConfSaver_h
