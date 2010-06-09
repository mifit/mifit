#ifndef mifit_ConfSaver_h
#define mifit_ConfSaver_h

#include <vector>
#include <list>

#include "GeomSaver.h"
#include "SaveAtom.h"

namespace chemlib
{

    class Residue;
    class MIMoleculeBase;

    typedef  std::list<SaveAtom>::iterator confAtomIter;

    class ConfSaver
    {
        Residue *_res;
        int _natoms;
        ConfSaver(const ConfSaver &rhs);
        ConfSaver&operator=(const ConfSaver &rhs);

        std::list<SaveAtom> atom_store;
        std::vector < std::vector < confAtomIter > > SaveSets;

    public:

        ConfSaver(Residue *res);
        virtual ~ConfSaver();
        void Save();
        void Restore(unsigned int set) const;
        void RestoreLast() const;
        int NumberSets() const;
        const Residue *GetResidue() const;
        void ConvertToGeomSaver(GeomSaver &gs, MIMoleculeBase *model);
    };
}
#endif // ifndef mifit_ConfSaver_h
