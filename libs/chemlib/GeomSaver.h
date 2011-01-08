#ifndef mifit_model_GeomSaver_h
#define mifit_model_GeomSaver_h

#include "SaveItem.h"
#include <ui/Logger.h>
#include "MIAtom_fwd.h"

namespace chemlib
{

    class Residue;
    class MIMoleculeBase;

/**
 * Saves geometry sets for undo'ing.  Top level of the undo system.
 */
    class GeomSaver
    {

        std::vector<SaveItem> SaveSets;
        unsigned int unique();

    public:

        GeomSaver();
        virtual ~GeomSaver();

        MIMoleculeBase *Model(unsigned int token);
        unsigned int Save(Residue *res, int nres, MIMoleculeBase *model);
        unsigned int Save(const MIAtomList&, MIMoleculeBase *model);
        int Restore(unsigned int set) const;
        bool RestoreColor(unsigned int token, unsigned int mask);
        unsigned int RestoreLast(MIMoleculeBase *model);
        int NumberSets() const
        {
            return SaveSets.size();
        }

        void Purge(MIMoleculeBase*);
        void Purge(MIAtom*);
        void Clear();
    };

}
#endif // ifndef mifit_model_GeomSaver_h
