#ifndef mifit_SaveItem_h
#define mifit_SaveItem_h

#include <string>
#include <vector>

#include "SaveAtom.h"

namespace chemlib
{

    class MIMoleculeBase;

/**
 * Save a state for undo function.
 */
    class SaveItem
    {
    public:
        SaveItem();
        SaveItem(MIMoleculeBase *node, std::string title);
        MIMoleculeBase *SaveMolecule;
        std::string Title;
        std::vector<SaveAtom> SaveSet;
    };

}
#endif // ifndef mifit_SaveItem_h
