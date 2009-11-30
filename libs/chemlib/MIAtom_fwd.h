#ifndef chemlib_MIAtom_fwd_h
#define chemlib_MIAtom_fwd_h

#include <vector>

namespace chemlib
{

    class MIAtom;

    typedef std::vector<MIAtom*> MIAtomList;
    typedef MIAtomList::const_iterator MIAtom_const_iter;
    typedef MIAtomList::iterator MIAtom_iter;

    const unsigned int MAXATOMNAME = 6;

}

#endif // ifndef chemlib_MIAtom_fwd_h
