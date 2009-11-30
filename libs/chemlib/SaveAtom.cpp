#include <nongui/nonguilib.h>
#include <vector>
#include "SaveAtom.h"
#include "model.h"

#include "MIAtom.h"


using namespace std;

namespace chemlib
{
int SaveAtom::Restore() const
{
    if (!MIAtom::isValid(from))
    {
        return 0;
    }
    if (from->x() == x && from->y() == y && from->z() == z)
    {
        return 0;
    }
    from->setPosition(x, y, z);
    return 1;
}

bool operator ==(SaveAtom a1, SaveAtom a2)
{
    if (a1.from != a2.from)
    {
        return false;
    }
    if (a1.x > a2.x+0.05 || a1.x < a2.x-0.05)
    {
        return false;
    }
    if (a1.y > a2.y+0.05 || a1.y < a2.y-0.05)
    {
        return false;
    }
    if (a1.z > a2.z+0.05 || a1.z < a2.z-0.05)
    {
        return false;
    }
    // color irrelevant for this comparison
    return true;
}

bool operator ==(vector<SaveAtom> s1, vector<SaveAtom> s2)
{
    if (s1.size() != s2.size())
    {
        return false;
    }
    for (unsigned int i = 0; i < s1.size(); i++)
    {
        if (!(s1[i] == s2[i]))
        {
            return false;
        }
    }
    return true;
}

SaveAtom::SaveAtom()
{
    x = y = z = 0.0f;
    color = 0;
    from = NULL;
}

SaveAtom::SaveAtom(MIAtom *a)
{
    x = a->x();
    y = a->y();
    z = a->z();
    color = a->color();
    from = a;
}

void SaveAtom::RestoreColor(unsigned int mask) const
{
    if (MIAtom::isValid(from))
    {
        from->removeType(mask);
    }
}

}
