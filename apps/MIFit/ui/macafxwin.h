#ifndef MACAFXWIN
#define MACAFXWIN

#include <stdio.h>

//@{
// Imitation of MFC CPoint for portability purposes.
//@}
class CPoint
{
public:
    CPoint()
    {
        x = y = 0;
    }

    CPoint(int ix, int iy)
    {
        x = ix;
        y = iy;
    }

    int x, y;
};

#endif //MACAFXWIN
