#ifndef MI_CRECT_H
#define MI_CRECT_H

#ifndef _WIN32
typedef struct tagRECT
{
    int left, top, right, bottom;
} RECT;
#else
#include <windows.h>
#endif

//@{
// Imitation of MFC CRect for portability purposes.
//@}
class CRect : public tagRECT
{
public:
    CRect()
    {
        top = left = bottom = right = 0;
    }

    CRect(int l, int t, int r, int b)
    {
        left = l;
        top = t;
        right = r;
        bottom = b;
    }

    int Width()
    {
        return right-left;
    }

    int Height()
    {
        return bottom-top;
    }

    bool Within(long x, long y)
    {
        return (x >= left && x <= right && y >= top && y <= bottom);
    }

};

#endif // ifndef MI_CRECT_H
