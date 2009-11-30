#ifndef MI_MATHLIB_H
#define MI_MATHLIB_H

#include <stdlib.h>

#include "APOINT.h"
#include "PLINE.h"
#include "LINE.h"
#include "mi_math.h"
#include "LSQMatrix.h"
#include "Matrices.h"
#include "crystmat.h"


#define DEG2RAD (0.017453292519943295769236907684886)       // PI / 180
#define RAD2DEG (57.295779513082320876798154814105)         // 180 / PI

#ifndef PI
#define PI (3.1415926535897932384626433832795028841972)
#endif

// overload ROUND instead of macro for performance reasons.  Unnecessary
// conversion to double (which would happen w/just a single macro version) is slower than necessary
inline int ROUND(float a)
{
    return (a > 0.0f ? (int)(a+0.5f) : (int)(a-0.5f));
}
inline int ROUND(double a)
{
    return (a > 0.0 ? (int)(a+0.5) : (int)(a-0.5));
}
inline int ROUND(int a)
{
    return (a > 0 ? (int)(a+0.5) : (int)(a-0.5));
}

#define SIGN(a)  ((a) >= 0 ? (1) : (-1))
#define fRAND_MAX ((float)RAND_MAX)
#define absolute(x) ((x) >= 0 ? (x) : -(x))

inline float frand2(float range)
{
    return ((2.0F*(float)rand()/fRAND_MAX*range) - range);
}

inline float grand(float range)
{
    return (frand2(range)*frand2(range)*frand2(range));
}

inline float frand()
{
    return (float)rand()/fRAND_MAX;
}

inline float frand(float range)
{
    return (float)rand()/fRAND_MAX*range;
}

inline int irand(int range)
{
    float r = frand();
    int i = (int)(r*(float)range);
    if (i < 0)
    {
        i = 0;
    }
    if (i >= range)
    {
        i = range -1;
    }
    return i;
}


inline bool MIIsNan(double x)
{
    return x != x;
}


#endif // ifndef MI_MATHLIB_H
