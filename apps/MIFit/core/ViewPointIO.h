#ifndef VIEWPOINTIO_H
#define VIEWPOINTIO_H

#include <cstdio>

class CArchive;
class ViewPoint;

class ViewPointIO
{
    ViewPointIO();
public:

    static void save(ViewPoint &vp, CArchive &ar);
    static void load(ViewPoint &vp, FILE *fp);

};

#endif // VIEWPOINTIO_H
