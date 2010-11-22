#include <cmath>
#include <cstdio>

#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>
#include <chemlib/MIMoleculeBase.h>

#include "ViewPoint.h"
#include "Cfiles.h"

using namespace chemlib;

namespace
{
    const float CRS = 100.0f;
}

ViewPoint::ViewPoint()
    : width_(200), height_(100),
      scale(20),
      frontClip_(5.0), backClip_(-5.0)

{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
            {
                viewmat[i][j] = 1.0;
            }
            else
            {
                viewmat[i][j] = 0;
            }
        }
    }
    center[0] = center[1] = center[2] = 0;
    zangle = (long)(1000.0*sin(perspective*M_PI/180.0));
}

void ViewPoint::UnDo()
{
    if (UndoList.size() > 0)
    {
        for (int i = 0; i < 3; i++)
        {
            center[i] = UndoList.back().center[i];
            for (int j = 0; j < 3; j++)
            {
                viewmat[i][j] = UndoList.back().viewmat[i][j];
            }
        }
        scale = UndoList.back().scale;
        frontClip_ = UndoList.back().frontClip;
        backClip_ = UndoList.back().backClip;
        //undoable = false;
        UndoList.pop_back();
    }
}

void ViewPoint::Do()
{
    ViewSave last;
    for (int i = 0; i < 3; i++)
    {
        last.center[i] = center[i];
        for (int j = 0; j < 3; j++)
        {
            last.viewmat[i][j] = viewmat[i][j];
        }
    }
    last.scale = scale;
    last.frontClip = frontClip_;
    last.backClip = backClip_;
    UndoList.push_back(last);
}

void ViewPoint::setmatrix(float mat[3][3])
{
    viewmat[0][0] = mat[0][0];
    viewmat[0][1] = mat[0][1];
    viewmat[0][2] = mat[0][2];
    viewmat[1][0] = mat[1][0];
    viewmat[1][1] = mat[1][1];
    viewmat[1][2] = mat[1][2];
    viewmat[2][0] = mat[2][0];
    viewmat[2][1] = mat[2][1];
    viewmat[2][2] = mat[2][2];
    orthomatrix(viewmat, viewmat);
}

void ViewPoint::copymatrix(float mat[3][3])
{
    //float msf = MSF;
    mat[0][0] = viewmat[0][0];
    mat[0][1] = viewmat[0][1];
    mat[0][2] = viewmat[0][2];
    mat[1][0] = viewmat[1][0];
    mat[1][1] = viewmat[1][1];
    mat[1][2] = viewmat[1][2];
    mat[2][0] = viewmat[2][0];
    mat[2][1] = viewmat[2][1];
    mat[2][2] = viewmat[2][2];
}

void ViewPoint::setSize(int width, int height)
{
    width_ = width;
    height_ = height;
}

void ViewPoint::rotate(float rx, float ry, float rz)
{
    incmatrix(rx, ry, rz, viewmat, viewmat);
    orthomatrix(viewmat, viewmat);
}

void ViewPoint::translate(int x, int y, int z)
{
    center[0] += x;
    center[1] += y;
    center[2] += z;
}

void ViewPoint::moveto(int x, int y, int z)
{
    center[0] = x;
    center[1] = y;
    center[2] = z;
}

void ViewPoint::moveto(float x, float y, float z)
{
    center[0] = ROUND(x*CRS);
    center[1] = ROUND(y*CRS);
    center[2] = ROUND(z*CRS);
}

void ViewPoint::setperspective(float p)
{
    perspective = p;
    if (p < 0.0)
    {
        perspective = 0.0;
    }
    zangle = (long)(1000.0*sin(perspective*M_PI/180.0));
}

QVector3D ViewPoint::transform(const QVector3D& point)
{
    return QVector3D(viewmat[0][0]*point.x() + viewmat[0][1]*point.y()
                     + viewmat[0][2]*point.z(),
                     viewmat[1][0]*point.x() + viewmat[1][1]*point.y()
                     + viewmat[1][2]*point.z(),
                     viewmat[2][0]*point.x() + viewmat[2][1]*point.y()
                     + viewmat[2][2]*point.z());
}


void ViewPoint::transform(float ix, float iy, float iz, int &xt, int &yt, int &zt)
{
    ix *= CRS;
    iy *= CRS;
    iz *= CRS;

    static const float cscale = 100.0;

    float x = ix - center[0];
    float y = iy - center[1];
    float z = iz - center[2];

    float rz = -(viewmat[2][0]*x+viewmat[2][1]*y+viewmat[2][2]*z);

    // z is reversed to keep system righthanded
    // if x is across from left to right and y is down
    // then z points into the screen
    float zp = rz*zangle;

    float rx = (viewmat[0][0]*x+viewmat[0][1]*y+viewmat[0][2]*z);
    float t = rx*zp/cscale/cscale/10.0;

    rx += t;
    rx = rx * scale;
    rx = (rx > 0) ? (rx+50.0) : (rx-50.0);
    rx /= cscale;
    xt = (int)rx;
    xt += width_/2;

    float ry = (viewmat[1][0]*x+viewmat[1][1]*y+viewmat[1][2]*z);
    t = ry*zp/cscale/cscale/10.0;
    ry += t;
    ry = ry * scale;
    ry = (ry > 0) ? (ry+50.0) : (ry-50.0);
    ry /= cscale;
    yt = (int)ry;
    yt += height_/2;

    rz *= scale;
    rz = (rz > 0) ? (rz+50.0) : (rz-50.0);
    rz /= cscale;
    zt = (int)rz;
}

void ViewPoint::Invert(int sx, int sy, int sz, float &x, float &y, float &z)
{
    float xdir, ydir, zdir, zp, t;
    float rx, ry, rz;
    static const float cscale = 100.0F;
    uinv(viewmat, umat);
    sx -= width_/2;
    sy -= height_/2;
    rx = sx*cscale;
    ry = sy*cscale;
    rz = sz*cscale;
    zp = rz*zangle;
    rx = rx/scale;
    ry = ry/scale;
    rz = -rz/scale;
    t = rx*zp/cscale/cscale/10.0F;
    rx -= t;
    t = ry*zp/cscale/cscale/10.0F;
    ry -= t;
    xdir = rx* umat[0][0] + ry*umat[0][1] + rz*umat[0][2];
    ydir = rx* umat[1][0] + ry*umat[1][1] + rz*umat[1][2];
    zdir = rx* umat[2][0] + ry*umat[2][1] + rz*umat[2][2];
    xdir += center[0];
    ydir += center[1];
    zdir += center[2];
    x = xdir/CRS;
    y = ydir/CRS;
    z = zdir/CRS;
}

void ViewPoint::Center(MIMoleculeBase *mol)
{
    if (!mol)
        return;

    float xsum = 0, ysum = 0, zsum = 0;
    float n = 0;
    ResidueListIterator res = mol->residuesBegin();
    for (; res != mol->residuesEnd(); ++res)
    {
        for (int i = 0; i < res->atomCount(); ++i)
        {
            xsum += res->atom(i)->x();
            ysum += res->atom(i)->y();
            zsum += res->atom(i)->z();
            n++;
        }
    }
    if (n == 0)
        return;

    center[0] = (long)(xsum/n * CRS);
    center[1] = (long)(ysum/n * CRS);
    center[2] = (long)(zsum/n * CRS);
}

void ViewPoint::zoom(float ds)
{
    scale *= ds;
}

void ViewPoint::getdirection(int xdrag, int ydrag, int &xdir, int &ydir, int &zdir)
{
    //Scroll view
    //determine direction relative to molecule axes
    //by inverting matrix
    //int xdir, ydir, zdir;

    uinv(viewmat, umat);
    xdir = ROUND((xdrag* umat[0][0]
                  +ydrag*umat[0][1])/scale*CRS);
    ydir = ROUND((xdrag* umat[1][0]
                  +ydrag*umat[1][1])/scale*CRS);
    zdir = ROUND((xdrag* umat[2][0]
                  +ydrag*umat[2][1])/scale*CRS);
}

void ViewPoint::scroll(int xdrag, int ydrag, int zdrag)
{
    //Scroll view
    //determine direction relative to molecule axes
    //by inverting matrix
    int xdir, ydir, zdir;

    uinv(viewmat, umat);
    xdir = ROUND((xdrag* umat[0][0]
                  +ydrag*umat[0][1]
                  +zdrag*umat[0][2])/scale*CRS);
    ydir = ROUND((xdrag* umat[1][0]
                  +ydrag*umat[1][1]
                  +zdrag*umat[1][2])/scale*CRS);
    zdir = ROUND((xdrag* umat[2][0]
                  +ydrag*umat[2][1]
                  +zdrag* umat[2][2])/scale*CRS);
    translate(-xdir, -ydir, -zdir);
}

void ViewPoint::slab(int s)
{
    if (s > 0 || abs(frontClip_-backClip_) > abs(s))
    {
        frontClip_ += s/2.0;
        backClip_ -= s/2.0;
    }
    else
    {
        backClip_ = (frontClip_ + backClip_)/2;
        frontClip_ = backClip_ + 1;
    }
    if (frontClip_ < backClip_)
    {
        frontClip_ = backClip_ + 1;
    }
}

void ViewPoint::slab(float f, float b)
{
    if (f < b)
    {
        std::swap(f, b);
    }
    frontClip_ = f;
    backClip_ = b;
}

void ViewPoint::PutVertical(float x2, float y2, float z2)
{
    // puts the vector x1 -> x2 vertical on the screen
    // first point is the center of the screen
    // second should be vertical (minimum sy)
    int sx2, sy2, sz2;
    int minsy;
    float bestx = 0, besty = 0, bestz = 0;
    float xrot, yrot, zrot;
    Do();  // this saves the starting viewpoint
    transform(x2, y2, z2, sx2, sy2, sz2);
    minsy = sy2;
    size_t i;
    for (i = 0; i < 300; i++)
    {
        xrot = frand2(180.0F);
        yrot = frand2(180.0F);
        zrot = frand2(180.0F);
        rotate(xrot, yrot, zrot);
        transform(x2, y2, z2, sx2, sy2, sz2);
        if (sy2 < minsy)
        {
            minsy = sy2;
            bestx = xrot;
            besty = yrot;
            bestz = zrot;
        }
        UnDo();
        Do();
    }
    UnDo();
    rotate(bestx, besty, bestz);
    Do();
    transform(x2, y2, z2, sx2, sy2, sz2);
    int mindx = abs(sx2 - width_/2)+ abs(sz2);
    int dx;
    bestx = besty = bestz = 0;
    float angle = 10.0F;
    for (i = 0; i < 50; i++)
    {
        xrot = frand2(angle)+bestx;
        yrot = frand2(angle)+besty;
        zrot = frand2(angle)+bestz;
        rotate(xrot, yrot, zrot);
        transform(x2, y2, z2, sx2, sy2, sz2);
        dx = abs(sx2 - width_/2) + abs(sz2);
        if (dx < mindx)
        {
            mindx = dx;
            bestx = xrot;
            besty = yrot;
            bestz = zrot;
            angle = angle * 0.80F;
        }
        UnDo();
        Do();
        if (dx == 0)
        {
            break;
        }
    }
    UnDo();
    rotate(bestx, besty, bestz);
}

void ViewPoint::PutOnLeft(float x2, float y2, float z2)
{
    // puts the point on the leftmost (min x) with only y axis rotation
    int sx2, sy2, sz2;
    int minsx;
    float bestx = 0, besty = 0, bestz = 0;
    float xrot, yrot, zrot;
    Do();  // this saves the starting viewpoint
    transform(x2, y2, z2, sx2, sy2, sz2);
    minsx = sx2;
    for (size_t i = 0; i < 100; i++)
    {
        xrot = 0;
        yrot = frand2(180.0F);
        zrot = 0;
        rotate(xrot, yrot, zrot);
        transform(x2, y2, z2, sx2, sy2, sz2);
        if (sx2 < minsx)
        {
            minsx = sx2;
            bestx = xrot;
            besty = yrot;
            bestz = zrot;
        }
        UnDo();
        Do();
    }
    UnDo();
    rotate(bestx, besty, bestz);
}

bool ViewPoint::CenterAtResidue(const Residue *res)
{
    if (!Monomer::isValid(res))
    {
        return false;
    }
    MIAtom *CA = atom_from_name("CA", *res);
    Do();
    if (CA)
    {
        MIAtom *a2 = atom_from_name("CB", *res);
        if (!a2)
        {
            a2 = atom_from_name("HA1", *res);
        }
        MIAtom *N = atom_from_name("N", *res);
        moveto(CA->x(), CA->y(), CA->z());
        if (a2)
        {
            PutVertical(a2->x(), a2->y(), a2->z());
            if (N)
            {
                PutOnLeft(N->x(), N->y(), N->z());
                PutVertical(a2->x(), a2->y(), a2->z());
            }
        }
    }
    else
    {
        float x = 0, y = 0, z = 0;
        for (int i = 0; i < res->atomCount(); i++)
        {
            x += res->atom(i)->x();
            y += res->atom(i)->y();
            z += res->atom(i)->z();
        }
        x /= res->atomCount();
        y /= res->atomCount();
        z /= res->atomCount();
        moveto(x, y, z);
    }
    return true;
}

void ViewPoint::setscale(qreal s)
{
    scale = s;
    if (scale < 1.1)
    {
        scale = 1.1;
    }
}

qreal ViewPoint::getscale()
{
    return scale;
}

qreal ViewPoint::frontClip()
{
    return frontClip_;
}

qreal ViewPoint::backClip()
{
    return backClip_;
}

void ViewPoint::setFrontClip(qreal f)
{
    frontClip_ = f;
}

void ViewPoint::setBackClip(qreal b)
{
    backClip_ = b;
}



float ViewPoint::getperspective()
{
    return perspective;
}

float ViewPoint::getcenter(int i)
{
    if (i >= 0 && i < 3)
    {
        return center[i]/CRS;
    }
    else
    {
        return 0.0;
    }
}

int ViewPoint::width()
{
    return width_;
}

int ViewPoint::height()
{
    return height_;
}


float ViewPoint::getwidth()
{
    return width_ / scale;
}

float ViewPoint::getheight()
{
    return height_ / scale;
}

bool ViewPoint::UnDoable()
{
    return UndoList.size() > 0;
}
