#include <cmath>
#include <cstdio>

#include <math/mathlib.h>
#include <chemlib/chemlib.h>
#include <chemlib/Monomer.h>
#include <chemlib/MIMoleculeBase.h>

#include "ViewPoint.h"
#include "Cfiles.h"
//#include "Viewpnt.h"
//#include "macafxwin.h"

using namespace chemlib;

const int ViewPoint::STICKS = 0;
const int ViewPoint::BALLANDSTICK = 1;
const int ViewPoint::CPK = 2;
const int ViewPoint::BALLANDCYLINDER = 3;

ViewPoint::ViewPoint()
    : lineThickness(1)
{
    cscale = 100;
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
    changed = 0;
    zangle = (long)(1000.0*sin(perspective*M_PI/180.0));
    scale = 200;
    frontclip = 500;
    backclip = -500;
    width_ = 200;
    xcenter = width_/2;
    height_ = 100;
    ycenter = height_/2;
    topview = 0;
    ballmode = 0;
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
        frontclip = UndoList.back().frontclip;
        backclip = UndoList.back().backclip;
        //undoable = false;
        UndoList.pop_back();
        clampScale();
        changed = 1;
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
    last.frontclip = frontclip;
    last.backclip = backclip;
    UndoList.push_back(last);
    changed = 1;
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
    orthomat();
    changed = 1;
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
    xcenter = width_/2;
    ycenter = height_/2;
    changed = 1;
}

void ViewPoint::clampScale()
{
    if (scale < 1)
    {
        scale = 1;
        changed = 1;
    }
}

void ViewPoint::rotate(float rx, float ry, float rz)
{
    incmatrix(rx, ry, rz, viewmat, viewmat);
    orthomatrix(viewmat, viewmat);
    changed = 1;
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
    changed = 1;
}

float ViewPoint::x(float x, float y, float z)
{
    return (viewmat[0][0]*x+viewmat[0][1]*y+viewmat[0][2]*z);
}

float ViewPoint::y(float x, float y, float z)
{
    return (viewmat[1][0]*x+viewmat[1][1]*y+viewmat[1][2]*z);
}

float ViewPoint::z(float x, float y, float z)
{
    return (-(viewmat[2][0]*x+viewmat[2][1]*y+viewmat[2][2]*z));
}

void ViewPoint::transform(float ix, float iy, float iz, int &xt, int &yt, int &zt)
{
    ix *= (float)CRS;
    iy *= (float)CRS;
    iz *= (float)CRS;

    static float rx, ry, rz, x, y, z, zp, t;
    static const float cscale = 100.0F;

    x = ix - center[0];
    y = iy - center[1];
    z = iz - center[2];

    rz = -(viewmat[2][0]*x+viewmat[2][1]*y+viewmat[2][2]*z);

    // z is reversed to keep system righthanded
    // if x is across from left to right and y is down
    // then z points into the screen
    zp = rz*zangle;

    rx = (viewmat[0][0]*x+viewmat[0][1]*y+viewmat[0][2]*z);
    t = rx*zp/cscale/cscale/10.0F;

    rx += t;
    rx = rx * scale/10.0F;
    rx = (rx > 0) ? (rx+50.0F) : (rx-50.0F);
    rx /= cscale;
    xt = (int)rx;
    xt += xcenter;

    ry = (viewmat[1][0]*x+viewmat[1][1]*y+viewmat[1][2]*z);
    t = ry*zp/cscale/cscale/10L;
    ry += t;
    ry = ry * scale/10L;
    ry = (ry > 0) ? (ry+50L) : (ry-50L);
    ry /= cscale;
    yt = (int)ry;
    yt += ycenter;

    rz *= scale/10L;
    rz = (rz > 0) ? (rz+50L) : (rz-50L);
    rz /= cscale;
    zt = (int)rz;
}

void ViewPoint::Invert(int sx, int sy, int sz, float &x, float &y, float &z)
{
    float xdir, ydir, zdir, zp, t;
    float rx, ry, rz;
    static float cscale = 100.0F;
    uinv(viewmat, umat);
    sx -= xcenter;
    sy -= ycenter;
    rx = (float)sx*cscale;
    ry = (float)sy*cscale;
    rz = (float)sz*cscale;
    zp = rz*zangle;
    rx = rx/scale*10.0F;
    ry = ry/scale*10.0F;
    rz = -rz/scale*10.0F;
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
    x = xdir/(float)CRS;
    y = ydir/(float)CRS;
    z = zdir/(float)CRS;
}

void ViewPoint::Center(MIMoleculeBase *mol)
{
    float xsum = 0, ysum = 0, zsum = 0;
    float n = 0;
    int i;

    if (!mol)
    {
        return;
    }
    ResidueListIterator res = mol->residuesBegin();
    for (; res != mol->residuesEnd(); ++res)
    {
        for (i = 0; i < res->atomCount(); i++)
        {
            xsum += res->atom(i)->x()*CRS;
            ysum += res->atom(i)->y()*CRS;
            zsum += res->atom(i)->z()*CRS;
            n++;
        }
    }
    if (n == 0)
    {
        return;
    }
    center[0] = (long)(xsum/n);
    center[1] = (long)(ysum/n);
    center[2] = (long)(zsum/n);
    changed = 1;
}

void ViewPoint::zoom(float ds)
{
    long s = scale;
    scale = (long)((float)scale*ds);
    // make sure that integer roundoff did not cause scale not to change
    if (ds < 1.0 && s == scale)
    {
        scale--;
    }
    if (ds > 1.0 && s == scale)
    {
        scale++;
    }
    clampScale();
    changed = 1;
}

void ViewPoint::getdirection(int xdrag, int ydrag, int &xdir, int &ydir, int &zdir)
{
    //Scroll view
    //determine direction relative to molecule axes
    //by inverting matrix
    //int xdir, ydir, zdir;

    uinv(viewmat, umat);
    xdir = ROUND((xdrag* umat[0][0]
                  +ydrag*umat[0][1])/(float)scale*MSF);
    ydir = ROUND((xdrag* umat[1][0]
                  +ydrag*umat[1][1])/(float)scale*MSF);
    zdir = ROUND((xdrag* umat[2][0]
                  +ydrag*umat[2][1])/(float)scale*MSF);
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
                  +zdrag*umat[0][2])/(float)scale*MSF);
    ydir = ROUND((xdrag* umat[1][0]
                  +ydrag*umat[1][1]
                  +zdrag*umat[1][2])/(float)scale*MSF);
    zdir = ROUND((xdrag* umat[2][0]
                  +ydrag*umat[2][1]
                  +zdrag* umat[2][2])/(float)scale*MSF);
    translate(-xdir, -ydir, -zdir);
    changed = 1;
}

void ViewPoint::slab(int s)
{
    if (s > 0 || abs(frontclip-backclip) > abs(s))
    {
        frontclip += s/2;
        backclip -= s/2;
    }
    else
    {
        backclip = (frontclip + backclip)/2;
        frontclip = backclip +1;
    }
    if (frontclip < backclip)
    {
        frontclip = backclip +1;
    }
    changed = 1;
}

void ViewPoint::slab(float f, float b)
{
    if (f < b)
    {
        float t;
        t = f;
        f = b;
        b = t;
    }
    f *= CRS;
    b *= CRS;
    frontclip = (int)f;
    backclip = (int)b;
    changed = 1;
}

void ViewPoint::Save(CArchive &ar)
{
    std::string s;
    if (ar.IsStoring())
    {
        s = format("translation %0.2f %0.2f %0.2f\n", getcenter(0),
                   getcenter(1), getcenter(2));
        ar.Write(s.c_str(), s.size());
        s = format("rotation %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f %0.4f\n",
                   (float)viewmat[0][0], (float)viewmat[0][1], (float)viewmat[0][2],
                   (float)viewmat[1][0], (float)viewmat[1][1], (float)viewmat[1][2],
                   (float)viewmat[2][0], (float)viewmat[2][1], (float)viewmat[2][2]);
        ar.Write(s.c_str(), s.size());
        s = format("zoom %0.3f\nperspective %0.3f\n", (float)scale/10.F, perspective);
        ar.Write(s.c_str(), s.size());
        s = format("frontclip %0.2f\nbackclip %0.2f\n", (float)frontclip/CRS, (float)backclip/CRS);
        ar.Write(s.c_str(), s.size());
        s = format("transform\n");
        ar.Write(s.c_str(), s.size());
    }
    else
    {
    }
}

void ViewPoint::Load(FILE *fp)
{
    float v1, v2, v3, v4, v5, v6, v7, v8, v9;
    std::string buf;

    std::auto_ptr<io> ioObj(io::defaultIo());
    io &file = *ioObj;
    file.attach(fp);
    while (file.readLine(buf) != 0)
    {
        std::transform(buf.begin(), buf.end(), buf.begin(), tolower);
        if (!strncmp(buf.c_str(), "rotation", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f%f%f%f%f%f%f%f%f", &v1, &v2, &v3, &v4, &v5, &v6, &v7, &v8, &v9) == 9)
            {
                viewmat[0][0] = v1;
                viewmat[0][1] = v2;
                viewmat[0][2] = v3;
                viewmat[1][0] = v4;
                viewmat[1][1] = v5;
                viewmat[1][2] = v6;
                viewmat[2][0] = v7;
                viewmat[2][1] = v8;
                viewmat[2][2] = v9;
                orthomatrix(viewmat, viewmat);
                changed = 1;
            }
        }
        else if (!strncmp(buf.c_str(), "zoom", 4))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                v1 *= 10.0F;
                scale = ROUND(v1);
                changed = 1;
            }
        }
        else if (!strncmp(buf.c_str(), "frontclip", 9))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                frontclip = ROUND(v1*CRS);
                changed = 1;
            }
        }
        else if (!strncmp(buf.c_str(), "backclip", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                backclip = ROUND(v1*CRS);
                changed = 1;
            }
        }
        else if (!strncmp(buf.c_str(), "perspect", 8))
        {
            if (sscanf(buf.c_str(), "%*s%f", &v1) == 1)
            {
                setperspective(v1);
                changed = 1;
            }
        }
        else if (!strncmp(buf.c_str(), "translation", 11))
        {
            if (sscanf(buf.c_str(), "%*s%f%f%f", &v1, &v2, &v3) == 3)
            {
                center[0] = ROUND(v1 * CRS);
                center[1] = ROUND(v2 * CRS);
                center[2] = ROUND(v3 * CRS);
                changed = 1;
            }
        }
    }
    if (changed)
    {
        clampScale();
    }
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
    int mindx = abs(sx2 - xcenter)+ abs(sz2);
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
        dx = abs(sx2 - xcenter) + abs(sz2);
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
        x /= (float)res->atomCount();
        y /= (float)res->atomCount();
        z /= (float)res->atomCount();
        moveto(x, y, z);
    }
    changed = 1;
    return true;
}

int ViewPoint::colordepth(int c)
{
    if (c < 0)
    {
        return 0;
    }
    else if (depthcue_color)
    {
        return (c);
    }
    else
    {
        return (0);
    }
}

void ViewPoint::setDepthCuedLineWidth(bool on)
{
    depthcue_width = on;
}

bool ViewPoint::isDepthCuedLineWidth()
{
    return depthcue_width;
}

void ViewPoint::setDepthCuedColors(bool on)
{
    depthcue_color = on;
}

bool ViewPoint::isDepthCuedColors()
{
    return depthcue_color;
}

void ViewPoint::setDimNonactiveModels(bool on)
{
    dimNonactiveModels = on;
}

bool ViewPoint::isDimNonactiveModels()
{
    return dimNonactiveModels;
}

void ViewPoint::setAmountToDimNonactiveModels(float percent)
{
    amountToDimNonactiveModels = percent;
}

float ViewPoint::getAmountToDimNonactiveModels()
{
    return amountToDimNonactiveModels;
}

void ViewPoint::SetBallandStick()
{
    ballandstick = BALLANDSTICK;
}

void ViewPoint::SetBallandStick(int s)
{
    ballandstick = s;
}

void ViewPoint::SetBallandCylinder()
{
    ballandstick = BALLANDCYLINDER;
}

void ViewPoint::SetSpaceFilling()
{
    ballandstick = CPK;
}

void ViewPoint::SetSticks()
{
    ballandstick = STICKS;
}

void ViewPoint::ToggleBallandStick()
{
    ballandstick = ballandstick ^ 1;
}

int ViewPoint::GetBallandStick()
{
    return ballandstick;
}

int ViewPoint::GetBallSize()
{
    return ballsize;
}

void ViewPoint::SetBallSize(int b)
{
    ballsize = b;
}

int ViewPoint::GetBallMode()
{
    return ballmode;
}

void ViewPoint::SetBallMode(int b)
{
    ballmode = b;
}

int ViewPoint::GetCylinderSize()
{
    return cylindersize;
}

void ViewPoint::SetCylinderSize(int b)
{
    cylindersize = b;
}

int ViewPoint::linewidth(int w)
{
    if (depthcue_width)
    {
        return (w);
    }
    else
    {
        return (lineThickness);
    }
}

void ViewPoint::setscale(long s)
{
    scale = s;
    if (scale < 11)
    {
        scale = 11;
    }
    clampScale();
    changed = 1;
}

long ViewPoint::getscale()
{
    return (scale);
}

float ViewPoint::getbackclip()
{
    return ((float)backclip/(float)CRS);
}

float ViewPoint::getfrontclip()
{
    return ((float)frontclip/(float)CRS);
}

int ViewPoint::getbackclipi()
{
    return (ROUND(backclip*scale/(float)CRS/10.0F));
}

int ViewPoint::getfrontclipi()
{
    return (ROUND(frontclip*scale/(float)CRS/10.0F));
}

void ViewPoint::setfrontclip(float f)
{
    if (f > 0.0)
    {
        frontclip = (int)(f*CRS);
    }
    else
    {
        frontclip = 0;
    }
    changed = 1;
}

void ViewPoint::setbackclip(float b)
{
    if (b < 0.0)
    {
        backclip = (int)(b*CRS);
    }
    else
    {
        backclip = 0;
    }
    changed = 1;
}

float ViewPoint::getperspective()
{
    return perspective;
}

float ViewPoint::getcenter(int i)
{
    if (i >= 0 && i < 3)
    {
        return ((float)center[i]/100.0F);
    }
    else
    {
        return (0.0F);
    }
}

int ViewPoint::getcenteri(int i)
{
    if (i >= 0 && 1 < 3)
    {
        return (center[i]);
    }
    else
    {
        return (0);
    }
}

void ViewPoint::orthomat()
{
    orthomatrix(viewmat, viewmat);
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
    return (float)((float)(width_)*(float)CRS/(float)scale/10.0F);
}

float ViewPoint::getheight()
{
    return (float)((float)(height_)*(float)CRS/(float)scale/10.0F);
}

void ViewPoint::LineThickness(int n)
{
    if (n < 0)
    {
        n = 0;
    }
    lineThickness = n;
}

int ViewPoint::GetLineThickness()
{
    return lineThickness;
}

void ViewPoint::SetLineThickness(int w)
{
    lineThickness = w;
}

bool ViewPoint::ischanged()
{
    return (changed);
}

void ViewPoint::clearchanged()
{
    changed = false;
}

void ViewPoint::setchanged()
{
    changed = true;
}

void ViewPoint::setTopView(int on)
{
    topview = on;
}

int ViewPoint::getTopView()
{
    return topview;
}

bool ViewPoint::UnDoable()
{
    return UndoList.size() > 0;
}
