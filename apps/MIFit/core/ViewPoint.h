#ifndef MIFIT_VIEWPOINT_H_
#define MIFIT_VIEWPOINT_H_

class ViewPoint;

#include <math/mathlib.h>
#include <chemlib/chemlib.h>

class CArchive;

class ViewSave
{
public:
    float viewmat[3][3];
    long center[3];
    int frontclip; // clip in model coords.
    int backclip;
    long scale; // multiply by this factor.
};

/**
 * The viewpoint contains the current center rotation and all other information about the viewpoint.
 */
const unsigned int NTYPES = 9;
const float CRS = 100.0f;

class ViewPoint
{
protected:
    float viewmat[3][3];
    float umat[3][3];
    long center[3];

    // coordinate bounding box - if outside
    // there is no need to transform - can't be on screen.
    long xubound, xlbound, yubound, ylbound, zubound, zlbound;

    long cscale; // coords are (int)Angstroms *cscale

    float perspective;
    long zangle; //sin of perspective angle * 100;
    long scale;   // multiply by this factor.
    float rotated;

    // clip in model coords.
    int frontclip;
    int backclip;

    std::vector<ViewSave> UndoList;

    bool changed;

    bool depthcue_color;
    bool depthcue_width;
    bool dimNonactiveModels;
    float amountToDimNonactiveModels;
    int ballandstick;
    int ballsize;
    int ballmode;
    int cylindersize;
    int topview;
    int imageScale;
    int lineThickness;

    void clampScale();

public:

    static const int COLORPRINT;
    static const int BWPRINT;

    static const int LEFT;
    static const int RIGHT;

    //definitions for ballandstick
    static const int STICKS;
    static const int BALLANDSTICK;
    static const int CPK;
    static const int BALLANDCYLINDER;

    bool atom_cross;

    int Yup;
    long xmax, ymax, xmin, ymin; // viewport size; screen coords.
    long sxmax, symax, sxmin, symin; // viewport clipbox with extra edges; screen coords.
    long xcenter, ycenter; // viewport center in screen coords
    int printstyle;
    float Bmin, Bmax;
    int iradius[NTYPES];

    ViewPoint();
    virtual ~ViewPoint()
    {
    }

    bool CenterAtResidue(const chemlib::RESIDUE *res);
    void PutOnLeft(float x, float y, float z);
    void PutVertical(float x, float y, float z);

    void SetRadii();

    void Save(CArchive&);
    void Load(FILE *fp);

    void setDepthCuedLineWidth(bool on);
    bool isDepthCuedLineWidth();
    void setDepthCuedColors(bool on);
    bool isDepthCuedColors();
    void setDimNonactiveModels(bool on);
    bool isDimNonactiveModels();
    void setAmountToDimNonactiveModels(float percent);
    float getAmountToDimNonactiveModels();
    void SetBallandStick();
    void SetBallandStick(int s);
    void SetBallandCylinder();
    void SetSpaceFilling();
    void SetSticks();
    void ToggleBallandStick();
    int GetBallandStick();
    int GetBallSize();
    void SetBallSize(int b);
    int GetBallMode();
    void SetBallMode(int b);
    int GetCylinderSize();
    void SetCylinderSize(int b);
    int linewidth(int w);
    int colordepth(int c);

    void rotx(float a);
    void roty(float a);
    void rotz(float a);

    void rotate(float rx, float ry, float rz);

    float *getmatrix();
    void copymatrix(float mat[3][3]);
    void setmatrix(float mat[3][3]);

    void transx(float);
    void transy(float);
    void transz(float);

    int transform(float x, float y, float z, int &xt, int &yt, int &zt);
    float x(float, float, float);
    float y(float, float, float);
    float z(float, float, float);
    void moveto(int, int, int);
    void moveto(float, float, float);
    void translate(int, int, int);
    void setscale(long s);
    long getscale();
    void zoom(float ds);
    void scroll(int, int);
    void scroll(int, int, int);
    void slab(int);
    void slab(float f, float b);
    float getbackclip();
    float getfrontclip();
    int getbackclipi();
    int getfrontclipi();
    int getsbackclip();
    int getsfrontclip();
    void setfrontclip(float f);
    void setbackclip(float b);

    float getperspective();
    void setperspective(float);
    long getzangle();
    void setzangle(long s);
    void SetCenter(int x, int y);
    void Center(chemlib::MIMoleculeBase*);
    int getcenterx();
    int getcentery();
    float getcenter(int i);
    int getcenteri(int i);
    void orthomat();
    void SetBounds();
    void SetWindow(int, int, int, int);
    float getwidth();
    float getheight();
    void LineThickness(int n);
    int GetLineThickness();
    void SetLineThickness(int w);
    float getrotated();
    void clearrotated();
    void setrotated(float r);
    bool ischanged();
    void clearchanged();
    void setchanged();
    void getdirection(int, int, int&, int&, int&);
    void setTopView(int on);
    int getTopView();

    void Do();
    void UnDo();
    bool UnDoable();

    int GetImageScale();
    void SetImageScale(int s);
    void Invert(int sx, int sy, int sz, float &x, float &y, float &z);
};

inline int ViewPoint::colordepth(int c)
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

inline void ViewPoint::setDepthCuedLineWidth(bool on)
{
    depthcue_width = on;
}

inline bool ViewPoint::isDepthCuedLineWidth()
{
    return depthcue_width;
}

inline void ViewPoint::setDepthCuedColors(bool on)
{
    depthcue_color = on;
}

inline bool ViewPoint::isDepthCuedColors()
{
    return depthcue_color;
}

inline void ViewPoint::setDimNonactiveModels(bool on)
{
    dimNonactiveModels = on;
}

inline bool ViewPoint::isDimNonactiveModels()
{
    return dimNonactiveModels;
}

inline void ViewPoint::setAmountToDimNonactiveModels(float percent)
{
    amountToDimNonactiveModels = percent;
}

inline float ViewPoint::getAmountToDimNonactiveModels()
{
    return amountToDimNonactiveModels;
}

inline void ViewPoint::SetRadii()
{
    for (int unsigned i = 0; i < NTYPES; i++)
    {
        iradius[i] = (int)(chemlib::MIAtom::MIAtomRadiusForType(i) * (float)scale / 10.0);
    }
}

inline void ViewPoint::SetBallandStick()
{
    ballandstick = BALLANDSTICK;
}

inline void ViewPoint::SetBallandStick(int s)
{
    ballandstick = s;
}

inline void ViewPoint::SetBallandCylinder()
{
    ballandstick = BALLANDCYLINDER;
}

inline void ViewPoint::SetSpaceFilling()
{
    ballandstick = CPK;
}

inline void ViewPoint::SetSticks()
{
    ballandstick = STICKS;
}

inline void ViewPoint::ToggleBallandStick()
{
    ballandstick = ballandstick ^ 1;
}

inline int ViewPoint::GetBallandStick()
{
    return ballandstick;
}

inline int ViewPoint::GetBallSize()
{
    return ballsize;
}

inline void ViewPoint::SetBallSize(int b)
{
    ballsize = b;
}

inline int ViewPoint::GetBallMode()
{
    return ballmode;
}

inline void ViewPoint::SetBallMode(int b)
{
    ballmode = b;
    SetRadii();
}

inline int ViewPoint::GetCylinderSize()
{
    return cylindersize;
}

inline void ViewPoint::SetCylinderSize(int b)
{
    cylindersize = b;
}

inline int ViewPoint::linewidth(int w)
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

inline void ViewPoint::rotx(float a)
{
    rotate(a, 0.0F, 0.0F);
}

inline void ViewPoint::roty(float a)
{
    rotate(0.0F, a, 0.0F);
}

inline void ViewPoint::rotz(float a)
{
    rotate(0.0F, 0.0F, a);
}

inline float*ViewPoint::getmatrix()
{
    return &(viewmat[0][0]);
}

inline void ViewPoint::setscale(long s)
{
    scale = s;
    if (scale < 11)
    {
        scale = 11;
    }
    SetBounds();
}

inline long ViewPoint::getscale()
{
    return (scale);
}

inline float ViewPoint::getbackclip()
{
    return ((float)backclip/(float)CRS);
}

inline float ViewPoint::getfrontclip()
{
    return ((float)frontclip/(float)CRS);
}

inline int ViewPoint::getbackclipi()
{
    return (ROUND(backclip*scale/(float)CRS/10.0F));
}

inline int ViewPoint::getfrontclipi()
{
    return (ROUND(frontclip*scale/(float)CRS/10.0F));
}

inline int ViewPoint::getsbackclip()
{
    return (int)(backclip*scale/(float)CRS/10.0F)+ycenter;
}

inline int ViewPoint::getsfrontclip()
{
    return (int)(frontclip*scale/(float)CRS/10.0F)+ycenter;
}

inline void ViewPoint::setfrontclip(float f)
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

inline void ViewPoint::setbackclip(float b)
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

inline float ViewPoint::getperspective()
{
    return perspective;
}

inline long ViewPoint::getzangle()
{
    return zangle;
}

inline void ViewPoint::setzangle(long s)
{
    zangle = s;
}

inline void ViewPoint::SetCenter(int x, int y)
{
    center[0] = x;
    center[1] = y;
    changed = 1;
}

inline int ViewPoint::getcenterx()
{
    return xcenter;
}

inline int ViewPoint::getcentery()
{
    return ycenter;
}

inline float ViewPoint::getcenter(int i)
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

inline int ViewPoint::getcenteri(int i)
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

inline void ViewPoint::orthomat()
{
    orthomatrix(viewmat, viewmat);
}

inline float ViewPoint::getwidth()
{
    return (float)((float)(xmax-xmin)*(float)CRS/(float)scale/10.0F);
}

inline float ViewPoint::getheight()
{
    return (float)((float)(ymax-ymin)*(float)CRS/(float)scale/10.0F);
}

inline void ViewPoint::LineThickness(int n)
{
    if (n < 0)
    {
        n = 0;
    }
    lineThickness = n;
}

inline int ViewPoint::GetLineThickness()
{
    return lineThickness;
}

inline void ViewPoint::SetLineThickness(int w)
{
    lineThickness = w;
}

inline float ViewPoint::getrotated()
{
    return (rotated);
}

inline void ViewPoint::clearrotated()
{
    rotated = 0.0F;
}

inline void ViewPoint::setrotated(float r)
{
    rotated = r;
}

inline bool ViewPoint::ischanged()
{
    return (changed);
}

inline void ViewPoint::clearchanged()
{
    changed = false;
}

inline void ViewPoint::setchanged()
{
    changed = true;
}

inline void ViewPoint::setTopView(int on)
{
    topview = on;
}

inline int ViewPoint::getTopView()
{
    return topview;
}

inline bool ViewPoint::UnDoable()
{
    return UndoList.size() > 0;
}

inline int ViewPoint::GetImageScale()
{
    return imageScale;
}

inline void ViewPoint::SetImageScale(int s)
{
    imageScale = s;
}

bool IsOnScreen(const chemlib::MIAtom *a, ViewPoint *vp);
bool IsOnScreen(const APOINT &p, ViewPoint *vp);




#endif // ifndef MIFIT_VIEWPOINT_H_
