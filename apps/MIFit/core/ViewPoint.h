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

    bool CenterAtResidue(const chemlib::Residue *res);
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

private:
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
};

bool IsOnScreen(const chemlib::MIAtom *a, ViewPoint *vp);
bool IsOnScreen(const APOINT &p, ViewPoint *vp);


#endif // ifndef MIFIT_VIEWPOINT_H_
