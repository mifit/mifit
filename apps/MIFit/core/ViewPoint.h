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
    qreal frontClip;
    qreal backClip;
    long scale; // multiply by this factor.
};

class QVector3D;

/**
 * The viewpoint contains the current center rotation and all other information about the viewpoint.
 */
const unsigned int NTYPES = 9;
const float CRS = 100.0f;

class ViewPoint
{
public:

    ViewPoint();

    void PutOnLeft(float x, float y, float z);
    void PutVertical(float x, float y, float z);

    void rotate(float rx, float ry, float rz);

    void copymatrix(float mat[3][3]);
    void setmatrix(float mat[3][3]);

    void transform(float x, float y, float z, int &xt, int &yt, int &zt);
    QVector3D transform(const QVector3D& point);
    void moveto(int, int, int);
    void moveto(float, float, float);
    void translate(int, int, int);
    void setscale(long s);
    long getscale();
    void zoom(float ds);
    void scroll(int, int, int);
    void slab(int);
    void slab(float f, float b);
    qreal frontClip();
    qreal backClip();
    void setFrontClip(qreal f);
    void setBackClip(qreal b);

    float getperspective();
    void setperspective(float);
    float getcenter(int i);
    int getcenteri(int i);
    void setSize(int, int);
    int width();  // in screen coordinates
    int height();  // in screen coordinates
    float getwidth();  // in ? coordinates
    float getheight();  // in ? coordinates
    void getdirection(int, int, int&, int&, int&);

    void Do();
    void UnDo();
    bool UnDoable();

    void Invert(int sx, int sy, int sz, float &x, float &y, float &z);

    bool CenterAtResidue(const chemlib::Residue *res);
    void Center(chemlib::MIMoleculeBase*);


private:
    float viewmat[3][3];
    float umat[3][3];
    long center[3];

    int width_, height_; // viewport size; screen coords.

    float perspective;
    long zangle; //sin of perspective angle * 100;
    long scale;   // multiply by this factor.

    // clip in model coords.
    qreal frontClip_;
    qreal backClip_;

    std::vector<ViewSave> UndoList;

    void clampScale();
};

#endif // ifndef MIFIT_VIEWPOINT_H_
