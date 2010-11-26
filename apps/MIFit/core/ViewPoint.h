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
    float center[3];
    qreal frontClip;
    qreal backClip;
    qreal scale;
};

class QVector3D;

/**
 * The viewpoint contains the current center rotation and all other information about the viewpoint.
 */
class ViewPoint
{
public:

    ViewPoint();

    void PutOnLeft(const QVector3D& point);
    void PutVertical(const QVector3D& point);

    void rotate(float rx, float ry, float rz);

    void copymatrix(float mat[3][3]);
    void setmatrix(float mat[3][3]);

    QVector3D transform(const QVector3D& point) const;
    void moveto(float, float, float);
    void setScale(qreal s);
    qreal scale() const;
    void zoom(float ds);
    void scroll(int, int, int);
    void changeSlab(float delta);
    void setSlab(float frontClip, float backClip);
    qreal frontClip() const;
    qreal backClip() const;
    void setFrontClip(qreal f);
    void setBackClip(qreal b);

    float perspective() const;
    void setPerspective(float);
    float getcenter(int i) const;
    void setSize(int, int);
    int width() const;  // in screen coordinates
    int height() const;  // in screen coordinates
    float getwidth() const;  // in ? coordinates
    float getheight() const;  // in ? coordinates
    void getdirection(int, int, float&, float&, float&);

    void Do();
    void UnDo();
    bool UnDoable() const;

    QVector3D Invert(int sx, int sy, int sz);

private:
    float viewmat[3][3];
    float umat[3][3];
    float center[3];

    int width_, height_; // viewport size; screen coords.

    float perspective_;
    float zangle; //sin of perspective angle;
    qreal scale_;   // multiply by this factor.

    // clip in model coords.
    qreal frontClip_;
    qreal backClip_;

    std::vector<ViewSave> UndoList;

    void transform(const QVector3D& point, int &xt, int &yt, int &zt) const;
};

#endif // ifndef MIFIT_VIEWPOINT_H_
