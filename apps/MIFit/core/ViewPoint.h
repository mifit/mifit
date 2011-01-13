#ifndef MIFIT_VIEWPOINT_H_
#define MIFIT_VIEWPOINT_H_

class ViewPoint;

#include <math/mathlib.h>
#include <math/Matrix4.h>
#include <math/Quaternion.h>
#include <math/Vector3.h>
#include <chemlib/chemlib.h>

class CArchive;

class ViewSave
{
public:
    float viewmat[3][3];
    mi::math::Quaternion<float> view;
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

    void rotate(float rx, float ry, float rz);

    void set(const mi::math::Quaternion<float> &q);
    void set(const mi::math::Matrix4<float> &mat);
    const mi::math::Quaternion<float> &orientation() const;

    void moveto(float, float, float);
    void moveBy(const mi::math::Vector3<float> &delta);

    void setScale(qreal s);
    qreal scale() const;
    void zoom(float ds);
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

    void Do();
    void UnDo();
    bool UnDoable() const;

private:
    mi::math::Quaternion<float> view;
    float center[3];

    int width_, height_; // viewport size; screen coords.

    float perspective_;
    float zangle; //sin of perspective angle;
    qreal scale_;   // multiply by this factor.

    // clip in model coords.
    qreal frontClip_;
    qreal backClip_;

    std::vector<ViewSave> UndoList;
};

#endif // ifndef MIFIT_VIEWPOINT_H_
