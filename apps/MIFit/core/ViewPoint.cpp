#include "ViewPoint.h"

#include <cmath>
#include <cstdio>
#include <math/Quaternion.h>
#include <math/Vector3.h>

using namespace mi::math;

ViewPoint::ViewPoint()
    : view(Vector3<float>(0, 0, 1), 0),
      width_(200), height_(100),
      perspective_(0),
      scale_(20),
      frontClip_(5.0), backClip_(-5.0)

{
    zangle = sin(perspective_ * M_PI/180.0);
}

void ViewPoint::UnDo()
{
    if (UndoList.size() > 0)
    {
        center_ = UndoList.back().center;
        view = UndoList.back().view;
        scale_ = UndoList.back().scale;
        frontClip_ = UndoList.back().frontClip;
        backClip_ = UndoList.back().backClip;
        UndoList.pop_back();
    }
}

void ViewPoint::Do()
{
    ViewSave last;
    last.center = center_;
    last.view = view;
    last.scale = scale_;
    last.frontClip = frontClip_;
    last.backClip = backClip_;
    UndoList.push_back(last);
}

void ViewPoint::setView(const Quaternion<float> &q)
{
    view = q;
}

void ViewPoint::setView(const Matrix4<float> &m)
{
    view.set(m);
    view.normalize();
}

const mi::math::Quaternion<float> &ViewPoint::orientation() const
{
    return view;
}


void ViewPoint::setSize(int width, int height)
{
    width_ = width;
    height_ = height;
}

void ViewPoint::rotate(float rx, float ry, float rz)
{
    Quaternion<float> q(rx/2 * DEG2RAD, ry/2 * DEG2RAD, rz/2 * DEG2RAD, 1.0);
    view = q * view;
}

const mi::math::Vector3<float> &ViewPoint::center() const
{
    return center_;
}

void ViewPoint::setCenter(const mi::math::Vector3<float> &pos)
{
    center_ = pos;
}

void ViewPoint::moveBy(const mi::math::Vector3<float> &delta)
{
    center_ += delta;
}


float ViewPoint::perspective() const
{
    return perspective_;
}

void ViewPoint::setPerspective(float p)
{
    perspective_ = std::max(0.0f, p);
    zangle = sin(perspective_*M_PI/180.0);
}

void ViewPoint::zoom(float ds)
{
    scale_ *= ds;
}

void ViewPoint::changeSlab(float delta)
{
    frontClip_ += delta/2.0;
    backClip_ -= delta/2.0;
    if (frontClip_ < backClip_)
        frontClip_ = backClip_ + 1;
}

void ViewPoint::setSlab(float frontClip, float backClip)
{
    frontClip_ = frontClip;
    backClip_ = backClip;
    if (frontClip_ < backClip_)
        std::swap(frontClip_, backClip_);
}

void ViewPoint::setScale(qreal s)
{
    scale_ = s;
    if (scale_ < 1.1)
        scale_ = 1.1;
}

qreal ViewPoint::scale() const
{
    return scale_;
}

qreal ViewPoint::frontClip() const
{
    return frontClip_;
}

qreal ViewPoint::backClip() const
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

int ViewPoint::width() const
{
    return width_;
}

int ViewPoint::height() const
{
    return height_;
}


float ViewPoint::getwidth() const
{
    return width_ / scale_;
}

float ViewPoint::getheight() const
{
    return height_ / scale_;
}

bool ViewPoint::UnDoable() const
{
    return UndoList.size() > 0;
}
