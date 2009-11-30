#include <opengl/Axis.h>

#include <opengl/QuatUtil.h>

using namespace mi::math;

namespace mi
{
namespace opengl
{

const Vector3<float> Axis::X(1.0f, 0.0f, 0.0f);

const Vector3<float> Axis::Y(0.0f, 1.0f, 0.0f);

const Vector3<float> Axis::Z(0.0f, 0.0f, 1.0f);

Axis::Axis()
    : x(1.0f, 0.0f, 0.0f),
      y(0.0f, 1.0f, 0.0f),
      z(0.0f, 0.0f, 1.0f)
{
}

void Axis::set(const Axis &otherAxis)
{
    x.set(otherAxis.x);
    y.set(otherAxis.y);
    z.set(otherAxis.z);
}

void Axis::rotate(Quaternion<float> &rotation)
{
    QuatUtil::rotateVector(rotation, x);
    QuatUtil::rotateVector(rotation, y);
    QuatUtil::rotateVector(rotation, z);
}

void Axis::translate(Vector3<float> &translation)
{
    x.add(translation);
    y.add(translation);
    z.add(translation);
}

void Axis::reset()
{
    x.set(X);
    y.set(Y);
    z.set(Z);
}

}
}
