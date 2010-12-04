#include <opengl/Viewpoint.h>

#include <opengl/QuatUtil.h>
#include <opengl/OpenGL.h>

using namespace mi::math;

namespace mi
{
namespace opengl
{

Viewpoint::Viewpoint()
    : position(0.0f, 0.0f, 0.0f),
      rotation(0.0f, 0.0f, 0.0f, 1.0f)
{
}

void Viewpoint::set(const Viewpoint &viewpoint)
{
    position.set(viewpoint.getPosition());
    rotation.set(viewpoint.getRotation());
}

Vector3<float> Viewpoint::getPosition() const
{
    return position;
}

void Viewpoint::setPosition(const Vector3<float> &position)
{
    this->position.set(position);
}

void Viewpoint::translate(Vector3<float> &translation)
{
    position.add(translation);
}

Quaternion<float> Viewpoint::getRotation() const
{
    return rotation;
}

void Viewpoint::setRotation(const Quaternion<float> &rotation)
{
    this->rotation.set(rotation);
}

void Viewpoint::rotate(Quaternion<float> &rotation)
{
    this->rotation.multiply(rotation);
}

void Viewpoint::render()
{
    static float *floatBuffer = new float[16];
    glTranslatef(position.x, position.y, position.z);
    QuatUtil::store(rotation, floatBuffer);
    glMultMatrixf(floatBuffer);
}

void Viewpoint::transformVector(Vector3<float> &vector)
{
    QuatUtil::rotateVector(rotation, vector);
}

void Viewpoint::untransformVector(Vector3<float> &vector)
{
    Quaternion<float> reverseRotation(rotation);
    reverseRotation.conjugate();
    QuatUtil::rotateVector(reverseRotation, vector);
}

Vector3<float> Viewpoint::getTargetVector(float focalLength)
{
    Vector3<float> targetVector(0.0f, 0.0f, -focalLength);
    transformVector(targetVector);
    return targetVector;
}

Vector3<float> Viewpoint::getTarget(float focalLength)
{
    Vector3<float> target(getTargetVector(focalLength));
    target.add(getPosition());
    return target;
}

Vector3<float> Viewpoint::getViewVector()
{
    return getTargetVector(1.0f);
}

Vector3<float> Viewpoint::getRightVector()
{
    Vector3<float> rightVector(1.0f, 0.0f, 0.0f);
    transformVector(rightVector);
    return rightVector;
}

Vector3<float> Viewpoint::getUpVector()
{
    Vector3<float> upVector(0.0f, 1.0f, 0.0f);
    transformVector(upVector);
    return upVector;
}

}
}
