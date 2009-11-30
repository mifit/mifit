#include <opengl/QuatUtil.h>

using namespace mi::math;

namespace mi
{
namespace opengl
{

Quaternion<float> QuatUtil::alignVectors(const Vector3<float> &source, const Vector3<float> &target)
{
    Vector3<float> sourceVector(source);
    sourceVector.normalize();
    Vector3<float> targetVector(target);
    targetVector.normalize();
    Vector3<float> cross;
    cross.cross(sourceVector, targetVector);
    if (cross.lengthSquared() <= 1.0e-7f)
    {
        if (sourceVector.dot(targetVector) <= -1.0f)
        {
            return Quaternion<float>(1.0f, 0.0f, 0.0f, 0.0f);
        }
        return Quaternion<float>(0.0f, 0.0f, 0.0f, 1.0f);
    }
    float dotPlus1 = 1.0f + sourceVector.dot(targetVector);
    float s = std::sqrt(0.5f * dotPlus1);
    cross.scale(0.5f / s);
    return Quaternion<float>(cross.x, cross.y, cross.z, s);
}

void QuatUtil::rotateVector(const Quaternion<float> &q, Vector3<float> &v)
{
    v = q.rotate(v);
    //    Quaternion<float> copy(q);
    //    Quaternion<float> qv(v.x, v.y, v.z, 0.0f);
    //    copy.multiply(qv);
    //    Quaternion<float> conj(q);
    //    conj.inverse();
    //    copy.multiply(conj);
    //	copy.normalize();
    //    copy.scale(v.length());
    //    v.set(copy.x, copy.y, copy.z);
}

void QuatUtil::rotateVector(Quaternion<float> q, Vector3<float> &source, Vector3<float> &destination)
{
    destination = q.rotate(source);
    //    Quaternion<float> copy(q);
    //    Quaternion<float> qv(source.x, source.y, source.z, 0.0f);
    //    copy.multiply(qv);
    //    Quaternion<float> conj(q);
    //    conj.inverse();
    //    copy.multiply(conj);
    //	copy.normalize();
    //    copy.scale(source.length());
    //    destination.set(copy.x, copy.y, copy.z);
}

void QuatUtil::store(Quaternion<float> q, float *buffer)
{
    static Matrix4<float> rotationMatrix;
    rotationMatrix.set(q);
    rotationMatrix.getInColumnMajorOrder(buffer);
}

}
}

