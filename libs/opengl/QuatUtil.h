#ifndef mi_opengl_QuatUtil_h
#define mi_opengl_QuatUtil_h

#include <math/Vector3.h>
#include <math/Quaternion.h>
#include <math/Matrix4.h>

namespace mi
{
    namespace opengl
    {

        class QuatUtil
        {
        public:
            static mi::math::Quaternion<float> alignVectors(const mi::math::Vector3<float> &source, const mi::math::Vector3<float> &target);

            static void rotateVector(const mi::math::Quaternion<float> &q, mi::math::Vector3<float> &v);

            static void rotateVector(mi::math::Quaternion<float> q, mi::math::Vector3<float> &source, mi::math::Vector3<float> &destination);

            static void store(mi::math::Quaternion<float> q, float *buffer);
        };

    }
}

#endif // ifndef mi_opengl_QuatUtil_h
