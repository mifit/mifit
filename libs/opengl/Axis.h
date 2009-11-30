#ifndef mi_opengl_Axis_h
#define mi_opengl_Axis_h

#include <math/Vector3.h>
#include <math/Quaternion.h>

namespace mi
{
    namespace opengl
    {

        class Axis
        {

        public:

            static const mi::math::Vector3<float> X;

            static const mi::math::Vector3<float> Y;

            static const mi::math::Vector3<float> Z;

            mi::math::Vector3<float> x;

            mi::math::Vector3<float> y;

            mi::math::Vector3<float> z;

            Axis();

            void set(const Axis &otherAxis);

            void rotate(mi::math::Quaternion<float> &rotation);

            void translate(mi::math::Vector3<float> &translation);

            void reset();

        };

    }
}

#endif // ifndef mi_opengl_Axis_h
