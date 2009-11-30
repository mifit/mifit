#ifndef mi_opengl_interact_SimpleMouseRotator_h
#define mi_opengl_interact_SimpleMouseRotator_h

#include <math/Vector3.h>
#include <math/Quaternion.h>

namespace mi
{
    namespace opengl
    {

        class Viewpoint;

        namespace interact
        {

            class SimpleMouseRotator
            {

                Viewpoint *viewpoint;

                int previousX;

                int previousY;

                float angle;

                mi::math::Vector3<float> verticalAxis;
                mi::math::Vector3<float> horizontalAxis;

            public:

                SimpleMouseRotator(Viewpoint *viewpoint, float angle = 1.0f);

                ~SimpleMouseRotator();

                float getAngle();

                void setAngle(float angle);

                mi::math::Vector3<float> getVerticalAxis();

                void setVerticalAxis(mi::math::Vector3<float> &verticalAxis);

                mi::math::Vector3<float> getHorizontalAxis();

                void setHorizontalAxis(mi::math::Vector3<float> &horizontalAxis);

                void beginRotate(int x, int y);

                void rotate(int x, int y);

                void endRotate();

            };

        }
    }
}

#endif // ifndef mi_opengl_interact_SimpleMouseRotator_h
