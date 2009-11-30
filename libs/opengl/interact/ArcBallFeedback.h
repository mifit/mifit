#ifndef ARCBALLFEEDBACK_H_
#define ARCBALLFEEDBACK_H_

#include <math/Point3.h>
#include <math/Vector3.h>

namespace mi
{
    namespace opengl
    {

        class Circle;
        class Arc;

        namespace interact
        {

            class ArcBallFeedback
            {
                mi::opengl::Circle *ball;

                mi::opengl::Arc *arc;

                mi::math::Point3<float> arcColor;

                mi::math::Point3<float> ballColor;

            public:

                ArcBallFeedback(float radius, float *color);

                void render();

                void setRadius(float radius);

                void setFrom(mi::math::Vector3<float> &from);

                void setTo(mi::math::Vector3<float> &to);

            };
        }
    }
}

#endif // ifndef ARCBALLFEEDBACK_H_
