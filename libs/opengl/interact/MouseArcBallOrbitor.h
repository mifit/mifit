#ifndef mi_opengl_interact_MouseArcBallOrbitor_h
#define mi_opengl_interact_MouseArcBallOrbitor_h

#include <opengl/interact/MouseArcBall.h>
#include <math/Point3.h>
#include <math/Quaternion.h>

namespace mi
{
    namespace opengl
    {

        class Viewpoint;

        namespace interact
        {

            class MouseArcBallOrbitor : public MouseArcBall
            {

                mi::math::Quaternion<float> startQuat;

                float distanceToTarget;

                virtual void initializeRotation();

                virtual void applyRotation();

            public:

                MouseArcBallOrbitor(mi::opengl::Viewpoint *viewpoint, float distanceToTarget);

                MouseArcBallOrbitor(mi::opengl::Viewpoint *viewpoint, mi::math::Point3<float> &center, float radius, float distanceToTarget);

                void setDistanceToTarget(float distanceToTarget);

            };

        }
    }
}

#endif // ifndef mi_opengl_interact_MouseArcBallOrbitor_h

