#ifndef mi_opengl_Frustum_h
#define mi_opengl_Frustum_h

#include <math/Vector3.h>

namespace mi
{
    namespace opengl
    {

        class Viewport;
        class Plane;

        class Frustum
        {

            float nearClipping;

            float farClipping;

            float fieldOfView;

            bool perspective;

            float frustumTop;

            float frustumBottom;

            float frustumLeft;

            float frustumRight;

            float frustumOffset;

            Plane **planes;

            float focalLength;

            float *projection;
            float *modelview;
            float *planesMatrix;

            bool picking;

            int pickPointX;

            int pickPointY;

            float getPerspectiveTop(float z);

        public:

            static const int TOP;

            static const int BOTTOM;

            static const int LEFT;

            static const int RIGHT;

            static const int NEARPLANE;

            static const int FARPLANE;


            Frustum();
            ~Frustum();

            void setFieldOfView(float fieldOfView);

            bool isPerspective();

            void setPerspective(bool perspective);

            void beginPicking(int x, int y);

            void endPicking();

            void render(Viewport &viewport);

            void render(Viewport &viewport, float offset);
            void setHeight(float height);
            void updateFrustum(Viewport &viewport);

            mi::math::Vector3<float> *getCorners();

            enum CullingResult
            {
                OUTSIDE, INTERSECT, INSIDE
            };


            float getFrustumWidth();
            float getFrustumHeight();

            float getNearClipping();

            float getFieldOfView();

            float getFrustumOffset();

            void setFrustumOffset(float frustumAsymmetry);

            void updatePlanes();

            CullingResult pointInFrustum(mi::math::Vector3<float> &position);

            CullingResult sphereInFrustum(mi::math::Vector3<float> &position, float radius);
            float getPerspectiveHeight(float z);

            float getFrustumHeight(float z);

            void setNearClipping(float nearClipping);
            float getFarClipping();

            void setFarClipping(float farClipping);

            float getFocalLength();

            void setFocalLength(float focalLength);

            float getHeight();

            float getWidth();

            Plane **getPlanes();

        };

    }
}


#endif // ifndef mi_opengl_Frustum_h
