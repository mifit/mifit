#ifndef mi_opengl_StereoView_h
#define mi_opengl_StereoView_h

#include <math/Vector3.h>

namespace mi
{
    namespace opengl
    {

        class Frustum;
        class Camera;
        class Viewport;
        class Scene;
        class Renderable;

        class StereoView
        {

            Frustum *frustum;

            Frustum *overlayFrustum;

            Camera *camera;

            Camera *stereoCamera;

            Camera *overlayCamera;

            Viewport *viewport;

            bool glStereo;

            bool stereo;

            bool hardwareStereo;

            Viewport *rightViewport;

            Viewport *leftViewport;

            float eyeSeparation;

            float imageSeparation;

            bool crossEyed;

            void resetView(int flags);

            Viewport *prepareRightEye(float cameraOffset, float frustumAsymmetry);

            Viewport *prepareLeftEye(float cameraOffset, float frustumAsymmetry);

            void calcStereoParameters(float *cameraOffset, float *frustumAsymmetry);

        public:

            StereoView(Frustum *frustum, Camera *camera);

            ~StereoView();

            bool isGLStereo();

            bool isStereo();

            void setStereo(bool stereo);

            bool isHardwareStereo();

            void setHardwareStereo(bool hardwareStereo);

            Viewport *getViewport();

            Viewport *getLeftViewport();

            Viewport *getRightViewport();

            void render(Scene &scene);

            void render(Renderable &renderable);

            Frustum *getFrustum();

            void setFrustum(Frustum *frustum);

            Camera *getCamera();

            void setCamera(Camera *camera);

            void setSize(int width, int height);

            bool isCrossEyed();

            void setCrossEyed(bool crossEyed);

            float getEyeSeparation();

            void setEyeSeparation(float eyeSeparation);

            float getImageSeparation();

            void setImageSeparation(float imageSeparation);

        };
    }
}

#endif // ifndef mi_opengl_StereoView_h
