#ifndef mi_opengl_Scene_h
#define mi_opengl_Scene_h

#include <opengl/Renderable.h>

namespace mi
{
    namespace opengl
    {

        class Viewport;

        class Scene : public Renderable
        {
        public:

            virtual void setViewport(Viewport *viewport) = 0;

            virtual void initializeForRender() = 0;

            virtual void preCameraRender() = 0;

            virtual void renderOverlay() = 0;

            virtual void render2DOverlay() = 0;

        };
    }
}

#endif // ifndef mi_opengl_Scene_h
