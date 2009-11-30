#ifndef MIFIT_UI_DICTEDITSCENE_H_
#define MIFIT_UI_DICTEDITSCENE_H_

#include <opengl/Scene.h>

#include <chemlib/chemlib.h>

class Displaylist;
class GLRenderer;
class Stack;

namespace mi
{
    namespace opengl
    {
        class Light;
        class Viewpoint;
        class Frustum;
        class Camera;
        namespace interact
        {
            class ArcBallFeedback;
            class MouseArcBallOrbitor;
            class MouseTranslator;
            class TargetFeedback;
        }
    }
}

class DictEditScene : public mi::opengl::Scene
{

    friend class DictEditCanvas;

    mi::opengl::Light *light;

    mi::opengl::Viewpoint *viewpoint;

    mi::opengl::Frustum *frustum;

    mi::opengl::Camera *camera;

    mi::opengl::interact::MouseArcBallOrbitor *cameraMouseOrbitor;

    mi::opengl::interact::MouseTranslator *mouseTranslator;

    mi::opengl::interact::ArcBallFeedback *arcBallFeedback;

    mi::opengl::interact::TargetFeedback *targetFeedback;

    float glUnitsPerPixel;

    bool showAtomLabels;

    GLRenderer *renderer;

    Displaylist *models;

    Stack *atomStack;

    chemlib::Bond *pickedBond;

    chemlib::ANGLE *pickedAngle;

    chemlib::PLANE *pickedPlane;

    mi::opengl::Viewport *viewport;

public:

    DictEditScene();
    virtual ~DictEditScene();

    virtual void setViewport(mi::opengl::Viewport *viewport);

    virtual void initializeForRender();

    virtual void preCameraRender();

    virtual void render();

    virtual void renderOverlay();

    virtual void render2DOverlay();

};

#endif // ifndef MIFIT_UI_DICTEDITSCENE_H_
