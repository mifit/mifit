#ifndef mifit_ui_CMolwViewAtomPickingRenderable_h
#define mifit_ui_CMolwViewAtomPickingRenderable_h

#include <opengl/Renderable.h>

namespace chemlib
{
    class MIAtom;
}
class Displaylist;
class GLRenderer;
class RenderStyle;

namespace mi
{
    namespace opengl
    {
        class Camera;
        class Frustum;
        class StereoView;
    }
}

class CMolwViewAtomRenderable;

class CMolwViewAtomPickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    GLRenderer *renderer;

    CMolwViewAtomRenderable *atomRenderable;

public:

    CMolwViewAtomPickingRenderable(mi::opengl::StereoView *stereoView);
    virtual ~CMolwViewAtomPickingRenderable();

    virtual void render();

    chemlib::MIAtom *getAtom(int id);

    void setModels(Displaylist *models);

    void setFrustum(mi::opengl::Frustum *frustum);

    void setCamera(mi::opengl::Camera *camera);

    void setRenderStyle(RenderStyle renderStyle);

};

class CMolwViewAtomRenderable : public mi::opengl::Renderable
{

    Displaylist *displaylist;

    GLRenderer *renderer;

public:

    CMolwViewAtomRenderable(GLRenderer *renderer);

    virtual void render();

    void setModels(Displaylist *models);

};

#endif // ifndef mifit_ui_CMolwViewAtomPickingRenderable_h
