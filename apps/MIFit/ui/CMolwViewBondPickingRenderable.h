#ifndef mifit_ui_CMolwViewBondPickingRenderable_h
#define mifit_ui_CMolwViewBondPickingRenderable_h

#include <opengl/Renderable.h>
#include <vector>

namespace chemlib
{
    class Bond;
}
class Molecule;
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

class CMolwViewBondRenderable;

class CMolwViewBondPickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    GLRenderer *renderer;

    CMolwViewBondRenderable *bondRenderable;

public:

    CMolwViewBondPickingRenderable(mi::opengl::StereoView *stereoView);
    virtual ~CMolwViewBondPickingRenderable();

    virtual void render();

    chemlib::Bond *getBond(int id);

    void setMolecule(Molecule *molecule);

    void setFrustum(mi::opengl::Frustum *frustum);

    void setCamera(mi::opengl::Camera *camera);

    void setRenderStyle(RenderStyle renderStyle);

};

class CMolwViewBondRenderable : public mi::opengl::Renderable
{

    friend class CMolwViewBondPickingRenderable;

    Molecule *molecule;

    GLRenderer *renderer;

public:

    CMolwViewBondRenderable(GLRenderer *renderer);

    virtual void render();

};

#endif // ifndef mifit_ui_CMolwViewBondPickingRenderable_h
