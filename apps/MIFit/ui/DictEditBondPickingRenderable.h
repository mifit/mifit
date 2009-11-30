#ifndef MIFIT_UI_DICTEDITBONDPICKINGRENDERABLE_H_
#define MIFIT_UI_DICTEDITBONDPICKINGRENDERABLE_H_

#include <opengl/Renderable.h>
#include <vector>

#include <chemlib/chemlib.h>
class GLRenderer;

namespace mi
{
    namespace opengl
    {
        class Frustum;
        class StereoView;
    }
}

class DictEditBondRenderable;

class DictEditBondPickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    GLRenderer *renderer;

    DictEditBondRenderable *bondRenderable;

public:

    DictEditBondPickingRenderable(mi::opengl::StereoView *stereoView, std::vector<chemlib::Bond> &bonds, mi::opengl::Frustum *frustum);
    virtual ~DictEditBondPickingRenderable();

    virtual void render();

    chemlib::Bond *getBond(int id);

};

class DictEditBondRenderable : public mi::opengl::Renderable
{

    std::vector<chemlib::Bond> &bonds;

    GLRenderer *renderer;

public:

    DictEditBondRenderable(GLRenderer *renderer, std::vector<chemlib::Bond> &bonds);

    virtual void render();

};

#endif // ifndef MIFIT_UI_DICTEDITBONDPICKINGRENDERABLE_H_
