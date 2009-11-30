#ifndef MIFIT_UI_DICTEDITANGLEPICKINGRENDERABLE_H_
#define MIFIT_UI_DICTEDITANGLEPICKINGRENDERABLE_H_

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

class DictEditAngleRenderable;

class DictEditAnglePickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    GLRenderer *renderer;

    DictEditAngleRenderable *angleRenderable;

public:

    DictEditAnglePickingRenderable(mi::opengl::StereoView *stereoView, std::vector<chemlib::ANGLE> &angles, mi::opengl::Frustum *frustum);
    virtual ~DictEditAnglePickingRenderable();

    virtual void render();

    chemlib::ANGLE *getAngle(int id);

};

class DictEditAngleRenderable : public mi::opengl::Renderable
{

    std::vector<chemlib::ANGLE> &angles;

    GLRenderer *renderer;

public:

    DictEditAngleRenderable(GLRenderer *renderer, std::vector<chemlib::ANGLE> &angles);

    virtual void render();




};

#endif // ifndef MIFIT_UI_DICTEDITANGLEPICKINGRENDERABLE_H_
