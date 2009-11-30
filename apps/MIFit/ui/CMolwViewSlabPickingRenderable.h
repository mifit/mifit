#ifndef mifit_ui_CMolwViewSlabPickingRenderable_h
#define mifit_ui_CMolwViewSlabPickingRenderable_h

#include <math/Vector3.h>
#include <opengl/Renderable.h>
#include <vector>

namespace mi
{
    namespace opengl
    {
        class Camera;
        class Frustum;
        class StereoView;
        class Viewpoint;
    }
}

class CMolwViewSlabRenderable;

class CMolwViewSlabPickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    CMolwViewSlabRenderable *slabRenderable;

public:

    CMolwViewSlabPickingRenderable(mi::opengl::StereoView *stereoView);
    virtual ~CMolwViewSlabPickingRenderable();

    virtual void render();

    void setCorners(mi::math::Vector3<float> *corners);

    void setViewpoint(const mi::opengl::Viewpoint &viewpoint);

};

class CMolwViewSlabRenderable : public mi::opengl::Renderable
{

    friend class CMolwViewSlabPickingRenderable;

    mi::opengl::Viewpoint *viewpoint;

    mi::math::Vector3<float> *corners;

public:

    CMolwViewSlabRenderable();
    virtual ~CMolwViewSlabRenderable();

    virtual void render();

    void setCorners(mi::math::Vector3<float> *corners);

    void setViewpoint(const mi::opengl::Viewpoint &viewpoint);

};

#endif // ifndef mifit_ui_CMolwViewSlabPickingRenderable_h
