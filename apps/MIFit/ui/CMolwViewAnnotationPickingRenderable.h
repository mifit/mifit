#ifndef mifig_ui_CMolwViewAnnotationPickingRenderable_h
#define mifig_ui_CMolwViewAnnotationPickingRenderable_h

#include <opengl/Renderable.h>
#include <list>

class Annotation;
class GLRenderer;
class Displaylist;
class MIGLWidget;

namespace mi
{
    namespace opengl
    {
        class Frustum;
        class StereoView;
    }
}

class CMolwViewAnnotationRenderable;

class CMolwViewAnnotationPickingRenderable : public mi::opengl::Renderable
{

    mi::opengl::StereoView *stereoView;

    GLRenderer *renderer;

    CMolwViewAnnotationRenderable *annotationRenderable;

public:

    CMolwViewAnnotationPickingRenderable(mi::opengl::StereoView *stereoView, mi::opengl::Frustum *frustum, MIGLWidget *canvas);
    virtual ~CMolwViewAnnotationPickingRenderable();

    virtual void render();

    Annotation *getAnnotation(int id);

    void updateTextScale(float glUnitsPerPixel);

    void setModels(Displaylist *displaylist);

    void setFontSize(int size);

};

class CMolwViewAnnotationRenderable : public mi::opengl::Renderable
{

    Displaylist *displaylist;

    GLRenderer *renderer;

public:

    CMolwViewAnnotationRenderable(GLRenderer *renderer);

    virtual void render();

    void setModels(Displaylist *displaylist);

};

#endif // ifndef mifig_ui_CMolwViewAnnotationPickingRenderable_h
