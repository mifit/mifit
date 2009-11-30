#include "CMolwViewAnnotationPickingRenderable.h"

#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

#include "GLRenderer.h"
#include "Displaylist.h"
#include "MIGLWidget.h"

using namespace mi::opengl;

CMolwViewAnnotationPickingRenderable::CMolwViewAnnotationPickingRenderable(StereoView *stereoView, Frustum *frustum, MIGLWidget *canvas)
    : stereoView(stereoView)
{

    renderer = new GLRenderer();
    renderer->setPickingEnabled(true);
    renderer->setFrustum(frustum);
    renderer->setRenderStyle(RenderStyle::getStick());
    renderer->setQGLWidget(canvas);

    annotationRenderable = new CMolwViewAnnotationRenderable(renderer);

}

CMolwViewAnnotationPickingRenderable::~CMolwViewAnnotationPickingRenderable()
{
    delete renderer;
    delete annotationRenderable;
}

void CMolwViewAnnotationPickingRenderable::render()
{
    renderer->initializeForRender();
    stereoView->render(*annotationRenderable);
}

Annotation*CMolwViewAnnotationPickingRenderable::getAnnotation(int id)
{
    Annotation *a = renderer->getAnnotation(id);
    renderer->clearPickNames();
    return a;
}

void CMolwViewAnnotationPickingRenderable::updateTextScale(float glUnitsPerPixel)
{
    renderer->updateTextScale(glUnitsPerPixel);
}

void CMolwViewAnnotationPickingRenderable::setModels(Displaylist *displaylist)
{
    annotationRenderable->setModels(displaylist);
}

void CMolwViewAnnotationPickingRenderable::setFontSize(int size)
{
    renderer->setFontSize(size);
}

CMolwViewAnnotationRenderable::CMolwViewAnnotationRenderable(GLRenderer *renderer)
    : displaylist(NULL),
      renderer(renderer)
{

}

void CMolwViewAnnotationRenderable::setModels(Displaylist *displaylist)
{
    this->displaylist = displaylist;
}

void CMolwViewAnnotationRenderable::render()
{
    if (displaylist == NULL)
    {
        return;
    }
    renderer->clearPickNames();
    renderer->DrawAnnotations(displaylist->getMolecules());
}

