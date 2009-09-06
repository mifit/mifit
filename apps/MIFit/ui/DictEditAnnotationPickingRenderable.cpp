#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

#include "DictEditAnnotationPickingRenderable.h"
#include "GLRenderer.h"
#include "Displaylist.h"
#include "DictEditCanvas.h"

using namespace mi::opengl;

DictEditAnnotationPickingRenderable::DictEditAnnotationPickingRenderable(StereoView* stereoView, Displaylist* displaylist, Frustum* frustum, DictEditCanvas* canvas)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  renderer->setFrustum(frustum);
  renderer->setRenderStyle(RenderStyle::getDefaultStick());
  renderer->setQGLWidget(canvas);

  annotationRenderable = new DictEditAnnotationRenderable(renderer, displaylist);

}

DictEditAnnotationPickingRenderable::~DictEditAnnotationPickingRenderable() {
  delete renderer;
  delete annotationRenderable;
}

void DictEditAnnotationPickingRenderable::render() {
  renderer->initializeForRender();
  stereoView->render(*annotationRenderable);
}

Annotation* DictEditAnnotationPickingRenderable::getAnnotation(int id) {
  Annotation* a = renderer->getAnnotation(id);
  renderer->clearPickNames();
  return a;
}

void DictEditAnnotationPickingRenderable::updateTextScale(float glUnitsPerPixel) {
  renderer->updateTextScale(glUnitsPerPixel);
}

DictEditAnnotationRenderable::DictEditAnnotationRenderable(GLRenderer* renderer, Displaylist* displaylist)
  : displaylist(displaylist), renderer(renderer) {

}

void DictEditAnnotationRenderable::render() {
  renderer->clearPickNames();
  renderer->DrawAnnotations(displaylist->getMolecules());
}

