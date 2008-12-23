#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/StereoView.h>
#include <list>

#include "chemlib.h"

#include "CMolwViewAtomPickingRenderable.h"
#include "GLRenderer.h"
#include "RenderStyle.h"
#include "Displaylist.h"


using namespace chemlib;
using namespace mi::opengl;

CMolwViewAtomPickingRenderable::CMolwViewAtomPickingRenderable(StereoView* stereoView)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  atomRenderable = new CMolwViewAtomRenderable(renderer);

}

CMolwViewAtomPickingRenderable::~CMolwViewAtomPickingRenderable() {
  delete renderer;
  delete atomRenderable;
}

void CMolwViewAtomPickingRenderable::render() {
  stereoView->render(*atomRenderable);
}

MIAtom* CMolwViewAtomPickingRenderable::getAtom(int id) {
  MIAtom* a = renderer->getAtom(id);
  renderer->clearPickNames();
  return a;

}

void CMolwViewAtomPickingRenderable::setModels(Displaylist* displaylist) {
  atomRenderable->setModels(displaylist);
}

void CMolwViewAtomPickingRenderable::setFrustum(Frustum* frustum) {
  renderer->setFrustum(frustum);
}

void CMolwViewAtomPickingRenderable::setCamera(Camera* camera) {
  renderer->setCamera(camera);
}

void CMolwViewAtomPickingRenderable::setRenderStyle(RenderStyle renderStyle) {
  renderer->setRenderStyle(renderStyle);
}

CMolwViewAtomRenderable::CMolwViewAtomRenderable(GLRenderer* renderer)
  : displaylist(NULL), renderer(renderer) {

}

void CMolwViewAtomRenderable::render() {
  if (displaylist == NULL) {
    return;
  }
  renderer->clearPickNames();
  std::list<Molecule*>::iterator node = displaylist->getMolecules().begin();
  while (node != displaylist->getMolecules().end()) {
    renderer->drawMolecule(*node);
    renderer->drawSymmetryMolecule(*node, false);
    node++;
  }
}

void CMolwViewAtomRenderable::setModels(Displaylist* displaylist) {
  this->displaylist = displaylist;
}

