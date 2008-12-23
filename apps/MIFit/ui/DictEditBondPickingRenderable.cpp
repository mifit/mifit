#ifdef USE_DICT_EDITOR

#include "DictEditBondPickingRenderable.h"

#include "GLRenderer.h"
#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

using namespace mi::opengl;
using namespace chemlib;

DictEditBondPickingRenderable::DictEditBondPickingRenderable(StereoView* stereoView, std::vector<Bond>& bonds, Frustum* frustum)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  renderer->setFrustum(frustum);
  renderer->setRenderStyle(RenderStyle::getDefaultStick());

  bondRenderable = new DictEditBondRenderable(renderer, bonds);

}

DictEditBondPickingRenderable::~DictEditBondPickingRenderable() {
  delete renderer;
  delete bondRenderable;
}

void DictEditBondPickingRenderable::render() {
  stereoView->render(*bondRenderable);
}

Bond* DictEditBondPickingRenderable::getBond(int id) {
  Bond* e = renderer->getBond(id);
  renderer->clearPickNames();
  return e;
}

DictEditBondRenderable::DictEditBondRenderable(GLRenderer* renderer, std::vector<Bond>& bonds)
  : bonds(bonds), renderer(renderer) {

}

void DictEditBondRenderable::render() {
  renderer->clearPickNames();
  renderer->drawBonds(bonds);
}


#endif // USE_DICT_EDITOR
