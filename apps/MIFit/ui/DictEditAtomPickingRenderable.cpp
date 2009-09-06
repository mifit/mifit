#include <chemlib/chemlib.h>

#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

#include "DictEditAtomPickingRenderable.h"
#include "GLRenderer.h"

using namespace mi::opengl;
using namespace chemlib;

DictEditAtomPickingRenderable::DictEditAtomPickingRenderable(StereoView* stereoView, Displaylist* displaylist, Frustum* frustum)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  renderer->setFrustum(frustum);
  RenderStyle style;
  style.set(RenderStyle::getDefaultBallAndLine());
  style.setBondLine(false);
  renderer->setRenderStyle(style);

  atomRenderable = new DictEditAtomRenderable(renderer, displaylist);

}

DictEditAtomPickingRenderable::~DictEditAtomPickingRenderable() {
  delete renderer;
  delete atomRenderable;
}

void DictEditAtomPickingRenderable::render() {
  stereoView->render(*atomRenderable);
}

MIAtom* DictEditAtomPickingRenderable::getAtom(int id) {
  MIAtom* a = renderer->getAtom(id);
  renderer->clearPickNames();
  return a;
}

DictEditAtomRenderable::DictEditAtomRenderable(GLRenderer* renderer, Displaylist* displaylist)
  : displaylist(displaylist), renderer(renderer) {

}

void DictEditAtomRenderable::render() {
  renderer->clearPickNames();
  renderer->Draw2(displaylist, false, false);
}

