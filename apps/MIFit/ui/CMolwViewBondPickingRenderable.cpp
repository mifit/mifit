#include "CMolwViewBondPickingRenderable.h"

#include <chemlib/chemlib.h>
#include "GLRenderer.h"
#include "core/RenderStyle.h"
#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

using namespace chemlib;
using namespace mi::opengl;

CMolwViewBondPickingRenderable::CMolwViewBondPickingRenderable(StereoView* stereoView)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  renderer->setJoinBondsOfSameColor(false);

  bondRenderable = new CMolwViewBondRenderable(renderer);

}

CMolwViewBondPickingRenderable::~CMolwViewBondPickingRenderable() {
  delete renderer;
  delete bondRenderable;
}

void CMolwViewBondPickingRenderable::render() {
  stereoView->render(*bondRenderable);
}

Bond* CMolwViewBondPickingRenderable::getBond(int id) {
  Bond* b = renderer->getBond(id);
  renderer->clearPickNames();
  return b;
}

void CMolwViewBondPickingRenderable::setMolecule(Molecule* molecule) {
  bondRenderable->molecule = molecule;
}

void CMolwViewBondPickingRenderable::setFrustum(Frustum* frustum) {
  renderer->setFrustum(frustum);
}

void CMolwViewBondPickingRenderable::setCamera(Camera* camera) {
  renderer->setCamera(camera);
}

void CMolwViewBondPickingRenderable::setRenderStyle(RenderStyle renderStyle) {
  renderer->setRenderStyle(renderStyle);
}

CMolwViewBondRenderable::CMolwViewBondRenderable(GLRenderer* renderer)
  : molecule(NULL), renderer(renderer) {

}

void CMolwViewBondRenderable::render() {
  if (molecule == NULL) {
    return;
  }
  renderer->clearPickNames();
  renderer->drawBonds(molecule->getBonds());
}

