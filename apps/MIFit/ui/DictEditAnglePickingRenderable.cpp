#ifdef USE_DICT_EDITOR

#include "DictEditAnglePickingRenderable.h"

#include "GLRenderer.h"
#include <opengl/Frustum.h>
#include <opengl/StereoView.h>

using namespace mi::opengl;
using namespace chemlib;

DictEditAnglePickingRenderable::DictEditAnglePickingRenderable(StereoView* stereoView, std::vector<ANGLE>& angles, Frustum* frustum)
  : stereoView(stereoView) {

  renderer = new GLRenderer();
  renderer->setPickingEnabled(true);
  renderer->setFrustum(frustum);
  renderer->setRenderStyle(RenderStyle::getDefaultStick());

  angleRenderable = new DictEditAngleRenderable(renderer, angles);

}

DictEditAnglePickingRenderable::~DictEditAnglePickingRenderable() {
  delete renderer;
  delete angleRenderable;
}

void DictEditAnglePickingRenderable::render() {
  stereoView->render(*angleRenderable);
}

ANGLE* DictEditAnglePickingRenderable::getAngle(int id) {
  ANGLE* a = renderer->getAngle(id);
  renderer->clearPickNames();
  return a;
}

DictEditAngleRenderable::DictEditAngleRenderable(GLRenderer* renderer, std::vector<ANGLE>& angles)
  : angles(angles), renderer(renderer) {

}

void DictEditAngleRenderable::render() {
  renderer->clearPickNames();
  renderer->drawAngles(angles);
}


#endif // USE_DICT_EDITOR
