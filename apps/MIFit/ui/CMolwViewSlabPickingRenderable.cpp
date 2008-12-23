#include "chemlib.h"
#include "CMolwViewSlabPickingRenderable.h"

#include <opengl/Camera.h>
#include <opengl/Frustum.h>
#include <opengl/OpenGL.h>
#include <opengl/QuatUtil.h>
#include <opengl/StereoView.h>
#include <opengl/Sphere.h>

using namespace mi::math;
using namespace mi::opengl;
using namespace chemlib;

CMolwViewSlabPickingRenderable::CMolwViewSlabPickingRenderable(StereoView* stereoView)
  : stereoView(stereoView) {

  slabRenderable = new CMolwViewSlabRenderable();

}

CMolwViewSlabPickingRenderable::~CMolwViewSlabPickingRenderable() {
  delete slabRenderable;
}

void CMolwViewSlabPickingRenderable::render() {
  stereoView->render(*slabRenderable);
}

void CMolwViewSlabPickingRenderable::setCorners(Vector3<float>* corners) {
  slabRenderable->setCorners(corners);
}

void CMolwViewSlabPickingRenderable::setViewpoint(const Viewpoint& viewpoint) {
  slabRenderable->setViewpoint(viewpoint);
}

CMolwViewSlabRenderable::CMolwViewSlabRenderable() {
  viewpoint = new Viewpoint();
  corners = new Vector3<float>[8];
}

CMolwViewSlabRenderable::~CMolwViewSlabRenderable() {
  delete viewpoint;
  delete corners;
}

void CMolwViewSlabRenderable::setCorners(Vector3<float>* corners) {
  for (int i = 0; i < 8; ++i) {
    this->corners[i].set(corners[i]);
  }
}

void CMolwViewSlabRenderable::setViewpoint(const Viewpoint& viewpoint) {
  this->viewpoint->set(viewpoint);
}

void CMolwViewSlabRenderable::render() {

  Quaternion<float> orientation = viewpoint->getRotation();
  Vector3<float> position = viewpoint->getPosition();

  for (int i = 0; i < 8; ++i) {
    Vector3<float>& corner = corners[i];
    corner.z *= -1.0f;
    QuatUtil::rotateVector(orientation, corner);
    corner.add(position);
  }

  glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glDisable(GL_DEPTH_TEST);

  glColor3f(1.0f, 1.0f, 1.0f);
  glLineWidth(2.0f);
  glPushName((GLuint) 1);
  glBegin(GL_LINES);
  glVertex3f(corners[0].x, corners[0].y, corners[0].z);
  glVertex3f(corners[1].x, corners[1].y, corners[1].z);
  glVertex3f(corners[2].x, corners[2].y, corners[2].z);
  glVertex3f(corners[3].x, corners[3].y, corners[3].z);
  glEnd();
  glPopName();

  glPushName((GLuint) 2);
  glBegin(GL_LINES);
  glVertex3f(corners[4].x, corners[4].y, corners[4].z);
  glVertex3f(corners[5].x, corners[5].y, corners[5].z);
  glVertex3f(corners[6].x, corners[6].y, corners[6].z);
  glVertex3f(corners[7].x, corners[7].y, corners[7].z);
  glEnd();
  glPopName();

  glPopAttrib();

}

