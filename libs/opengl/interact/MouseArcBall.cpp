#include <opengl/interact/MouseArcBall.h>

#include <opengl/Viewpoint.h>

using namespace mi::math;
using namespace mi::opengl;

namespace mi {
namespace opengl {
namespace interact {

MouseArcBall::MouseArcBall(Viewpoint* viewpoint, Point3<float>& center, float radius)
  : historicalQuat(0.0f, 0.0f, 0.0f, 1.0f), currentQuat(0.0f, 0.0f, 0.0f, 1.0f),
  currentVector(0.0f, 0.0f, 0.0f), relativeToViewpoint(NULL), relativeStartQuat(0.0f, 0.0f, 0.0f, 1.0f),
  startVector(NULL), invertRotation(false) {

  this->viewpoint = viewpoint;
  setPlace(center, radius);
}

MouseArcBall::MouseArcBall(Viewpoint* viewpoint)
  : historicalQuat(0.0f, 0.0f, 0.0f, 1.0f), currentQuat(0.0f, 0.0f, 0.0f, 1.0f),
  currentVector(0.0f, 0.0f, 0.0f), relativeToViewpoint(NULL), relativeStartQuat(0.0f, 0.0f, 0.0f, 1.0f),
  startVector(NULL) {

  this->viewpoint = viewpoint;
  setPlace(Point3<float>(0.0f, 0.0f, 0.0f), 1.0f);
}

void MouseArcBall::setRelativeToViewpoint(Viewpoint* relativeToViewpoint) {
  this->relativeToViewpoint = relativeToViewpoint;
}

void MouseArcBall::setPlace(const Point3<float>& center, float radius) {
  this->center.set(center);
  this->radius = radius;
}

Point3<float>& MouseArcBall::getCenter() {
  return center;
}

float MouseArcBall::getRadius() {
  return radius;
}

void MouseArcBall::setRadius(float radius) {
  this->radius = radius;
}

bool MouseArcBall::isInvertRotation() {
  return invertRotation;
}

void MouseArcBall::setInvertRotation(bool invertRotation) {
  this->invertRotation = invertRotation;
}

Vector3<float>* MouseArcBall::getStartVector() {
  return startVector;
}

Vector3<float>* MouseArcBall::getCurrentVector() {
  return &currentVector;
}

void MouseArcBall::windowToSphereCoordinates(Vector3<float>& vect) {
  vect.subtract(center);
  vect.scale(1.0f / radius);
  float magnitude = vect.lengthSquared();
  if (magnitude > 1.0f) {
    float scale = (float) (1.0f / std::sqrt(magnitude));
    vect.scale(scale);
  } else {
    vect.z = (float) std::sqrt(1.0 - magnitude);
  }
}

void MouseArcBall::setCurrentVector(int x, int y) {
  float adjustedX = (2.0f * x - controlSizeX) / controlSizeX * controlAspectX;
  float adjustedY = (controlSizeY - 2.0f * y) / controlSizeY * controlAspectY;
  currentVector.set(adjustedX, adjustedY, 0.0f);
  // Convert window coordinates to sphere coordinates
  windowToSphereCoordinates(currentVector);
}

void MouseArcBall::beginRotate(int x, int y) {
  initializeRotation();
  setCurrentVector(x, y);
  startVector = new Vector3<float>(currentVector);
}

void MouseArcBall::rotate(int x, int y) {
  if (startVector == NULL) {
    return;
  }
  setCurrentVector(x, y);

  // Calculate new rotation
  Vector3<float> cross;
  cross.cross(*startVector, currentVector);
  float dot = startVector->dot(currentVector);
  currentQuat.set(cross.x, cross.y, cross.z, dot);
  if (invertRotation) {
    currentQuat.inverse();
  }

  // Add old rotation
  currentQuat.multiply(historicalQuat);

  applyRotation();
}

void MouseArcBall::endRotate() {
  if (startVector == NULL) {
    return;
  }

  delete startVector;
  startVector = NULL;
  applyRotation();
  historicalQuat.set(currentQuat);
}

void MouseArcBall::initializeRotation() {
  if (relativeToViewpoint == NULL) {
    historicalQuat.set(viewpoint->getRotation());
  } else {
    relativeStartQuat.set(relativeToViewpoint->getRotation());
    historicalQuat.inverse(relativeStartQuat);
    historicalQuat.multiply(viewpoint->getRotation());
  }
}

void MouseArcBall::applyRotation() {
  if (relativeToViewpoint == NULL) {
    viewpoint->setRotation(currentQuat);
  } else {
    Quaternion<float> q(relativeStartQuat);
    q.multiply(currentQuat);
    viewpoint->setRotation(q);
  }
}

void MouseArcBall::setControlSize(int x, int y) {
  controlSizeX = x;
  controlSizeY = y;
  if (controlSizeX > controlSizeY) {
    controlAspectX = (float) controlSizeX / controlSizeY;
    controlAspectY = 1.0f;
  } else {
    controlAspectX = 1.0f;
    controlAspectY = (float) controlSizeY / controlSizeX;
  }
}

void MouseArcBall::getControlSize(int& x, int& y) {
  x = controlSizeX;
  y = controlSizeY;
}

}
}
}
