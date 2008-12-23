#include <opengl/interact/SimpleMouseRotator.h>

#include <opengl/Axis.h>
#include <opengl/Viewpoint.h>
#include <opengl/QuatUtil.h>

using namespace mi::math;

namespace mi {
namespace opengl {
namespace interact {

SimpleMouseRotator::SimpleMouseRotator(Viewpoint* viewpoint, float angle) {
  this->viewpoint = viewpoint;
  this->angle = angle;
  verticalAxis.set(Axis::X);
  horizontalAxis.set(Axis::Y);
}

SimpleMouseRotator::~SimpleMouseRotator() {
}

float SimpleMouseRotator::getAngle() {
  return angle;
}

void SimpleMouseRotator::setAngle(float angle) {
  this->angle = angle;
}

Vector3<float> SimpleMouseRotator::getVerticalAxis() {
  return verticalAxis;
}

void SimpleMouseRotator::setVerticalAxis(Vector3<float>& verticalAxis) {
  this->verticalAxis = verticalAxis;
}

Vector3<float> SimpleMouseRotator::getHorizontalAxis() {
  return horizontalAxis;
}

void SimpleMouseRotator::setHorizontalAxis(Vector3<float>& horizontalAxis) {
  this->horizontalAxis = horizontalAxis;
}

void SimpleMouseRotator::beginRotate(int x, int y) {
  previousX = x;
  previousY = y;
}

void SimpleMouseRotator::endRotate() {
}

void SimpleMouseRotator::rotate(int x, int y) {
  // Account for change of original axes
  Quaternion<float> previousRotation(viewpoint->getRotation());
  previousRotation.inverse();
  Vector3<float> xAxis;
  Vector3<float> yAxis;
  QuatUtil::rotateVector(previousRotation, verticalAxis, xAxis);
  QuatUtil::rotateVector(previousRotation, horizontalAxis, yAxis);

  Quaternion<float> rotationChange(0.0f, 0.0f, 0.0f, 1.0f);
  Quaternion<float> xQuat;
  xQuat.set(yAxis, (float)toRadians(angle));
  if (previousX > x) {
    xQuat.inverse();
  }
  xQuat.normalize();
  int xDelta = std::abs(previousX - x);
  for (int i = 0; i < xDelta; ++i) {
    rotationChange.multiply(xQuat);
  }
  Quaternion<float> yQuat;
  yQuat.set(xAxis, (float)toRadians(angle));
  if (previousY > y) {
    yQuat.inverse();
  }
  yQuat.normalize();
  int yDelta = std::abs(previousY - y);
  for (int j = 0; j < yDelta; ++j) {
    rotationChange.multiply(yQuat);
  }
  rotationChange.normalize();
  previousX = x;
  previousY = y;
  Quaternion<float> r(viewpoint->getRotation());
  r.multiply(rotationChange);
  r.normalize();
  viewpoint->setRotation(r);
}

}
}
}

