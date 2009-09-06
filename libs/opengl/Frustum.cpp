#include <opengl/Frustum.h>

#include <opengl/Viewport.h>
#include <opengl/OpenGL.h>

using namespace mi::math;

namespace mi {
namespace opengl {

class Plane {

public:

  Vector3<float> normal;

  Vector3<float> point;

  float d;

  Plane();

  Plane(Vector3<float>& v1, Vector3<float>& v2, Vector3<float>& v3);

  void set(Vector3<float>& v1, Vector3<float>& v2, Vector3<float>& v3);

  void set(Vector3<float>& normal, Vector3<float>& point);

  void setCoefficients(float a, float b, float c, float d);
  float distance(Vector3<float>& p);

};

Plane::Plane() {
}

Plane::Plane(Vector3<float>& v1, Vector3<float>& v2, Vector3<float>& v3) {
  set(v1, v2, v3);
}

void Plane::set(Vector3<float>& v1, Vector3<float>& v2, Vector3<float>& v3) {
  Vector3<float> aux1(v1);
  aux1.subtract(v2);
  Vector3<float> aux2(v3);
  aux2.subtract(v2);

  normal.cross(aux2, aux1);

  normal.normalize();
  point.set(v2);
  d = -normal.dot(point);
}

void Plane::set(Vector3<float>& normal, Vector3<float>& point) {
  this->normal.set(normal);
  this->normal.normalize();
  d = -(this->normal.dot(point));
}

void Plane::setCoefficients(float a, float b, float c, float d) {
  normal.set(a, b, c);
  float l = normal.length();
  normal.set(a / l, b / l, c / l);
  this->d = d / l;
}

float Plane::distance(Vector3<float>& p) {
  return (d + normal.dot(p));
}

const int Frustum::TOP = 0;

const int Frustum::BOTTOM = 1;

const int Frustum::LEFT = 2;

const int Frustum::RIGHT = 3;

const int Frustum::NEARPLANE = 4;

const int Frustum::FARPLANE = 5;

Frustum::Frustum() : nearClipping(0.01f), farClipping(500.0f), fieldOfView(35.0f),
  perspective(true), frustumOffset(0.0f), picking(false) {

  planes = new Plane*[6];
  for (int i = 0; i < 6; ++i) {
    planes[i] = new Plane;
  }
  projection = new float[16];
  modelview = new float[16];
  planesMatrix = new float[16];

}

Frustum::~Frustum() {
  delete[] planesMatrix;
  delete[] modelview;
  delete[] projection;
  for (int i = 0; i < 6; ++i) {
    delete planes[i];
  }
  delete[] planes;
}

void Frustum::setFieldOfView(float fieldOfView) {
  this->fieldOfView = fieldOfView;
}

bool Frustum::isPerspective() {
  return perspective;
}

void Frustum::setPerspective(bool perspective) {
  this->perspective = perspective;
}

void Frustum::beginPicking(int x, int y) {
  picking = true;
  pickPointX = x;
  pickPointY = y;
}

void Frustum::endPicking() {
  picking = false;
}

void Frustum::render(Viewport& viewport) {
  render(viewport, frustumOffset);
}

void Frustum::render(Viewport& viewport, float offset) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  if (picking) {
    GLint vp[4];
    int* vpi = viewport.toArray();
    for (unsigned int i = 0; i < 4; ++i) {
      vp[i] = vpi[i];
    }
    gluPickMatrix(pickPointX, viewport.getHeight() - pickPointY, 5.0f, 5.0f, vp);
  }

  updateFrustum(viewport);

  if (perspective) {
    glFrustum(frustumLeft + offset, frustumRight + offset, frustumBottom, frustumTop, nearClipping, farClipping);
  } else {
    glOrtho(frustumLeft + offset, frustumRight + offset, frustumBottom, frustumTop, nearClipping, farClipping);
  }

  glMatrixMode(GL_MODELVIEW);
}

void Frustum::setHeight(float height) {
  fieldOfView = (float) toDegrees(2.0 * std::atan(height * 0.5 / focalLength));
}

void Frustum::updateFrustum(Viewport& viewport) {
  float aspect = (float) viewport.getWidth() / (float) viewport.getHeight();

  float z;
  if (perspective) {
    z = nearClipping;
  } else {
    z = focalLength;
  }
  frustumTop = getPerspectiveTop(z);
  frustumBottom = -frustumTop;
  frustumLeft = aspect * frustumBottom;
  frustumRight = aspect * frustumTop;

}

Vector3<float>* Frustum::getCorners() {

  Vector3<float>* corners = new Vector3<float>[8];
  corners[0] = Vector3<float>(frustumLeft, frustumBottom, nearClipping);
  corners[1] = Vector3<float>(frustumRight, frustumBottom, nearClipping);
  corners[2] = Vector3<float>(frustumRight, frustumTop, nearClipping);
  corners[3] = Vector3<float>(frustumLeft, frustumTop, nearClipping);

  if (perspective) {
    float aspect = (frustumRight - frustumLeft) / (frustumTop - frustumBottom);
    float top = getPerspectiveTop(farClipping);
    corners[4] = Vector3<float>(-top * aspect, -top, farClipping);
    corners[5] = Vector3<float>(top * aspect, -top, farClipping);
    corners[6] = Vector3<float>(top * aspect, top, farClipping);
    corners[7] = Vector3<float>(-top * aspect, top, farClipping);
  } else {
    corners[4] = Vector3<float>(frustumLeft, frustumBottom, farClipping);
    corners[5] = Vector3<float>(frustumRight, frustumBottom, farClipping);
    corners[6] = Vector3<float>(frustumRight, frustumTop, farClipping);
    corners[7] = Vector3<float>(frustumLeft, frustumTop, farClipping);
  }
  return corners;
}

float Frustum::getFrustumWidth() {
  return frustumRight - frustumLeft;
}

float Frustum::getFrustumHeight() {
  return frustumTop - frustumBottom;
}

float Frustum::getNearClipping() {
  return nearClipping;
}

float Frustum::getFieldOfView() {
  return fieldOfView;
}

float Frustum::getFrustumOffset() {
  return frustumOffset;
}

void Frustum::setFrustumOffset(float frustumAsymmetry) {
  this->frustumOffset = frustumAsymmetry;
}

static float getMatrixValue(float* result, int row, int col) {
  return result[col * 4 + row - 5];
}

void Frustum::updatePlanes() {

  glGetFloatv(GL_PROJECTION_MATRIX, projection);
  glGetFloatv(GL_MODELVIEW_MATRIX, modelview);

  glPushMatrix();
  glLoadMatrixf(projection);
  glMultMatrixf(modelview);
  glGetFloatv(GL_MODELVIEW_MATRIX, planesMatrix);
  glPopMatrix();

  planes[NEARPLANE]->setCoefficients(getMatrixValue(planesMatrix, 3, 1) + getMatrixValue(planesMatrix, 4, 1),
    getMatrixValue(planesMatrix, 3, 2) + getMatrixValue(planesMatrix, 4, 2), getMatrixValue(planesMatrix, 3, 3)
    + getMatrixValue(planesMatrix, 4, 3), getMatrixValue(planesMatrix, 3, 4)
    + getMatrixValue(planesMatrix, 4, 4));
  planes[FARPLANE]->setCoefficients(-getMatrixValue(planesMatrix, 3, 1) + getMatrixValue(planesMatrix, 4, 1),
    -getMatrixValue(planesMatrix, 3, 2) + getMatrixValue(planesMatrix, 4, 2), -getMatrixValue(planesMatrix, 3, 3)
    + getMatrixValue(planesMatrix, 4, 3), -getMatrixValue(planesMatrix, 3, 4)
    + getMatrixValue(planesMatrix, 4, 4));
  planes[BOTTOM]->setCoefficients(getMatrixValue(planesMatrix, 2, 1) + getMatrixValue(planesMatrix, 4, 1),
    getMatrixValue(planesMatrix, 2, 2) + getMatrixValue(planesMatrix, 4, 2), getMatrixValue(planesMatrix, 2, 3)
    + getMatrixValue(planesMatrix, 4, 3), getMatrixValue(planesMatrix, 2, 4)
    + getMatrixValue(planesMatrix, 4, 4));
  planes[TOP]->setCoefficients(-getMatrixValue(planesMatrix, 2, 1) + getMatrixValue(planesMatrix, 4, 1),
    -getMatrixValue(planesMatrix, 2, 2) + getMatrixValue(planesMatrix, 4, 2), -getMatrixValue(planesMatrix, 2, 3)
    + getMatrixValue(planesMatrix, 4, 3), -getMatrixValue(planesMatrix, 2, 4)
    + getMatrixValue(planesMatrix, 4, 4));
  planes[LEFT]->setCoefficients(getMatrixValue(planesMatrix, 1, 1) + getMatrixValue(planesMatrix, 4, 1),
    getMatrixValue(planesMatrix, 1, 2) + getMatrixValue(planesMatrix, 4, 2), getMatrixValue(planesMatrix, 1, 3)
    + getMatrixValue(planesMatrix, 4, 3), getMatrixValue(planesMatrix, 1, 4)
    + getMatrixValue(planesMatrix, 4, 4));
  planes[RIGHT]->setCoefficients(-getMatrixValue(planesMatrix, 1, 1) + getMatrixValue(planesMatrix, 4, 1),
    -getMatrixValue(planesMatrix, 1, 2) + getMatrixValue(planesMatrix, 4, 2), -getMatrixValue(planesMatrix, 1, 3)
    + getMatrixValue(planesMatrix, 4, 3), -getMatrixValue(planesMatrix, 1, 4)
    + getMatrixValue(planesMatrix, 4, 4));

}

Frustum::CullingResult Frustum::pointInFrustum(Vector3<float>& position) {

  CullingResult result = INSIDE;
  for (int i = 0; i < 6; i++) {

    if (planes[i]->distance(position) < 0) {
      return OUTSIDE;
    }
  }
  return result;

}

Frustum::CullingResult Frustum::sphereInFrustum(Vector3<float>& position, float radius) {

  CullingResult result = INSIDE;
  float distance;

  for (int i = 0; i < 6; i++) {
    distance = planes[i]->distance(position);
    if (distance < -radius) {
      return OUTSIDE;
    } else if (distance < radius) {
      result = INTERSECT;
    }
  }
  return result;

}

float Frustum::getPerspectiveHeight(float z) {
  return 2.0f * getPerspectiveTop(z);
}

float Frustum::getPerspectiveTop(float z) {
  return z * (float) std::tan(toRadians(fieldOfView) * 0.5);
}

float Frustum::getFrustumHeight(float z) {
  float heightAtZ;
  if (perspective) {
    heightAtZ = getPerspectiveHeight(z);
  } else {
    heightAtZ = frustumTop - frustumBottom;
  }
  return heightAtZ;
}

void Frustum::setNearClipping(float nearClipping) {
  this->nearClipping = nearClipping;
}

float Frustum::getFarClipping() {
  return farClipping;
}

void Frustum::setFarClipping(float farClipping) {
  this->farClipping = farClipping;
}

float Frustum::getFocalLength() {
  return focalLength;
}

void Frustum::setFocalLength(float focalLength) {
  this->focalLength = focalLength;
}

float Frustum::getHeight() {
  return frustumTop - frustumBottom;
}

float Frustum::getWidth() {
  return frustumRight - frustumLeft;
}

Plane** Frustum::getPlanes() {
  return planes;
}

}
}
