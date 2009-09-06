#include <opengl/interact/ArcBallFeedback.h>

#include <opengl/Arc.h>
#include <opengl/Circle.h>
#include <opengl/OpenGL.h>

using namespace mi::math;
using namespace mi::opengl;

namespace mi {
namespace opengl {
namespace interact {

ArcBallFeedback::ArcBallFeedback(float radius, float* color) {
  arc = new Arc();
  arcColor.set(color[0], color[1], color[2]);
  ballColor.set(arcColor);
  // Dim the ball circle
  ballColor.scale(0.8f);
  ball = new Circle(radius);
}

void ArcBallFeedback::render() {

  glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT | GL_COLOR_BUFFER_BIT);
  glShadeModel(GL_FLAT);
  glDisable(GL_FOG);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);

  glColor3f(ballColor.x, ballColor.y, ballColor.z);
  ball->render();
  glColor3f(arcColor.x, arcColor.y, arcColor.z);
  arc->render();

  glPopAttrib();

}

void ArcBallFeedback::setRadius(float radius) {
  ball->setRadius(radius);
}

void ArcBallFeedback::setFrom(Vector3<float>& from) {
  Vector3<float> v(from);
  v.scale(ball->getRadius());
  arc->setFrom(v);
}

void ArcBallFeedback::setTo(Vector3<float>& to) {
  Vector3<float> v(to);
  v.scale(ball->getRadius());
  arc->setTo(v);
}

}
}
}

