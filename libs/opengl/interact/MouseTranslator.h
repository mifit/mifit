#ifndef mi_opengl_interact_MouseTranslator_h
#define mi_opengl_interact_MouseTranslator_h

#include <math/Vector3.h>

namespace mi {
namespace opengl {

class Viewpoint;
class Axis;

namespace interact {

class MouseTranslator {

  Viewpoint* viewpoint;

  Viewpoint* relativeToViewpoint;

  int previousX;

  int previousY;

  mi::math::Vector3<float>* startPosition;

  mi::math::Vector3<float>* currentPosition;

  float scaling;

  Axis* axes;

  Axis* currentAxes;

public:

  MouseTranslator(Viewpoint* viewpoint, Viewpoint* relativeToViewpoint, float scaling);

  ~MouseTranslator();

  float getScaling();

  void setScaling(float angle);

  mi::math::Vector3<float>* getStartPosition();

  mi::math::Vector3<float>* getCurrentPosition();

  void setAxes(const Axis& axes);

  void beginTranslate(int x, int y);

  void endTranslate();

  void translate(int x, int y);

};

}
}
}

#endif
