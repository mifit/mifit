#ifndef mi_opengl_Circle_h
#define mi_opengl_Circle_h

namespace mi {
namespace opengl {

class Circle {

  float radius;

  int detail;

  float angle;

public:

  Circle(float radius);

  Circle(float radius, float angle);

  float getRadius();

  void setRadius(float radius);

  int getDetail();

  void setDetail(int detail);

  float getAngle();

  void setAngle(float angle);

  void render();

};

}
}

#endif
