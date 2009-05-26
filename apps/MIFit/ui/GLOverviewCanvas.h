#ifndef MIFIT_UI_GLOVERVIEWCANVAS_H_
#define MIFIT_UI_GLOVERVIEWCANVAS_H_

#include <math/Vector3.h>
#include <vector>

namespace chemlib {
class Bond;
}

namespace mi {
namespace opengl {
class Camera;
class Frustum;
class Viewport;
}
}

#include <QGLWidget>

class GLOverviewCanvas : public QGLWidget {
private:

  /**
   * Controls the size of the plus marking the center of the viewpoint.
   */
  float viewpointCenterMarkerRadius;

  /**
   * Controls the size of the circle marking the selected residue.
   */
  float selectedResidueMarkerRadius;

  /**
   * Controls the resolution of the circle marking the selected residue.
   */
  float selectedResidueMarkerResolution;


  /**
   * Draws a circle which marks the selected residue at the given position.
   * @param x the x position of the marker
   * @param y the y position of the marker
   * @param z the z position of the marker
   */
  void drawSelectedResidueMarker(float x, float y, float z);

  /**
   * Adjusts color for depth.
   */
  float* getColor(int colorIndex);

  /**
   * Sets a rotation matrix (float[16]) to the identity matrix.
   */
  void setIdentity(float rotation[]);

  /**
   * Sets the GLCanvas to the current OpenGL context.
   */
  void setCurrent();

  std::vector<chemlib::Bond> bonds;

  mi::opengl::Viewport* viewport;

  mi::opengl::Frustum* frustum;

  mi::opengl::Camera* camera;

  mi::math::Vector3<float> min;

  mi::math::Vector3<float> max;

  mi::math::Vector3<float> center;

  mi::math::Vector3<float> cameraOffset;

  void createAlphaCarbonTrace();


public:
  GLOverviewCanvas(QWidget *parent);
  virtual ~GLOverviewCanvas();

  void initializeGL();
  void resizeGL(int width, int height);
  void paintGL();

  virtual QSize sizeHint() const;

  /**
   * Processes mouse events by forwarding them on to the view.
   */
  void handleMouseEvent(QMouseEvent *e);
  void mousePressEvent(QMouseEvent *e);
  void mouseReleaseEvent(QMouseEvent *e);
  void mouseMoveEvent(QMouseEvent *e);
  void mouseDoubleClickEvent(QMouseEvent *e);

  /**
   * Processes key events by forwarding them on to the view.
   */
  void keyPressEvent(QKeyEvent *e);
};

#endif /*MIFIT_UI_GLNAVIGATORCANVAS_H_*/
