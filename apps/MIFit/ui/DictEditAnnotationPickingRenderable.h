#ifndef MIFIT_UI_DICTEDITANNOTATIONPICKINGRENDERABLE_H_
#define MIFIT_UI_DICTEDITANNOTATIONPICKINGRENDERABLE_H_

#include <opengl/Renderable.h>
#include <list>


class Annotation;
class GLRenderer;
class Displaylist;
class DictEditCanvas;

namespace mi {
namespace opengl {
class Frustum;
class StereoView;
}
}

class DictEditAnnotationRenderable;

class DictEditAnnotationPickingRenderable : public mi::opengl::Renderable {

  mi::opengl::StereoView* stereoView;

  GLRenderer* renderer;

  DictEditAnnotationRenderable* annotationRenderable;

public:

  DictEditAnnotationPickingRenderable(mi::opengl::StereoView* stereoView, Displaylist* displaylist, mi::opengl::Frustum* frustum, DictEditCanvas* canvas);
  virtual ~DictEditAnnotationPickingRenderable();

  virtual void render();

  Annotation* getAnnotation(int id);

  void updateTextScale(float glUnitsPerPixel);

};

class DictEditAnnotationRenderable : public mi::opengl::Renderable {

  Displaylist* displaylist;

  GLRenderer* renderer;

public:

  DictEditAnnotationRenderable(GLRenderer* renderer, Displaylist* displaylist);

  virtual void render();

};

#endif
