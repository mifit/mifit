#ifndef MIFIT_UI_DICTEDITATOMPICKINGRENDERABLE_H_
#define MIFIT_UI_DICTEDITATOMPICKINGRENDERABLE_H_

#include <opengl/Renderable.h>

#include <chemlib/chemlib.h>

class Displaylist;
class GLRenderer;

namespace mi {
namespace opengl {
class Frustum;
class StereoView;
}
}

class DictEditAtomRenderable;

class DictEditAtomPickingRenderable : public mi::opengl::Renderable {

  mi::opengl::StereoView* stereoView;

  GLRenderer* renderer;

  DictEditAtomRenderable* atomRenderable;

public:

  DictEditAtomPickingRenderable(mi::opengl::StereoView* stereoView, Displaylist* displaylist, mi::opengl::Frustum* frustum);
  virtual ~DictEditAtomPickingRenderable();

  virtual void render();

  chemlib::MIAtom* getAtom(int id);

};

class DictEditAtomRenderable : public mi::opengl::Renderable {

  Displaylist* displaylist;

  GLRenderer* renderer;

public:

  DictEditAtomRenderable(GLRenderer* renderer, Displaylist* displaylist);

  virtual void render();

};

#endif
