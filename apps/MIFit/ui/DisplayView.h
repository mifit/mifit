#ifndef mifit_ui_DisplayView_h
#define mifit_ui_DisplayView_h

#include "ViewSyncedPanel.h"

class MIGLWidget;

class DisplayView : public ViewSyncedPanel {

protected:

  virtual QWidget* createPanelForView(MIGLWidget* view, QWidget* parent);
  virtual void destroyContentsForView(MIGLWidget* view, QWidget* panel);

public:

  DisplayView(QWidget* parent);
  virtual ~DisplayView();

};

#endif
