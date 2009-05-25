#include "ViewSyncedPanel.h"

#include <boost/bind.hpp>
#include <boost/signal.hpp>

#include <QVBoxLayout>

#include "MIGLWidget.h"


ViewSyncedPanel::ViewSyncedPanel(QWidget* parent) : QWidget(parent), currentPanel(NULL) {

  layout = new QVBoxLayout(this);
  setLayout(layout);

  createDefaultNullPanel();
  setViewPanel(NULL);

  MIGLWidget::viewActivated.connect(boost::bind(&ViewSyncedPanel::viewActivated, this, _1));
  MIGLWidget::viewDeactivated.connect(boost::bind(&ViewSyncedPanel::viewDeactivated, this, _1));

}

ViewSyncedPanel::~ViewSyncedPanel() {
}

void ViewSyncedPanel::createDefaultNullPanel() {
  QWidget* panel = new QWidget(this);
  viewToPanel[NULL] = panel;
  layout->addWidget(panel);
}

void ViewSyncedPanel::createViewPanel(MIGLWidget* view) {
  QWidget* panel = createPanelForView(view, this);

  viewToPanel[view] = panel;
  layout->addWidget(panel);
}

void ViewSyncedPanel::viewActivated(MIGLWidget* view) {
  setViewPanel(view);
}

void ViewSyncedPanel::viewDeactivated(MIGLWidget* view) {
  setViewPanel(NULL);
  if (panelForViewExists(view)) {
    QWidget* viewPanel = viewToPanel[view];
    destroyContentsForView(view, viewPanel);
    viewToPanel.erase(view);
    layout->removeWidget(viewPanel);
    delete viewPanel;
  }
}

void ViewSyncedPanel::setViewPanel(MIGLWidget* view) {
  if (currentPanel != NULL) {
    currentPanel->setVisible(false);
  }
  if (!panelForViewExists(view)) {
    createViewPanel(view);
  }
  currentPanel = viewToPanel[view];
  currentPanel->setVisible(true);
  //layout->update();
}

bool ViewSyncedPanel::panelForViewExists(MIGLWidget* view) {
  return viewToPanel.find(view) != viewToPanel.end();
}

